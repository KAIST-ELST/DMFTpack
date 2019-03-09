//#include <iostream>
//#include <stdio.h>
//#include <stdlib.h>
//#include <time.h>
#include "mpi.h"
#include <fstream>
#include "ImgTimeFtn.h"
#include "tight_common.h"
#include "model.h"
#include "TB.h"
#include <Eigen/Eigenvalues>

void get_preNAOs(Eigen::MatrixXcd weightMatrix,
                 Eigen::MatrixXcd & principal_number_tfm, Eigen::MatrixXcd & weightMatrix_preNAOs);
void SpreadFtn( int knum,
                std::vector<Eigen::MatrixXcd> S_overlap, std::vector<Eigen::MatrixXcd> &  transformMatrix_k, std::vector<int> & accumulated_Num_SpinOrbital   );

void direct_projection( Eigen::MatrixXcd & Hk, Eigen::MatrixXcd & Sk, Eigen::MatrixXcd  & evec,  Eigen::VectorXd & eval  ) ;
void lowdin_symmetric_orthogonalization( Eigen::MatrixXcd & Hk, Eigen::MatrixXcd & Sk, Eigen::MatrixXcd  & evec,  Eigen::VectorXd & eval, Eigen::MatrixXcd & transformMatrix    );
//void naturalAtomicOrbitals_population_weighted_symmetric_orthogonalization_k(Eigen::MatrixXcd & Hk, Eigen::MatrixXcd & Sk, Eigen::MatrixXcd  & evec,  Eigen::VectorXd & eval,int NumOrbit, double muDFT, Eigen::MatrixXcd & transformMatrix) ;
void naturalAtomicOrbitals_population_weighted_symmetric_orthogonalization_r(
    Eigen::MatrixXcd & Hk, Eigen::MatrixXcd & Sk, int kpoint,  std::vector<int> & accumulated_Num_SpinOrbital,
    Eigen::MatrixXcd  & evec,  Eigen::VectorXd & eval, Eigen::MatrixXcd weightMatrix, Eigen::MatrixXcd & transformMatri, Eigen::MatrixXcd principal_number_tfm);



Eigen::MatrixXcd getweight_dual(  std::vector<Eigen::MatrixXcd> S_overlap, std::vector<Eigen::MatrixXcd> dual_DM_direct ) ;
Eigen::MatrixXcd getweight_hyb(  std::vector<Eigen::MatrixXcd> S_overlap, std::vector<Eigen::MatrixXcd> dual_DM_direct ) ;
Eigen::MatrixXcd getweight_direct(std::vector<Eigen::MatrixXcd> S_overlap,std::vector<Eigen::MatrixXcd> dual_DM_direct ) ;
//Eigen::MatrixXcd getweight_sc(    std::vector<Eigen::MatrixXcd> S_overlap, std::vector<Eigen::MatrixXcd>  dual_DM_direct,  Eigen::MatrixXcd  weightMatrix, int * accumulated_Num_SpinOrbital   ) ;

int read_OverlapMat(Eigen::MatrixXi &S_overlap_Rindex, Eigen::VectorXcd  &S_overlap_RMatrix, const std::string &hamiltonian, std::vector<int> &  accumulated_Num_SpinOrbital);
int read_DensityMat(Eigen::MatrixXi &NumMat_Rindex,         Eigen::VectorXcd  &NumMat_RMatrix, const std::string &hamiltonian, std::vector<int> & accumulated_Num_SpinOrbital) ;

//void read_HR(Eigen::MatrixXi &H_Rindex, Eigen::VectorXcd  &H_RMatrix, const std::string &hamiltonian, int * accumulated_Num_SpinOrbital, int NumOrbit);



int Construct_Hk_Sk(
    int knum, int knum_mpiGlobal,   Eigen::MatrixXi  H_Rindex, Eigen::VectorXcd H_RMatrix, double ** kmesh, std::vector<int> & accumulated_Num_SpinOrbital,
    std::vector<Eigen::MatrixXcd> & H_k_inModelSpace,  std::vector<Eigen::MatrixXcd> & S_overlap,
    std::vector<Eigen::MatrixXcd> & KS_eigenVectors_orthoBasis, Eigen::VectorXd  * KS_eigenEnergy

)
{


    int temp;
    /*Read Overlap matrix*/
    Eigen::MatrixXi S_overlap_Rindex    ;
    Eigen::VectorXcd S_overlap_RMatrix  ;


    int    overlap_exist = read_OverlapMat(S_overlap_Rindex, S_overlap_RMatrix, std::string("OverlapMatrix.HWR"), accumulated_Num_SpinOrbital);


////////////////////////////////////////////////////////////
//Construct low energy Hamiltonian from  Hk, overlapmatrix, S(k)*/
////////////////////////////////////////////////////////////
    /* Calculate H(k) and S(k)   */
    /* Convention : |R>=1/sqrt(N_k) \sum_k |k> exp(-ikR) => |k>  = 1/sqrt(N_k) \sum_R |R> exp(ikR) */
    /* Convention : H(R) = 1/N_k sum_k H(k) exp(ikR)*/
    /* Convention : H(k) = sum_R H(R) exp(-ikR)*/
    /* Convention : S(k) = sum_R S(R) exp(-ikR)*/
    std::cout.precision(5);

    /*DFT_KS hamiltonian in k-space*/
    FromValToKS.resize(knum);
    for(int k = 0;  k < knum; k++) {
//        std::cout <<k<<"aa\n";
        H_k_inModelSpace[k].setZero(NumOrbit, NumOrbit);
        S_overlap[k].setIdentity(NumOrbit, NumOrbit);
        for(int indx=0; indx<H_RMatrix.size(); indx++) {
            H_k_inModelSpace[k](H_Rindex(indx,3),H_Rindex(indx,4))
            += H_RMatrix(indx)* exp ( -I*( (kmesh[k][0]*ax*H_Rindex(indx,0))+(kmesh[k][1]*ay*H_Rindex(indx,1))+(kmesh[k][2]*az*H_Rindex(indx,2))) )  ;

//            std::cout << H_Rindex(indx,0) <<" " << H_Rindex(indx, 1)  <<" " << H_Rindex(indx, 2)<<"\n";
//            std::cout << H_Rindex(indx,3) <<" " << H_Rindex(indx, 4) <<"\n";
//            std::cout << H_RMatrix(indx) <<"\n";
        }
//        std::cout <<k<<"bb\n";
//std::cout << H_k_inModelSpace[k] <<"\n";
        /*S(k),   H*S |\psi> = E |\psi>      */
        if (overlap_exist==1) {
            S_overlap[k].setZero(NumOrbit,NumOrbit);
            for(int indx=0; indx<S_overlap_RMatrix.size(); indx++) {
                S_overlap[k](S_overlap_Rindex(indx,3)*2, S_overlap_Rindex(indx,4)*2) +=
                    S_overlap_RMatrix(indx)* exp ( -I*( (kmesh[k][0]*ax*S_overlap_Rindex(indx,0))+
                                                        (kmesh[k][1]*ay*S_overlap_Rindex(indx,1))+
                                                        (kmesh[k][2]*az*S_overlap_Rindex(indx,2))) )  ;


            }
//            std::cout <<k<<"cc\n";
            int Norbital_spatial = NumOrbit/Spin_DegreeOfFreedom;
            for (int i0=0; i0<Norbital_spatial; i0++) {
                for (int m0=0; m0<Norbital_spatial; m0++) {
                    S_overlap[k]((2*i0+1), (2*m0+1)) = S_overlap[k](2*i0,2*m0);
                }
            }
        }//if OverlapMatrix
    }//k
    if(mpi_rank==0 )   std::cout << "<TB> Construct Hk from HR :\n";


    for(int k = 0;  k < knum; k++) {

        Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXcd> ces1(NumOrbit);
        ces1.compute( H_k_inModelSpace[k], S_overlap[k] );       //  HS \psi = S  \psi E

        KS_eigenEnergy[k] = ces1.eigenvalues();
        KS_eigenVectors_orthoBasis[k] = ces1.eigenvectors();

    }

    return overlap_exist;
}






Eigen::MatrixXcd  weightMatrix;

void  ConstructModelHamiltonian
(
    int knum, int knum_mpiGlobal, double ** kmesh, std::vector<int> & accumulated_Num_SpinOrbital,
    std::vector<Eigen::MatrixXcd> & H_k_inModelSpace,  std::vector<Eigen::MatrixXcd> & S_overlap,
    std::vector<Eigen::MatrixXcd> & KS_eigenVectors_orthoBasis, std::vector<Eigen::MatrixXcd> & transformMatrix_k, Eigen::VectorXd  * KS_eigenEnergy,  int overlap_exist,double muDFT
)
{



    std::vector<Eigen::MatrixXcd>  dual_DM_direct;
    dual_DM_direct.resize(knum);
    Eigen::MatrixXi   NumMat_Rindex;
    Eigen::VectorXcd  NumMat_RMatrix;

//    /*
    //read density Matrix
    int temp =          read_DensityMat(NumMat_Rindex, NumMat_RMatrix, std::string("dual_DM_direct.HWR"), accumulated_Num_SpinOrbital);
    if (temp ==1 and (
                localOrbitalType.find(std::string("recip_F")) != std::string::npos or
                localOrbitalType.find(std::string("direct_F")) != std::string::npos  or
                localOrbitalType.find(std::string("hyb_F")) != std::string::npos
            )
       ) {
        double diffnorm=0;
        double diffnorm_diag=0;
        for(int k = 0;  k < knum; k++) {
            dual_DM_direct[k].setZero(NumOrbit, NumOrbit);
            for(int indx=0; indx<NumMat_RMatrix.size(); indx++) {
                int n=  NumMat_Rindex(indx,0) ;
                int l=  NumMat_Rindex(indx,1) ;
                int m=  NumMat_Rindex(indx,2) ;
                int i0= NumMat_Rindex(indx,3) ;
                int m0= NumMat_Rindex(indx,4) ;
                dual_DM_direct[k](i0,m0) += NumMat_RMatrix(indx) * exp ( -I*( (kmesh[k][0]*ax*n)+(kmesh[k][1]*ay*l)+(kmesh[k][2]*az*m)) )  ;
            }
        }
    }
    else {
        for(int k = 0;  k < knum; k++) {
            dual_DM_direct[k].setZero(NumOrbit, NumOrbit);

            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> ces1(S_overlap[k]);
            Eigen::MatrixXcd S_invsq =  (  ces1.operatorInverseSqrt() );
            Eigen::MatrixXcd S_sq = S_invsq.inverse();
            Eigen::MatrixXcd Hk = S_invsq * H_k_inModelSpace[k] * S_invsq;

            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> ces(Hk);
            ces.compute(Hk);

            for (int n=0; n<NumOrbit; n++) {
                dual_DM_direct[k](n,n) =  1./(1+std::exp(beta*(ces.eigenvalues()(n)-muDFT))) ;
            }

            dual_DM_direct[k] = (ces.eigenvectors() * dual_DM_direct[k] * (ces.eigenvectors()).adjoint()).eval();
            dual_DM_direct[k] =  (S_invsq * dual_DM_direct[k] * S_sq).eval();


        }//k
    }
    ifroot std::cout <<"Now, we have occupation matrix\n";
    //    */



    if (overlap_exist==1) {

//        /*
        Eigen::MatrixXcd occ_temp(NumOrbit, NumOrbit);
        Eigen::MatrixXcd      occ(NumOrbit, NumOrbit);

        occ_temp.setZero(NumOrbit,NumOrbit);
        occ.setZero(NumOrbit,NumOrbit);


        for(int k = 0;  k < knum; k++) {
            Eigen::MatrixXcd temp = 0.5*(dual_DM_direct[k]);
            occ_temp +=  ((temp + temp.adjoint()))/knum_mpiGlobal;
        }//k
        MPI_Allreduce(occ_temp.data(), occ.data(), occ.size(), MPI_DOUBLE_COMPLEX, MPI_SUM,  MPI_COMM_WORLD);
        if(mpi_rank==0 )   std::cout << "<TB> Mulliken populations:\n";
        ifroot{
            double sum =0;
            for(int at=0; at<NumAtom; at++) {
                for(int h1=accumulated_Num_SpinOrbital[at]; h1<accumulated_Num_SpinOrbital[at+1]; h1++) {
                    if(isOrbitalHartrDFT[h1]) std::cout << at+1<< " " << h1-accumulated_Num_SpinOrbital[at]+1 <<" " <<  occ(h1,h1) <<"\n";
                    sum+=real(occ(h1,h1));
                }
            }
            std::cout << "sum_rule (Mulliken):" <<sum<<"\n";
        }




//orthogonalization
        if ( localOrbitalType.find(std::string("nao")) != std::string::npos ) {
            for(int k = 0;  k < knum; k++) {
                for (int i=0; i<NumOrbit; i+=2) {
                    for (int j=0; j<NumOrbit; j+=2) {
                        cmplx temp = ( dual_DM_direct[k](i,j) + dual_DM_direct[k](i+1,j+1) )/2.0;
                        dual_DM_direct[k](i,j   ) = temp;
                        dual_DM_direct[k](i+1,j+1) = temp;
                        dual_DM_direct[k](i+1,j   ) = 0.0;
                        dual_DM_direct[k](i,j+1   ) = 0.0;
                    }
                }//n
            }//k

            static bool we_have_weightMatrix = false;
            if(we_have_weightMatrix==false or
                    localOrbitalType.find(std::string("recip_F")) != std::string::npos or
                    localOrbitalType.find(std::string("direct_F")) != std::string::npos  or
                    localOrbitalType.find(std::string("hyb_F")) != std::string::npos) {


                //get weightMatrix
                if ( localOrbitalType.find(std::string("nao_recip")) != std::string::npos or localOrbitalType.find(std::string("nao_compl")) != std::string::npos) { //complement
                    weightMatrix = getweight_dual(S_overlap, dual_DM_direct);
                }
                else if (localOrbitalType.find(std::string("nao_hyb")) != std::string::npos) {
                    weightMatrix = getweight_hyb(S_overlap,dual_DM_direct);
                }
                else if (localOrbitalType.find(std::string("nao_direct")) != std::string::npos) {
                    weightMatrix = getweight_direct(S_overlap,dual_DM_direct);
                }
//                else if ( localOrbitalType.find(std::string("nao_sc")) != std::string::npos ) {
//                    //initial condition for  _sc
//                    if(we_have_weightMatrix ==false) {
//                        weightMatrix.setIdentity(NumOrbit, NumOrbit);
//                        weightMatrix  *= NumberOfElectron/NumOrbit;
//                    }
//                    weightMatrix = getweight_sc( S_overlap, dual_DM_direct, weightMatrix, accumulated_Num_SpinOrbital);
//                }
                else {
                    std::cout << "Please check input file, nao_OOO\n";
                    exit(1);
                }
                we_have_weightMatrix = true;

                ifroot{
                    std::cout << "<TB> weightMatrix:\n";
                    for(int at=0; at<NumAtom; at++) {
                        for(int h1=accumulated_Num_SpinOrbital[at]; h1<accumulated_Num_SpinOrbital[at+1]; h1++) {
                            if(isOrbitalHartrDFT[h1]) std::cout << at+1<< " " << h1-accumulated_Num_SpinOrbital[at]+1 <<" " <<  weightMatrix(h1,h1) <<"\n";
                        }
                    }
                }
            }//weightMatrix cal

            Eigen::MatrixXcd weightMatrix_preNAOs;
            Eigen::MatrixXcd  principal_number_tfm;
            get_preNAOs(weightMatrix, principal_number_tfm, weightMatrix_preNAOs);
            ifroot std::cout << "we have pre-NAOs\n";

            for(int k = 0;  k < knum; k++) {
                naturalAtomicOrbitals_population_weighted_symmetric_orthogonalization_r(
                    H_k_inModelSpace[k],S_overlap[k], k, accumulated_Num_SpinOrbital,
                    KS_eigenVectors_orthoBasis[k], KS_eigenEnergy[k],
                    weightMatrix_preNAOs, transformMatrix_k[k], principal_number_tfm   );
            }
        }//nao
        else if (localOrbitalType == std::string("lowdin")  ) {
            for(int k = 0;  k < knum; k++) {
                lowdin_symmetric_orthogonalization(H_k_inModelSpace[k], S_overlap[k], KS_eigenVectors_orthoBasis[k], KS_eigenEnergy[k], transformMatrix_k[k]);
            }
        }
        else if (localOrbitalType == std::string("direct")  ) {
            for(int k = 0;  k < knum; k++) {
                direct_projection(H_k_inModelSpace[k], S_overlap[k], KS_eigenVectors_orthoBasis[k], KS_eigenEnergy[k]);
            }
        }
//        SpreadFtn( knum,  S_overlap, transformMatrix_k, accumulated_Num_SpinOrbital);
    }//overlap
//    else {
//        for(int k = 0;  k < knum; k++) {
//            /* eigen problem  For DFT Hamiltonian,   */
//            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> ces(NumOrbit);
//            ces.compute(H_k_inModelSpace[k]);
//            KS_eigenVectors_orthoBasis[k] = ces.eigenvectors();
//            for(int i =0; i<NumOrbit; i++) {
//                KS_eigenEnergy[k][i] = (ces.eigenvalues()[i]) ;  //(V[0][0], V[1][0], V[2][0],..V[n][0]) = 0th eigenvector
//            }
//        }
//    }
    /*info:spread function*/
    ifroot    std::cout << "Hk was constructed..\n"  ;




    /*Find NBAND and FromValToKS with given KS_eigenEnergy and muDFT*/
    ///////////////////////////////////////////////////////////////////
    double  tempd2=0;

    Eigen::VectorXd  HartreWeightInWindows_local;
    Eigen::VectorXd  HartreWeightInWindows_global;
    HartreWeightInWindows_local.setZero(N_peratom_HartrOrbit);
    HartreWeightInWindows_global.setZero(N_peratom_HartrOrbit);

    for(int k=0; k< knum ; k++) {
        if (downfolding==1) {
            /*Check d-orbital character in  energy window*/
            for (int m=0; m<N_peratom_HartrOrbit; m++) {
                int m0 = HartrIndex_inDFT[m];
                for (int i0=0; i0<NumOrbit; i0++) {
                    tempd2 += std::pow(std::abs(KS_eigenVectors_orthoBasis[k](m0,i0)),2);
                    if( (KS_eigenEnergy[k][i0]-muDFT) < lower_model_window or KS_eigenEnergy[k][i0]-muDFT > upper_model_window  ) {
                        HartreWeightInWindows_local[m]+=std::pow(std::abs(KS_eigenVectors_orthoBasis[k](m0,i0)),2);
                    }
                }
            }
        }
    }//k
    double tempd2GL=0;
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allreduce(&tempd2, &tempd2GL, 1, MPI_DOUBLE, MPI_SUM,  MPI_COMM_WORLD);
    ifroot std::cout << "d-orbital normalization       : "  << tempd2GL/(knum_mpiGlobal*N_peratom_HartrOrbit)<<"\n";
    MPI_Allreduce(HartreWeightInWindows_local.data(), HartreWeightInWindows_global.data(), HartreWeightInWindows_global.size(), MPI_DOUBLE, MPI_SUM,  MPI_COMM_WORLD);
    for (int m=0; m<N_peratom_HartrOrbit; m++) {
        ifroot std::cout << "d-orbital in rest energy space: "  << HartreWeightInWindows_global[m] /(knum_mpiGlobal) <<"\n";
    }

}//ConstructModel








/*
void naturalAtomicOrbitals_population_weighted_symmetric_orthogonalization_k(Eigen::MatrixXcd & Hk, Eigen::MatrixXcd & Sk, Eigen::MatrixXcd  & evec,  Eigen::VectorXd & eval,int NumOrbit, double muDFT , Eigen::MatrixXcd &transformMatrix ) {

    Eigen::MatrixXcd KS_evec_direct, weightMatrix;
    weightMatrix.setZero(NumOrbit, NumOrbit);
    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXcd> ces1(NumOrbit);
    ces1.compute( Hk, Sk );       //  H  \psi = E S  \psi
    evec = ces1.eigenvectors();
    eval = ces1.eigenvalues();

    KS_evec_direct = (Sk * ces1.eigenvectors()); //full,     KS_evec_direct[k][i][n] = <\psi_i | kn>
    for (int i=0; i<NumOrbit; i++) {
        for (int n=0; n<NumOrbit; n++) {
            weightMatrix(i,i) += 1./(1+std::exp(beta*(eval(n)-muDFT)))  * std::pow(std::abs(KS_evec_direct(i,n)),2);
        }
    }

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> ces2( (weightMatrix*Sk*weightMatrix) );
    transformMatrix =     weightMatrix * ces2.operatorInverseSqrt();
    Hk = transformMatrix.adjoint() * Hk * transformMatrix;

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> ces3( Hk );
    ces3.compute(Hk);
    evec = ces3.eigenvectors();
    eval = ces3.eigenvalues();
}
*/






