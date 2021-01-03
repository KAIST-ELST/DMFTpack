//#include <iostream>
//#include <stdio.h>
//#include <stdlib.h>
#include "mpi.h"
#include "ImgTimeFtn.h"
#include "tight_common.h"
#include "model.h"
#include "TB.h"
#include <Eigen/Eigenvalues>
void  downfolding_ftn
(
    int knum,  int knum_mpiGlobal,
    std::vector<int> & NBAND,   std::vector<Eigen::MatrixXcd> &   H_k_inModelSpace,
    std::vector<Eigen::MatrixXcd> & KS_eigenVectors_orthoBasis, Eigen::VectorXd  * KS_eigenEnergy,
    double muDFT
) {
    ifroot std::cout <<"<downfolding>\n";




    if (downfolding==1) {
        for(int k = 0;  k < knum; k++) {
            /*project out rest space ;*/
            /*d-orbital space extraction*/
            Eigen::MatrixXcd  S_overlap_dSpace;
            DF_CorrBase[k].setZero( N_peratom_HartrOrbit*NumCorrAtom, NBAND[k]); // At the final stage, DF_CorrBase[k](NBAND, NBAND) = <downfolded orbitals | nk >
            int p00=0;
            for(int cl1=0; cl1<NumCluster; cl1++) {
//                for(int p0=HartrRange_DFT[at1][0] ; p0<HartrRange_DFT[at1][1] ; p0++) {}
                for(int p0_=0; p0_<NumHartrOrbit_per_cluster ; p0_++) {
                    int p0 = HartrIndex_inDFT[cl1*NumHartrOrbit_per_cluster+p0_];
                    for(int i1=0; i1<NBAND[k]; i1++) {
                        /*DF_CorrBase = <\phi_\alpha in model | nk > ( size = dim(CorrSpace) * NBAND ) */
                        DF_CorrBase[k](p00,i1) =  KS_eigenVectors_orthoBasis[k](p0,FromValToKS.at(k).at(i1)) ;
                    }
                    p00++;
                }
            }
            S_overlap_dSpace = DF_CorrBase[k] * DF_CorrBase[k].adjoint();    //dim=Cor*Cor    <Pc|Pc'>. Here c= correlated orbital, P= low-energy projector

            /*Orthogonalize d-space*/
            Eigen::MatrixXcd    temp;
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> ces(S_overlap_dSpace);
            temp =  ( DF_CorrBase[k].adjoint() * ces.operatorInverseSqrt() ).adjoint();    //(size=dim(CorrSpace)*NBAND)
            DF_CorrBase[k] = temp;   /*  <\phi_\alpha in model | nk >   with orthogonalized |Pc> */
            if(mpi_rank==0 and k==0) std::cout << "We have model projected-d orbital\n";

            /*QR decomposition; find d+p space in orthonomal base*/
            temp = DF_CorrBase[k];           //(size=dim(CorrSpace)*NBAND)

            DF_CorrBase[k] =((DF_CorrBase[k].adjoint()).householderQr().householderQ()).adjoint();       //dim = recovered to NBAND * NBAND, QR-decompose
            if(mpi_rank==0 and k==0) std::cout << "qr-decompose done\n";
            for (int i=0; i<NBAND[k]; i++) {
                for (int j=0; j<NBAND[k]; j++) {
                    if(i<N_peratom_HartrOrbit*NumCorrAtom) DF_CorrBase[k](i,j) = temp(i,j);  //Do not alter the projected-d orbitals,dim(DF_CorrBase) = NBAND*NBAND
                }
            }

            if(mpi_rank==0 and k==0) std::cout << "We have model projected-ligand orbital\n";
            H_k_inModelSpace[k].setZero(NBAND[k],NBAND[k]);
            for (int i=0; i<NBAND[k]; i++) H_k_inModelSpace[k](i,i) = KS_eigenEnergy[k][FromValToKS[k][i]];
            H_k_inModelSpace[k]          = DF_CorrBase[k] * H_k_inModelSpace[k] * DF_CorrBase[k].adjoint();

            if(mpi_rank==0 and k==0)   std::cout << "ConstructModel: downfolding done\n";
        }//k
        std::vector<int>  HartreeOrbital_idx(NumCorrAtom * 2);
        for(int i=0; i<NumCorrAtom; i++) {
            HartreeOrbital_idx[i*2+0] = N_peratom_HartrOrbit*i;
            HartreeOrbital_idx[i*2+1] = N_peratom_HartrOrbit*(i+1);
        }
        for (int i0=0; i0<NumOrbit; i0++) {
            if (i0 < N_peratom_HartrOrbit*NumCorrAtom) isOrbitalHartr[i0] = 1;
            else                                       isOrbitalHartr[i0] = 0;
        }
        // HartrRange was modified here.
        setCorrelatedSpaceIndex(HartreeOrbital_idx, NumCorrAtom);

    }/*downfolding*/

//    for(int k = 0;  k < knum; k++) {
////        for (int i=0; i<NBAND[k]; i+=2) {
////            H_k_inModelSpace[k](i+0, i+0) += Zeeman_field_spin(0,0);
////            H_k_inModelSpace[k](i+0, i+1) += Zeeman_field_spin(0,1);
////            H_k_inModelSpace[k](i+1, i+0) += Zeeman_field_spin(1,0);
////            H_k_inModelSpace[k](i+1, i+1) += Zeeman_field_spin(1,1);
////        }
//        for (int at=0; at<NumCorrAtom; at++) {
//            for(int  i=at*N_peratom_HartrOrbit; i<at*N_peratom_HartrOrbit+N_peratom_HartrOrbit; i+=2) {
//                H_k_inModelSpace[k](HartrIndex[i+0], HartrIndex[i+0]) += Zeeman_field_spin_corratom[at](0,0);
//                H_k_inModelSpace[k](HartrIndex[i+0], HartrIndex[i+1]) += Zeeman_field_spin_corratom[at](0,1);
//                H_k_inModelSpace[k](HartrIndex[i+1], HartrIndex[i+0]) += Zeeman_field_spin_corratom[at](1,0);
//                H_k_inModelSpace[k](HartrIndex[i+1], HartrIndex[i+1]) += Zeeman_field_spin_corratom[at](1,1);
//
//            }
//        }
//    }//k




    //(Re-)construct local energy level
    ifroot        std::cout <<"Himp, on-site:";
    std::vector<Eigen::MatrixXcd>  HRtemp(NumCluster);
    std::vector<Eigen::MatrixXcd>      HR(NumCluster);
    for(int cl=0; cl<NumCluster; cl++) {
        HRtemp[cl].setZero(NumHartrOrbit_per_cluster, NumHartrOrbit_per_cluster);
        HR[cl].setZero(NumHartrOrbit_per_cluster, NumHartrOrbit_per_cluster);
        for(int h1=0; h1<NumHartrOrbit_per_cluster; h1++) {
//            std::cout <<  mpi_rank << " " << HartrIndex[at*N_peratom_HartrOrbit+h1] <<"\n";
            for(int h2=0; h2<NumHartrOrbit_per_cluster; h2++) {
                int h1F = HartrIndex[cl*NumHartrOrbit_per_cluster+h1];
                int h2F = HartrIndex[cl*NumHartrOrbit_per_cluster+h2];
                HRtemp[cl](h1,h2)=0;
                HR[cl](h1,h2)=0;
                for(int k=0 ; k < knum; k++) {
                    HRtemp[cl](h1,h2) +=  (H_k_inModelSpace[k](h1F,h2F)  ) ;
                }
                HRtemp[cl](h1,h2) /= knum_mpiGlobal;
            }
        }
        MPI_Allreduce(HRtemp[cl].data(), HR[cl].data(), HRtemp[cl].size(), MPI_DOUBLE_COMPLEX, MPI_SUM,  MPI_COMM_WORLD);
    }


    ifroot        std::cout <<" (decomp) \n";
    ifroot std::cout << "We have " << NumCluster <<" clusters with "<< NumHartrOrbit_per_cluster <<" orbitals for each cluster\n";
    impurity_site_Hamiltonian.setZero( NumCluster*NumHartrOrbit_per_cluster, NumCluster*NumHartrOrbit_per_cluster);
    for(int cl=0; cl<NumCluster; cl++) {
        for(int i=0; i<NumHartrOrbit_per_cluster; i++) {
            int i0=cl*NumHartrOrbit_per_cluster+i;
            for(int j=0; j<NumHartrOrbit_per_cluster; j++) {
                int j0=cl*NumHartrOrbit_per_cluster+j;
                impurity_site_Hamiltonian( i0,j0) = HR[cl](i,j);
            }
            ifroot std::cout << std::fixed << std::setprecision(4)<< real(impurity_site_Hamiltonian(i0,i0))  <<" ";
        }
        ifroot        std::cout <<"\n";
    }
//    ifroot        std::cout <<"Himp, Matrix:\n";
//    for(int cl=0; cl<NumCluster; cl++) {
//        ifroot std::cout << std::fixed << std::setprecision(4)<<
//                         (impurity_site_Hamiltonian.block(cl*NumHartrOrbit_per_cluster,cl*NumHartrOrbit_per_cluster,
//                                 NumHartrOrbit_per_cluster, NumHartrOrbit_per_cluster) )  <<"\n";
//    }
}







void low_energy_subspace_in_KS_basis(
    int knum,  int knum_mpiGlobal,
    std::vector<int> & NBAND,  std::vector<std::vector<int> >  & FromValToKS,double muDFT,
    std::vector<Eigen::MatrixXcd> & KS_eigenVectors_orthoBasis, Eigen::VectorXd  * KS_eigenEnergy
) {
    /*Find NBAND and FromValToKS with given KS_eigenEnergy and muDFT*/
    ///////////////////////////////////////////////////////////////////
//    double  tempd2=0;
//    Eigen::VectorXd  HartreWeightInWindows_local;
//    HartreWeightInWindows_local.setZero(N_peratom_HartrOrbit);
//    Eigen::VectorXd  HartreWeightInWindows_global;
//    HartreWeightInWindows_global.setZero(N_peratom_HartrOrbit);

    for(int k=0; k< knum ; k++) {

        if (downfolding==1) {
            NBAND[k] = 0;
            double Emin=upper_model_window;
            double Emax=lower_model_window;
            for(int i=0; i<NumOrbit; i++) {
                if(  lower_model_window < (KS_eigenEnergy[k][i]-muDFT) and (KS_eigenEnergy[k][i]-muDFT) <  upper_model_window  ) {
                    NBAND[k]++;
                    if (k==0 and mpi_rank==0 and Emin>KS_eigenEnergy[k][i]-muDFT) Emin=KS_eigenEnergy[k][i]-muDFT;
                    if (k==0 and mpi_rank==0 and Emax<KS_eigenEnergy[k][i]-muDFT) Emax=KS_eigenEnergy[k][i]-muDFT;
                }

            }

            /*index for correlated (d-p) space */
            FromValToKS.at(k).resize(NBAND[k]);
            int i1=0;
            for(int i=0; i<NumOrbit; i++) {
                if(  lower_model_window < (KS_eigenEnergy[k][i]-muDFT) and (KS_eigenEnergy[k][i]-muDFT) <  upper_model_window  ) {
                    FromValToKS.at(k).at(i1) = i;
                    i1++;
                }
            }
            assert ( NBAND[k] != 0) ;
            assert ( i1 == NBAND[k]);
            if(mpi_rank==0 and k==0) std::cout << "<downfolding> NBAND at k=0 is " <<  NBAND[0] <<"\n";
            if(mpi_rank==0 and k==0) std::cout << "<downfolding> E_max at k=0 is " <<  Emax <<"\n";
            if(mpi_rank==0 and k==0) std::cout << "<downfolding> E_min at k=0 is " <<  Emin <<"\n";
            if(NBAND[k] < N_peratom_HartrOrbit*NumCorrAtom) {
                std::cout << "NBAND at rank " <<mpi_rank <<"k-point: " << k <<" = " << NBAND[k] <<"\n";
                exit(1);
            }
        }
        else {
            /* No downfolding, */
            FromValToKS.at(k).resize(NumOrbit);
            if(mpi_rank==0 and k==0 ) std::cout << "no downfolding\n";
            NBAND[k] = NumOrbit;
            for(int i=0; i<NumOrbit; i++) {
                FromValToKS[k][i] = i;
            }
        }
    }//k
}
