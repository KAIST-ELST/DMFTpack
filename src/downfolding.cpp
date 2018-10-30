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




    if (downfolding==1) {
        for(int k = 0;  k < knum; k++) {
            /*project out rest space ;*/
            /*d-orbital space extraction*/
            Eigen::MatrixXcd  S_overlap_dSpace;
            DF_CorrBase[k].setZero( N_peratom_HartrOrbit*NumCorrAtom , NBAND[k]);
            int p00=0;
            for(int at1=0; at1<NumCorrAtom; at1++) {
                for(int p0=HartrRange_DFT[at1][0] ; p0<HartrRange_DFT[at1][1] ; p0++) {
                    for(int i1=0; i1<NBAND[k]; i1++) {
                        /*DF_CorrBase = <\phi_\alpha in model | nk > ( size = dim(CorrSpace) * NBAND ) */
                        DF_CorrBase[k](p00,i1) =  KS_eigenVectors_orthoBasis[k](p0,FromValToKS.at(k).at(i1));
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
            //Following  is the definition of the DF_CorrBase;::  <\alpha_ortho|nk in W>
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





    //Construct  effective Hamiltonian
//    for(int k = 0;  k < knum; k++) {
//        for(int at1=0; at1<NumCorrAtom; at1++) {
//            for(int p0=HartrRange[at1][0] ; p0<HartrRange[at1][1] ; p0++) {
//                for(int q0=HartrRange[at1][0] ; q0<HartrRange[at1][1] ; q0++) {
//                    if(doublecounting == 1)         H_k_inModelSpace[k](p0,q0)  -=  Sw_doublecounting(LongRangeOrder[p0],LongRangeOrder[q0]);
//                    H_k_inModelSpace[k](p0,q0)  +=  Sw_Hartree(LongRangeOrder[p0],LongRangeOrder[q0]);
//                }
//            }
//        }
//        if (mpi_rank ==0 and k==0  )          std::cout << "Effective Hamiltonian in model space = H_DFT - doublecounting + Sw_Hartree\n";
//    }


    //(Re-)construct local energy level
    ifroot        std::cout <<"Himp, on-site:";
    std::vector<Eigen::MatrixXcd>  HRtemp(NumCorrAtom);
    std::vector<Eigen::MatrixXcd>      HR(NumCorrAtom);
    for(int at=0; at<NumCorrAtom; at++) {
        HRtemp[at].setZero(N_peratom_HartrOrbit, N_peratom_HartrOrbit);
        HR[at].setZero(N_peratom_HartrOrbit, N_peratom_HartrOrbit);
        for(int h1=0; h1<N_peratom_HartrOrbit; h1++) {
//            std::cout <<  mpi_rank << " " << HartrIndex[at*N_peratom_HartrOrbit+h1] <<"\n";
            for(int h2=0; h2<N_peratom_HartrOrbit; h2++) {
                int h1F = HartrIndex[at*N_peratom_HartrOrbit+h1];
                int h2F = HartrIndex[at*N_peratom_HartrOrbit+h2];
                HRtemp[at](h1,h2)=0;
                HR[at](h1,h2)=0;
                for(int k=0 ; k < knum; k++) {
                    HRtemp[at](h1,h2) +=  (H_k_inModelSpace[k](h1F,h2F)  ) ;
                }
                HRtemp[at](h1,h2) /= knum_mpiGlobal;
            }
        }
        MPI_Allreduce(HRtemp[at].data(), HR[at].data(), HRtemp[at].size(), MPI_DOUBLE_COMPLEX, MPI_SUM,  MPI_COMM_WORLD);
    }


    ifroot        std::cout <<" (decomp) \n";
    for(int at=0; at<NumCorrAtom; at++) {
        impurity_site_Hamiltonian[at].setZero(N_peratom_HartrOrbit,N_peratom_HartrOrbit);
        for(int i=0; i<N_peratom_HartrOrbit; i++) {
            for(int j=0; j<N_peratom_HartrOrbit; j++) {
                impurity_site_Hamiltonian[at](i,j)= HR[at](i,j);
            }
            ifroot std::cout << std::fixed << std::setprecision(4)<< real(impurity_site_Hamiltonian[at](i,i))  <<" ";
        }
        ifroot        std::cout <<"\n";
    }
    ifroot        std::cout <<"Himp, Matrix:\n";
    for(int at=0; at<NumCorrAtom; at++) {
        ifroot std::cout << std::fixed << std::setprecision(4)<< (impurity_site_Hamiltonian[at])  <<"\n";
    }
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
            /*Check d-orbital character in  energy window*/
//            for (int m=0; m<N_peratom_HartrOrbit; m++) {
//                int m0 = HartrIndex_inDFT[m];
//                for (int i0=0; i0<NumOrbit; i0++) {
//                    tempd2 += std::pow(std::abs(KS_eigenVectors_orthoBasis[k](m0,i0)),2);
//                    if( (KS_eigenEnergy[k][i0]-muDFT) < lower_model_window or KS_eigenEnergy[k][i0]-muDFT > upper_model_window  ) {
//                        HartreWeightInWindows_local[m]+=std::pow(std::abs(KS_eigenVectors_orthoBasis[k](m0,i0)),2);
//                    }
//                }
//            }
            NBAND[k] = 0;
            for(int i=0; i<NumOrbit; i++) {
                if(  lower_model_window < (KS_eigenEnergy[k][i]-muDFT) and (KS_eigenEnergy[k][i]-muDFT) <  upper_model_window  ) NBAND[k]++;
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
            if(mpi_rank==0 and k==0) std::cout << "<down-folding> NBAND at k=0 is " <<  NBAND[0] <<"\n";
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
//    double tempd2GL=0;
//    MPI_Allreduce(&tempd2, &tempd2GL, 1, MPI_DOUBLE, MPI_SUM,  MPI_COMM_WORLD);
//    ifroot std::cout << "d-orbital normalization       : "  << tempd2GL/(knum_mpiGlobal*N_peratom_HartrOrbit)<<"\n";
//    MPI_Allreduce(HartreWeightInWindows_local.data() , HartreWeightInWindows_global.data(), HartreWeightInWindows_global.size(), MPI_DOUBLE, MPI_SUM,  MPI_COMM_WORLD);
//    for (int m=0; m<N_peratom_HartrOrbit; m++) {
//        ifroot std::cout << "d-orbital in rest energy space: "  << HartreWeightInWindows_global[m] /(knum_mpiGlobal) <<"\n";
//    }

    ///////////////////////////////////////////////////////////////////
}
