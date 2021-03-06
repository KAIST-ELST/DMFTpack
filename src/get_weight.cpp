#include "mpi.h"
#include "ImgTimeFtn.h"
#include "tight_common.h"
#include "model.h"
#include "TB.h"
#include "pulay.h"
#include <Eigen/Eigenvalues>
#include <algorithm>
void  sort_eigen_vectors( Eigen::MatrixXcd & evec, Eigen::VectorXd & eval, int dim ) ;


void get_preNAOs (Eigen::MatrixXcd weightMatrix,
                  Eigen::MatrixXcd & principal_number_tfm, Eigen::MatrixXcd & weightMatrix_preNAOs) {

//The principal quantum number is redefined to diagonalize the local density matrix.

///*subshell averaging*/
///*
    weightMatrix_preNAOs.setZero(NumOrbit, NumOrbit);
    for(int nl=0; nl<num_subshell; nl++) {
        for(int nll=0; nll<num_subshell; nll++) {
            int size_nl  =  subshell(nl+1) - subshell(nl);
            int size_nll =  subshell(nll+1) - subshell(nll);
            if(rot_sym(nl) == rot_sym(nll)) {
                if(size_nl !=  size_nll) {
                    std::cout << "WARNING: Check rot_sym " << nl <<"," << nll << "\n";
                    exit(1);
                }
                cmplx av_occ =         (weightMatrix.block( subshell(nl), subshell(nll), size_nl, size_nll ) ).trace() / size_nl;
                for(int i=0;  i<size_nl;  i++) {
                    weightMatrix_preNAOs(subshell(nl)+i, subshell(nll)+i)  = av_occ;
                }
            }
        }
    }
//            */
//N, diagonalize Al subspace. We now redefine principal quantum number!.
// /*
    Eigen::MatrixXcd densityMat_Al_block;
    densityMat_Al_block.setZero(num_subshell, num_subshell);
    principal_number_tfm.setZero(weightMatrix_preNAOs.rows(), weightMatrix_preNAOs.cols());

    for(int nl=0; nl<num_subshell; nl++) {
        for(int nll=0; nll<num_subshell; nll++) {
            int size_nl  =  subshell(nl+1) - subshell(nl);
            int size_nll =  subshell(nll+1) - subshell(nll);
            if(rot_sym(nl) == rot_sym(nll)) {
                densityMat_Al_block(nl, nll)  =  (weightMatrix_preNAOs.block( subshell(nl), subshell(nll), size_nl, size_nll ) ).trace();
            }
        }
    }
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> ces ( num_subshell );
    ces.compute(densityMat_Al_block);
    Eigen::VectorXd eval;
    Eigen::MatrixXcd evec;
    eval = ces.eigenvalues();
    evec = ces.eigenvectors();
    sort_eigen_vectors(  evec,  eval,num_subshell) ;


    for(int nl=0; nl<num_subshell; nl++) {
        for(int nll=0; nll<num_subshell; nll++) {

            int size_nl  =  subshell(nl+1) - subshell(nl);
            if(rot_sym(nl) == rot_sym(nll)) {
                for(int i = 0; i<size_nl; i++) {
                    principal_number_tfm( subshell(nl)+i, subshell(nll)+i   ) = evec(nl,nll);
                }//for,i
            }//if
        }
    }
    weightMatrix_preNAOs =    (principal_number_tfm.adjoint() * weightMatrix_preNAOs * principal_number_tfm).eval();
// */
}




//////////////////////////////////////////////////////////////////////////////
void get_NAO_transfrom( Eigen::MatrixXcd & Sk, Eigen::MatrixXcd &  KS_evec_k,
                        Eigen::MatrixXcd weightMatrix_preNAOs, Eigen::MatrixXcd & transformMatrix,
                        int kpoint,std::vector<int> & accumulated_Num_SpinOrbital,
                        Eigen::MatrixXcd principal_number_tfm,  bool WSW   ) {

//preNAO
    Eigen::MatrixXcd KS_evec_preNAOs = Sk * KS_evec_k;
    KS_evec_preNAOs = principal_number_tfm.adjoint()* KS_evec_preNAOs ;

    Eigen::MatrixXcd Sk_preNAOs;
    Sk_preNAOs = (principal_number_tfm.adjoint()* Sk *principal_number_tfm).eval();

    Eigen::MatrixXcd Os, Ow,   principal_number_tfm_rediag;

    principal_number_tfm_rediag.setIdentity(weightMatrix_preNAOs.rows(), weightMatrix_preNAOs.cols());
    Os.setIdentity(weightMatrix_preNAOs.rows(), weightMatrix_preNAOs.cols());
    Ow.setIdentity(weightMatrix_preNAOs.rows(), weightMatrix_preNAOs.cols());

    if(WSW==true) {
//S, Gram-Schmidt transformation for Rydberg set, to remove large overlap between minimal set and rydberg set.
//    /*
        int dim_minimal=0;
        int dim_rydberg=0;
        for(int nl=0; nl<num_subshell; nl++) {
            for(int m1= subshell(nl); m1< subshell(nl+1); m1++) {
                if(Rydberg_set[nl] == 1) {
                    dim_rydberg++;
                }
                else dim_minimal++;
            }
        }
        Eigen::MatrixXcd T_sub;
        Eigen::MatrixXcd S_sub[2][2];
        S_sub[0][0].setZero(dim_minimal, dim_minimal);
        S_sub[0][1].setZero(dim_minimal, dim_rydberg);
        S_sub[1][0].setZero(dim_rydberg, dim_minimal);
        S_sub[1][1].setZero(dim_rydberg, dim_rydberg);

        T_sub.setZero(dim_minimal, dim_rydberg);

        int ii[2], jj[2];
        ii[0]=0;
        ii[1]=0;
        for(int nl=0; nl<num_subshell; nl++) {
            for(int i = subshell(nl);  i<subshell(nl +1); i++) {
                jj[0]=0;
                jj[1]=0;
                for(int nll=0; nll<num_subshell; nll++) {
                    for(int j = subshell(nll); j<subshell(nll+1); j++) {

                        S_sub[Rydberg_set[nl]][Rydberg_set[nll]](ii[Rydberg_set[nl]], jj[Rydberg_set[nll]])
                            = Sk_preNAOs(i,j);

                        jj[Rydberg_set[nll]] ++;
                    }//jj
                }//nll
                ii[Rydberg_set[nl]] ++;
            }//i
        }
        T_sub =  -S_sub[0][0].inverse() * S_sub[0][1];


        int iii=0;
        for(int nl=0; nl<num_subshell; nl++) {
            if(Rydberg_set[nl]==0 ) {
                for(int i = subshell(nl);  i<subshell(nl +1); i++) {

                    int jjj=0;
                    for(int nll=0; nll<num_subshell; nll++) {
                        if( Rydberg_set[nll]==1) {
                            for(int j = subshell(nll); j<subshell(nll+1); j++) {
                                Os( i, j )  = T_sub(iii, jjj);
                                jjj++;
                            }//for, j
                        }
                    }//for nll

                    iii++;
                }//for, i
            }//if_nl
        }


        Sk_preNAOs =              (Os.adjoint() * Sk_preNAOs           * Os).eval();  //adjoint? inverse?
        KS_evec_preNAOs =         (Os.adjoint() *KS_evec_preNAOs   ).eval();
        weightMatrix_preNAOs =    (Os.adjoint() * weightMatrix_preNAOs * Os).eval();
        double check_Os=0;
        for(int nl=0; nl<num_subshell; nl++) {
            for(int nll=nl; nll<num_subshell; nll++) {
                if( Rydberg_set[nl]==1 and Rydberg_set[nll]==0) {
                    for(int i = subshell(nl); i<subshell(nl+1); i++) {
                        for(int j = subshell(nll); j<subshell(nll+1); j++) {
                            if(check_Os < std::abs(Sk_preNAOs( i, j )))  check_Os = std::abs(Sk_preNAOs(i,j)) ;
                        }//for, j
                    }
                }
            }//for nl
        }
        if(check_Os > 1e-5) {
            std::cout << "Overlap between minimal set and the Rydberg set\n";
        }
    }


////Ow
//occupancy-weighted symmetric orthogonalization for minimal set and SO for Rydberg set,
//since we don't interested in Rydberg set
    Eigen::MatrixXd weightMatrix_diag =   (weightMatrix_preNAOs.real().diagonal().asDiagonal());
//    for(int nl=0; nl<num_subshell; nl++) {
//        if( Rydberg_set[nl]==1) {
//            for(int j = subshell(nl); j<subshell(nl+1); j++) {
//                weightMatrix_diag( j , j )  =  1.0;
//            }//for, j
//        }
//    }//for nl

// /*
//input:  weightMatrix_diag, KS_evec_preNAOs    , NBAND_k,   Sk_preNAOs,    FromValToKS_k, accumulated_Num_SpinOrbital
//output:  weightMatrix_tilde

    Eigen::MatrixXcd weightMatrix_tilde = weightMatrix_diag;
//    bool constrainedNAOs   =  true;
    bool constrainedNAOs   =  false;
    if(constrainedNAOs == true) {
        Eigen::MatrixXd Pwbar, Pc;
        Eigen::MatrixXcd c_Pwbar_cinv;
        Pwbar.setIdentity(NumOrbit, NumOrbit);
        for (int i=0; i<NBAND[kpoint]; i++) {
            Pwbar(FromValToKS[kpoint][i], FromValToKS[kpoint][i]) =  0.0;
        }
        c_Pwbar_cinv  = Sk_preNAOs.inverse() * KS_evec_preNAOs * Pwbar * KS_evec_preNAOs.adjoint();

        Pc.setZero(NumOrbit, NumOrbit);
        for(int at=0; at<NumAtom; at++) {
            for(int h1=accumulated_Num_SpinOrbital[at]; h1<accumulated_Num_SpinOrbital[at+1]; h1++) {
                if(isOrbitalHartrDFT[h1])  Pc(h1,h1) = 1.0;
            }
        }

        double mixing =0.1;
        int iter = 0;
        double resid, resid_prev=0;

        Eigen::MatrixXcd H = (weightMatrix_tilde.adjoint()*Sk_preNAOs*weightMatrix_tilde );
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> ces( H );
        H= ces.operatorInverseSqrt();


        Eigen::MatrixXcd H_diag_inv;
        H_diag_inv.setIdentity(NumOrbit, NumOrbit);
        for(int at=0; at<NumAtom; at++) {
            for(int h1=accumulated_Num_SpinOrbital[at]; h1<accumulated_Num_SpinOrbital[at+1]; h1++) {
                for(int h2=accumulated_Num_SpinOrbital[at]; h2<accumulated_Num_SpinOrbital[at+1]; h2++) {
                    if(isOrbitalHartrDFT[h1] and isOrbitalHartrDFT[h2])  H_diag_inv(h1,h2) = H(h1,h2);
                }
            }
        }
        H_diag_inv = H_diag_inv.inverse();

        Eigen::MatrixXcd H_pHp_inv = H*(Pc*H_diag_inv*Pc);
        Eigen::MatrixXcd  H_pHp_inv_out ;
        pulayMixing mixing_weight_matrix(3, 20, 1, NumOrbit, NumOrbit, true );
        do {

//
//            for (int m=0; m<N_peratom_HartrOrbit; m++) {
//                double norm=0;
//                double norm1=0;
//                double norm2=0;
//                double norm3=0;
//                int m0 = HartrIndex_inDFT[m];
//                for (int i0=0; i0<NBAND[kpoint]; i0++) {
//                    Eigen::MatrixXcd temp = H*weightMatrix_tilde.adjoint()*KS_evec_preNAOs;
//                    norm+=std::pow(std::abs(  temp(m0,FromValToKS[kpoint][i0])),2);
//                    norm1+=std::pow(std::abs(  KS_evec_preNAOs(m0,FromValToKS[kpoint][i0])),2);
//                    norm3+=std::pow(std::abs(  (Sk*KS_evec_k)(m0,FromValToKS[kpoint][i0])),2);
//                }
//                for (int i0=0; i0<NumOrbit; i0++) {
//                    norm2+=std::pow(std::abs(  KS_evec_preNAOs(m0,i0) ),2);
//
//                }
//                std::cout <<m0<<": "<< norm << " " << norm1<<" " <<norm3<<" "  << norm2  <<"\n";
//            }
//

            Eigen::MatrixXcd id;
            id.setIdentity( NumOrbit,NumOrbit);
            weightMatrix_tilde = weightMatrix_diag - c_Pwbar_cinv *weightMatrix_diag*H_pHp_inv  ;
//
//            H = (Sk_preNAOs );
//            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> ces11( NumOrbit );
//            ces11.compute(H);
//            std::cout << std::setprecision(15) << "min_eval1: "<<ces11.eigenvalues()[0] <<"\n";
//
//            H = (weightMatrix_tilde.adjoint()*weightMatrix_tilde );
//            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> ces2( NumOrbit );
//            ces2.compute(H);
//            std::cout << std::setprecision(15) << "min_eval2: "<<ces2.eigenvalues()[0] <<"\n";



            H = (weightMatrix_tilde.adjoint()*Sk_preNAOs*weightMatrix_tilde );

            for (int n=0; n<NumOrbit; n++) {
                for (int m=0; m<NumOrbit; m++) {
                    if( ( std::isnan(std::abs(H(n,m)) ) or   std::isinf(std::abs(H(n,m)))      ))   {
                        std::cout <<"nan1\n";
                    }
                }
            }
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> ces( H );
            H= ces.operatorInverseSqrt();

            for (int n=0; n<NumOrbit; n++) {
                for (int m=0; m<NumOrbit; m++) {
                    if( ( std::isnan(std::abs(H(n,m)) ) or   std::isinf(std::abs(H(n,m)))      ))   {




                        H = (weightMatrix_tilde.adjoint()*Sk_preNAOs*weightMatrix_tilde );
                        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> ces33( NumOrbit );
                        ces33.compute(H);
                        std::cout << std::setprecision(15) << "min_eval3: "<<ces33.eigenvalues()[0] <<"\n";


                        std::cout <<"nan2\n";
                    }
                }
            }







            H_diag_inv.setIdentity(NumOrbit, NumOrbit);
            for(int at=0; at<NumAtom; at++) {
                for(int h1=accumulated_Num_SpinOrbital[at]; h1<accumulated_Num_SpinOrbital[at+1]; h1++) {
                    for(int h2=accumulated_Num_SpinOrbital[at]; h2<accumulated_Num_SpinOrbital[at+1]; h2++) {
                        if(isOrbitalHartrDFT[h1] and isOrbitalHartrDFT[h2])  H_diag_inv(h1,h2) = H(h1,h2);
                    }
                }
            }
            H_diag_inv = H_diag_inv.inverse();



            H_pHp_inv_out =  H*(Pc*H_diag_inv*Pc);

            resid = (H_pHp_inv - H_pHp_inv_out).norm() ;
            if(iter ==0 ) resid_prev = resid;
            else if((iter)%10 ==0) {
                mixing_checker(resid, resid_prev, mixing, 0.1, 1e-4);
                resid_prev = resid;
//                std::cout << "rank"<< mpi_rank <<  iter<< ": "<<resid << " " << mixing  <<"\n";

            }

            mixing_weight_matrix.mixing( &H_pHp_inv, &H_pHp_inv_out, mixing, iter, 1);
            H_pHp_inv = H_pHp_inv_out;
            iter++;
        } while(resid>1e-4);
//        std::cout <<"done, rank:" << mpi_rank <<"\n";
        weightMatrix_tilde = weightMatrix_diag - c_Pwbar_cinv *weightMatrix_diag*  H_pHp_inv  ;
    }
// */


    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> ces2( (weightMatrix_tilde.adjoint() *  Sk_preNAOs* weightMatrix_tilde) );
    Ow = weightMatrix_tilde * ces2.operatorInverseSqrt();

    Sk_preNAOs =              ( Ow.adjoint() * Sk_preNAOs           * Ow  ).eval();
    weightMatrix_preNAOs =    ( Ow.adjoint() * weightMatrix_preNAOs * Ow  ).eval();


    /*
    if(WSW==true) {
    //N_Restoration
           Eigen::MatrixXcd densityMat_Al_block;
           densityMat_Al_block.setZero(num_subshell, num_subshell);
           principal_number_tfm_rediag.setZero(weightMatrix.rows(), weightMatrix.cols());


           for(int nl=0; nl<num_subshell; nl++) {
               for(int nll=0; nll<num_subshell; nll++) {
                   int size_nl  =  subshell(nl+1)  - subshell(nl) ;
                   int size_nll =  subshell(nll+1) - subshell(nll);
                   if(rot_sym(nl) == rot_sym(nll)) {
                       densityMat_Al_block(nl, nll)  =
                           (weightMatrix_preNAOs.block( subshell(nl), subshell(nll)  , size_nl  , size_nll ) ).trace();
                   }
               }
           }
           Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> ces ( num_subshell );
           ces.compute(densityMat_Al_block);
           Eigen::VectorXd eval       =ces.eigenvalues() ;
           Eigen::MatrixXcd evec      =ces.eigenvectors();


           sort_eigen_vectors(  evec,  eval,num_subshell) ;

           for(int nl=0; nl<num_subshell; nl++) {
               for(int nll=0; nll<num_subshell; nll++) {

                   int size_nl  =  subshell(nl+1) - subshell(nl);
                   if(rot_sym(nl) == rot_sym(nll)) {
                       for(int i = 0; i<size_nl; i++) {
                           principal_number_tfm_rediag( subshell(nl)+i, subshell(nll)+i   ) = evec(nl,nll);
                       }//for,i
                   }//if
               }
           }

           Sk_preNAOs =              (principal_number_tfm_rediag.adjoint() * Sk_preNAOs           * principal_number_tfm_rediag).eval();
           weightMatrix_preNAOs =    (principal_number_tfm_rediag.adjoint() * weightMatrix_preNAOs * principal_number_tfm_rediag).eval();
    }
    // */
    transformMatrix =     principal_number_tfm * Os * Ow * principal_number_tfm_rediag ;  //<old|new>
}



Eigen::MatrixXcd getweight_dual( std::vector<Eigen::MatrixXcd> S_overlap, std::vector<Eigen::MatrixXcd> dual_DM_direct  ) {
    Eigen::MatrixXcd  weightMatrix;
    weightMatrix.setZero(NumOrbit, NumOrbit);

    Eigen::MatrixXcd  weightMatrix_locl;
    weightMatrix_locl.setZero(NumOrbit, NumOrbit);
    for(int k = 0;  k < knum; k++) {
        Eigen::MatrixXcd dual_DM_dual =  dual_DM_direct[k] * S_overlap[k].inverse();
        weightMatrix_locl += dual_DM_dual ;
    }
    MPI_Allreduce(weightMatrix_locl.data(), weightMatrix.data(), weightMatrix.size(), MPI_DOUBLE_COMPLEX, MPI_SUM,  MPI_COMM_WORLD);
    weightMatrix /= knum_mpiGlobal;
    return weightMatrix;
}//weight_dual

Eigen::MatrixXcd getweight_hyb( std::vector<Eigen::MatrixXcd> S_overlap, std::vector<Eigen::MatrixXcd> dual_DM_direct  ) {
    Eigen::MatrixXcd  weightMatrix;
    weightMatrix.setZero(NumOrbit, NumOrbit);

    Eigen::MatrixXcd  weightMatrix_locl;
    weightMatrix_locl.setZero(NumOrbit, NumOrbit);
    for(int k = 0;  k < knum; k++) {
        Eigen::MatrixXcd  dual_DM_direct_hyb;
        dual_DM_direct_hyb  = 0.5* dual_DM_direct[k];
        dual_DM_direct_hyb  = (dual_DM_direct_hyb + dual_DM_direct_hyb.adjoint()).eval();
        weightMatrix_locl += dual_DM_direct_hyb;
    }
    MPI_Allreduce(weightMatrix_locl.data(), weightMatrix.data(), weightMatrix.size(), MPI_DOUBLE_COMPLEX, MPI_SUM,  MPI_COMM_WORLD);
    weightMatrix /= knum_mpiGlobal;
    return weightMatrix;
}//weight_dual


Eigen::MatrixXcd getweight_direct(std::vector<Eigen::MatrixXcd> S_overlap,  std::vector<Eigen::MatrixXcd> dual_DM_direct  ) {
    Eigen::MatrixXcd  weightMatrix;
    Eigen::MatrixXcd  weightMatrix_locl;

    weightMatrix.setZero(NumOrbit, NumOrbit);
    weightMatrix_locl.setZero(NumOrbit, NumOrbit);

    for(int k = 0;  k < knum; k++) {
        Eigen::MatrixXcd temp = (S_overlap[k] * dual_DM_direct[k]);
        weightMatrix_locl += temp;
    }
    MPI_Allreduce(weightMatrix_locl.data(), weightMatrix.data(), weightMatrix.size(), MPI_DOUBLE_COMPLEX, MPI_SUM,  MPI_COMM_WORLD);
    weightMatrix /= knum_mpiGlobal;
    return weightMatrix;
}//weight_direct

//Eigen::MatrixXcd getweight_sc(  std::vector<Eigen::MatrixXcd>  S_overlap, std::vector<Eigen::MatrixXcd> dual_DM_direct , Eigen::MatrixXcd  weightMatrix, int *accumulated_Num_SpinOrbital ) {
//
//
//    double     mixing = 0.9;
//    double resid, resid_prev=0;
//    Eigen::MatrixXcd id;
//    Eigen::MatrixXcd direct_DM_direct[knum];
//    for(int k = 0;  k < knum; k++) {
//        direct_DM_direct[k] = (S_overlap[k] * dual_DM_direct[k]);
//    }
//    id.setIdentity(NumOrbit,NumOrbit);
//
//    Eigen::MatrixXcd    weightMatrix_out;
//    weightMatrix_out = weightMatrix;
//    int iter = 0;
//    pulayMixing Mixing_NAOsWeight(3,100, 1, NumOrbit, NumOrbit);
//    do {
//
//
//        Eigen::MatrixXcd  weightMatrix_locl;
//        weightMatrix_locl.setZero(NumOrbit, NumOrbit);
//        for(int k = 0;  k < knum; k++) {
//            Eigen::MatrixXcd transformMatrix;  // \ket{new_i} = \sum_j \ket{old_j} T_{ji}
//            get_NAO_transfrom( S_overlap[k], weightMatrix , transformMatrix, false);
//
//            Eigen::MatrixXcd temp;
//            temp = (transformMatrix.adjoint() * direct_DM_direct[k] * transformMatrix).eval();
//            weightMatrix_locl += temp;
//
//        }//k
//        weightMatrix_out.setZero(NumOrbit, NumOrbit);
//        MPI_Allreduce(weightMatrix_locl.data() , weightMatrix_out.data(), weightMatrix_out.size(), MPI_DOUBLE_COMPLEX, MPI_SUM,  MPI_COMM_WORLD);
//        weightMatrix_out /= knum_mpiGlobal;
//
//
//        resid = (weightMatrix-weightMatrix_out).norm() ;
//        if(iter ==0 ) resid_prev = resid;
//        else if(iter%10 ==0) {
//            mixing_checker(resid, resid_prev, mixing, 0.9, 1e-3);
//            resid_prev = resid;
//        }
//
//        ifroot std::cout <<"weight RD:" << resid  <<"  mixing:" << mixing << " for iter:" << iter <<  "\n";
//        Mixing_NAOsWeight.mixing( &weightMatrix, &weightMatrix_out, mixing,  iter, 2);
//        weightMatrix  = weightMatrix_out;
//        iter++;
//    } while(resid > 1e-4);
//
//
//    return weightMatrix_out;
//}//weight_sc









void lowdin_symmetric_orthogonalization( Eigen::MatrixXcd & Hk, Eigen::MatrixXcd & Sk, Eigen::MatrixXcd  & evec,  Eigen::VectorXd & eval, Eigen::MatrixXcd & transformMatrix   ) {
//    Eigen::MatrixXcd    transformMatrix;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> ces1(Sk);
    transformMatrix =  (  ces1.operatorInverseSqrt() );
    Hk = transformMatrix * Hk *transformMatrix;     //\ket{j_ortho} = \ket{DFT_AO} T_{ij}

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> ces(Hk);
    ces.compute(Hk);
    eval = ces.eigenvalues();
    evec = ces.eigenvectors();
}


void naturalAtomicOrbitals_population_weighted_symmetric_orthogonalization_r(
    Eigen::MatrixXcd & Hk, Eigen::MatrixXcd & Sk, int kpoint,  std::vector<int> & accumulated_Num_SpinOrbital,
    Eigen::MatrixXcd & evec,  Eigen::VectorXd & eval,
    Eigen::MatrixXcd weightMatrix_preNAOs,
    Eigen::MatrixXcd & transformMatrix,
    Eigen::MatrixXcd  principal_number_tfm) {

    /*subshell averaging*/

//    Eigen::MatrixXcd  Sk_preNAOs = preNAOs.adjoint() * (Sk) *preNAOs;    // <i|j> = <i|a> S^-1 <b|c> S^-1 <d|j> = <i|a> S^-1 <d|j>, plz check
    Eigen::MatrixXcd  KS_evec_k  =  evec;
    get_NAO_transfrom(Sk, KS_evec_k, weightMatrix_preNAOs, transformMatrix, kpoint, accumulated_Num_SpinOrbital,
                      principal_number_tfm, true);   //  \ket{j_ortho} = \ket{DFT_AO} T_{ij}

    Hk  = (transformMatrix.adjoint() * Hk * transformMatrix).eval();

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> ces3( Hk );
    ces3.compute(Hk);
    evec = ces3.eigenvectors();
    eval = ces3.eigenvalues();
}


void direct_projection( Eigen::MatrixXcd & Hk, Eigen::MatrixXcd & Sk, Eigen::MatrixXcd  & evec,  Eigen::VectorXd & eval   ) {
//    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXcd> ces1(NumOrbit);
//    ces1.compute( Hk, Sk );       //  HS \psi = S  \psi E
//
//    eval = ces1.eigenvalues();
    evec = (Sk*evec).eval();
}






void  sort_eigen_vectors( Eigen::MatrixXcd & evec, Eigen::VectorXd & eval, int dim ) {
    Eigen::VectorXd eval_temp;
    Eigen::MatrixXcd evec_temp;
    eval_temp.setZero(dim);
    evec_temp.setZero(dim,dim);
    double RD;


    std::vector<int> occup(dim,0);
    for(int i=0; i<dim; i++) {
        eval_temp[i] =  eval[i];
        for(int j=0; j<dim; j++) {
            evec_temp(i,j) = evec(i,j);
        }
    }

//for given n, evec_sort[i,n]=evec[i,m],
//where m=max{ abs(evec[n,m']) | m'}
    for(int n=0; n<dim; n++) {
        int m_max_for_n=-1;
        double RD_max=0;
        for(int m=0; m<dim; m++) {
            RD = std::abs(evec_temp(n,m));
            if(RD>RD_max) {
                m_max_for_n = m;   //will be placed at jmax
                RD_max = RD;
            }
        }
        eval(n) = eval_temp(m_max_for_n);
        for(int k=0; k<dim; k++) {
            evec(k,n) = evec_temp(k,m_max_for_n);
        }
    }

    //int jmax=0;
    //for(int i=0; i<dim; i++) {
    //    ith-eigenvectors is assigned to jmax
    //    RD=0.0;
    //    jmax=0;
    //    for(int j=0; j<dim; j++) {
    //        RD = std::abs(evec_temp(j,i));
    //        if(RD>RD_prev and occup.at(j)==0) {
    //            jmax = j;   //will be placed at jmax
    //        }
    //        RD_prev = RD;
    //    }
    //    occup.at(jmax) =1;
    //    eval(jmax) = eval_temp(i);
    //    for(int k=0; k<dim; k++) {
    //        evec(k,jmax) = evec_temp(k,i);
    //    }
    //}



    //    do {
    //
    //        RD=0.0;
    //        for(int i=0; i<dim; i++) {
    //            RD +=  real(evec_temp(i, permu[i]));   //max Re <\phi_i | v_p[i] >
    //        }
    //
    //        if( RD> RD_prev) {
    //            for(int i=0; i<dim; i++) {
    //                permu_max[i] = permu[i];
    //            }
    //        }
    //        RD_prev = RD;
    //
    //    } while(next_permutation(permu.begin(),permu.end()));
    //    for(int i=0; i<dim; i++) {
    //        eval[i] =  eval_temp[permu_max[i]];
    //        for(int j=0; j<dim; j++) {
    //            evec(i,j) = evec_temp(permu_max[i],permu_max[j]);
    //        }
    //    }
}
