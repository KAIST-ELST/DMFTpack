#include "tight_common.h"
#include "model.h"
#include "TB.h"
#include <Eigen/Dense>

//void Construct_H0_local(Eigen::MatrixXcd * Heff_loc, double mu, ImgFreqFtn & SelfE_w, int atom);

void Construct_Gloc(int impurityDim, std::vector<int> impurityOrbit,
                    ImgFreqFtn & SelfE_w, std::vector<Eigen::MatrixXcd>   Gw, double mu,
                    ImgFreqFtn &  weiss_fieldTB, int atom,     Eigen::MatrixXcd & SolverBasis  ) {



    int   h2, h1F, h2F, h1H,h2H;
    cmplx iw;

    std::vector<Eigen::MatrixXcd>     projGw(N_freq);
    for(int n=0; n<N_freq; n++) {
        projGw[n].setZero(impurityDim, impurityDim);
    }

    for(int h1=0; h1< impurityDim; h1++) {
        for(int h2=0; h2< impurityDim; h2++) {
            int h1F = impurityOrbit.at(h1);
            int h2F = impurityOrbit.at(h2);
            for(int n=0; n<N_freq; n++) {
                projGw[n](h1,h2) =  Gw[n](h1F,h2F);
            }
        }
    }



    for(int n=0; n<N_freq; n++) {
        iw=I*pi*(2.*n+1.)/beta;
        for(int h1=0; h1< impurityDim ; h1++) {
            for(int h2=0; h2< impurityDim ; h2++) {
                weiss_fieldTB.setValueSubMat(n, atom,  h1, h2,
                                             projGw[n](h1,h2) )  ;
            }//h2
        }//h1
    }//n

}


void Construct_hyb_delta(int impurityDim, std::vector<int> impurityOrbit,
                         ImgFreqFtn & SelfE_w, std::vector<Eigen::MatrixXcd>   Gw, double mu,
                         ImgFreqFtn &  weiss_fieldTB, int atom,     Eigen::MatrixXcd & SolverBasis  ) {
    int   h2, h1F, h2F, h1H,h2H;
    cmplx iw;
    Eigen::MatrixXcd   projimpurity_site_Hamiltonian;
//Eigen::MatrixXcd projSolverBasis ;
//    int start_hyb_index = atom*impurityDim;

    std::vector<Eigen::MatrixXcd>     projSw(N_freq);
    std::vector<Eigen::MatrixXcd>     projGw(N_freq);
    for(int n=0; n<N_freq; n++) {
        projSw[n].setZero(impurityDim,impurityDim);
        projGw[n].setZero(impurityDim, impurityDim);
    }
    projimpurity_site_Hamiltonian.setZero(impurityDim, impurityDim);
//    projSolverBasis.setZero(impurityDim, impurityDim);

    for(int h1=0; h1< impurityDim; h1++) {
        for(int h2=0; h2< impurityDim; h2++) {
            int h1F = impurityOrbit.at(h1);
            int h2F = impurityOrbit.at(h2);
            for(int n=0; n<N_freq; n++) {
                projSw[n](h1,h2) =  SelfE_w.getValue(n,h1F,h2F);
                projGw[n](h1,h2) =  Gw[n](h1F,h2F);
            }
            projimpurity_site_Hamiltonian(h1,h2) = impurity_site_Hamiltonian(h1F,h2F);
//            projSolverBasis(h1,h2) = SolverBasis(h1F,h2F);
        }
    }



//      if ( SOLVERtype.find(std::string("SEG")) != std::string::npos ) {//diagonal part dmft
//
//  //        Eigen::MatrixXcd  impurity_site_Hamiltonian_bestBasis = projSolverBasis.adjoint() * projimpurity_site_Hamiltonian * projSolverBasis;
//          for(int n=0; n<N_freq; n++) {
//  //            projGw[n] =     (projSolverBasis).adjoint() * projGw[n]  *   projSolverBasis;
//  //            projSw[n] =     (projSolverBasis).adjoint() * projSw[n]  *   projSolverBasis;
//              for(int h1=0; h1<impurityDim; h1++) {
//  //            int h1F = impurityOrbit.at(h1);
//                  iw=I*pi*(2.*n+1.)/beta;
//                  weiss_fieldTB.setValueSubMat( n, atom,  h1, h1,
//                                                iw+mu
//                                                - projimpurity_site_Hamiltonian(h1,h1)-   projSw[n](h1,h1)
//                                                - 1.0/projGw[n](h1,h1));
//                  //NOTE : (1/G_{ii}) \neq G^{-1}_{ii}
//              }//n
//          }//h1
//      }//approxLevel==0
//      else {//single-site DMFT with off-diagonal part hybridization
    for(int n=0; n<N_freq; n++) {
        iw=I*pi*(2.*n+1.)/beta;
        if(impurityDim!=0)     projGw[n] = (projGw[n].inverse());
        for(int h1=0; h1< impurityDim ; h1++) {
//                int h1at = start_hyb_index  +h1;
//            int h1F = impurityOrbit.at(h1);
            for(int h2=0; h2< impurityDim ; h2++) {
//                    int h2at = start_hyb_index +h2;
//            int h2F = impurityOrbit.at(h1);
                weiss_fieldTB.setValueSubMat(n, atom,  h1, h2,
                                             -projimpurity_site_Hamiltonian(h1,h2) -    projSw[n](h1,h2) -projGw[n](h1,h2) )  ;
            }//h2
            weiss_fieldTB.setValueSubMat(n,atom, h1,h1,
                                         weiss_fieldTB.getValueSubMat(n,atom,h1,h1)+iw +mu )  ;
        }//h1
    }//n
//    }//approxLevel==1




}//Construct_hyb_delta





void write_hyb_t ( ImgFreqFtn & weiss_field, Eigen::MatrixXcd * delta_t, int atom) {
    Eigen::MatrixXcd * delta_w  = new Eigen::MatrixXcd [N_freq];
    for(int n=0; n<N_freq; n++) {
        delta_w[n].setZero(NSpinOrbit_per_atom,NSpinOrbit_per_atom);
        for(int h1=0; h1<NSpinOrbit_per_atom; h1++) {
            for(int h2=0; h2<NSpinOrbit_per_atom; h2++) {
                delta_w[n](h1,h2) = weiss_field.getValue(n,h1,h2);
            }
        }//h1
    }
    std::string HYB = std::string("delta_w in write_hyb_t()");
    FourierTransform (delta_w, delta_t, false, HYB,1 ) ;
    delete [] delta_w;
}
