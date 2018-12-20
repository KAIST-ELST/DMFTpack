#include "tight_common.h"
#include "model.h"
#include "TB.h"
#include <Eigen/Dense>

void Construct_H0_local(Eigen::MatrixXcd * Heff_loc, double mu, ImgFreqFtn & SelfE_w, int atom);

void Construct_hyb_delta(int impurityDim, std::vector<int> impurityOrbit,
                         ImgFreqFtn & SelfE_w, std::vector<Eigen::MatrixXcd>   Gw, double mu,
                         ImgFreqFtn &  weiss_fieldTB, int atom,     Eigen::MatrixXcd & SolverBasis  ) {
    int   h2, h1F, h2F, h1H,h2H;
    cmplx iw;
    Eigen::MatrixXcd   projimpurity_site_Hamiltonian, projSolverBasis ;
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



    if ( SOLVERtype.find(std::string("SEG")) != std::string::npos ) {//diagonal part dmft

//        Eigen::MatrixXcd  impurity_site_Hamiltonian_bestBasis = projSolverBasis.adjoint() * projimpurity_site_Hamiltonian * projSolverBasis;
        for(int n=0; n<N_freq; n++) {
//            projGw[n] =     (projSolverBasis).adjoint() * projGw[n]  *   projSolverBasis;
//            projSw[n] =     (projSolverBasis).adjoint() * projSw[n]  *   projSolverBasis;
            for(int h1=0; h1<impurityDim; h1++) {
//            int h1F = impurityOrbit.at(h1);
                iw=I*pi*(2.*n+1.)/beta;
                weiss_fieldTB.setValueSubMat( n, atom,  h1, h1,
                                              iw+mu
                                              - impurity_site_Hamiltonian(h1,h1)-   projSw[n](h1,h1)
                                              - 1.0/projGw[n](h1,h1));
                //NOTE : (1/G_{ii}) \neq G^{-1}_{ii}
            }//n
        }//h1
    }//approxLevel==0
    else {//single-site DMFT with off-diagonal part hybridization
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
    }//approxLevel==1



//    for(int i0=0; i0<impurityDim; i0++) {
//          for(int l0=0; l0<impurityDim; l0++) {
//            int i0A = atom*impurityDim+i0;
//            int l0A = atom*impurityDim+l0;
//            double r3, r2, r1,w3,w2,w1;
//            double rr3, rr2, rr1,ww3,ww2,ww1;
//
//            r3  = -(imag(weiss_fieldTB.getValue(N_freq-3,i0A,l0A)))  / std::pow(weiss_fieldTB.getValue(N_freq-3),1);
//            r2  = -(imag(weiss_fieldTB.getValue(N_freq-2,i0A,l0A)))  / std::pow(weiss_fieldTB.getValue(N_freq-2),1);
//            r1  = -(imag(weiss_fieldTB.getValue(N_freq-1,i0A,l0A)))  / std::pow(weiss_fieldTB.getValue(N_freq-1),1);
//            w3  = 1./std::pow(weiss_fieldTB.getValue(N_freq-3),2);
//            w2  = 1./std::pow(weiss_fieldTB.getValue(N_freq-2),2);
//            w1  = 1./std::pow(weiss_fieldTB.getValue(N_freq-1),2);
//
//            rr3 = -(real(weiss_fieldTB.getValue(N_freq-3,i0A,l0A)))  / std::pow(weiss_fieldTB.getValue(N_freq-3),2);
//            rr2 = -(real(weiss_fieldTB.getValue(N_freq-2,i0A,l0A)))  / std::pow(weiss_fieldTB.getValue(N_freq-2),2);
//            rr1 = -(real(weiss_fieldTB.getValue(N_freq-1,i0A,l0A)))  / std::pow(weiss_fieldTB.getValue(N_freq-1),2);
//            ww3 = 1/std::pow(weiss_fieldTB.getValue(N_freq-3),4);
//            ww2 = 1/std::pow(weiss_fieldTB.getValue(N_freq-2),4);
//            ww1 = 1/std::pow(weiss_fieldTB.getValue(N_freq-1),4);
//
//            weiss_fieldTB.setValue(N_freq+0, i0A, l0A,   0       ) ;
//            weiss_fieldTB.setValue(N_freq+1, i0A, l0A,   (r1+r2+r3)/(w1+w2+w3)       ) ;
//            weiss_fieldTB.setValue(N_freq+2, i0A, l0A,   (rr1+rr2+rr3)/(ww1+ww2+ww3) ) ;
//            weiss_fieldTB.setValue(N_freq+3, i0A, l0A, 0) ;
//            weiss_fieldTB.setValue(N_freq+4, i0A, l0A, 0) ;
//        }//l0
//    }//i0

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
    FourierTransform (delta_w, delta_t, HYB,1 ) ;
    delete [] delta_w;
}
