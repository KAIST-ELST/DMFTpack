#include "tight_common.h"
#include "model.h"
#include "TB.h"
#include <Eigen/Dense>

void Construct_H0_local(Eigen::MatrixXcd * Heff_loc, double mu, ImgFreqFtn & SelfE_z, int atom);

void Construct_hyb_delta(int impurityDim, std::vector<int> impurityOrbit,
                         ImgFreqFtn & SelfE_z, std::vector<Eigen::MatrixXcd>   Gw, double mu,
                         ImgFreqFtn &  weiss_fieldTB, int atom, Eigen::MatrixXcd & SolverBasis,     int segmentsolver  ) {
    int   h2, h1F, h2F, h1H,h2H;
    cmplx iw;
    Eigen::MatrixXcd   projimpurity_site_Hamiltonian;

    std::vector<Eigen::MatrixXcd>     projSw(N_freq);
    std::vector<Eigen::MatrixXcd>     projGw(N_freq);
    for(int n=0; n<N_freq; n++) {
        projSw[n].setZero(impurityDim,impurityDim);
        projGw[n].setZero(impurityDim, impurityDim);
    }
    projimpurity_site_Hamiltonian.setZero(impurityDim, impurityDim);



    Eigen::MatrixXcd temp;
    temp=  SolverBasis.adjoint() * impurity_site_Hamiltonian * SolverBasis;

    for(int h1=0; h1< impurityDim; h1++) {
        for(int h2=0; h2< impurityDim; h2++) {
            int h1F = impurityOrbit.at(h1);
            int h2F = impurityOrbit.at(h2);
            projimpurity_site_Hamiltonian(h1,h2) = temp(h1F,h2F);
        }
    }

    for(int n=0; n<N_freq; n++) {
        Eigen::MatrixXcd temp1, temp2;
        temp1 =     ( SolverBasis.adjoint() * Gw[n]  *   SolverBasis ).eval();
        temp2 =     ( SolverBasis.adjoint() * SelfE_z.getMatrix(n) *   SolverBasis ).eval();
        for(int h1=0; h1< impurityDim; h1++) {
            for(int h2=0; h2< impurityDim; h2++) {
                int h1F = impurityOrbit.at(h1);
                int h2F = impurityOrbit.at(h2);
                projGw[n](h1,h2) =  temp1(h1F,h2F);
                projSw[n](h1,h2) =  temp2(h1F,h2F);
            }
        }
    }






    if (segmentsolver==1 ) {
        //NOTE : (1/G_{ii}) \neq G^{-1}_{ii}
        Eigen::MatrixXcd temp;
        temp =  projimpurity_site_Hamiltonian.diagonal().asDiagonal();
        projimpurity_site_Hamiltonian =  temp;
        for(int n=0; n<N_freq; n++) {
            Eigen::MatrixXcd temp1, temp2;
            temp1 =     (projGw[n].diagonal().asDiagonal());
            temp2 =     (projSw[n].diagonal().asDiagonal());
            projGw[n] =     temp1;
            projSw[n] =     temp2;
        }//n
    }//approxLevel==0
    for(int n=0; n<N_freq; n++) {
        iw=I*pi*(2.*n+1.)/beta;
        if(impurityDim!=0)     projGw[n] = (projGw[n].inverse());
        for(int h1=0; h1< impurityDim ; h1++) {
            for(int h2=0; h2< impurityDim ; h2++) {
                weiss_fieldTB.setValueSubMat(n, atom,  h1, h2,
                                             -projimpurity_site_Hamiltonian(h1,h2) -    projSw[n](h1,h2) -projGw[n](h1,h2) )  ;
            }//h2
            weiss_fieldTB.setValueSubMat(n,atom, h1,h1,
                                         weiss_fieldTB.getValueSubMat(n,atom,h1,h1)+iw +mu )  ;
        }//h1
//ifroot std::cout << "jhs3 \n"<<  weiss_fieldTB.getMatrix(n)<<"\n";
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
    FourierTransform (delta_w, delta_t, false, HYB,0 ) ;
    delete [] delta_w;
}
