#include "pulay.h"
#include <Eigen/QR>
#include <Eigen/Eigenvalues>
void mixing_checker(double resid, double resid_prev, double & mixing,  double mixing_max,  double mixing_min) {
    if       ( resid_prev > resid ) mixing = std::min (1.07    *mixing ,  mixing_max);
    else if  ( resid_prev < resid ) mixing = std::max (1./1.09 *mixing ,  mixing_min);
}



pulayMixing::pulayMixing(int mixing_history_, int start_mixing_,int dim_i_, int dim_j_, int dim_k_, bool parallel_pulay_) {
    mixing_history = mixing_history_;
    start_mixing = start_mixing_;
    dim_i = dim_i_;
    dim_j = dim_j_;
    dim_k = dim_k_;
    parallel_pulay = parallel_pulay_;
    if(parallel_pulay or mpi_rank==0) {
        input_history = new  Eigen::MatrixXcd * [mixing_history];
        res_history   = new  Eigen::MatrixXcd * [mixing_history];
        for(int k=0; k<mixing_history; k++) {
            input_history[k] = new Eigen::MatrixXcd [dim_i];
            res_history[k]   = new  Eigen::MatrixXcd [dim_i];
            for(int i=0; i<dim_i; i++) {
                input_history[k][i].setZero(dim_j, dim_k);
                res_history[k][i].setZero(dim_j, dim_k);
            }
        }
        OptimalMatrix.setZero(mixing_history+1, mixing_history+1);
        for(int k=0; k<mixing_history; k++) {
            OptimalMatrix(mixing_history,k)= -1;
            OptimalMatrix(k,mixing_history)= -1;
            OptimalMatrix(mixing_history-1, mixing_history-1)=-100.;
        }
        OptimalMatrix(mixing_history, mixing_history)=0.;
    }
}
pulayMixing::~pulayMixing() {
    if(parallel_pulay or mpi_rank==0) {
        for(int k=0; k<mixing_history; k++) {
            delete [] input_history[k];
            delete [] res_history[k];
        }
        delete [] input_history;
        delete [] res_history;
    }
}
void pulayMixing::mixing(Eigen::MatrixXcd  * inputDensity_n, Eigen::MatrixXcd * outputDensity_n, double  mixingPulay, int SCGFloop, int mixingStep) {
    Eigen::MatrixXcd temp(dim_i, dim_j*dim_k);
    temp.setZero(dim_i, dim_j*dim_k);
    if(parallel_pulay or mpi_rank==0) {
//        double maxValue=0;
        int removeComp=mixing_history-1;
//        //find maximum res components
//        if(SCGFloop > mixing_history  and OptimalMatrix(mixing_history-1, mixing_history-1) > 0) {
//            for(int k=0; k < mixing_history; k++) {
//                if (maxValue < OptimalMatrix(k,k) ) {
//                    maxValue = OptimalMatrix(k,k);
//                    removeComp= k;
//                }
//            }
//        }
        //shift history
        for(int k=removeComp; k>0; k--) {
            for(int i=0; i<dim_i; i++) {
                input_history[k][i]=input_history[k-1][i];
                res_history[k][i]=  res_history[k-1][i];
            }
        }
        for(int i=0; i<dim_i; i++) {
            input_history[0][i]=inputDensity_n[i];
            res_history[0][i]=outputDensity_n[i] - inputDensity_n[i];
        }

        if(SCGFloop>start_mixing) {
            //TODO : efficiency
            //shift optimal matrix
            Eigen::MatrixXd temp(mixing_history-1, mixing_history-1);
            for(int k=0; k < mixing_history-1; k++) {
                int indxk=k;
                if (k>=removeComp)  indxk++;
                for(int l=0; l < mixing_history-1; l++) {
                    int indxl=l;
                    if (l>=removeComp)  indxl++;
                    temp(k,l)= OptimalMatrix(indxk,indxl);
                }
            }
            for(int k=1; k < mixing_history; k++) {
                for(int l=1; l < mixing_history; l++) {
                    OptimalMatrix(k,l) =       temp(k-1,l-1);
                }
            }

            //new elements for optimal matrix
            for(int k=0; k<mixing_history; k++) {
                OptimalMatrix(0,k) = 0;
                for(int i=0; i<dim_i; i++) {
                    OptimalMatrix(0,k) += std::real(  (res_history[0][i].adjoint() * res_history[k][i]).trace()   );
                }
                OptimalMatrix(k,0) =  OptimalMatrix(0,k);
            }
        }

        std::vector<Eigen::MatrixXcd> optimalDensity(dim_i);
        std::vector<Eigen::MatrixXcd> optimalres(dim_i);
        for(int i=0; i<dim_i; i++) {
            optimalDensity[i].setZero(dim_j, dim_k);
            optimalres[i].setZero(dim_j, dim_k);
        }

        //Linear solver to find optimal history mixing values
//        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ces( OptimalMatrix );
        Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(OptimalMatrix);
        bool lin_dep_check = lu_decomp.isInvertible();
        if(SCGFloop>mixing_history + start_mixing and SCGFloop%mixingStep==0 and lin_dep_check   ) {
            Eigen::VectorXd mixing_beta(mixing_history+1);
            Eigen::VectorXd zeroVec(mixing_history+1);
            zeroVec.setZero(mixing_history+1);
            zeroVec(mixing_history) = -1;
            mixing_beta = OptimalMatrix.colPivHouseholderQr().solve(zeroVec);

            for(int k=0; k<mixing_history; k++) {
                for(int i=0; i<dim_i; i++) {
                    optimalDensity[i] += mixing_beta[k] * input_history[k][i];
                    optimalres[i]     += mixing_beta[k] * res_history[k][i];
                }
            }
        }
        else {
            if( SCGFloop>mixing_history + start_mixing and SCGFloop%mixingStep==0 and   !(lin_dep_check)   ) std::cout << "Warning: PulayMixing,linear dependence arise\n";
            for(int i=0; i<dim_i; i++) {
                optimalDensity[i] = inputDensity_n[i];
                optimalres[i] =  res_history[0][i];
            }
        }
        //double minValue=9999;
        //int minComp=0;
        //for(int k=0; k < mixing_history; k++) {
        //    if (minValue > OptimalMatrix(k,k) ) {
        //        minValue = OptimalMatrix(k,k);
        //        minComp = k;
        //    }
        //}
        for(int i=0; i<dim_i; i++) {
            outputDensity_n[i] =  optimalDensity[i] + mixingPulay * optimalres[i];
//            outputDensity_n[i] =  (1-mixingPulay) * input_history[minComp][i] + mixingPulay  * (  optimalDensity[i] + mixingPulay * optimalres[i] );
//            outputDensity_n[i] =  (1-mixingPulay) * input_history[0][i] + mixingPulay  * (  optimalDensity[i] + mixingPulay * optimalres[i] );
        }
        for(int i=0; i<dim_i; i++) {
            for(int j=0; j<dim_j; j++) {
                for(int k=0; k<dim_k; k++) {
                    temp(i,j*dim_k + k) =  outputDensity_n[i](j,k);
                }
            }
        }
    }//ifroot
    if( parallel_pulay == false ) MPI_Bcast(temp.data(),  temp.size(), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

    for(int i=0; i<dim_i; i++) {
        for(int j=0; j<dim_j; j++) {
            for(int k=0; k<dim_k; k++) {
                outputDensity_n[i](j,k) =        temp(i,j*dim_k + k) ;
            }
        }
    }
}
//void pulayMixing::mixing( std::vector<Eigen::MatrixXcd> & inputDensity_n, std::vector<Eigen::MatrixXcd>  & outputDensity_n, double  mixingPulay, int SCGFloop, int mixingStep) {
void pulayMixing::mixing( ImgFreqFtn & inputDensity_n, ImgFreqFtn  & outputDensity_n, double  mixingPulay, int SCGFloop, int mixingStep) {
    Eigen::MatrixXcd *  input = new Eigen::MatrixXcd [dim_i];
    Eigen::MatrixXcd * output = new Eigen::MatrixXcd [dim_i];
    for(int i=0; i<dim_i; i++) {
        input[i] =   inputDensity_n.getMatrix(i);
        output[i] = outputDensity_n.getMatrix(i);
    }
    mixing( input, output, mixingPulay, SCGFloop, mixingStep);
    for(int i=0; i<dim_i; i++) {
        inputDensity_n.setMatrix(i, input[i] );
        outputDensity_n.setMatrix(i, output[i]);
    }
    delete [] input;
    delete [] output;


}


void pulayMixing::mixing(Eigen::VectorXd  * inputDensity_n, Eigen::VectorXd * outputDensity_n, double  mixingPulay, int SCGFloop, int mixingStep) {
    assert(dim_i==1);
    assert(dim_j==1);
    Eigen::MatrixXcd input(dim_j,dim_k);
    Eigen::MatrixXcd output(dim_j,dim_k);

    input.row(0) =  inputDensity_n[0];
    output.row(0) = outputDensity_n[0];
    mixing(&input, &output, mixingPulay, SCGFloop, mixingStep);



    inputDensity_n[0] = (input.row(0)  ).real();
    outputDensity_n[0]= (output.row(0) ).real();
}




//void pulayMixing::my_mixing(Eigen::MatrixXcd  * inputDensity_n, Eigen::MatrixXcd * outputDensity_n, double  mixingPulay, int SCGFloop, int mixingStep) {
//    Eigen::MatrixXcd temp(dim_i, dim_j*dim_k);
//    temp.setZero(dim_i, dim_j*dim_k);
////        double maxValue=0;
//    int removeComp=mixing_history-1;
////        //find maximum res components
////        if(SCGFloop > mixing_history  and OptimalMatrix(mixing_history-1, mixing_history-1) > 0) {
////            for(int k=0; k < mixing_history; k++) {
////                if (maxValue < OptimalMatrix(k,k) ) {
////                    maxValue = OptimalMatrix(k,k);
////                    removeComp= k;
////                }
////            }
////        }
//    //shift history
//    for(int k=removeComp; k>0; k--) {
//        for(int i=0; i<dim_i; i++) {
//            input_history[k][i]=input_history[k-1][i];
//            res_history[k][i]=  res_history[k-1][i];
//        }
//    }
//    for(int i=0; i<dim_i; i++) {
//        input_history[0][i]=inputDensity_n[i];
//        res_history[0][i]=outputDensity_n[i] - inputDensity_n[i];
//    }
//
//    if(SCGFloop>start_mixing) {
//        //TODO : efficiency
//        //shift optimal matrix
//        Eigen::MatrixXd temp(mixing_history-1, mixing_history-1);
//        for(int k=0; k < mixing_history-1; k++) {
//            int indxk=k;
//            if (k>=removeComp)  indxk++;
//            for(int l=0; l < mixing_history-1; l++) {
//                int indxl=l;
//                if (l>=removeComp)  indxl++;
//                temp(k,l)= OptimalMatrix(indxk,indxl);
//            }
//        }
//        for(int k=1; k < mixing_history; k++) {
//            for(int l=1; l < mixing_history; l++) {
//                OptimalMatrix(k,l) =       temp(k-1,l-1);
//            }
//        }
//
//        //new elements for optimal matrix
//        for(int k=0; k<mixing_history; k++) {
//            OptimalMatrix(0,k) = 0;
//            for(int i=0; i<dim_i; i++) {
//                OptimalMatrix(0,k) += std::real(  (res_history[0][i].adjoint() * res_history[k][i]).trace()   );
//            }
//            OptimalMatrix(k,0) =  OptimalMatrix(0,k);
//        }
//    }
//
//    std::vector<Eigen::MatrixXcd> optimalDensity(dim_i);
//    std::vector<Eigen::MatrixXcd> optimalres(dim_i);
//    for(int i=0; i<dim_i; i++) {
//        optimalDensity[i].setZero(dim_j, dim_k);
//        optimalres[i].setZero(dim_j, dim_k);
//    }
//
//    //Linear solver to find optimal history mixing values
//    if(SCGFloop>mixing_history + start_mixing and SCGFloop%mixingStep==0  ) {
//        Eigen::VectorXd mixing_beta(mixing_history+1);
//        Eigen::VectorXd zeroVec(mixing_history+1);
//        zeroVec.setZero(mixing_history+1);
//        zeroVec(mixing_history) = -1;
//        mixing_beta = OptimalMatrix.colPivHouseholderQr().solve(zeroVec);
//
//        for(int k=0; k<mixing_history; k++) {
//            for(int i=0; i<dim_i; i++) {
//                optimalDensity[i] += mixing_beta[k] * input_history[k][i];
//                optimalres[i]     += mixing_beta[k] * res_history[k][i];
//            }
//        }
//    }
//    else {
//        for(int i=0; i<dim_i; i++) {
//            optimalDensity[i] = inputDensity_n[i];
//            optimalres[i] =  res_history[0][i];
//        }
//    }
//    //double minValue=9999;
//    //int minComp=0;
//    //for(int k=0; k < mixing_history; k++) {
//    //    if (minValue > OptimalMatrix(k,k) ) {
//    //        minValue = OptimalMatrix(k,k);
//    //        minComp = k;
//    //    }
//    //}
//    for(int i=0; i<dim_i; i++) {
//        outputDensity_n[i] =  optimalDensity[i] + mixingPulay * optimalres[i];
////            outputDensity_n[i] =  (1-mixingPulay) * input_history[minComp][i] + mixingPulay  * (  optimalDensity[i] + mixingPulay * optimalres[i] );
////            outputDensity_n[i] =  (1-mixingPulay) * input_history[0][i] + mixingPulay  * (  optimalDensity[i] + mixingPulay * optimalres[i] );
//    }
//    for(int i=0; i<dim_i; i++) {
//        for(int j=0; j<dim_j; j++) {
//            for(int k=0; k<dim_k; k++) {
//                temp(i,j*dim_k + k) =  outputDensity_n[i](j,k);
//            }
//        }
//    }
////    MPI_Bcast(temp.data(),  temp.size(), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
//    for(int i=0; i<dim_i; i++) {
//        for(int j=0; j<dim_j; j++) {
//            for(int k=0; k<dim_k; k++) {
//                outputDensity_n[i](j,k) =        temp(i,j*dim_k + k) ;
//            }
//        }
//    }
//}
