#include <Eigen/Eigenvalues>
#include "tight_common.h"
#include "model.h"
#include "TB.h"



void dos(Eigen::VectorXd * KS_eigenEnergy,std::vector<Eigen::MatrixXcd> & KS_eigenVectors_orthoBasis, double E_window, double & muDFT  ) {
    int  index ;
    double de = real(dE);
    ifroot printf("********************************\n");
    ifroot printf("                DOS\n");
    ifroot printf("********************************\n");
    ifroot std::cout << "Chemical potential:" << muDFT<<"\n";

    /*initialize*/
    Eigen::VectorXd      TdosData;
    Eigen::VectorXd  TdosData_RDC;
    Eigen::MatrixXd  PdosData_RDC;
    Eigen::MatrixXd      PdosData;
    double E;
    TdosData.setZero(Spectral_EnergyGrid);
    TdosData_RDC.setZero(Spectral_EnergyGrid);
    PdosData_RDC.setZero(Spectral_EnergyGrid,NumCorrAtom*N_peratom_HartrOrbit);
    PdosData.setZero(Spectral_EnergyGrid,NumCorrAtom*N_peratom_HartrOrbit);

    /*calculate dos*/
    for(int k=0; k<knum; k++) {
        for(int band=0; band<NumOrbit; band++) {
            for(int index=0; index<Spectral_EnergyGrid; index++) {
                double energy_i      =  -E_window+index*de  ;
                double energy_i_high =     energy_i+de/2.   ;
                double energy_i_low  =     energy_i-de/2.   ;
                /*Lorentzian broadening :\integral 1/x**2 +1 = atan(x)*/
                if (    std::abs( energy_i- KS_eigenEnergy[k][band] +muDFT ) <=  50.* infinitesimal  ) {
                    double Aw = (   std::atan( (energy_i_high - KS_eigenEnergy[k][band]+muDFT)/infinitesimal )
                                    -std::atan( (energy_i_low - KS_eigenEnergy[k][band]+muDFT)/infinitesimal )  ) / pi;
                    TdosData[index] += Aw ;
                    assert(Aw>=0);
                    for(int at1=0; at1<NumCorrAtom; at1++) {
                        for(int m=0; m< N_peratom_HartrOrbit ; m++) {
//                            int mDFT = m+HartrRange_DFT[at1][0];
                            int mDFT = HartrIndex_inDFT[at1*N_peratom_HartrOrbit+m];
                            PdosData(index,at1*N_peratom_HartrOrbit+m) +=  Aw * pow(std::abs(KS_eigenVectors_orthoBasis[k](mDFT,band)),2)    ;
                        }
                    }
                }//if
                else if ( energy_i-KS_eigenEnergy[k][band]+muDFT > 50.*infinitesimal) break;
            }/*dos data index*/
        }/*band*/
    }//k
    MPI_Allreduce(PdosData.data(), PdosData_RDC.data(), PdosData.size(), MPI_DOUBLE, MPI_SUM,  MPI_COMM_WORLD);
    MPI_Allreduce(TdosData.data(), TdosData_RDC.data(), TdosData.size(), MPI_DOUBLE, MPI_SUM,  MPI_COMM_WORLD);



    /*Normalization*/

    TdosData_RDC/=(knum_mpiGlobal*NumCorrAtom*N_peratom_HartrOrbit);
    PdosData_RDC/=knum_mpiGlobal;
    double dosnormalization;
    for(int at1=0; at1<NumCorrAtom; at1++) {
        for(int m=0; m< N_peratom_HartrOrbit ; m++) {
            dosnormalization=0;
            for(int i=0; i<Spectral_EnergyGrid; i++) dosnormalization += PdosData_RDC(i,at1*N_peratom_HartrOrbit+m);
            ifroot std::cout << "normalization for atom" << at1 <<" and orbital"<< m <<":" << dosnormalization<<"\n";
        }
    }


    /*DOS distribution*/
    ifroot{
        double  aver[NumCorrAtom*N_peratom_HartrOrbit];
        double aver2[NumCorrAtom*N_peratom_HartrOrbit];
        for(int at1=0; at1<NumCorrAtom; at1++) {
            for(int m1=0; m1< N_peratom_HartrOrbit ; m1++) {
                aver[m1]=0;
                aver2[m1]=0;
                double   E=  -E_window-(de)-muDFT;
                for(int i=0; i<Spectral_EnergyGrid; i++) {
                    E+=de;
                    aver[m1]  +=PdosData_RDC(i,at1*N_peratom_HartrOrbit+m1) * E;
                    aver2[m1] +=PdosData_RDC(i,at1*N_peratom_HartrOrbit+m1) * std::pow(E,2);
                }
                if(at1==0) std::cout << "D2,orbital" << m1<<": " << (aver2[m1] - std::pow(aver[m1],2)) <<"\n";
            }
        }
    }


    /*number operator*/
    cmplx    occupation_mat[NumCorrAtom][NumCorrAtom][N_peratom_HartrOrbit][N_peratom_HartrOrbit];
    cmplx occupation_matLoc[NumCorrAtom][NumCorrAtom][N_peratom_HartrOrbit][N_peratom_HartrOrbit];
    for(int at1=0; at1<NumCorrAtom; at1++) {
        for(int at2=0; at2<NumCorrAtom; at2++) {
            for(int m1=0; m1<N_peratom_HartrOrbit; m1++) {
                for(int m2=0; m2<N_peratom_HartrOrbit; m2++) {
                    occupation_mat[at1][at2][m1][m2] =0.;
                    occupation_matLoc[at1][at2][m1][m2] =0.;
                }
            }
        }
    }
    /* numMat_ij = 1/Nk \sum_{kn}  <kn|i><j|kn> */
    for(int k=0; k<knum; k++) {
        for(int at1=0; at1<NumCorrAtom; at1++) {
            for(int at2=0; at2<NumCorrAtom; at2++) {
                for(int m1=0; m1<N_peratom_HartrOrbit; m1++) {
                    for(int m2=0; m2<N_peratom_HartrOrbit; m2++) {
//                        int m1DFT = m1+HartrRange_DFT[at1][0];
//                        int m2DFT = m2+HartrRange_DFT[at2][0];
                        int m1DFT = HartrIndex_inDFT[at1*N_peratom_HartrOrbit+m1];
                        int m2DFT = HartrIndex_inDFT[at2*N_peratom_HartrOrbit+m2];
                        for(int band=0; band<NumOrbit; band++) {
                            occupation_matLoc[at1][at2][m1][m2]
                            += conj(KS_eigenVectors_orthoBasis[k](m1DFT,band)) * KS_eigenVectors_orthoBasis[k](m2DFT,band) * 1./(1.+std::exp(beta*(KS_eigenEnergy[k][band]-muDFT)));   //   |<band|orbital>|^2
                        }
                    }//m2
                }//m1
            }//at2
        }//at1
    }
    MPI_Allreduce(occupation_matLoc, occupation_mat, NumCorrAtom*NumCorrAtom*N_peratom_HartrOrbit*N_peratom_HartrOrbit, MPI_DOUBLE_COMPLEX, MPI_SUM,  MPI_COMM_WORLD);
    for(int at1=0; at1<NumCorrAtom; at1++) {
        ifroot std::cout <<  "Corr ATOM=" << at1<<"\n";
        for(int m1=0; m1<N_peratom_HartrOrbit; m1++) {
            for(int m2=0; m2<N_peratom_HartrOrbit; m2++) {
                occupation_mat[at1][at1][m1][m2] /=knum_mpiGlobal;
                ifroot std::cout << occupation_mat[at1][at1][m1][m2] <<" ";
            }
            ifroot std::cout << "\n";
        }
        ifroot std::cout << "\n";
    }





    /*Write result*/
    ifroot {
        for(int at1=0; at1<NumCorrAtom; at1++) {
            Eigen::MatrixXcd occupation_mat_corr(N_peratom_HartrOrbit,N_peratom_HartrOrbit);
            for(int i=0; i < N_peratom_HartrOrbit; i++) {
                for(int j=0; j < N_peratom_HartrOrbit; j++) {
                    occupation_mat_corr(i,j) = occupation_mat[at1][at1][i][j] ;
                }
            }
//            std::cout <<occupation_mat_corr <<"\n";

            /* rotation 2D
                       Eigen::MatrixXcd my2Drotation_Mat(2,2), temp(2,2);
            //            my2Drotation_Mat(0,0)=  std::cos(11.72/180.0*3.141592) ;
            //            my2Drotation_Mat(0,1)=  std::sin(11.72/180.0*3.141592) ;
            //            my2Drotation_Mat(1,0)= -std::sin(11.72/180.0*3.141592) ;
            //            my2Drotation_Mat(1,1)=  std::cos(11.72/180.0*3.141592);

                       my2Drotation_Mat(0,0)=  std::cos(0/180.0*3.141592) ;
                       my2Drotation_Mat(0,1)=  std::sin(0/180.0*3.141592) ;
                       my2Drotation_Mat(1,0)= -std::sin(0/180.0*3.141592) ;
                       my2Drotation_Mat(1,1)=  std::cos(0/180.0*3.141592);
                       for(int spin=0; spin<2; spin++) {
                           for(int i=0; i<2; i++) {
                               for(int j=0; j<2; j++) {
                                   temp(i,j)=0;
                                   for(int m=0; m<2; m++) {
                                       for(int l=0; l<2; l++) {
                                           temp(i,j) +=  my2Drotation_Mat(m,i) * occupation_mat_corr(2*(1+m)+spin, 2*(1+l)+spin ) * my2Drotation_Mat(l,j);
                                       }
                                   }

                               }
                           }

                           for(int m=0; m<2; m++) {
                               for(int l=0; l<2; l++) {
                                   occupation_mat_corr(2*(1+m)+spin, 2*(1+l)+spin ) = temp(m,l);
                               }
                           }
                       }//spin
            // */

            double entropy=0, traceN=0, ev;
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> ces(N_peratom_HartrOrbit);
            ces.compute(occupation_mat_corr);
            for(int i=0; i < N_peratom_HartrOrbit; i++) {
                ev = ces.eigenvalues()[i];
                if(std::abs(ev-0.)>1e-10 and std::abs(ev-1.)>1e-10) {
                    entropy -= ev*std::log(ev) + (1-ev)*std::log((1-ev));
                }
                traceN += ev;
                ifroot std::cout <<std::setprecision(4)
                                 <<"For index" <<i <<  " Occ:" << real(occupation_mat_corr(i,i))<<":"
                                 <<  " ev:"  <<ev <<":"
                                 <<  " S:"   << (-ev*std::log(ev) - (1-ev)*std::log((1-ev)))/std::log(2.) << "\n";
            }//i

            ifroot std::cout << "eigen vector of local occ:\n" << ces.eigenvectors() <<"\n";
            /*
            // jeff projection
                      Eigen::VectorXcd  jeff_half_u(6);
                      Eigen::VectorXcd  jeff_half_d(6);
                      Eigen::VectorXcd  jeff_three_half_mj_p3(6);
                      Eigen::VectorXcd  jeff_three_half_mj_m3(6);
                      Eigen::VectorXcd  jeff_three_half_mj_p1(6);
                      Eigen::VectorXcd  jeff_three_half_mj_m1(6);

                      cmplx a = std::sqrt(1./3.);
                      cmplx b = std::sqrt(1./2.);
                      cmplx c = std::sqrt(2./3.);
                      cmplx d = std::sqrt(1./6.);

            ///////////////////////////////////////////////////////////////////
            ////                               xyu xyd; yzu d;   xzu   d
            //j1
                      jeff_half_u           << a ,0. ,0.,a    ,0.   , -I*a;
                      jeff_half_d           << 0.,a  ,-a,0.   ,-I*a , 0.;
            //j2
                      jeff_three_half_mj_p3 << 0.,0.  ,-b,0.   ,I*b,0.;
                      jeff_three_half_mj_m3 << 0.,0.  ,0.,b    ,0. ,I*b;
            //j3
                      jeff_three_half_mj_p1 << c ,0. ,0,-d   ,0. ,I*d;
                      jeff_three_half_mj_m1 << 0.,c  ,d,0    ,I*d,0.;
            ///////////////////////////////////////////////////////////////////



                      for(int n=0; n<N_peratom_HartrOrbit; n++) { //occupation eigvec index
                          cmplx j1p=0, j2p=0, j3p=0;
                          cmplx j1m=0, j2m=0, j3m=0;
                          std::cout << "occupation eigenvector, |" <<n<<">\n";
                          for(int i=0; i<N_peratom_HartrOrbit; i++) {
                              j1p += conj(jeff_half_u(i)) * ces.eigenvectors()(i,n);
                              j1m += conj(jeff_half_d(i)) * ces.eigenvectors()(i,n);

                              j2p += conj(jeff_three_half_mj_p3(i)) * ces.eigenvectors()(i,n);
                              j2m += conj(jeff_three_half_mj_m3(i)) * ces.eigenvectors()(i,n);

                              j3p += conj(jeff_three_half_mj_p1(i)) * ces.eigenvectors()(i,n);
                              j3m += conj(jeff_three_half_mj_m1(i)) * ces.eigenvectors()(i,n);
                          }
                          std::cout << "<j1|-abs:"<< std::abs(std::abs(j1p)) << ", " << std::abs(std::abs(j1m)) << "; " << (std::abs(j1p)+std::abs((j1m)))/2.  <<  "\n";
                          std::cout << "<j2|-abs;"<< std::abs(std::abs(j2p)) << ", " << std::abs(std::abs(j2m)) << "; " << (std::abs(j2p)+std::abs((j2m)))/2.  <<  "\n";
                          std::cout << "<j3|-abs;"<< std::abs(std::abs(j3p)) << ", " << std::abs(std::abs(j3m)) << "; " << (std::abs(j3p)+std::abs((j3m)))/2.  <<  "\n\n";

                          std::cout << "j1-sq:"<< std::pow(std::abs(j1p),2) + std::pow(std::abs(j1m),2) << ";" << std::sqrt(std::pow(std::abs(j1p),2) + std::pow(std::abs(j1m),2)) <<"\n";
                          std::cout << "j2-sq;"<< std::pow(std::abs(j2p),2) + std::pow(std::abs(j2m),2) << ";" << std::sqrt(std::pow(std::abs(j2p),2) + std::pow(std::abs(j2m),2)) <<"\n";
                          std::cout << "j3-sq;"<< std::pow(std::abs(j3p),2) + std::pow(std::abs(j3m),2) << ";" << std::sqrt(std::pow(std::abs(j3p),2) + std::pow(std::abs(j3m),2)) <<"\n\n";
                      }

            // */


            entropy /= log(2.);
            ifroot std::cout <<  "entropy:=" << entropy<<"\n";
            ifroot std::cout <<  "TotalN :=" << traceN<<"\n";

            FILE *datap2;   //dos.dat
//            datap2 = fopen( (std::string("dos.dat")+std::to_string(  (at1+1) )).c_str()       , "w");
            datap2 = fopen( (std::string("dos.dat")+ intToString(  (at1+1) )).c_str(), "w");
            fprintf(datap2, "E    Totaldos");
            for(int k=0; k<N_peratom_HartrOrbit; k++) {
                fprintf(datap2, "  Pdos%d",k);
            }
            fprintf(datap2, "\n");

            E=-E_window-(de);
            for(int n=0; n<Spectral_EnergyGrid; n++) {
                E+=de;
                fprintf(datap2, "%0.5f    %0.8f", E,TdosData_RDC[n]/de);
                for(int k=0; k<N_peratom_HartrOrbit; k++) {
                    fprintf(datap2, "  %0.8f",PdosData_RDC(n,at1*N_peratom_HartrOrbit+k)/de);
                }
                fprintf(datap2, "\n");
            }
            fclose(datap2);


            double IDOS[N_peratom_HartrOrbit+1];
            datap2 = fopen( (std::string("Integdos.dat")+ intToString(at1+1)).c_str(), "w");
            fprintf(datap2, "E    Totaldos");
            IDOS[0] = 0.0;
            for(int k=0; k<N_peratom_HartrOrbit; k++) {
                IDOS[k+1] =0.0;
                fprintf(datap2, "  Pdos%d",k);
            }
            fprintf(datap2, "\n");

            E=-E_window-(de);
            for(int n=0; n<Spectral_EnergyGrid; n++) {
                E+=de;
                IDOS[0] += TdosData_RDC[n];
                fprintf(datap2, "%0.5f    %0.8f", E,IDOS[0]);
                for(int k=0; k<N_peratom_HartrOrbit; k++) {
                    IDOS[k+1] += PdosData_RDC(n,at1*N_peratom_HartrOrbit+k);
                    fprintf(datap2, "  %0.8f",IDOS[k+1]);
                }
                fprintf(datap2, "\n");
            }
            fclose(datap2);


        }//atom

    }//ifroot
    double TNumEle=0;
    double TNumEle_local=0;
    for(int k=0; k<knum; k++) {
        for(int band=0; band<NumOrbit; band++) {
            if (KS_eigenEnergy[k][band]-muDFT < -E_window  ) TNumEle_local +=1.0;
        }
    }
    MPI_Allreduce(&(TNumEle_local), &(TNumEle), 1, MPI_DOUBLE, MPI_SUM,  MPI_COMM_WORLD);
    TNumEle/=knum_mpiGlobal;

    E=-E_window-(de);
    for(int n=0; n<Spectral_EnergyGrid; n++) {
        E+=de;
        if( E<0)    TNumEle += (TdosData_RDC[n]*NumCorrAtom*N_peratom_HartrOrbit);
    }
    ifroot std::cout << "Total Num of electron from dos = " << TNumEle <<"\n";
}




































double FromHkToNele (double muDFT, std::vector<Eigen::VectorXd>  & KS_eigenEnergy) {

    double TNumEle=0;
    double TNumEle_local=0;
    for(int k=0; k<knum; k++) {
        for(int band=0; band<NumOrbit; band++) {
            TNumEle_local +=  1./(1+std::exp( beta*(KS_eigenEnergy[k][band]-muDFT) ));
        }
    }
    MPI_Allreduce(&(TNumEle_local), &(TNumEle), 1, MPI_DOUBLE, MPI_SUM,  MPI_COMM_WORLD);
    TNumEle/=knum_mpiGlobal;
    return TNumEle;
}

double Nele_non_Inter(
    int knum, int knum_mpiGlobal,
    std::vector<Eigen::MatrixXcd> & H_k_inModelSpace,  std::vector<Eigen::MatrixXcd> & S_overlap
)
{
//find chemical potentiol for non-interacting band

    std::vector<Eigen::VectorXd> KS_eigenEnergy(knum);

    for(int k=0; k< knum; k++) {
        Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXcd> ces1(NumOrbit);
        ces1.compute( H_k_inModelSpace[k], S_overlap[k] );       //  HS \psi = S  \psi E
        KS_eigenEnergy[k] = ces1.eigenvalues();
    }

    double Nele[2];
    double dmu = 1.0;
    double mu = 1.0;


    for (int i=0; i<2; i++)  Nele[i]=-1;
    Nele[0] =  FromHkToNele(mu, KS_eigenEnergy);

    /*Find chemical potential*/
    double muUB=0, muLB=0;
    double elecNumGoal = NumberOfElectron;
    int nearAns=0;
    if(Nele[0] > NumberOfElectron) mu-=(dmu-dmu);
    else if(Nele[0] < NumberOfElectron) mu-=(dmu+dmu);
    int step=0;
    while (  (Nele[0] != elecNumGoal)   and  std::abs(dmu)>1e-6 ) {
        mu += dmu;
        for(int i=2-1; i>0 ; i--) Nele[i]=Nele[i-1];
        Nele[0] = FromHkToNele(mu, KS_eigenEnergy);

        if (nearAns ==0) {
            if ( Nele[0] < elecNumGoal) {
                muLB = mu;
                if (Nele[1]>elecNumGoal && Nele[1] >0 ) nearAns = 1;
                else dmu = fabs(dmu);
            }
            if ( Nele[0] > elecNumGoal) {
                muUB = mu;
                if (Nele[1] < elecNumGoal && Nele[1] > 0)  nearAns =1;
                else dmu = -1*fabs(dmu);
            }
            if ( Nele[0] == elecNumGoal) break;
            assert( std::abs(mu) <  600)   ;
        }
        else if (nearAns ==1) {
            if (Nele[0] > elecNumGoal) muUB = mu;
            else if(Nele[0] < elecNumGoal) muLB =mu;
            double mu_next;
            mu_next = (muUB+muLB)/2.;
            dmu = mu_next - mu;
        }
    }//while
    Nele[0] = FromHkToNele(mu, KS_eigenEnergy);

    return mu;
}
