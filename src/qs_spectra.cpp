#include "mpi.h"
//#include <fstream>
#include "ImgTimeFtn.h"
#include "tight_common.h"
#include "model.h"
#include "TB.h"
#include <Eigen/Eigenvalues>
//void qs_spectra(ImgFreqFtn & SelfE_w ){
void qs_spectra(ImgFreqFtn & SelfE_w,std::vector<Eigen::MatrixXcd>  KS_eigenVectors_orthoBasis,Eigen::VectorXd *KS_eigenEnergy, double  mu,
                std::vector<Eigen::MatrixXcd> H_k_inModelSpace, Eigen::MatrixXcd & SolverBasis) {
    ///*set Swinf used in dos and band*/
    //Eigen::MatrixXcd  Swinf[NumCorrAtom];
    //for(int at1=0; at1<NumCorrAtom; at1++) {
    //    Swinf[at1].setZero(N_peratom_HartrOrbit, N_peratom_HartrOrbit);
    //}
    //for(int at1=0; at1<NumCorrAtom; at1++) {
    //    for(int idx=0;  idx <  Utensor.size(); idx++) {
    //           if ( isOrbitalCorrinHart[Uindex(idx,0)] and
    //                isOrbitalCorrinHart[Uindex(idx,1)] and
    //                isOrbitalCorrinHart[Uindex(idx,2)] and
    //                isOrbitalCorrinHart[Uindex(idx,3)] ) {
    //            int aH = at1 * N_peratom_HartrOrbit + Uindex(idx,0);
    //            int bH = at1 * N_peratom_HartrOrbit + Uindex(idx,1);
    //            int cH = at1 * N_peratom_HartrOrbit + Uindex(idx,2);
    //            int dH = at1 * N_peratom_HartrOrbit + Uindex(idx,3);
    //            Swinf[at1](Uindex(idx,1),Uindex(idx,3)) +=   Utensor[idx] * NumMatrix(aH,cH);
    //            Swinf[at1](Uindex(idx,1),Uindex(idx,2)) -=   Utensor[idx] * NumMatrix(aH,dH);
    //        }
    //    }
    //}



    ifroot std::cout << "start to  real freq calculation..."<<mpi_rank<<"\n";
    if (mode == std::string("dos"))    dos(KS_eigenEnergy,   SolverBasis, std::abs(real(E0)), muDFT);
    if (mode == std::string("optical") ) { //optical cond
//            optical_conductivity(H_k_inModelSpace, SelfE_w, mu, knum, H_R, kmesh,0,0);  //sjh : in-plain
//            optical_conductivity(H_k_inModelSpace, SelfE_w, mu, knum, H_R, kmesh,2,2);  //sjh : out-of-plain
    }//mode 10
    if (mode == std::string("optical") or mode == std::string("qsdos") ) { //qsdos
        ifroot printf("********************************\n");
        ifroot printf("         quasi-ptl dos  \n");
        ifroot printf("********************************\n");
        ifroot std::cout <<"kpoints:"<< knum_mpiGlobal << " pi:" << pi <<"\n";
        ifroot std::cout << "qsdos:Chemical potential : " << mu <<"\n";
        Eigen::MatrixXcd  retGkw, retGkw_full, retGkw0;
        double dosnorm[NumHartrOrbit_per_cluster];
        for(int i0=0; i0<NumHartrOrbit_per_cluster; i0++) {
            dosnorm[i0]  = 0;
        }





        int itssta, itsend;
        Eigen::MatrixXd            spectralWeight(Spectral_EnergyGrid, NumCorrAtom+1);   //pdos_d_atom1; pdos_d_atom2; ... ; tdos
        Eigen::MatrixXd  spectralWeight_mpiGlobal(Spectral_EnergyGrid, NumCorrAtom+1);
//        Eigen::VectorXd            spectralWeight0(NumCorrAtom+1);   //pdos_d_atom1; pdos_d_atom2; ... ; tdos
//        Eigen::VectorXd  spectralWeight_mpiGlobal0(NumCorrAtom+1);
        Eigen::MatrixXd            spectralWeight_full(Spectral_EnergyGrid,num_subshell);   //pdos_d_atom1; pdos_d_atom2; ... ; tdos
        Eigen::MatrixXd  spectralWeight_mpiGlobal_full(Spectral_EnergyGrid,num_subshell);

        for(int itsRank=0 ; itsRank<mpi_numprocs; itsRank++) {   // Bcast w(itsRand) -> w(global)
            ifroot std::cout << itsRank <<"/" << mpi_numprocs <<"\n";
            para_range(0,Spectral_EnergyGrid-1,mpi_numprocs,itsRank, &itssta, &itsend);
            for(int n=itssta; n<=itsend; n++) {
                for(int i0=0; i0<NumCorrAtom+1; i0++) {
                    spectralWeight(n,i0)=0;
                    spectralWeight_mpiGlobal(n,i0)=0;
                }
                for(int lm=0; lm < num_subshell; lm++) {
                    spectralWeight_full(n,lm)=0;
                    spectralWeight_mpiGlobal_full(n,lm)=0;
                }


                for (int k=0; k< knum; k++) {
                    retarded_GreenFtn2(retGkw_full, retGkw, retGkw0, H_k_inModelSpace, SelfE_w, mu, k,n);   //itsSE
                    for(int i0=0; i0<NBAND[k]; i0++)         spectralWeight(n, NumCorrAtom)    -= imag(retGkw(i0,i0));  //diag1

                    Eigen::MatrixXcd  retGkw_d   =  retGkw.block(0,0, NumHartrOrbit,NumHartrOrbit);
//                    Eigen::MatrixXcd  retGkw0_d  = retGkw0.block(0,0, NumHartrOrbit,NumHartrOrbit);
                    retGkw_full =   ( KS_eigenVectors_orthoBasis[k]  * retGkw_full * KS_eigenVectors_orthoBasis[k].adjoint()   ).eval();

                    for(int at=0; at<NumCorrAtom; at++) {
                        for(int i0=0; i0<N_peratom_HartrOrbit; i0++) {
                            int i = at*N_peratom_HartrOrbit+i0;
                            spectralWeight(n,at) -= imag(retGkw_d(i,i));  //diag2
                        }
                    }
                    for(int lm=0; lm < num_subshell; lm++) {
                        for(int i = subshell[lm]; i < subshell[lm+1]; i++) {
                            spectralWeight_full(n,lm) -= imag(retGkw_full(i,i));  //diag2
                        }
                    }
                }//k
                for(int at=0; at< NumCorrAtom+1; at++)  spectralWeight(n,at)/=(pi*knum_mpiGlobal);
                for(int lm=0; lm<  num_subshell; lm++)  spectralWeight_full(n,lm)/=(pi*knum_mpiGlobal);
            } //n
        }
        MPI_Reduce(spectralWeight.data(),           spectralWeight_mpiGlobal.data(), spectralWeight.size(),      MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(spectralWeight_full.data(), spectralWeight_mpiGlobal_full.data(), spectralWeight_full.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

        FILE *datap1;
        FILE *datap3;
        ifroot datap1 = fopen("qsdos.dat", "w");
        ifroot datap3 = fopen("qsdos_lm.dat", "w");
        for(int n= 0 ; n <  Spectral_EnergyGrid ; n++) {
            ifroot {
                fprintf(datap1, "%+0.3f", real(lower_spectrum_window+n*dE)   );          //Pdos1, 2, 3,...Tdos,
                fprintf(datap3, "%+0.3f", real(lower_spectrum_window+n*dE)   );          //Pdos1, 2, 3,...Tdos,
                for(int i0=0; i0<NumCorrAtom; i0++) {
                    fprintf(datap1, "    %0.8f", spectralWeight_mpiGlobal(n,i0)    );
                    dosnorm[i0] += spectralWeight_mpiGlobal(n,i0);
                }
                for(int i0=0; i0< num_subshell; i0++) {
                    fprintf(datap3, "    %0.8f", spectralWeight_mpiGlobal_full(n,i0)    );
                }
                fprintf(datap1, "      %0.8f", spectralWeight_mpiGlobal( n, NumCorrAtom  ) );

                fprintf(datap1, "\n" );
                fprintf(datap3, "\n" );
            }//ifroot
        } //n
        ifroot fflush(datap1);
        ifroot fflush(datap3);
        ifroot std::cout << "FILEOUT : qsdos.dat\n";
        for(int i0=0; i0<NumHartrOrbit_per_cluster; i0++) {
            ifroot std::cout << "dosnorm : " << dosnorm[i0]*real(dE) << "\n";
        }
        ifroot fclose(datap1);
        ifroot fclose(datap3);
        ifroot{
            FILE * GNU = fopen("qsdos.gnuplot", "w");
            fprintf(GNU,"set term x11 dashed\n p \\");
            fprintf(GNU,"\n\"./qsdos.dat\" u 1:($%d) w l  lw 2 lc rgb  \"black\"   lt 1     title \"%d\"   ,\\",  NumCorrAtom+2, 0);
            fprintf(GNU,"\n\"./qsdos0.dat\" u 1:($%d) w l  lw 2 lc rgb  \"black\"   lt 1     title \"%d\"   ,\\", NumCorrAtom+2, 0);

            fprintf(GNU,"\n\npause -1 ");
            fclose(GNU);
        }
    }//mode 100
    if (mode == std::string("band")  )      band(KS_eigenEnergy,muDFT, knum );
    if (mode == std::string("hartreeband")  )     band(H_k_inModelSpace, muTB,   knum);
    if (mode == std::string("qsband") ) {
        band(KS_eigenEnergy, muDFT, knum);
//        band(H_k_inModelSpace, muTB,knum);
        Eigen::MatrixXcd  retGkw_full;
        Eigen::MatrixXcd  retGkw, retGkw0;
        ifroot printf("********************************\n");
        ifroot printf("         quasi-ptl Band  \n");
        ifroot printf("********************************\n");
        ifroot std::cout << "qsband:Chemical potential : " << mu <<"\n";
//        double         spectralWeight[knum][Spectral_EnergyGrid];
//        double partial_spectralWeight[knum][Spectral_EnergyGrid][N_peratom_HartrOrbit];
        std::vector<std::vector<double> >                       spectralWeight(knum);
        std::vector<std::vector<std::vector<double> > > partial_spectralWeight(knum);
        for (int k=0; k< knum; k++) {
            spectralWeight[k].resize(Spectral_EnergyGrid);
            partial_spectralWeight[k].resize(Spectral_EnergyGrid);
            for(int n=0; n<Spectral_EnergyGrid; n++) {
                partial_spectralWeight[k][n].resize(N_peratom_HartrOrbit*NumCorrAtom);
            }
        }
        FILE *datap1;
        int itssta, itsend;
        for(int itsRank=0 ; itsRank<mpi_numprocs; itsRank++) {
            para_range(0,Spectral_EnergyGrid-1,mpi_numprocs,itsRank, &itssta, &itsend);
//            ImgFreqFtn itsSE(0);
//            itsSE.realFreq( E0,  real(dE), (itsend-itssta+1), N_peratom_HartrOrbit, NumCorrAtom, 0, itssta );
//            if(mpi_rank == itsRank) itsSE.update( SelfE_w, 1, 0, 0);
//            itsSE.mpiBcast(itsRank, MPI_COMM_WORLD);
            for (int k=0; k< knum; k++) {
                for(int n=itssta; n<=itsend; n++) {
                    retarded_GreenFtn2(retGkw_full, retGkw, retGkw0,  H_k_inModelSpace, SelfE_w,  mu, k,n);  //itsSE

                    spectralWeight[k][n]=0;
                    for(int i0=0; i0<NumOrbit; i0++) {
                        spectralWeight[k][n] -= imag( retGkw_full(i0,i0) );
                    }
                    spectralWeight[k][n]/=pi;

                    for(int at=0; at<NumCorrAtom; at++) {
                        for(int i0=0; i0<N_peratom_HartrOrbit; i0++) {
                            int i = at*N_peratom_HartrOrbit+i0;
                            partial_spectralWeight[k][n][i] = -imag(retGkw(HartrIndex[i],HartrIndex[i])) /pi;
                        }
                    }
                }//n
            }//k
        }//its

        ifroot    datap1 = fopen("qsband.dat", "w");
        ifroot    fclose(datap1);
        for(int itsRank=0 ; itsRank<mpi_numprocs; itsRank++) {
            if(mpi_rank==itsRank) {
                datap1 = fopen("qsband.dat", "a");
                for (int k=0; k< knum; k++) {
                    for(int n=0; n<Spectral_EnergyGrid; n++) {
                        fprintf(datap1, " %d   %0.8f   %+0.3f    %0.8f", k+myksta, kdist_band[k+myksta], real(E0+n*dE),  spectralWeight[k][n]  );    //for k, n
                        for(int orb=0; orb<N_peratom_HartrOrbit*NumCorrAtom; orb++)  fprintf(datap1, "    %+0.8f", partial_spectralWeight[k][n][orb]);
                        fprintf(datap1, "\n" );
                    }//n
                    fprintf(datap1, "\n" );
                }//k
                fclose(datap1);
            }
            MPI_Barrier(MPI_COMM_WORLD);
            sleep(1);
        }
    }//mode11
}
