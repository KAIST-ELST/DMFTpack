#include "mpi.h"
//#include <fstream>
#include "ImgTimeFtn.h"
#include "tight_common.h"
#include "model.h"
#include "TB.h"
#include <Eigen/Eigenvalues>
//void qs_spectra(ImgFreqFtn & SelfE_w ){
void qs_spectra(ImgFreqFtn & SelfE_w,std::vector<Eigen::MatrixXcd>  KS_eigenVectors_orthoBasis,Eigen::VectorXd *KS_eigenEnergy, double  mu,
                std::vector<Eigen::MatrixXcd> H_k_inModelSpace ) {
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
    if (mode == std::string("dos"))    dos(KS_eigenEnergy,   KS_eigenVectors_orthoBasis, std::abs(real(E0)), muDFT);
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
        //version Note1 : "NumOrbit" should be "N_peratom_HartrOrbit"
        //version Note2 : we   dont need to upfolding?
//        double            spectralWeight[N_peratom_HartrOrbit+2];
//        double  spectralWeight_mpiGlobal[N_peratom_HartrOrbit+2];
        Eigen::MatrixXcd  retGkw, retGkw_full;
        double dosnorm[NumHartrOrbit_per_cluster];
        for(int i0=0; i0<NumHartrOrbit_per_cluster; i0++) {
            dosnorm[i0]  = 0;
        }
        FILE *datap1;
        ifroot datap1 = fopen("qsdos.dat", "w");
        ifroot std::cout << "FILEOUT : qsdos.dat\n";


// dos, out of energy window
        Eigen::VectorXd      TdosData;
        Eigen::VectorXd  TdosData_RDC;
        TdosData.setZero(Spectral_EnergyGrid);
        TdosData_RDC.setZero(Spectral_EnergyGrid);
        for(int k=0; k<knum; k++) {
            for(int band=0; band<NumOrbit; band++) {
                if ( (KS_eigenEnergy[k][band]-muDFT) > upper_model_window or  (KS_eigenEnergy[k][band]-muDFT) < lower_model_window) {
//                if ( true) {}
                    for(int index=0; index<Spectral_EnergyGrid; index++) {
                        double energy_i      =  lower_spectrum_window+index*real(dE)  ;
                        double energy_i_high =     energy_i+real(dE)/2.   ;
                        double energy_i_low  =     energy_i-real(dE)/2.   ;
                        /*Lorentzian broadening :\integral 1/x**2 +1 = atan(x)*/
                        if (    std::abs( energy_i- KS_eigenEnergy[k][band] +mu ) <=  50.* infinitesimal  ) {
                            double Aw = (    std::atan( (energy_i_high - KS_eigenEnergy[k][band]+mu)/infinitesimal )
                                             -std::atan( (energy_i_low  - KS_eigenEnergy[k][band]+mu)/infinitesimal )  ) / pi;
                            TdosData[index] += Aw ;
                        }//if
                        else if ( energy_i-KS_eigenEnergy[k][band]+mu > 50.*infinitesimal) break;
                    }/*dos data index*/
                }
            }/*band*/
        }//k
        MPI_Allreduce(TdosData.data(), TdosData_RDC.data(), TdosData.size(), MPI_DOUBLE, MPI_SUM,  MPI_COMM_WORLD);
        TdosData_RDC /= (knum_mpiGlobal * real(dE));



        int itssta, itsend;
        Eigen::VectorXd            spectralWeight(NumHartrOrbit_per_cluster+2);
        Eigen::VectorXd  spectralWeight_mpiGlobal(NumHartrOrbit_per_cluster+2);
//        Eigen::MatrixXcd   retGw(NumHartrOrbit_per_cluster, NumHartrOrbit_per_cluster);
//        Eigen::MatrixXcd   retGw_mpiGlobal(NumHartrOrbit_per_cluster, NumHartrOrbit_per_cluster);
//        retGw.setZero(NumHartrOrbit_per_cluster, NumHartrOrbit_per_cluster);
        for(int itsRank=0 ; itsRank<mpi_numprocs; itsRank++) {
            para_range(0,Spectral_EnergyGrid-1,mpi_numprocs,itsRank, &itssta, &itsend);
            ImgFreqFtn itsSE(0);
            itsSE.realFreq( lower_spectrum_window,  real(dE), (itsend-itssta+1), NumHartrOrbit_per_cluster, NumCorrAtom,0, itssta );
            if(mpi_rank == itsRank) itsSE.update( SelfE_w, 1, 0, 0);
            itsSE.mpiBcast(itsRank, MPI_COMM_WORLD);
            for(int n=itssta; n<=itsend; n++) {
                for(int i0=0; i0<NumHartrOrbit_per_cluster+2; i0++) {
                    spectralWeight(i0)=0;
                    spectralWeight_mpiGlobal(i0)=0;
                }

                Eigen::MatrixXcd Pdos_basis;
                if(impurityBasisSwitch) {
                    Eigen::MatrixXcd projimpurity_site_Ham = impurity_site_Hamiltonian.block(0,0, NumHartrOrbit_per_cluster,NumHartrOrbit_per_cluster);
                    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> ces( NumHartrOrbit_per_cluster);
                    ces.compute(    projimpurity_site_Ham   );
                    Pdos_basis = ces.eigenvectors();
                }
                else {
                    Pdos_basis.setIdentity(NumHartrOrbit_per_cluster, NumHartrOrbit_per_cluster);
                }
                for (int k=0; k< knum; k++) {
                    retarded_GreenFtn2(retGkw_full, retGkw, H_k_inModelSpace, itsSE, mu, k,n);
                    Eigen::MatrixXcd  retGkw_d  = retGkw.block(0,0, NumHartrOrbit_per_cluster,NumHartrOrbit_per_cluster);
                    for(int i0=0; i0<NBAND[k]; i0++)         spectralWeight(NumHartrOrbit_per_cluster+1)    -= imag(retGkw(i0,i0));  //diag1

                    retGkw_d = (Pdos_basis.adjoint() * retGkw_d * Pdos_basis).eval();

                    for(int i0=0; i0<NumHartrOrbit_per_cluster; i0++) {
                        spectralWeight(i0) -= imag(retGkw_d(i0,i0));  //diag2
                        spectralWeight(NumHartrOrbit_per_cluster) -= imag(retGkw_d(i0,i0));  //diag2
                    }



//                    for(int i0=0; i0<NBAND[k]; i0++) retGw += retGkw.block(0,0, NumHartrOrbit_per_cluster,NumHartrOrbit_per_cluster);
                }//k
                MPI_Reduce(spectralWeight.data(), spectralWeight_mpiGlobal.data(), NumHartrOrbit_per_cluster+2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);


//                MPI_Reduce(retGw.data(), retGw_mpiGlobal.data(), retGw.size(), MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD);
//                for(int i0=0; i0<NBAND[k]; i0++)   spectralWeight_mpiGlobal(NumHartrOrbit_per_cluster+1)    -= imag(retGw_mpiGlobal(i0,i0));  //diag1
//                for(int i0=HartrRange[0][0]; i0<HartrRange[0][1]; i0++) {
//                    spectralWeight_mpiGlobal(i0)                    = -imag(retGw_mpiGlobal(i0,i0));  //diag2
//                    spectralWeight_mpiGlobal(NumHartrOrbit_per_cluster) -= imag(retGw_mpiGlobal(i0,i0));  //diag2
//                }

                for(int i0=0; i0<NumHartrOrbit_per_cluster+2; i0++)  spectralWeight_mpiGlobal(i0)/=(pi*knum_mpiGlobal);

                ifroot {
                    fprintf(datap1, "%+0.3f", real(lower_spectrum_window+n*dE)   ); //Pdos1, 2, 3,...Tdos,
                    for(int i0=HartrRange[0][0]; i0<HartrRange[0][1]; i0++) {
                        fprintf(datap1, "    %0.8f", spectralWeight_mpiGlobal(i0)    );
                        dosnorm[i0] += spectralWeight_mpiGlobal(i0);
                    }
                    fprintf(datap1, "      %0.8f", spectralWeight_mpiGlobal(NumHartrOrbit_per_cluster  ) );
                    fprintf(datap1, "      %0.8f", spectralWeight_mpiGlobal(NumHartrOrbit_per_cluster+1)+TdosData_RDC[n] );
                    fprintf(datap1, "\n" );
                }//ifroot
            } //n
            ifroot fflush(datap1);
        }
        for(int i0=0; i0<NumHartrOrbit_per_cluster; i0++) {
            ifroot std::cout << "dosnorm : " << dosnorm[i0]*real(dE) << "\n";
        }
        ifroot fclose(datap1);
        ifroot{
            FILE * GNU = fopen("qsdos.gnuplot", "w");
            fprintf(GNU,"set term x11 dashed\n p \\");
            fprintf(GNU,"\n\"./qsdos.dat\" u 1:($%d) w l  lw 2 lc rgb  \"black\"   lt 1     title \"%d\"   ,\\", NumHartrOrbit_per_cluster+2, 0);
            fprintf(GNU,"\n\"./qsdos.dat\" u 1:($%d) w l  lw 2 lc rgb  \"gray\"   lt 1      title \"%d\"   ,\\", NumHartrOrbit_per_cluster+3, 0);
            fprintf(GNU,"\n\"dos.dat\" u 1:2     w l title \"total\" lc rgb \"black\"");

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
        Eigen::MatrixXcd  retGkw;
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
            ImgFreqFtn itsSE(0);
            itsSE.realFreq( E0,  real(dE), (itsend-itssta+1), N_peratom_HartrOrbit, NumCorrAtom, 0, itssta );
            if(mpi_rank == itsRank) itsSE.update( SelfE_w, 1, 0, 0);
            itsSE.mpiBcast(itsRank, MPI_COMM_WORLD);
            for (int k=0; k< knum; k++) {
                for(int n=itssta; n<=itsend; n++) {
                    retarded_GreenFtn2(retGkw_full, retGkw,  H_k_inModelSpace, itsSE,  mu, k,n);

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
