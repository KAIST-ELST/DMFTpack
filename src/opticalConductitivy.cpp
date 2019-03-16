//Peierls approximation, calculate sigma_a, sigma_b, sigma_c,  where a,b and c are the latt const direction
//Ref : arXiv1203.5711 and other materials for linear response, single band optical conductiviy in DMFT.
//151209 Jae-Hoon Sim
//160203 Jae-Hoon Sim
//160409 Jae-Hoon Sim


//2019Mar09 Jae-Hoon Sim


//
//Options:
//SOLVERTYPE = 0
//mode = 10
//Spectral_EnergyWindow =  max freq in realFreq_Sw.dat
//Spectral_EnergyGrid   =  realFreq_Sw.dat line number -1
//Don't worry about orbital=Sw_inf=0;

//#include "./eigen-eigen-b30b87236a1b/Eigen/Eigenvalues"
#include "./eigen-eigen-dc6cfdf9bcec/Eigen/Core"
#include "./eigen-eigen-dc6cfdf9bcec/Eigen/Eigenvalues"
#include "tight_common.h"
#include "mpi.h"
#include "TB.h"
#include <time.h>



void optical_conductivity(cmplx ***H_k_inDualBaseSet, ImgFreqFtn Sw, double mu, int knum, cmplx ***** H_R, double **kmesh,int alp_direc, int beta_direc) {
    int Norbital_spatial = NumOrbit/Spin_DegreeOfFreedom;
    assert(alp_direc<dimension and  beta_direc<dimension);
    assert(alp_direc==beta_direc); // Hall res. cannot be obtained
    ifroot std::cout << "********************************\n";
    ifroot std::cout << "      Optical Conductivity " << alp_direc<<beta_direc<<"\n"     ;
    ifroot std::cout << "********************************\n";
    double Re_sigma[N_freq], DosDos[N_freq],kprojDos[N_freq];
    //note  w=n*dw for \sigma1,
    // w = 0 ~ 2*w0  ==> n< N_freq

    double Re_sigma_mpiloc[N_freq], DosDos_mpiloc[N_freq];
    double fermi_dist_Ftn[N_freq];
    double w,  OpticalCond_t2g, OpticalCond_eg, OpticalCond_tot, OpticalCond_yzxz;
    int nPlus, nMinus, k,j, m1,m2,m10,m20 , R1,R2,R3, ntemp, ltemp, mtemp ;
    int subDim = NumOrbit / NumSubSpace;
    int submat,tempint1, tempint2;
    time_t timeStart, timeEnd;
    int mysta, myend;
    cmplx tempc, tempc2;
    Eigen::MatrixXcd ** OpticalCond_transMat = new Eigen::MatrixXcd * [N_freq];
    for(int n=0; n<N_freq; n++) {
        OpticalCond_transMat[n]= new Eigen::MatrixXcd  [NumSubSpace];
        for(int submat=0; submat < NumSubSpace; submat++) {
            OpticalCond_transMat[n][submat]= Eigen::MatrixXcd::Zero(subDim, subDim);
        }
    }

    if(imag(w0) <  3*real(dw) ) {
        ifroot std::cout <<  "###############################################\n";
        ifroot std::cout <<  "Warning: Use the infinitesimal larger than dw\n";
        ifroot std::cout <<  "###############################################\n";
    }
    timeStart = clock();
    ifroot std::cout << "Chemical potential : " << mu << "\n";

//read atom position.
//you can take   please note that the atomic position should be ginven as Cartesian form

    double  LatticeUnitLength[dimension] ;
    double siteDistance;
    int atomicCoordi =0 ;  //0=FRAC, 1=Ang
    double LatticeVector[dimension][dimension];


    cmplx retGkw[NumOrbit*NumOrbit];
    double atomicPosition[NumOrbit/Spin_DegreeOfFreedom][dimension];
    double atomicPosition_direct[NumOrbit/Spin_DegreeOfFreedom][dimension];
    /***********************************************************
    Read lattice information
    ************************************************************/
    FILE * dataIN  = fopen("lattice.dat", "r");
    if(dataIN == NULL) {
        std::cout << "Cannot read 'lattice.dat' file...\n";
        assert(0);
        exit(1);
    }
    fscanf(dataIN,"%d\n",&atomicCoordi );
    ifroot {
        if (atomicCoordi==0)       std::cout << "Fractional  coordinate\n";
        else if (atomicCoordi==1)  std::cout << "Direct(ang) coordinate\n";
        else exit(1);
        assert(atomicCoordi==0 or atomicCoordi==1);
    }
    for (int k=0; k<dimension; k++) {
        fscanf(dataIN,"%lf  %lf %lf\n",&LatticeVector[k][0],&LatticeVector[k][1],&LatticeVector[k][2] );
        ifroot std::cout << LatticeVector[k][0] <<"  "<<LatticeVector[k][1] <<"  "<<LatticeVector[k][2] <<"\n" ;
    }
    LatticeUnitLength[0]=std::sqrt(std::pow(LatticeVector[0][0],2)+std::pow(LatticeVector[0][1],2)+std::pow(LatticeVector[0][2],2));
    LatticeUnitLength[1]=std::sqrt(std::pow(LatticeVector[1][0],2)+std::pow(LatticeVector[1][1],2)+std::pow(LatticeVector[1][2],2));
    LatticeUnitLength[2]=std::sqrt(std::pow(LatticeVector[2][0],2)+std::pow(LatticeVector[2][1],2)+std::pow(LatticeVector[2][2],2));

    for (int k=0; k<NumOrbit/Spin_DegreeOfFreedom; k++) {
        fscanf(dataIN,"%lf  %lf %lf\n",&atomicPosition[k][0],&atomicPosition[k][1],&atomicPosition[k][2] );
        ifroot std::cout << atomicPosition[k][0] <<"  "<<atomicPosition[k][1] <<"  "<<atomicPosition[k][2] <<"\n" ;
    }
    for (int k=0; k<NumOrbit/Spin_DegreeOfFreedom; k++) {
        for (int i=0; i<dimension; i++) {
            if(atomicCoordi==0)
                atomicPosition_direct[k][i] = atomicPosition[k][0]*LatticeVector[0][i]+
                                              atomicPosition[k][1]*LatticeVector[1][i]+
                                              atomicPosition[k][2]*LatticeVector[2][i];
            else if(atomicCoordi==1)  atomicPosition_direct[k][i] = atomicPosition[k][i];
        }//i
        for (int l=0; l<k; l++) {
            siteDistance =0;
            for (int i=0; i<dimension; i++) {
                siteDistance+= std::pow((atomicPosition_direct[k][i] - atomicPosition_direct[l][i]) ,2);
            }
            siteDistance=std::sqrt(siteDistance);
            if(mpi_rank==0 and   FromOrbitalToAtom[2*k]!=FromOrbitalToAtom[2*l]    and  (siteDistance<1. or siteDistance > 10. )) {
                std::cout << "****************************************************\n";
                std::cout << "Warning! Too small distance, "<<siteDistance<<", between\n"
                          << FromOrbitalToAtom[2*k]+1<<"th atom at "
                          <<       atomicPosition_direct[k][0]
                          <<", "<< atomicPosition_direct[k][1]
                          <<", "<< atomicPosition_direct[k][2] <<" and\n"
                          << FromOrbitalToAtom[2*l]+1<<"th atom at "
                          << atomicPosition_direct[l][0]
                          <<", "<< atomicPosition_direct[l][1]
                          <<", "<< atomicPosition_direct[l][2] <<"\n";
                std::cout << "****************************************************\n";
            }
        }//l
    }//k

    /***********************************************************
    Group velocity:
    ************************************************************/
    Eigen::MatrixXcd ** groupVel_W = new Eigen::MatrixXcd *  [dimension];
    for (j=0; j<dimension; j++) {
        groupVel_W[j] = new Eigen::MatrixXcd   [NumSubSpace];
        for(submat=0; submat < NumSubSpace; submat++) {
            groupVel_W[j][submat] =  Eigen::MatrixXcd::Zero(subDim,subDim);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);




//Initialize Egkw, groupVelocity and
    Eigen::MatrixXcd ** Akw = new Eigen::MatrixXcd *  [N_freq];
    for(int n=0; n<N_freq; n++) {
        kprojDos[n] = 0;
        Akw[n] = new  Eigen::MatrixXcd   [NumSubSpace];
        for(submat=0; submat < NumSubSpace; submat++) {
            Akw[n][submat] =  Eigen::MatrixXcd::Zero(subDim,subDim);
        }
    }
    for(int n=0; n<N_freq; n++) {
        Re_sigma[n] =0;
        Re_sigma_mpiloc[n] = 0;
        DosDos[n] =0;
        DosDos_mpiloc[n] = 0;
    }
    for(int n=0; n<N_freq; n++) {
        w= real(w0+n*dw);
        fermi_dist_Ftn[n] = 1.0/(1.+exp(beta*w)) ;
    }


    ifroot std::cout << "NumSubSpace:"<<NumSubSpace<<" (should be =2)\n";
    ifroot std::cout << "Start Main loop\n";




    /*overlap matrix with position operator: beyond Peierls*/


    cmplx ** S_overlap_position_R ;
    int beyond_peierls=1;
    if(beyond_peierls) {
        FILE *Overlap;
        if (alp_direc==0) Overlap= fopen("OverlapMatrix_x.HWR","r");
        if (alp_direc==1) Overlap= fopen("OverlapMatrix_y.HWR","r");
        if (alp_direc==2) Overlap= fopen("OverlapMatrix_z.HWR","r");
        S_overlap_position_R = new cmplx * [Norbital_spatial*Norbital_spatial];
        for (int i=0; i< Norbital_spatial*Norbital_spatial; i++) {
            S_overlap_position_R[i] = new cmplx [NumLattxyz];
        }
        for (int i=0; i< Norbital_spatial; i++) {
            for (int j=0; j<Norbital_spatial; j++) {
                for (int n=0; n < 2*BraLatt_x+1; n++) {
                    for (int l=0; l < 2*BraLatt_y+1; l++) {
                        for (int m=0; m < 2*BraLatt_z+1; m++) {
                            S_overlap_position_R[i*Norbital_spatial+j][(n)+ (l)*NumLattx  + (m) * NumLattxy ]= 0 ;//0+I*0;
                        }
                    }
                }
            }
        }
        int i,j,Atom1,Atom2, n,l,m;
        double  HopRe,HopIm;
        while(!feof(Overlap)) {
            fscanf(Overlap,"%d %d %d   %d %d    %d %d     %lf  %lf\n",&n,&l,&m,  &Atom1, &Atom2,    &i,&j,   &HopRe, &HopIm);
            if (n >= -BraLatt_x && l >= -BraLatt_y && m >= -BraLatt_z && n <= BraLatt_x && l <= BraLatt_y && m <= BraLatt_z  ) {
                i0=(accumulated_Num_SpinOrbital[Atom1-1]/Spin_DegreeOfFreedom) +(i-1);
                m0=(accumulated_Num_SpinOrbital[Atom2-1]/Spin_DegreeOfFreedom) +(j-1);
                S_overlap_position_R[(i0)*Norbital_spatial+(m0)][(n+BraLatt_x)+(l+BraLatt_y)*NumLattx +(m+BraLatt_z)*NumLattxy] = (HopRe + I*HopIm);
            }
        }
        if (mpi_rank ==0)  cout << "We have Overlap matrix"  <<"\n";
        fclose(Overlap);
    }




    /***********************************************************
    optical conductivity main loop
    ************************************************************/
    ifroot timeStart = clock();
    for(int k=0; k<knum; k++) {
//construct groupVel_W

        for(submat=0; submat < NumSubSpace; submat++) {
            /*Peierls approximation*/
            for (int m10=0; m10<subDim; m10++) {
                for (int m20=0; m20<subDim; m20++) {
                    groupVel_W[alp_direc][submat](m10,m20) =0;
                }
            }
            assert((groupVel_W[alp_direc][submat].adjoint() - groupVel_W[alp_direc][submat]).norm() < 10e-10); //hermiticity
            /*additional contribution*/
            if( beyond_peierls) {
                cmplx  * H_k_inHybridSet = new cmplx [subDim*subDim];
                cmplx  * S_overlap_position_k = new cmplx [NumOrbit*NumOrbit];
                for (int m10=0; m10<subDim; m10++) {
                    int m1=(NumSubSpace*m10+submat);
                    for (int m20=0; m20<subDim; m20++) {
                        int m2=(NumSubSpace*m20+submat);
                        H_k_inHybridSet[m10*subDim+m20]=0;
                        for (int m30=0; m30<subDim; m30++) {
                            int m3=(NumSubSpace*m30+submat);
                            H_k_inHybridSet[m10*subDim+m20] += invS_overlap_k[k][m1*NumOrbit+m3] * H_k_inDirect[k][m3][m2];
                        }
                    }
                }

                /*overlap with position in (k,dual)-space */
                for (int i0=0; i0<NumOrbit; i0++) {
                    for (int m0=0; m0<NumOrbit; m0++) {
                        S_overlap_position_k[i0*NumOrbit+m0]=0;
                    }
                }

                for (int m10=0; m10<Norbital_spatial; m10++) {
                    for (int m20=0; m20<Norbital_spatial; m20++) {

                        S_overlap_position_k[(2*m10)*NumOrbit+(2*m20)]=0;


                        for (int m30=0; m30<Norbital_spatial; m30++) {
                            for (int m40=0; m40<Norbital_spatial; m40++) {
                                cmplx temp=0;   //=S_overlap_position_k[m30_m40] in Direct space
                                {
                                    for (int n=0; n < NumLattx; n++) {
                                        ntemp = n-BraLatt_x;
                                        for (int l= 0; l < NumLatty; l++) {
                                            int ltemp = l-BraLatt_y;
                                            for(int m=0; m<NumLattz; m++) {
                                                int mtemp = m-BraLatt_z;
                                                temp +=
                                                    S_overlap_position_R[(m30)*Norbital_spatial+(m40)][(n)+(l)*NumLattx+(m)*NumLattxy]
                                                    * exp ( -I*((kmesh[k][0]*ax*ntemp)+(kmesh[k][1]*ay*ltemp)+ (kmesh[k][2]*az*mtemp)) )  ;
                                            }
                                        }
                                    }
                                }
                                S_overlap_position_k[(2*m10)*NumOrbit+(2*m20)]
                                += invS_overlap_k[k][2*m10*NumOrbit+2*m30] * temp *invS_overlap_k[k][2*m40*NumOrbit+2*m20];

                            }//m40
                        }//m3

                        S_overlap_position_k[(2*m10+1)*NumOrbit+(2*m20+1)]= S_overlap_position_k[(2*m10)*NumOrbit+(2*m20)];
                    }
                }


                for (int m10=0; m10<subDim; m10++) {
                    m1=(NumSubSpace*m10+submat);
                    for (int m20=0; m20<subDim; m20++) {
                        m2=(NumSubSpace*m20+submat);
                        for (int m30=0; m30<subDim; m30++) {
                            int m3=(NumSubSpace*m30+submat);
                            groupVel_W[alp_direc][submat](m10,m20) += H_k_inHybridSet[m10*subDim+m30] * S_overlap_position_k[m3*NumOrbit+m2];
                        }
                        groupVel_W[alp_direc][submat](m10,m20) *= I;
                    }
                }
                for (int m10=0; m10<subDim; m10++) {
                    for (int m20=0; m20<subDim; m20++) {
                        groupVel_W[alp_direc][submat](m10,m20) += conj(groupVel_W[alp_direc][submat](m20,m10));
                    }
                }
                delete [] H_k_inHybridSet;
                delete [] S_overlap_position_k;
            }//beyond Peierls


           /*Peierls part*/
            for (int m10=0; m10<subDim; m10++) {
                m1=(NumSubSpace*m10+submat);
                for (int m20=0; m20<subDim; m20++) {
                    m2=(NumSubSpace*m20+submat);

                    for (int R1=0; R1<2*BraLatt_x+1; R1++) {
                        ntemp = R1-BraLatt_x;
                        for (int R2=0; R2<2*BraLatt_y+1; R2++) {
                            ltemp = R2-BraLatt_y;
                            for (int R3=0; R3<2*BraLatt_z+1; R3++) {
                                mtemp = R3-BraLatt_z;
                                groupVel_W[alp_direc][submat](m10,m20) +=
                                    ( I*
                                      ( (-atomicPosition_direct[m1/Spin_DegreeOfFreedom][alp_direc]
                                         +atomicPosition_direct[m2/Spin_DegreeOfFreedom][alp_direc])
                                        -( ntemp*LatticeVector[0][alp_direc]
                                           +ltemp*LatticeVector[1][alp_direc]
                                           +mtemp*LatticeVector[2][alp_direc])
                                      )* H_R[m1][m2][R1][R2][R3]
                                      *exp ( -I*((kmesh[k][0]*ax*ntemp)+(kmesh[k][1]*ay*ltemp)+(kmesh[k][2]*az*mtemp)))
                                    );
                            }//R3
                        }//R2
                    }//R1

                }//m2
            }//m1



        }//submat
        if(mpi_rank==0 and k==0) {
            timeEnd=clock();
            std::cout << (((timeEnd-timeStart))/(CLOCKS_PER_SEC))/3600 <<":"
                      << ((((timeEnd-timeStart))/(CLOCKS_PER_SEC))%3600)/60 << ":"
                      << ((((timeEnd-timeStart))/(CLOCKS_PER_SEC))%3600)%60
                      <<"for group_velocity construction\n";
            timeStart = clock();
        }

//construct Akw
        for(int n=0; n<N_freq; n++) {
            retarded_GreenFtn(retGkw, H_k_inDirect, Sw, mu, k,n, NumSubSpace);
            for(submat=0; submat < NumSubSpace; submat++) {
                for (m10=0; m10<subDim; m10++) {
                    m1=(NumSubSpace*m10+submat);
                    for (m20=0; m20<subDim; m20++) {
                        m2=(NumSubSpace*m20+submat);
                        Akw[n][submat](m10,m20) =   I * (retGkw[m1*NumOrbit+m2] -std::conj(retGkw[m2*NumOrbit+m1])) /(2* pi);  //-Im(G) = I*(G-G*) /2
                    }
                }//m10
                assert((Akw[n][submat].adjoint() - Akw[n][submat]).norm() < 10e-10);
            }//submat
        }//n


        for(int n=0; n<N_freq; n++) {
            for(int submat=0; submat < NumSubSpace; submat++) {
                OpticalCond_transMat[n][submat]     = (Akw[n][submat].selfadjointView<Eigen::Upper>()*groupVel_W[alp_direc][submat]);
            }//submat
        }
        if(mpi_rank==0 and k==0) {
            timeEnd=clock();
            std::cout << (((timeEnd-timeStart))/(CLOCKS_PER_SEC))/3600 <<":"
                      << ((((timeEnd-timeStart))/(CLOCKS_PER_SEC))%3600)/60 << ":"
                      << ((((timeEnd-timeStart))/(CLOCKS_PER_SEC))%3600)%60
                      <<"for transition Matrix\n";
            timeStart = clock();
        }
        for(int n=0; n<N_freq; n++) {
            //w= real(n*dw) = 0~2*w(n_max);   ==>      w' + w = (-w0+n0*dw) +  (n*dw) =  -w0+  ( n0 + n   ) * dw
            //
            //w' + w = -w0 + (n0+n)*dw > -kT  ==>  n0 > N_freq/2 -(kT/dw) - n
            //w'     = -w0 +  n0   *dw  < +kT ==>  n0 < N_freq/2 +(kT/dw)
            tempint1 = std::max( (int)(N_freq/2. -  (5./beta)/real(dw) -n),0);
            tempint2 = std::min( (int)(N_freq/2. +  (5./beta)/real(dw)), N_freq-n);
            for(n0=tempint1; n0<tempint2; n0++) { //sum_wprime   from -w(n_max)+w/2  to w(n_max) - w/2
                nMinus =  n0 ;     //wp
                nPlus =   n0+n ;   //wp + w
                assert(fermi_dist_Ftn[nMinus] >=  fermi_dist_Ftn[nPlus]);
                OpticalCond_tot =0;
                for(int submat=0; submat < NumSubSpace; submat++) {
                    OpticalCond_tot  = real((OpticalCond_transMat[nPlus][submat] * OpticalCond_transMat[nMinus][submat]).trace());
                    Re_sigma_mpiloc[n] +=   OpticalCond_tot *(fermi_dist_Ftn[nMinus] - fermi_dist_Ftn[nPlus]) ;
                }//submat
            }//n0  = wprim
        }//n
        if(mpi_rank==0 and k==0) {
            timeEnd=clock();
            std::cout << (((timeEnd-timeStart))/(CLOCKS_PER_SEC))/3600 <<":"
                      << ((((timeEnd-timeStart))/(CLOCKS_PER_SEC))%3600)/60 << ":"
                      << ((((timeEnd-timeStart))/(CLOCKS_PER_SEC))%3600)%60
                      <<"for Re[sigma]\n";
            std::cout << "we have "<<knum<<" k-points for each processors\n";
            std::cout << "..then.."<<(timeEnd-timeStart)/(CLOCKS_PER_SEC)<<"*"<<knum<<"="<< (timeEnd-timeStart)/(CLOCKS_PER_SEC)*knum <<"\n";
            timeStart = clock();
        }
    } //k
//////////////////////////////////////////////////////////////////////////////////////////////////
    std::cout << mpi_rank <<"OpticalCond main loop finished.\n";
    for(int n=0; n<N_freq; n++) {
        w= real(n*dw);
        MPI_Reduce(Re_sigma_mpiloc+n, Re_sigma+n, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Bcast (Re_sigma+n,        1           ,MPI_DOUBLE, 0, MPI_COMM_WORLD);

        Re_sigma[n]  *= (2*pi) / ((k_pointx*LatticeUnitLength[0]*k_pointy*LatticeUnitLength[1]*k_pointz*LatticeUnitLength[2]) *w)
                        * conductivity_fromUnit_to_invOhmcm ;
    }//n
    timeEnd = clock(); //time
    ifroot std::cout << "end of optical conductivity cal " << (timeEnd - timeStart)/(CLOCKS_PER_SEC) << "\n";



    if(mpi_rank ==0) {
        FILE *datap1;
        char chr[100];
        sprintf(chr, "%d%d", alp_direc,beta_direc);

///////////////////////////////////////////////////////
        std::string filename = std::string("opticalConductivity") + chr+".dat";
        datap1 = fopen(filename.c_str()       , "w");
        for( int n=0; n<N_freq; n++) {
            fprintf(datap1, "%+0.3f    %0.8f\n" , real(n*dw) ,  Re_sigma[n]    );
        } //n
        fclose(datap1);
    }



    for(int n=0; n<N_freq; n++) {
        delete [] OpticalCond_transMat[n];
    }
    delete [] OpticalCond_transMat;

    for (j=0; j<dimension; j++) {
        delete [] groupVel_W[j];
    }
    for (j=0; j<N_freq; j++) {
        delete [] Akw[j];
    }
    delete [] groupVel_W;
    delete [] Akw;
}

