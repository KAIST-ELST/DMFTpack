#include <fstream>
#include   "tight_common.h"
#include   "model.h"
#include   "TB.h"
#include <Eigen/Eigenvalues>
//
#define   history                  2
#define  highFreq2 2
#define valenceHAMILTONIAN(n,k)                                                                                                                                           \
{                                                                                                                                                                         \
	Eigen::MatrixXcd SWLOCAL;\
	SWLOCAL.setZero(NumCorrAtom*N_peratom_HartrOrbit, NumCorrAtom*N_peratom_HartrOrbit);\
	if(n<N_freq*highFreq2){\
		cmplx   iw=I*pi*(2.*n+1.)/beta;\
		if(n<N_freq){\
			for(int a=0; a<N_peratom_HartrOrbit*NumCorrAtom; a++){\
				for(int b=0; b<N_peratom_HartrOrbit*NumCorrAtom; b++){\
					SWLOCAL(a,b) = SelfE_w.getValue(n,a,b);\
				}\
			}\
		}\
		else{\
			for(int a=0; a<N_peratom_HartrOrbit*NumCorrAtom; a++){\
				for(int b=0; b<N_peratom_HartrOrbit*NumCorrAtom; b++){\
					SWLOCAL(a,b)   = moments[0](a,b);\
					SWLOCAL(a,b)  += moments[1](a,b)/(iw);\
					SWLOCAL(a,b)  += moments[2](a,b)/std::pow(iw,2);\
					SWLOCAL(a,b)  += moments[3](a,b)/std::pow(iw,3);\
				}\
			}\
		}\
		for(int cl1=0; cl1<NumCluster; cl1++) {\
			for(int p=cl1*NumHartrOrbit_per_cluster ; p < (cl1+1)*NumHartrOrbit_per_cluster ; p++){\
			for(int q=cl1*NumHartrOrbit_per_cluster ; q < (cl1+1)*NumHartrOrbit_per_cluster ; q++){\
					model_qp_hamiltonian[HartrIndex[p]*NBAND[k]+HartrIndex[q]] += \
					SWLOCAL(p,q);\
				}\
			}\
		}\
	}\
	else{\
		for(int a=0; a<N_peratom_HartrOrbit*NumCorrAtom; a++){\
			for(int b=0; b<N_peratom_HartrOrbit*NumCorrAtom; b++){\
				SWLOCAL(a, b) = moments[0](a,b);\
			}\
		}\
		for(int cl1=0; cl1<NumCluster; cl1++) {\
			for(int p=cl1*NumHartrOrbit_per_cluster ; p < (cl1+1)*NumHartrOrbit_per_cluster ; p++){\
			for(int q=cl1*NumHartrOrbit_per_cluster ; q < (cl1+1)*NumHartrOrbit_per_cluster ; q++){\
					model_qp_hamiltonian[HartrIndex[p]*NBAND[k]+HartrIndex[q]] += \
					SWLOCAL(p,q);\
				}\
			}\
		}\
	}\
}while(0);




double GiwToNele (double muTB, ImgFreqFtn & SelfE_w, cmplx *** eigenVal, int preComputation,std::vector<Eigen::MatrixXcd> &H_k_inModelSpace  );
//double traceNorm_full( cmplx *A, int n,int dim ) ;
//double traceNorm_diag( cmplx *A, int n,int dim ) ;
//void traceTest(cmplx **Gw) ;

double GreenFtn_w_adjust_mu( std::vector<Eigen::MatrixXcd> &   H_k_inModelSpace, ImgFreqFtn & SelfE_w, double mu,  double dmu, int mu_adjustLOCAL     ) {
    time_t timeStartGreen, timeEndGreen;
    timeStartGreen = clock();
    ifroot printf(" Adjusting Chemical Potential...\n");
//    int i,j,k,l,m,n,tau;
    cmplx *** eigenVal = new cmplx ** [N_freq*highFreq2+1];
    for (int n=0; n<N_freq*highFreq2+1; n++) {
        eigenVal[n] = new cmplx * [knum];
        for (int k=0; k<knum; k++) {
            eigenVal[n][k] = new cmplx  [NBAND[k]+1];
            eigenVal[n][k][NBAND[k]] = (cmplx) magicNumber;
        }
    }
    ifroot   printf("Num of electrons(ideal) : %2.8f \n", NumberOfElectron);

    std::vector<double> Nele(history);
    double muInput = mu;
    for (int i=0; i<history; i++)  Nele[i]=-1;
    Nele[0] =  GiwToNele(mu,SelfE_w,eigenVal,0, H_k_inModelSpace); ;
    assert( not( std::isnan(std::abs(Nele[0]))  or   std::isinf(std::abs(Nele[0]))      ));

    /*Find chemical potential*/
    //        double muRD;
    double mu_next;
    double muUB=0, muLB=0;
    double elecNumGoal = NumberOfElectron;
    int nearAns=0;
    if(Nele[0] > NumberOfElectron) mu-=(dmu-dmu);
    else if(Nele[0] < NumberOfElectron) mu-=(dmu+dmu);
    int step=0;
    while (  (Nele[0] != elecNumGoal)   and  std::abs(dmu)>1e-5 ) {
        mu += dmu;
        for(int i=history-1; i>0 ; i--) Nele[i]=Nele[i-1];
        Nele[0] = GiwToNele(mu, SelfE_w,  eigenVal, 1, H_k_inModelSpace  );

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
            mu_next = (muUB+muLB)/2.;
            dmu = mu_next - mu;
        }
    }//while
    timeEndGreen = clock(); //time
    if( mpi_rank==0 )  printf("Adjusted Chemical potential : From %2.5f to %2.5f  (RD= %+2.5f )\n",muInput, mu,  mu-muInput) ;
    Nele[0] = GiwToNele( mu, SelfE_w,  eigenVal, 1, H_k_inModelSpace  );

    for (int n=0; n<N_freq*highFreq2+1; n++) {
        for (int k=0; k<knum; k++) {
            assert ( eigenVal[n][k][NBAND[k]] == (cmplx) magicNumber);
        }
    }
    for (int n=0; n<N_freq*highFreq2+1; n++) {
        for (int k=0; k<knum; k++) {
            delete [] eigenVal[n][k];
        }
        delete [] eigenVal[n];
    }
    delete [] eigenVal;
    return mu;
}


void GreenFtn_w(  int NumCluster, int NumHartrOrbit_per_cluster, std::vector<Eigen::MatrixXcd> &  H_k_inModelSpace,
                  ImgFreqFtn & SelfE_w, std::vector<Eigen::MatrixXcd>  & Gw, double mu,                  std::vector<Eigen::MatrixXcd> & densityMatDFT
               ) {
    ifroot  printf("*GreenFtn calculation\n");





    Gw.resize(N_freq);
    for(int n=0; n < N_freq; n++) {
        Gw[n].setZero( NumCluster*NumHartrOrbit_per_cluster, NumCluster*NumHartrOrbit_per_cluster);
    }
    int COL = NumHartrOrbit_per_cluster * NumHartrOrbit_per_cluster;

    Eigen::MatrixXcd Gw_loc_mpiLoc (NumCluster, N_freq*COL);
    Eigen::MatrixXcd Gw_loc        (NumCluster, N_freq*COL);
    Gw_loc_mpiLoc.setZero(NumCluster, N_freq*COL);
    Gw_loc.setZero(NumCluster, N_freq*COL);

    for(int k=0; k<knum; k++) {
        densityMatDFT[k].setZero(NBAND[k],NBAND[k]);
    }

    Eigen::MatrixXcd NumMat_loc (NumCluster*NumHartrOrbit_per_cluster, NumCluster*NumHartrOrbit_per_cluster);
    Eigen::MatrixXcd NumMat_loc0 (NumCluster*NumHartrOrbit_per_cluster, NumCluster*NumHartrOrbit_per_cluster);
    NumMat_loc.setZero(NumCluster*NumHartrOrbit_per_cluster, NumCluster*NumHartrOrbit_per_cluster);
    NumMat_loc0.setZero(NumCluster*NumHartrOrbit_per_cluster, NumCluster*NumHartrOrbit_per_cluster);
    NumMatrix.setZero(NumCluster*NumHartrOrbit_per_cluster, NumCluster*NumHartrOrbit_per_cluster);

    /* Construct Gw_loc : main loop*/
    time_t timeStart, timeEnd;
    ifroot std::cout << "To invert low-evergy model matrix, "<<knum*N_freq*highFreq2*NumSubSpace << "times...:\n";
    timeStart=clock();
    Eigen::MatrixXcd * staticHamiltonian_Model= new Eigen::MatrixXcd [knum];
    Eigen::VectorXcd * static_ev              = new Eigen::VectorXcd [knum];
    /*Gkw0*/


    Eigen::MatrixXcd * SelfEnergy_Eig = new Eigen::MatrixXcd [N_freq];
    std::vector<Eigen::MatrixXcd >  Swmoments;
    for(int n=0; n < N_freq; n++) {
        SelfEnergy_Eig[n] = SelfE_w.getMatrix(n);
    }
    getAsymto_moments(Swmoments, SelfEnergy_Eig);
    delete [] SelfEnergy_Eig;

    for(int k=0; k<knum; k++) {
        static_ev[k].setZero(NBAND[k]);
        staticHamiltonian_Model[k] = H_k_inModelSpace[k];
        for(int i0=0; i0<NBAND[k]; i0++) {
            for(int m0=0; m0<NBAND[k]; m0++) {
                if(isOrbitalHartr[i0] and isOrbitalHartr[m0] and isSameAtom(i0,m0) ) {
                    staticHamiltonian_Model[k](i0,m0)  +=  Swmoments[0](i0,m0);
                }
            }
        }
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> ces( NBAND[k] );
        ces.compute(staticHamiltonian_Model[k]);

        for(int i0=0; i0<NBAND[k]; i0++) {
            static_ev[k](i0) = ces.eigenvalues()[i0];
        }
        staticHamiltonian_Model[k] = ces.eigenvectors();
    }


    for(int n=0; n<N_freq*highFreq2; n++) {
        cmplx   iw=I*pi*(2.*n+1.)/beta;
        //Sw for given freq. n
        Eigen::MatrixXcd SWLOCAL;
        SWLOCAL.setZero(NumHartrOrbit_per_cluster*NumCluster,
                        NumHartrOrbit_per_cluster*NumCluster);
        if(n<N_freq) {
            for(int a=0; a<NumCluster*NumHartrOrbit_per_cluster; a++) {
                for(int b=0; b<NumCluster*NumHartrOrbit_per_cluster; b++) {
                    SWLOCAL(a, b) = SelfE_w.getValue(n,a,b);
                }
            }
        }
        else {
            for(int a=0; a<NumCluster*NumHartrOrbit_per_cluster; a++) {
                for(int b=0; b<NumCluster*NumHartrOrbit_per_cluster; b++) {
                    SWLOCAL(a,b)   = Swmoments[0](a,b);
                    SWLOCAL(a,b)  += Swmoments[1](a,b)/(iw);
                    SWLOCAL(a,b)  += Swmoments[2](a,b)/std::pow(iw,2);
                    SWLOCAL(a,b)  += Swmoments[3](a,b)/std::pow(iw,3);
                }
            }
        }

        for(int k=0; k<knum; k++) {
            Eigen::MatrixXcd Gkw;  // =(iwn+mu -H - S(iwn))^{-1}
            Eigen::MatrixXcd Gkw0; //= (iwn +mu - H - S_HF)^{-1}
            Gkw = -H_k_inModelSpace[k];
            Gkw0.setZero(NBAND[k], NBAND[k]);
            Gkw0 =  -1*(static_ev[k].asDiagonal());
            //            Gkw0 *= -1.;
            /*Gkw for given submat;*/
            for(int i0=0; i0<NBAND[k]; i0++) {
                for(int m0=0; m0<NBAND[k]; m0++) {
                    if(isOrbitalHartr[i0] and isOrbitalHartr[m0] and isSameAtom(i0,m0) ) {
                        Gkw(i0,m0)  -= SWLOCAL(KS2Hartr[i0],KS2Hartr[m0]) ;
                    }
                }//msub
                Gkw(i0,i0)  += (iw + mu) ;
                Gkw0(i0,i0) += (iw + mu )  ;
                Gkw0(i0,i0) = 1./Gkw0(i0,i0);
            }//isub
            Gkw = Gkw.inverse();
            Gkw0 =  staticHamiltonian_Model[k]     * Gkw0 * staticHamiltonian_Model[k].adjoint();
            /*Gw_loc in the correlated space*/
            for(int cl=0 ; cl < NumCluster; cl++) {
                if(n<N_freq) {
//                                        for(int i0=0; i0<N_peratom_HartrOrbit; i0+=1) {}
//                                            for( int m0=0; m0<N_peratom_HartrOrbit; m0+=1) {}
//                    for(int p=HartrRange[at][0] ; p <HartrRange[at][1] ; p++) {}
//                        for(int q=HartrRange[at][0] ; q <HartrRange[at][1] ; q++) {}
                    for(int i0=0; i0<NumHartrOrbit_per_cluster; i0+=1) {
                        for( int m0=0; m0<NumHartrOrbit_per_cluster; m0+=1) {
//                            int p0 =  KS2Hartr[p]-at*N_peratom_HartrOrbit;
//                            int q0 =  KS2Hartr[q]-at*N_peratom_HartrOrbit;
                            int p=HartrIndex[cl*NumHartrOrbit_per_cluster + i0];
                            int q=HartrIndex[cl*NumHartrOrbit_per_cluster + m0];
                            Gw_loc_mpiLoc(cl, n*COL+(i0*NumHartrOrbit_per_cluster+m0)) +=Gkw( p, q );
                        }
                    }
                }
            }//cl
            /*NumMat*/
            for(int l1=0; l1<NumCluster*NumHartrOrbit_per_cluster; l1++) {
                for(int l2=0; l2<NumCluster*NumHartrOrbit_per_cluster; l2++) {
                    int ll1=HartrIndex[l1];
                    int ll2=HartrIndex[l2];
                    NumMat_loc(l1,l2)  += ( Gkw(ll1,ll2) +std::conj(Gkw(ll2,ll1)) );   // Gkw(-w) = [Gkw(w)]^\dagger
                    NumMat_loc0(l1,l2) += ( Gkw0(ll1,ll2)+std::conj(Gkw0(ll2,ll1)) );
                }/*l2*/
            }/*l1*/
            for(int l1=0; l1<NBAND[k]; l1++) {
                for(int l2=0; l2<NBAND[k]; l2++) {
                    densityMatDFT[k](l1,l2) += (  Gkw(l1,l2)+std::conj( Gkw(l2,l1)) );
                    densityMatDFT[k](l1,l2) -= ( Gkw0(l1,l2)+std::conj(Gkw0(l2,l1)) );
                }/*l2*/
            }/*l1*/
        }///k
        for(int ATOM=0 ; ATOM < NumCluster; ATOM++) {
            if(n<N_freq) {
                for(int h1=0; h1<COL; h1++) {
                    Gw_loc_mpiLoc(ATOM,n*COL+(h1)) /= knum_mpiGlobal;
                }
            }
        }
        timeEnd= clock();
        if(mpi_rank==0 and n==0)  std::cout << "We will spend " << ((timeEnd - timeStart)* N_freq*highFreq2)/(CLOCKS_PER_SEC)  << " sec\n";
    }/*n*/

    /*upfolding_number operator*/
    for(int k=0; k<knum; k++) {
        densityMatDFT[k] /= beta;
        for(int i0=0; i0<NBAND[k]; i0++) {
            if        (real(beta*(static_ev[k](i0) - mu)) >  100.0 ) static_ev[k](i0) = 0.0;
            else if   (real(beta*(static_ev[k](i0) - mu)) < -100.0 ) static_ev[k](i0) = 1.0;
            else      static_ev[k](i0)  = 1./( 1.+std::exp(  beta*(static_ev[k](i0)-mu) ));
        }
        densityMatDFT[k] += (   (staticHamiltonian_Model[k]) * static_ev[k].asDiagonal() ) * (staticHamiltonian_Model[k].adjoint()) ;
    }


    /*NumMat high freq. tails.. Caution: Order is important!!.*/
    for(int l0=0; l0<NumCluster*NumHartrOrbit_per_cluster; l0++) {
        for(int m0=0; m0<NumCluster*NumHartrOrbit_per_cluster; m0++) {
            NumMat_loc(l0,m0) /=  beta ;
            NumMat_loc0(l0,m0) /=  beta ;
        }
        NumMat_loc(l0,l0)   +=  0.5*knum;
        NumMat_loc0(l0,l0)  +=  0.5*knum;
    }
    for(int k=0; k<knum; k++) {
        Eigen::MatrixXcd  occtemp = (staticHamiltonian_Model[k] * static_ev[k].asDiagonal() ) * staticHamiltonian_Model[k].adjoint();
        for(int l0=0; l0<NumCluster*NumHartrOrbit_per_cluster; l0++) {
            for(int m0=0; m0<NumCluster*NumHartrOrbit_per_cluster; m0++) {
                int ll0=HartrIndex[l0];
                int mm0=HartrIndex[m0];
                NumMat_loc(l0,m0) += occtemp(ll0,mm0);
            }
        }
    }//k
    delete [] staticHamiltonian_Model;
    delete [] static_ev;




    for(int l0=0; l0<NumCluster*NumHartrOrbit_per_cluster; l0++) {
        for(int m0=0; m0<NumCluster*NumHartrOrbit_per_cluster; m0++) {
            NumMat_loc(l0,m0) /=  knum_mpiGlobal ;
            NumMat_loc0(l0,m0) /=  knum_mpiGlobal ;
        }
    }

    NumMat_loc -= NumMat_loc0;
    /*NumMatrix, mpireduce for k space sum*/
    MPI_Allreduce(NumMat_loc.data(), NumMatrix.data(), NumMatrix.size(), MPI_DOUBLE_COMPLEX, MPI_SUM,  MPI_COMM_WORLD);

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> ces (NumHartrOrbit_per_cluster*NumCluster);
    ces.compute(NumMatrix);
    for(int i=0; i < NumHartrOrbit_per_cluster*NumCluster; i++)  {
        if (ces.eigenvalues()[i] < 0 or ces.eigenvalues()[i] >1+1e-10)  {
            ifroot std::cout << "[GreenFtn]: NumMatrix is unphysical, " << ces.eigenvalues()[i]  <<"\n";
            //            exit(1);
        }
    }

    /*Gw, mpireduce for k space sum*/
    MPI_Allreduce(Gw_loc_mpiLoc.data(), Gw_loc.data(), Gw_loc.size(), MPI_DOUBLE_COMPLEX, MPI_SUM,  MPI_COMM_WORLD);

    for(int ATOM=0 ; ATOM < NumCluster; ATOM++) {
        for(int n=0; n<N_freq; n++) {
            for (int h1=0; h1< NumHartrOrbit_per_cluster; h1++) {
                for (int h2=0; h2< NumHartrOrbit_per_cluster; h2++) {
                    int h1F=ATOM*NumHartrOrbit_per_cluster + h1;
                    int h2F=ATOM*NumHartrOrbit_per_cluster + h2;
                    int h12=h1*NumHartrOrbit_per_cluster + h2;
                    Gw[n](h1F,h2F)=  Gw_loc(ATOM,n*COL+h12);
                }
            }
        }
    }



    ifroot {
        for(int ATOM=0 ; ATOM < NumCluster; ATOM++) {
            std::stringstream ss;
            ss << ATOM;
            FILE *datap2= fopen( (std::string("Gw_loc.dat"     )+ ss.str()).c_str(),"w");
            FILE *datap3= fopen( (std::string("Gw_loc.full.dat")+ ss.str()).c_str(),"w");
            for(int n=0; n<N_freq; n++) {
                fprintf(datap2, "%0.5f", imag(I*pi*(2*n+1)/beta));
                for(int h1=0; h1<NumHartrOrbit_per_cluster; h1++) {
                    int h1F = ATOM * NumHartrOrbit_per_cluster +h1;
                    fprintf(datap2, "    %e    %e", real(Gw[n](h1F,h1F)),  imag(Gw[n](h1F,h1F))); //BUG!
                    for(int h2=0; h2<NumHartrOrbit_per_cluster; h2++) {
                        int h2F = ATOM * NumHartrOrbit_per_cluster +h2;
                        if( std::abs(Gw[n](h1F,h2F)) >1e-8)
                            fprintf(datap3, "%d %d %d %e    %e\n",n,h1,h2, real(Gw[n](h1F,h2F)),  imag(Gw[n](h1F,h2F)));
                    }//h2
                }//h1
                fprintf(datap2, "\n");
            }//n
            fclose(datap2);
            fclose(datap3);
            std::cout << "FILEOUT:Gw_loc.dat" << ATOM <<"\n";
            std::cout << "FILEOUT:Gw_loc.full.dat" <<ATOM <<"\n";
        }//ATOM
    }//ifroot
}/*GreenFtn_w*/


void retarded_GreenFtn2( Eigen::MatrixXcd &retGkw_full,    Eigen::MatrixXcd & retGkw,
                         std::vector<Eigen::MatrixXcd> &  H_k_inModelSpace, ImgFreqFtn &  SE,
                         double mu, int k, int n) {
    Eigen::MatrixXcd retG0kw;
    retGkw.setZero(NBAND[k], NBAND[k]);
    retG0kw.setZero(NBAND[k], NBAND[k]);
    retGkw = -H_k_inModelSpace[k];
    retG0kw = -H_k_inModelSpace[k];
    for (int i0=0; i0<NBAND[k]; i0++) {
        for (int m0=0; m0<NBAND[k]; m0++) {
            if(isOrbitalHartr[i0] and isOrbitalHartr[m0] and isSameAtom(i0,m0) ) {
                retGkw(i0,m0) -= SE.getValue( n, (KS2Hartr[i0]), (KS2Hartr[m0]) );
            }
        }//m
        retGkw(i0,i0) += (   E0+n*dE   + mu + I*infinitesimal  );
        retG0kw(i0,i0) += (   E0+n*dE   + mu + I*infinitesimal  );
    }//i
    retGkw = retGkw.inverse();
    retG0kw = retG0kw.inverse();

    //from model basis to KS band base index
    Eigen::MatrixXcd retGkw_band =  (DF_CorrBase[k].adjoint() * retGkw * DF_CorrBase[k]).eval();
    retG0kw = (DF_CorrBase[k].adjoint() * retG0kw * DF_CorrBase[k]).eval();

    //upfolding to full DFT basis
    retGkw_full.setZero(NumOrbit,NumOrbit);
    for( int b1=0; b1<NumOrbit; b1++) {
        retGkw_full(b1,b1) = 1.0/((E0+n*dE)  - KS_eigenEnergy[k][b1]  +mu+ I*infinitesimal);
    }
    for( int band1=0; band1<NBAND[k]; band1++) {
        for( int band2=0; band2<NBAND[k]; band2++) {
            retGkw_full(FromValToKS[k][band1], FromValToKS[k][band2]) +=   (retGkw_band(band1,band2) - retG0kw(band1,band2));
        }
    }

#if DEBERG
    for (i0=0; i0<NBAND[k]; i0++) {
        assert(imag(retGkw(i0,i0)) <= 0 );
    }
#endif
}/*retarded_GreenFtn2*/


//double traceNorm( cmplx **A, int n ) {
//    int i,j;
//    double rst=0;
//    for(i=0; i<n; i++) {
//        for(j=0; j<n; j++) {
//            rst += pow(std::norm(A[i][j]),2);
//        }
//    }
//    return sqrt(rst);
//}
//double traceNorm_diag( cmplx *A, int n,int dim ) {
//    int i,j;
//    double rst=0;
//    for(i=0; i<n; i++) {
//        rst += pow(std::norm(A[i*dim+i]),2);
//    }
//    return sqrt(rst);
//}
//double traceNorm_full( cmplx *A, int n, int dim ) {
//    int i,j;
//    double rst=0;
//    for(i=0; i<n; i++) {
//        for(j=0; j<n; j++) {
//            rst += pow(std::norm(A[i*dim+j]),2);
//        }
//    }
//    return sqrt(rst);
//}

// Green ftn off-diagoal component?
//void traceTest(cmplx **Gw) {
//    double rst = 10e+10, temp;
//    for(int n=0; n<N_freq; n++) {
//        temp = traceNorm_diag(Gw[n],N_peratom_HartrOrbit, N_peratom_HartrOrbit) / traceNorm_full(Gw[n],N_peratom_HartrOrbit, N_peratom_HartrOrbit );
//        if(rst > temp)  rst = temp;
//    }
//    ifroot std:: cout <<"Gloc diag/full-diag = " << rst <<"\n";
//    //if (rst < 0.5 and SOLVERtype==1 and N_peratom_HartrOrbit!=0 ) {
//    //    ifroot std::cout << "###########################################\n";
//    //    ifroot std::cout << "      Warning: Please use hybmat,\n  we have large off diagonal hyb.components\n"    ;
//    //    ifroot std::cout << "###########################################\n";
//    //}
//}

double GiwToNele (double muTB, ImgFreqFtn & SelfE_w, cmplx *** eigenVal, int preComputation, std::vector<Eigen::MatrixXcd> &H_k_inModelSpace  ) {
    cmplx w;

    if(preComputation==0) {

        Eigen::MatrixXcd * SelfEnergy_Eig = new Eigen::MatrixXcd [N_freq];
        std::vector<Eigen::MatrixXcd >  moments;
        for(int n=0; n < N_freq; n++) {
            SelfEnergy_Eig[n] = SelfE_w.getMatrix(n);
        }
        getAsymto_moments(moments, SelfEnergy_Eig);
        ifroot std::cout <<"Sw_HF    :\n" <<  moments[0] <<"\n";
        ifroot std::cout <<"1/(iwn)  :\n" <<  moments[1] <<"\n";
        ifroot std::cout <<"1/(iwn)^2:\n" <<  moments[2] <<"\n";
        ifroot std::cout <<"1/(iwn)^3:\n" <<  moments[3] <<"\n";




        time_t timeStart, timeEnd;
        timeStart=clock();
        for(int k=0; k<knum; k++) {
            for(int n=0; n<N_freq*highFreq2; n++) {
                std::vector<cmplx > model_qp_hamiltonian(NBAND[k] * NBAND[k] );
                for(int i1=0; i1<NBAND[k]; i1++) {
                    for(int m1=0; m1<NBAND[k]; m1++) {
                        model_qp_hamiltonian[i1* NBAND[k] +    m1] = 0;
                    }
                }
                if(NumCorrAtom!=0 ) {
                    valenceHAMILTONIAN(n,k)
                }
                for(int i=0; i<NBAND[k]; i++) {
                    for(int j=0; j<NBAND[k]; j++) {
                        model_qp_hamiltonian[i*NBAND[k]+j] +=  H_k_inModelSpace[k](i,j);
                    }
                }
                Eigen::MatrixXcd model_qp_H(NBAND[k],NBAND[k]);
                for(int m1=0; m1<NBAND[k]; m1++) {
                    for(int m2=0; m2<NBAND[k]; m2++) {
                        model_qp_H(m1,m2) = model_qp_hamiltonian[m1 * NBAND[k] + m2];
                    }
                }
                Eigen::ComplexEigenSolver<Eigen::MatrixXcd> ces(NumOrbit);
                ces.compute(model_qp_H,0);
                for(int i=0; i<NBAND[k]; i++) {
                    eigenVal[n][k][i] = ces.eigenvalues()[i];
                }
                //#endif
            }//n
        }//k
        ifroot    std::cout << "And we have freq. dep. eigen-energy.\n";
        /*highFreq*/
        for(int k=0; k<knum; k++) {
            int n=N_freq*highFreq2;
            std::vector<cmplx > model_qp_hamiltonian(  NBAND[k] * NBAND[k] );
            for(int i1=0; i1<NBAND[k]; i1++) {
                for(int m1=0; m1<NBAND[k]; m1++) {
                    model_qp_hamiltonian[i1* NBAND[k] +    m1] = 0;
                }
            }
            if(NumCorrAtom !=0 ) {
                valenceHAMILTONIAN(n,k);
            }
            for(int i=0; i<NBAND[k]; i++) {
                for(int j=0; j<NBAND[k]; j++) {
                    model_qp_hamiltonian[i*NBAND[k]+j] +=  H_k_inModelSpace[k](i,j);
                }
            }
            Eigen::MatrixXcd model_qp_H(NBAND[k],NBAND[k]);
            for(int m1=0; m1<NBAND[k]; m1++) {
                for(int m2=0; m2<NBAND[k]; m2++) {
                    model_qp_H(m1,m2) = model_qp_hamiltonian[m1 * NBAND[k] + m2];
                }
            }
            Eigen::ComplexEigenSolver<Eigen::MatrixXcd> ces(NumOrbit);
            ces.compute(model_qp_H,0);
            for(int i=0; i<NBAND[k]; i++) {
                eigenVal[n][k][i] = ces.eigenvalues()[i];
            }
        }//k
        timeEnd= clock();
        ifroot  std::cout << "H+S Diagonalization to get energy spectrum:" << (timeEnd - timeStart)/(CLOCKS_PER_SEC) << " sec\n";
        delete []  SelfEnergy_Eig ;
    }//precomp


    /*num. of electron in the valence space*/
    /*simple sum*/
    double numOfElecLoc=0;
    double numOfElecGlob=0;
    for(int n=0; n<N_freq*highFreq2; n++) {
        w=w0+n*dw;
        for(int k=0; k<knum; k++) {
            for(int l0=0; l0<NBAND[k]; l0++) {
                numOfElecLoc += real(  (1./(w-eigenVal[n][k][l0]+muTB))
                                       -(1./(w-eigenVal[N_freq*highFreq2][k][l0]+muTB))
                                    );
            }
        }
    }
    numOfElecLoc*=2./beta; //2=negative mat.freq.
    /*imaginary part tail*/
    for(int k=0; k<knum; k++) {
        for(int l0=0; l0<NBAND[k]; l0++) {
            numOfElecLoc += 1./(1+std::exp( beta * (real(eigenVal[N_freq*highFreq2][k][l0])-muTB)));  //Note ::  sum 1/iw exp(-iw\beta) =  - sum 1/iw exp(-iw0+)  = 0.5
        }
    }

    /*core level*/
    for(int k=0; k<knum; k++) {
        for(int l=0; l<NumOrbit; l++) {
            numOfElecLoc +=  1./(1+std::exp( beta * (KS_eigenEnergy[k][l]-muTB)));
        }
        for(int l0=0; l0<NBAND[k]; l0++) {
            int   l=FromValToKS[k][l0];
            numOfElecLoc -=  1./(1+std::exp( beta * (KS_eigenEnergy[k][l]-muTB)));
        }
    }
    numOfElecLoc/=knum_mpiGlobal;
    MPI_Allreduce(&(numOfElecLoc), &(numOfElecGlob), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return numOfElecGlob;
}/*GiwToNele */
