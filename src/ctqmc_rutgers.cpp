//ED solver JHSIM  171225 KAIST

#include <iostream>
#include "ImgTimeFtn.h"
#include <Eigen/Sparse>

#include "tight_common.h"
#include <time.h>
#include <Eigen/Eigenvalues>
//#include <vector>

/*
 state = 0~2^(n)-1 , 10 digit index of many-body basis set, where n= NSpinOrbit_per_atom + Num.of state for imp.

 alpha = single ptl state =  0,1,...n-1:
 single ptl state for each alpha, we may assign occupation.
 Then, we have many-body basis, state(|1100111000011.... 0>)   = 1*2^(n-1) + 0*2^(n-2) + ... +  0*2^(0) = sum_{alpha} occ_alpha * 2^(alpha)

 state = SectorToGlobalIndx[k][nCk], where, sectionIndx denote the congerved quantum number, for example number of ptls (k here), and nCk is the dimension of this subspace.

*/




int occupation_Bits(int i, unsigned long long alp) ;
int F_dagger_sign(int i, unsigned long long alp) ;
int countSetBits(unsigned  long long  n);
double getTotEng_seg( Eigen::MatrixXcd Local_Hamiltonian_ED, double muTB,  unsigned long long   alp);

unsigned long choose(int n, int k) ;   //nCk
void SectorToGlobalIndxConstruction (std::vector<std::vector<int> > & SectorToGlobalIndx, Eigen::VectorXi & GlobalIndx2Sector, Eigen::VectorXi & Global2SectorIndex,
                                     Eigen::VectorXi  & dim_for_sector, Eigen::MatrixXi & SuperState_quantumNum,  int NSector,  int Nstate) ;

int F_element(int alpha, unsigned long long state1, unsigned long long state2) ;
int F_element(int alpha, unsigned long long state1) ;


void HamiltonianConstruc (  Eigen::MatrixXcd & Himp, int Nstate, int Sector, std::vector<std::vector<int> > & SectorToGlobalIndx,
                            Eigen::MatrixXcd * f_annMat,Eigen::VectorXi  & dim_for_sector, Eigen::VectorXi & GlobalIndx2Sector, Eigen::VectorXi & Global2SectorIndex,
                            Eigen::MatrixXcd Local_Hamiltonian_ED, std::vector<cmplx > Utensor, std::vector<Eigen::VectorXi> Uindex,
                            double muTB, unsigned long long dim  ) ;
void write_PARMS() ;

void write_cix_file (int NSector, int maxDim, int Nstate, Eigen::VectorXi dim_for_sector,  std::vector<std::vector<double> > & ImpTotalEnergy, Eigen::MatrixXcd ** f_annMat,
                     std::vector<cmplx > Utensor, std::vector<Eigen::VectorXi> Uindex,
                     Eigen::MatrixXi  SuperState_quantumNum, Eigen::MatrixXi f_dagger_ann_SuperState   ) ;


void write_Delta(   ImgFreqFtn & weiss_field );
void write_Delta_diag(   ImgFreqFtn & weiss_field_diag) ;




void ctqmc_rutgers(  Eigen::MatrixXcd Local_Hamiltonian_ED, double muTB,
                     ImgFreqFtn & weiss_field, std::vector<cmplx > Utensor, std::vector<Eigen::VectorXi> Uindex, int solverDim  ) {
    ifroot "write rutgers_input file\n";
    int   Nstate = solverDim;


    int   NSector;
    if(magnetism ==0 or magnetism == 1) {
        int Norbital = Nstate/2;
        NSector = 2 * ((Norbital+1)*(Norbital+2))/2 - (Norbital+1)     ;              // superstate = (n\up, n\down) and nu+nd = N = 0,1,...Nstate    (0,0); (1,0), (0,1); (2,0),(1,1),(0,2); ...
    }
    else                {
        NSector = Nstate+1;                                   //Total number of the election in the system , 0,1,2,... Nstate
    }



    Eigen::MatrixXi   SuperState_quantumNum ;
    Eigen::VectorXi   GlobalIndx2Sector(   (int)std::pow(2,Nstate) );
    Eigen::VectorXi   Global2SectorIndex(  (int)std::pow(2,Nstate) );
    Eigen::VectorXi dim_for_sector(NSector);
    Eigen::MatrixXcd **   f_annMat = new Eigen::MatrixXcd * [NSector];                //fermion annihilation operator in matrix form  <N-1|f|N>
    std::vector<std::vector<int> > SectorToGlobalIndx(NSector);
    for (int N=0; N < NSector; N++) {
        f_annMat[N] = new Eigen::MatrixXcd [Nstate];
        dim_for_sector(N)  = 0;
    }
    SuperState_quantumNum.setZero(NSector, 2) ;
    SectorToGlobalIndxConstruction (SectorToGlobalIndx, GlobalIndx2Sector, Global2SectorIndex,
                                    dim_for_sector, SuperState_quantumNum, NSector, Nstate); //SectorToGlobalIndx[Sector][stateinSector] = state in 2^Nstate space. represented in the number basis.

    Eigen::MatrixXcd Himp;              //Hamiltonian block

    Eigen::MatrixXcd * evec_prev = new Eigen::MatrixXcd [NSector];

    std::vector<std::vector<double> > ImpTotalEnergy;

    //////////
    int maxDim=1;
    for (int Sector=0 ; Sector< NSector; Sector++) {
        std::vector<double> temp(dim_for_sector[Sector]);
        ImpTotalEnergy.push_back(temp);
        if(maxDim < dim_for_sector[Sector]) maxDim = dim_for_sector[Sector];

    }



    ////
    Eigen::MatrixXi f_ann_SuperState;
    Eigen::MatrixXi f_dagger_ann_SuperState;
    f_ann_SuperState.setZero(NSector, Nstate);
    f_dagger_ann_SuperState.setZero(NSector, Nstate);
    for (int Sector=0 ; Sector< NSector; Sector++) {
        for (int state=0 ; state < Nstate ; state++) {
            f_dagger_ann_SuperState(Sector,state) = -1;
            f_ann_SuperState(Sector,state) = -1;
        }
    }
    for (int Sector=0 ; Sector< NSector; Sector++) {
        for (int state=0 ; state < Nstate ; state++) {
            if( magnetism == 2 ) {
                f_ann_SuperState(Sector, state) =  Sector-1 ;
                if(Sector>0) f_dagger_ann_SuperState(  Sector-1   , state  ) = Sector;
            }
            else {
                for (int SectorBra=0 ; SectorBra< NSector; SectorBra++) {
                    if(state % 2 == 0 and
                            SuperState_quantumNum(SectorBra,0) == SuperState_quantumNum(Sector, 0) - 1   and
                            SuperState_quantumNum(SectorBra,1) == SuperState_quantumNum(Sector,1)  -1               ) {
                        f_ann_SuperState(Sector, state) = SectorBra;
                        f_dagger_ann_SuperState(  SectorBra   , state  ) = Sector;
                        break;
                    }
                    else if(state % 2 == 1 and
                            SuperState_quantumNum(SectorBra,0) == SuperState_quantumNum(Sector, 0) - 1   and
                            SuperState_quantumNum(SectorBra,1) == SuperState_quantumNum(Sector,1)  +1               ) {
                        f_ann_SuperState(Sector, state) = SectorBra;
                        f_dagger_ann_SuperState(  SectorBra   , state  ) = Sector;
                        break;
                    }//else if
                }
            }
        }
    }






    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    //SOLVER main
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
//for sector 0
    evec_prev[0].setZero(1,1);
    evec_prev[0](0,0)=1;
    ImpTotalEnergy[0][0] = 0.0;

    for (int Sector=1 ; Sector< NSector; Sector++) {
        ifroot       std::cout << "\nSector: " << Sector
                               << "\nNumber of electron in the total system: (" << SuperState_quantumNum(Sector,0) << "," << SuperState_quantumNum(Sector,1)
                               <<") in "<<  Nstate <<" states\n"
                               << "Hilbert space dimension:"<<dim_for_sector[Sector]<<"\n";


        for (int state=0 ; state < Nstate ; state++) {
            if(  f_ann_SuperState(Sector,state) >= 0 )       f_annMat[Sector][state].setZero(  dim_for_sector[f_ann_SuperState(Sector,state)],  dim_for_sector[Sector] );
            for(int  i=0; i< dim_for_sector[Sector]; i ++) {
                if(occupation_Bits(state, SectorToGlobalIndx[Sector][i]) == 1) {
                    int temp= SectorToGlobalIndx[Sector][i] - std::pow(2,state);
                    assert( GlobalIndx2Sector[temp] == f_ann_SuperState(Sector,state) );
                    f_annMat[Sector][state](   Global2SectorIndex[  temp ]    ,   i  )
                        = F_element(state, temp, SectorToGlobalIndx[Sector][i])  ;
                }
            }
        }


        time_t timeStartHam, timeEndHam;
        timeStartHam = clock();
        //Construct Hamiltonian
        ifroot std::cout << "HamiltonianConstruc\n";
        Himp.setZero( dim_for_sector[Sector], dim_for_sector[Sector] );
        HamiltonianConstruc ( Himp,   Nstate, Sector,SectorToGlobalIndx,
                              f_annMat[Sector], dim_for_sector,  GlobalIndx2Sector, Global2SectorIndex,
                              Local_Hamiltonian_ED, Utensor, Uindex,
                              muTB, dim_for_sector[Sector] );
        timeEndHam =clock();
        ifroot std::cout << "construct:" << (timeEndHam - timeStartHam)/(CLOCKS_PER_SEC) <<"\n";

        //diagonalization
        Eigen::MatrixXd Himp_r = Himp.real();
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ces(Himp_r);
        ces.compute(Himp_r);
        for(int alp=0; alp<Nstate; alp++) {
            if(  f_ann_SuperState(Sector,alp) >= 0 )
                f_annMat[Sector][alp] =  evec_prev[f_ann_SuperState(Sector,alp)].adjoint()  *  f_annMat[Sector][alp] * ces.eigenvectors();  //f_annMat in the many-body eigen state basis
        }
        for (int i=0; i<dim_for_sector[Sector] ; i++) {
            (ImpTotalEnergy.at(Sector)).at(i) = ces.eigenvalues()[i];
        }

        evec_prev[Sector] = ces.eigenvectors();
    }//section

    ifroot std::cout << "\nwrite_cix_file\n";
    ifroot     write_cix_file (NSector,  maxDim, Nstate, dim_for_sector, ImpTotalEnergy, f_annMat,
                               Utensor, Uindex,
                               SuperState_quantumNum, f_dagger_ann_SuperState);
    write_Delta(    weiss_field);
    ifroot write_PARMS();



    for (int N=0; N < NSector; N++) {
        delete []        f_annMat[N];
    }
    delete [] f_annMat;
} //ct_qmc_rut

void ctqmc_rutgers_seg(  Eigen::MatrixXcd Local_Hamiltonian_ED, double muTB,
                         ImgFreqFtn & weiss_field,  std::vector<cmplx > Utensor, std::vector<Eigen::VectorXi> Uindex) {
    ifroot{

        FILE * cix = fopen("rutgers_input.cix", "w");
        fprintf(cix, "# Cix file for cluster DMFT with CTQMC\n");
        fprintf(cix, "# cluster_size, number of states, number of baths, maximum matrix size\n");
        fprintf(cix, "1 %d %d 1 \n",((int)std::pow(2,NSpinOrbit_per_atom)), NSpinOrbit_per_atom);
        fprintf(cix, "# baths, dimension, symmetry, global flip\n");
        for(int alp=0; alp<NSpinOrbit_per_atom; alp++) { //0,2,4...NSpinOrbit_per_atom-2
            fprintf(cix, "%d 1 %d %d\n",alp, alp, alp/2);
        }
        fprintf(cix, "# cluster energies for unique baths, eps[k]\n");
        for(int i=0; i<NSpinOrbit_per_atom; i++) {
            fprintf(cix, " 0");
        }
        fprintf(cix, "\n#   N   K   Sz size F^{+,0}, F^{+,1}, F^{+,2}, F^{+,3}, F^{+,4}, F^{+,5}, {Ea,Eb..} ;  {Sa, Sb..}\n");
        for(int alp=0; alp<std::pow(2,NSpinOrbit_per_atom); alp++) {
            fprintf(cix,      "%d  %d   0  0  1  ",
            alp+1, countSetBits(alp));
            for(int i=0; i<NSpinOrbit_per_atom; i++) {   //alp = \sum_i  n_i * 2^i; i=0,1,... NSpinOrbit_per_atom-1 <=> |n_{NSpinOrbit_per_atom-1}, .., n_1, n_0>
                int operationState = std::pow(2,i);   // (0,0,...,0,1,0,...0)
                fprintf(cix,      "%d ", (1-occupation_Bits(i,alp))*(alp+1+operationState));
            }

            fprintf(cix,      "%+0.8f", getTotEng_seg(Local_Hamiltonian_ED, muTB,alp));
            fprintf(cix,      " 0 \n");
        }//alp



        fprintf(cix, "# matrix elements\n");
        for(int alp=0; alp<std::pow(2,NSpinOrbit_per_atom); alp++) {
            for(int i=0; i<NSpinOrbit_per_atom; i++) {
                int operationState = std::pow(2,i);   // (0,0,...,0,1,0,...0)
                fprintf(cix,      "%d  %d   1  %d  ",
                        alp+1, (1-occupation_Bits(i,alp))*(alp+1+operationState), (1-occupation_Bits(i,alp))  );
                if(  (occupation_Bits(i,alp)) ==0 )  fprintf(cix,      " %+0.8f",             ((double)F_dagger_sign(i,alp))   );
                fprintf(cix,      "\n");
            }
        }
        fprintf(cix, "HB2                # Hubbard-I is used to determine high-frequency\n");

        fprintf(cix, "# UCoulomb : (m1,s1) (m2,s2) (m3,s2) (m4,s1)  Uc[m1,m2,m3,m4]\n");
        int length = (int) Utensor.size();
        int nH=0;
        while(nH < length) {
//            if( isOrbitalCorrinHart[Uindex[nH](0) ] and
//                    isOrbitalCorrinHart[Uindex[nH](1) ] and
//                    isOrbitalCorrinHart[Uindex[nH](2) ] and
//                    isOrbitalCorrinHart[Uindex[nH](3) ] ) {}
            fprintf(cix, " %d %d %d %d %+0.8f\n",
                    Uindex[nH](0),
                    Uindex[nH](1),
                    Uindex[nH](3),
                    Uindex[nH](2),  real(Utensor[nH]) );
            assert(imag(Utensor[nH]) < 1e-5);
//            {}
            nH++;
        }


        fprintf(cix, "# number of operators needed\n0");

        fclose(cix);


        ifroot std::cout << "\nwrite_cix_file\n";
        write_Delta_diag(    weiss_field );
        ifroot write_PARMS();
    }//ifroot
}//ctqmc_rutgers_seg

void write_PARMS() {
    FILE *fp;
    if ((fp = fopen("PARAMS","r")) == NULL) {
        FILE * PARMS_file = fopen("PARAMS", "w");
        fprintf(PARMS_file, "nom %d  # number of Matsubara frequencies\n", N_freq                 );
        fprintf(PARMS_file, "svd_lmax 30 # number of SVD functions to project the solution\n"                 );
        fprintf(PARMS_file, "svd_L 15 # To compute the SVD decomposition of the kernel for analytic continuation, we need to choose the cutoff on the real axis.(default: 10).\n"                 );
        fprintf(PARMS_file, "tsample  50     # how often to record the measurements\n"                        );
        fprintf(PARMS_file, "aom      1      # number of frequency points to determin high frequency tail\n"  );
        fprintf(PARMS_file, "M %llu     # Number of Monte Carlo steps\n", Num_MC_steps                                       );
        fprintf(PARMS_file, "beta %0.5f        # Inverse temperature\n", beta                                 );
        fprintf(PARMS_file, "U 0        # Coulomb repulsion (F0), This information is wrtten in cix file\n"   );
        fprintf(PARMS_file, "GlobalFlip 500000  # how often to perform global flip\n"                         );
        fprintf(PARMS_file, "exe %s  # Path to executable\n", SOLVERexe.c_str()                 );
        fprintf(PARMS_file, "mu  0   # Chemical potential\n"                                                  );
        fprintf(PARMS_file, "mode SM   # S stands for self-energy sampling, M stands for high frequency moment tail\n");
        fprintf(PARMS_file, "Delta delta_w_rutgers.dat  # Input bath function hybridization\n");
        fprintf(PARMS_file,  "cix rutgers_input.cix # Input file with atomic state\n"                          );
        fclose(PARMS_file);
        sleep(3);
        std::cout <<"Write ct-qmc input, PARMS\n";
    }
}


void write_Delta(   ImgFreqFtn & weiss_field) {

    int spindim=1;
    if(magnetism==0 or magnetism==1) spindim=2;

    if(mpi_rank==0) {
        std::string filename = std::string("delta_w_rutgers.dat");
        FILE *datap4 = fopen(filename.c_str(), "w");
        for(int w=0; w< weiss_field.getNFreq(); w++) {
            fprintf(datap4, "%0.8f", weiss_field.getValue(w) );
            for(int spin = 0 ; spin< spindim; spin++) {
                for(int n=spin; n< NSpinOrbit_per_atom; n+=spindim) {
                    for(int m=n; m< NSpinOrbit_per_atom; m+=spindim) {
                        fprintf(datap4, "     %0.10f  %0.10f", real( weiss_field.getValue(w, n,m)  ),
                                imag( weiss_field.getValue(w, n,m)));
                    }
                }
            }
            fprintf(datap4,"\n");
        }
        fclose(datap4);
        std::cout << "FILEOUT:" <<  filename <<"\n";
    }
}

void write_Delta_diag(   ImgFreqFtn & weiss_field) {
    if(mpi_rank==0) {
        std::string filename = std::string("delta_w_rutgers.dat");
        FILE *datap4 = fopen(filename.c_str(), "w");
        for(int w=0; w< weiss_field.getNFreq(); w++) {
            fprintf(datap4, "%0.8f", weiss_field.getValue(w) );
            for(int n=0; n< NSpinOrbit_per_atom; n++) {
                fprintf(datap4, "     %0.10f  %0.10f", real( weiss_field.getValue(w, n, n)  ),
                        imag( weiss_field.getValue(w,  n, n)));
            }
            fprintf(datap4,"\n");
        }
        fclose(datap4);
        std::cout << "FILEOUT:" <<  filename <<"\n";
    }
}


void write_cix_file (int NSector, int maxDim, int Nstate, Eigen::VectorXi dim_for_sector,  std::vector<std::vector<double> > & ImpTotalEnergy, Eigen::MatrixXcd ** f_annMat,
                     std::vector<cmplx > Utensor, std::vector<Eigen::VectorXi> Uindex,
                     Eigen::MatrixXi  SuperState_quantumNum, Eigen::MatrixXi f_dagger_ann_SuperState   ) {


    int spindim=1;
    if(magnetism==0 or magnetism==1) spindim=2;

    FILE * cix = fopen("rutgers_input.cix", "w");
    fprintf(cix, "# Cix file for cluster DMFT with CTQMC\n");
    fprintf(cix, "# cluster_size, number of states, number of baths, maximum matrix size\n");
    fprintf(cix, "1 %d %d %d\n",NSector, spindim, maxDim);   //custersize, numberof SuperState, number_of_bath, maxDim;
    fprintf(cix, "# baths, dimension, symmetry, global flip\n");
    int upper_triangle=0;
    int deltaForRutgersCTQMS[Nstate][Nstate];
    for(int spin = 0 ; spin< spindim; spin++) {
        for(int n=spin; n< NSpinOrbit_per_atom; n+=spindim) {
            for(int m=n; m< NSpinOrbit_per_atom; m+=spindim) {
                deltaForRutgersCTQMS[n][m] = upper_triangle;
                upper_triangle++;
            }
        }
    }

//    for(int alp=0; alp<Nstate; alp++) {
//        for(int bet=0; bet<Nstate; bet++) {
//            if (alp==bet) {
//                deltaForRutgersCTQMS[alp][bet] = upper_triangle;
//                upper_triangle++;
//            }
//            else if (alp <bet) {
//                deltaForRutgersCTQMS[alp][bet] = upper_triangle;
//                upper_triangle++;
//            }
//
//        }
//    }



    for(int spin = 0 ; spin< spindim; spin++) {
        fprintf(cix, "%d  %d\n", spin, Nstate/spindim );
        for(int alp=spin; alp<Nstate; alp+=spindim) {
            for(int bet=spin; bet<Nstate; bet+=spindim) {
                if (alp<=bet) {
                    fprintf(cix, " %d", deltaForRutgersCTQMS[alp][bet]);
                }
                else if (alp > bet) {
                    fprintf(cix, " %d*", deltaForRutgersCTQMS[bet][alp]);
                }
            }
            fprintf(cix, "\n");
        }
        fprintf(cix, "0\n");
    }



    fprintf(cix, "# cluster energies for unique baths, eps[k]\n");
    for(int alp=0; alp<Nstate; alp++) {
        fprintf(cix, " 0");
    }
    fprintf(cix, "\n");



    fprintf(cix, "#   N   K   Sz size F^{+,0}, F^{+,1}, F^{+,2}, F^{+,3}, F^{+,4}, F^{+,5}, {Ea,Eb..} ;  {Sa, Sb..}\n");
    for (int Sector=0; Sector < NSector; Sector++) {
        int sectorIndex = Sector+1;
        fprintf(cix, "%d %d 0 %d %d   ", sectorIndex, SuperState_quantumNum(Sector,0), SuperState_quantumNum(Sector,1 ),  (int)dim_for_sector[Sector]);

        for(int spin = 0 ; spin< spindim; spin++) {
            for (int alp=spin; alp < Nstate; alp+=spindim) {
                int nextSector = f_dagger_ann_SuperState(Sector, alp) + 1;
                fprintf(cix, " %d",  nextSector );
            }
        }
        for (int N=0; N < dim_for_sector[Sector]; N++) {
            fprintf(cix, " %e",ImpTotalEnergy[Sector][N]);

        }
        for (int N=0; N < dim_for_sector[Sector]; N++) {
            fprintf(cix, " 0");
        }
        fprintf(cix, "\n");
    }

    fprintf(cix, "# matrix elements\n");
//    for (int ketSectorIndex=1; ketSectorIndex <= NSector; ketSectorIndex++) {
//        int braSectorIndex = ketSectorIndex+1;

    for (int ketSector=0; ketSector < NSector; ketSector++) {

        int ketSectorIndex = ketSector + 1;
        int ketDim = dim_for_sector[ketSector];



        for(int spin = 0 ; spin< spindim; spin++)
            for(int alp=spin; alp<Nstate; alp+=spindim) {

                int braSector = f_dagger_ann_SuperState(ketSector,alp);
                int braSectorIndex = braSector + 1;
                int braDim = ( (braSector==-1) ? 0:  dim_for_sector[braSector] );
//            braSectorIndex = ( (braSectorIndex>NSector) ? 0:  braSectorIndex );


                fprintf(cix, "%d %d %d %d",ketSectorIndex, braSectorIndex, ketDim, braDim);
                Eigen::MatrixXcd  f_dagger;
                if(braSector != -1 )   f_dagger= f_annMat[braSector][alp].adjoint(); // <N+1|f^+ |N>  = <N| f | N+1>.adjoint()
                for (int ketN=0; ketN < ketDim; ketN++) {
                    for (int braN=0; braN < braDim; braN++) {
                        fprintf(cix, " %e",real(f_dagger(braN,ketN)) );
                        if( std::abs(imag(f_dagger(braN,ketN))) >  1e-5) {
                            std::cout.precision(10) ;
                            std::cout << "braDim:"<< braDim << "ketDim:"<< ketDim << "Im(f_dagger):" <<  imag(f_dagger(braN,ketN)) <<"\n" ;
                            std::cout.precision(4) ;
                            exit(1);
                        }
                    }
                    fprintf(cix, "   ");
                }
                fprintf(cix, "\n");
            }//alp
    }//ketSector
    fprintf(cix, "HB2                # Hubbard-I is used to determine high-frequency\n");

    fprintf(cix, "# UCoulomb : (m1,s1) (m2,s2) (m3,s2) (m4,s1)  Uc[m1,m2,m3,m4]\n");
    int length = (int) Utensor.size();
    int nH=0;
    while(nH < length) {
        int i = Uindex[nH](0)/2 +  (Nstate/2) * (Uindex[nH](0)%2);
        int j = Uindex[nH](1)/2 +  (Nstate/2) * (Uindex[nH](1)%2);
        int k = Uindex[nH](2)/2 +  (Nstate/2) * (Uindex[nH](2)%2);
        int l = Uindex[nH](3)/2 +  (Nstate/2) * (Uindex[nH](3)%2);
        fprintf(cix, " %d %d %d %d %+0.8f\n",
                i,
                j,
                l,
                k,  real(Utensor[nH]) );

        if ( imag(Utensor[nH]) > 1e-5) {
            std::cout << "imag part of Utensor" << Uindex[nH](0) << Uindex[nH](1) <<  Uindex[nH](2) <<  Uindex[nH](3) <<Utensor[nH] <<"\n";
            exit(1);
        }
        nH++;
    }
////////////////////////////////////////////////////////
//    nH=0;
//    while(nH < length) {
//        fprintf(cix, " %d %d %d %d %+0.8f\n",
//                Uindex[nH](0),
//                Uindex[nH](1),
//                Uindex[nH](3),
//                Uindex[nH](2),  real(Utensor[nH]) );
//        if ( imag(Utensor[nH]) > 1e-5) {
//            std::cout << "imag part of Utensor" << Uindex[nH](0) << Uindex[nH](1) <<  Uindex[nH](2) <<  Uindex[nH](3) <<Utensor[nH] <<"\n";
//            exit(1);
//        }
//        nH++;
//    }
//////////////////////////////////////////////////////////


    fprintf(cix, "# number of operators needed\n0");

    fclose(cix);

}//write_cix_file








void HamiltonianConstruc (  Eigen::MatrixXcd & Himp, int Nstate, int Sector, std::vector<std::vector<int> > & SectorToGlobalIndx,
                            Eigen::MatrixXcd * f_annMat,Eigen::VectorXi  & dim_for_sector, Eigen::VectorXi & GlobalIndx2Sector, Eigen::VectorXi & Global2SectorIndex,
                            Eigen::MatrixXcd Local_Hamiltonian_ED, std::vector<cmplx > Utensor, std::vector<Eigen::VectorXi> Uindex,
                            double muTB, unsigned long long dim  ) {


    Himp.setZero(dim,dim);
    for(int alp=0; alp<Nstate; alp++) {
        for(int bet=0; bet<Nstate; bet++) {
            for(int  j=0; j< dim_for_sector[Sector]; j++) {
                int global_i  =       SectorToGlobalIndx[Sector][j] - std::pow(2,bet) + std::pow(2,alp);
                if(     occupation_Bits(bet, SectorToGlobalIndx[Sector][j]) == 1  and
                        (occupation_Bits(alp, SectorToGlobalIndx[Sector][j]) == 0 or alp==bet)  and
                        GlobalIndx2Sector[global_i] == Sector    ) {
                    int i  =   Global2SectorIndex[ global_i]      ;
                    Himp(i,j) +=  Local_Hamiltonian_ED(alp,bet) * F_element(alp, global_i) * F_element( bet, SectorToGlobalIndx[Sector][j] );

                    if(alp==bet) Himp(i,j) -= muTB * F_element(alp, global_i) * F_element( bet, SectorToGlobalIndx[Sector][j] );
                }
            }
        }//bet
    }//alp

    if ( (Himp.imag()).norm() > 1e-5 and mpi_rank == 0 ) {
        std::cout << "imag part of H0:" <<   (Himp.imag()).norm() <<"\n";
        exit(1);
    }


    for(int ll=0; ll<Utensor.size(); ll++) {
        int s1 = Uindex[ll](0);
        int s2 = Uindex[ll](1);
        int s3 = Uindex[ll](2);
        int s4 = Uindex[ll](3);
//        Himp += 0.5* Utensor[ll] * f_annMat[s1].adjoint() * f_annMat[s2].adjoint() * f_annMat[s4] *f_annMat[s3];
//  <==>  Himp += 0.5* Utensor[ll]  (f_annMat[s1].adjoint()*f_annMat[s3] - f_annMat[s1].adjoint()*f_annMat[s4]  * f_annMat[s2].adjoint()  *f_annMat[s3]);

//        Himp -= 0.5* Utensor[ll] * f_annMat[s1].adjoint()  * f_annMat[s4]  * f_annMat[s2].adjoint()  *f_annMat[s3];



        for(int  j=0; j< dim_for_sector[Sector]; j++) {
            if(     s1 != s2 and s3 != s4 and
                    occupation_Bits(s3, SectorToGlobalIndx[Sector][j]) == 1  and
                    occupation_Bits(s4, SectorToGlobalIndx[Sector][j]) == 1  and
                    (occupation_Bits(s1, SectorToGlobalIndx[Sector][j]) == 0 or s1==s3 or s1==s4) and
                    (occupation_Bits(s2, SectorToGlobalIndx[Sector][j]) == 0 or s2==s3 or s2==s4) and
                    GlobalIndx2Sector[  SectorToGlobalIndx[Sector][j] - std::pow(2,s3) - std::pow(2,s4)    + std::pow(2,s2) + std::pow(2,s1)  ] == Sector
              ) {

                int global_i  =   SectorToGlobalIndx[Sector][j] - std::pow(2,s3) - std::pow(2,s4)    + std::pow(2,s2) + std::pow(2,s1)      ;
                int i  =   Global2SectorIndex[ global_i]      ;
                int global_k  =   SectorToGlobalIndx[Sector][j] - std::pow(2,s3)  ;
                int global_l  =   global_k - std::pow(2,s4) + std::pow(2,s2)   ;

                Himp(i,j) += 0.5*Utensor[ll] *
                             F_element(s1, global_i) *  F_element(s2, global_l) * F_element(s4, global_k) * F_element(s3, SectorToGlobalIndx[Sector][j]);

            }
        }
    }

    if ( (Himp.imag()).norm() > 1e-5  and mpi_rank == 0  ) {
        std::cout << "imag part of Himp:" <<   (Himp.imag()).norm() <<"\n";
        exit(1);
    } //// The matrix elements of the f_annMat in the many-body eigen state basis should be real
}



//int F_element(int alpha, unsigned long long state1, unsigned long long state2) {
//    if ( state1 == state2-pow(2,alpha) ) {
//        int sign,i=0,N=0;
//        for(i=0; i<alpha; i++) {
//            N+=occupation_Bits(i,state2);
//        }
//        sign = (int)std::pow(-1,N);
//        return sign ;
//    }
//    else return 0;
//}

void SectorToGlobalIndxConstruction (std::vector<std::vector<int> > & SectorToGlobalIndx, Eigen::VectorXi & GlobalIndx2Sector, Eigen::VectorXi & Global2SectorIndex,
                                     Eigen::VectorXi  & dim_for_sector, Eigen::MatrixXi & SuperState_quantumNum,  int NSector,  int Nstate) {
    unsigned long long   S[NSector];
    for (int N=0; N< NSector ; N++) {
        S[N]=0;  //state index in the sector
    }
    for ( int globalindex =0; globalindex < pow(2,Nstate); globalindex++) {//sjh
        int Ntot = 0;
        int Ndn  = 0;

        for (int i =0; i<Nstate; i++) Ntot += occupation_Bits ( i, globalindex);
        for (int i =1; i<Nstate; i+=2)
            Ndn += occupation_Bits ( i, globalindex);

        int N;  //N= sector index;
        if( magnetism ==2) {
            N= Ntot;
            SuperState_quantumNum(N,0) = Ntot;
            SuperState_quantumNum(N,1) = 0;
        }
        else {
//int Norbital = Nstate/2;
            int Ntot_hole = Nstate - Ntot;
            int Ndn_hole  = Nstate/2 - Ndn;
            if(Ntot <=  Nstate/2)        N = (Ntot*(Ntot+1))/2 + Ndn ; //(0,0); (1,0), (0,1); (2,0),(1,1),(0,2); ...
            else N =  (NSector - 1)   -  ( (Ntot_hole*(Ntot_hole+1))/2 + Ndn_hole   );
//            if(N>=NSector) std::cout << NSector<<" " << N <<" " << Ntot<<" " << Ndn <<" " << Nhole<<"\n";

            SuperState_quantumNum(N,0) = Ntot;
            SuperState_quantumNum(N,1) = Ntot- 2*Ndn;
        }
        assert(N< NSector);
        SectorToGlobalIndx[N].push_back(globalindex);
//        SectorToGlobalIndx(N,S[N]) = globalindex;
        GlobalIndx2Sector(globalindex) = N;
        Global2SectorIndex(globalindex) = S[N];
        dim_for_sector(N) += 1;

        S[N]+=1;
    }
}//SectorToGlobalIndxConstruction



unsigned long choose(int n, int k) {
    //nCk
    if (k > n) {
        return 0;
    }
    unsigned long long r = 1;
    int dmax;
    if (k>=n/2) {
        dmax = n-k;
    }
    else if (k<n/2) {
        dmax=k;
    }
    for (int d = 1; d <= dmax; ++d) {
        r *= n--;
        r /= d;
        //    std::cout <<n <<" "<<k<<" "<<d<<" "<<r<<"\n";
    }
    return r;
}


//int occupation(int alpha, unsigned long long state) {
//    return ( (state/((int)pow(2,alpha))) %2 );
//}


/* Function to get no of set bits in binary
 *    representation of positive integer n */
int countSetBits(unsigned  long long  n)
{
    int count = 0;
    while (n)
    {
        count += n & 1;
        n >>= 1;
    }
    return count;
}
int occupation_Bits(int i, unsigned long long alp) {
    return ( (alp >> i) & 1);   //
}
//int occupation_tot (unsigned long long alp, int Nstate) {
//    int Ntot =0;
//    for (int i =0; i< Nstate ; i++) {
//        Ntot+= occupation_Bits ( i, globalindex);
//    }
//    return Ntot;
//}
//int occupation_spin (int spin, unsigned long long alp, int Nstate ) {
//    int Ndn =0;
//    for (int i = spin; i<Nstate; i+=2) {
//        Ndn += occupation_Bits ( i, globalindex);
//    }
//    return Ndn;
//}
//
//
int F_dagger_sign(int i, unsigned long long alp) {
    int count = 0;
    for(int j=0; j<i; j++)
    {
        count += ((alp >> j) & 1);
    }
    return  std::pow(-1, count);
}

int F_element(int alpha, unsigned long long state1, unsigned long long state2) {
    if ( state1 == state2-pow(2,alpha) ) {
        return F_dagger_sign  ( alpha, state1);
    }
    else return 0;
}
int F_element(int alpha, unsigned long long state1) {
    return F_dagger_sign  ( alpha, state1);
}



double getTotEng_seg( Eigen::MatrixXcd Local_Hamiltonian_ED, double muTB,  unsigned long long   alp) {

    double TotEng=0;
    for(int i=0; i<NSpinOrbit_per_atom; i++) {
        TotEng +=  ( real(Local_Hamiltonian_ED(i,i))-muTB) *  occupation_Bits(i,alp);
    }//alp
    for(int i=0; i<NSpinOrbit_per_atom; i+=2) { //0,2,4...NSpinOrbit_per_atom-2
        TotEng += UHubb   * occupation_Bits(i,alp) *occupation_Bits(i+1,alp) ;
        for(int j=0; j<i; j+=2) {
            TotEng +=  Uprime         *  occupation_Bits(i,alp)  * occupation_Bits(j+1,alp);   //H_interaction
            TotEng +=  Uprime         *  occupation_Bits(i+1,alp)  * occupation_Bits(j,alp);     //H_interaction
            TotEng += (Uprime-JHund)  *  occupation_Bits(i,alp)  * occupation_Bits(j,alp);     //H_interaction
            TotEng += (Uprime-JHund)  *  occupation_Bits(i+1,alp)  * occupation_Bits(j+1,alp); //H_interaction
        }
    }//alp
    return TotEng;
}
