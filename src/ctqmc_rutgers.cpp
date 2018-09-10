//ED solver JHSIM  171225 KAIST

#include <iostream>
#include "ImgTimeFtn.h"
//#include "mpi.h"
//#include <complex>
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



int occupation(int alpha, unsigned long long state) ;

int occupation_Bits(int i, unsigned long long alp) ;
int F_dagger_sign(int i, unsigned long long alp) ;
int countSetBits(unsigned  long long  n);
double getTotEng( Eigen::MatrixXcd Local_Hamiltonian_ED, double muTB,  unsigned long long   alp);

unsigned long choose(int n, int k) ;   //nCk
void SectorToGlobalIndxConstruction(unsigned long long ** SectorToGlobalIndx , int Nstate) ;
int occupation(int alpha, unsigned long long state) ;
int F_element(int alpha, unsigned long long state1, unsigned long long state2) ;
void HamiltonianConstruc(  Eigen::MatrixXcd & Himp,  Eigen::MatrixXi **f_ann, int Nstate , unsigned long long ** SectorToGlobalIndx, int Sector, Eigen::MatrixXcd * f_annMat, Eigen::MatrixXcd Local_Hamiltonian_ED, double muTB, unsigned long long dim, unsigned long long dim_prev   ) ;


void write_PARMS() ;
void write_cix_file(int NSector, int maxDim, int Nstate, unsigned long long * dim_for_sector,  std::vector<std::vector<double> > & ImpTotalEnergy, Eigen::MatrixXcd ** f_annMat) ;
void write_Delta(   ImgFreqFtn & weiss_field );
void write_Delta_diag(   ImgFreqFtn & weiss_field_diag) ;




void ctqmc_rutgers(  Eigen::MatrixXcd Local_Hamiltonian_ED, double muTB,  ImgFreqFtn & weiss_field) {
    int   Nstate = NSpinOrbit_per_atom;
//    ifroot std::cout << "#############################\n"
//                     << "Start Hub-I solver:#of state="<< Nstate << "   U:"<<UHubb<<" U':"<<Uprime <<" J:" <<JHund
//                     << "\nNstate (imp+bath) should be smaller than 36  (int) or  70 (unsigned long long)\n";
//    if(!(UHubb == 2*JHund+Uprime) and mpi_rank==0) std::cout << "Warning: U!=2*J+U'\n";

    int   NSector;
    if(magnetism ==0) {
        NSector = (Nstate+1)*(Nstate+2)/2     ;              // superstate = (n\up, n\down) and nu+nd = N = 0,1,...Nstate    (0,0); (1,0), (0,1); (2,0),(1,1),(0,2); ...
    }
    else                {
        NSector = Nstate+1;                                   //Total number of the election in the system , 0,1,2,... Nstate
    }


    NSector = Nstate+1;                                   //Total number of the election in the system , 0,1,2,... Nstate



    unsigned long long  dim_for_sector[NSector];
    unsigned long long ** SectorToGlobalIndx;

    Eigen::MatrixXcd ** f_annMat;//fermion annihilation operator in matrix form  <N-1|f|N>
    SectorToGlobalIndx = new unsigned long long  * [NSector];
    f_annMat = new Eigen::MatrixXcd * [NSector];
    for (int N=0; N < NSector; N++) {
        SectorToGlobalIndx[N] = new unsigned long long [choose(Nstate,N)];
        f_annMat[N] = new Eigen::MatrixXcd [Nstate];
    }
    SectorToGlobalIndxConstruction(SectorToGlobalIndx, Nstate); //SectorToGlobalIndx[Sector][stateinSector] = state in 2^Nstate space. represented in the number basis.

    Eigen::MatrixXi ** f_ann;           //fermion annihilation operator information
    Eigen::MatrixXcd Himp;              //Hamiltonian block

    Eigen::MatrixXcd evec_prev(1,1);

    dim_for_sector[0]=1;
    evec_prev(0,0)=1;
    f_ann =  new Eigen::MatrixXi * [Nstate+1];
    std::vector<std::vector<double> > ImpTotalEnergy;

    int maxDim=1;
    for (int Sector=1 ; Sector<=Nstate; Sector++) {
        f_ann[Sector] =  new Eigen::MatrixXi [Nstate];
        dim_for_sector[Sector] = choose(Nstate,Sector);

        if(maxDim < dim_for_sector[Sector]) maxDim = dim_for_sector[Sector];
    }
    for (int Sector=0 ; Sector<=Nstate; Sector++) {
        std::vector<double> temp(dim_for_sector[Sector]);
        ImpTotalEnergy.push_back(temp);
    }
    ImpTotalEnergy[0][0] = 0.0;




    unsigned long long rankf_ann  ;
    for (int Sector=1 ; Sector<=Nstate; Sector++) {
        rankf_ann  =  choose(Nstate-1,Sector-1); //# of number state with occupied alp
        for (int alp=0; alp< Nstate; alp++) {
            f_ann[Sector][alp].setZero(rankf_ann,3) ;
        }//alp
        /*To get annihilation matrix*/
        ifroot std::cout << "Here we start to get annihilation matrix\n";
        for (int alp=0; alp< Nstate; alp++) {
            unsigned long  k=0;
            unsigned long  j=0;
            for (int i=0; i< dim_for_sector[Sector-1]; i++) {
                if(occupation(alp, SectorToGlobalIndx[Sector-1][i]) == 0 )
                    for ( ; j< dim_for_sector[Sector]; j++) {
                        if ( F_element(alp,SectorToGlobalIndx[Sector-1][i],SectorToGlobalIndx[Sector][j]) != 0 ) {
                            f_ann[Sector][alp](k,0) = j;
                            f_ann[Sector][alp](k,1) = i;
                            f_ann[Sector][alp](k,2) = F_element(alp,SectorToGlobalIndx[Sector-1][i],SectorToGlobalIndx[Sector][j]);
                            k++;
                            break;
                        }
                    }
            }//i
            assert(k==rankf_ann);
        }//alp
    }


/////////////////////////////////////////////////////////////////////////////////////////////////////////
    /*SOLVER main*/
/////////////////////////////////////////////////////////////////////////////////////////////////////////
    for (int Sector=1 ; Sector<=Nstate; Sector++) {
        ifroot       std::cout <<"\nNumber of electron in the total system:"<<Sector<<" in "<<  Nstate <<" states\n"
                               << "Hilbert space dimension:"<<dim_for_sector[Sector]<<"\n";

        time_t timeStartHam, timeEndHam;
        timeStartHam = clock();
        /*Construct Hamiltonian*/
        ifroot std::cout << "HamiltonianConstruc\n";
        HamiltonianConstruc( Himp, f_ann,   Nstate, SectorToGlobalIndx, Sector, f_annMat[Sector] , Local_Hamiltonian_ED, muTB , dim_for_sector[Sector], dim_for_sector[Sector-1]);
        timeEndHam =clock();
        ifroot std::cout << "construct:" << (timeEndHam - timeStartHam)/(CLOCKS_PER_SEC) <<"\n";

        /*diagonalization*/
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> ces(Himp);
        ces.compute(Himp);
        for(int alp=0; alp<Nstate; alp++) {
            f_annMat[Sector][alp] =  evec_prev.adjoint()  *  f_annMat[Sector][alp] * ces.eigenvectors();  //f_annMat in the many-body eigen state basis
        }
        for (int i=0; i<dim_for_sector[Sector] ; i++) {
            (ImpTotalEnergy.at(Sector)).at(i) = ces.eigenvalues()[i];
        }

        evec_prev = ces.eigenvectors();
    }//section

    ifroot std::cout << "write_cix_file\n";
    ifroot     write_cix_file(NSector,  maxDim, Nstate, dim_for_sector, ImpTotalEnergy, f_annMat);
    write_Delta(    weiss_field);
    ifroot write_PARMS();



    for (int N=0; N < NSector; N++) {
        delete []        SectorToGlobalIndx[N] ;
        delete []        f_annMat[N];
    }
    delete [] SectorToGlobalIndx;
    delete [] f_annMat;
    for (int Sector=1 ; Sector<=Nstate; Sector++) {
        delete []    f_ann[Sector] ;
    }
    delete [] f_ann;
} //ct_qmc_rut

void ctqmc_rutgers_seg(  Eigen::MatrixXcd Local_Hamiltonian_ED, double muTB,  ImgFreqFtn & weiss_field) {
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

            fprintf(cix,      "%+0.8f", getTotEng(Local_Hamiltonian_ED, muTB,alp));
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
            if( isOrbitalCorrinHart[Uindex[nH](0) ] and
                    isOrbitalCorrinHart[Uindex[nH](1) ] and
                    isOrbitalCorrinHart[Uindex[nH](2) ] and
                    isOrbitalCorrinHart[Uindex[nH](3) ] ) {
                fprintf(cix, " %d %d %d %d %+0.8f\n",
                        Hart2Corr[Uindex[nH](0)],
                        Hart2Corr[Uindex[nH](1)],
                        Hart2Corr[Uindex[nH](3)],
                        Hart2Corr[Uindex[nH](2)],  real(Utensor[nH]) );
                assert(imag(Utensor[nH]) < 1e-5);
            }
            nH++;
        }


        fprintf(cix, "# number of operators needed\n0");

        fclose(cix);


        ifroot std::cout << "write_cix_file\n";
        write_Delta_diag(    weiss_field );
        ifroot write_PARMS();
    }//ifroot
}

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
    if(mpi_rank==0) {
        std::string filename = std::string("delta_w_rutgers.dat");
        FILE *datap4 = fopen(filename.c_str(), "w");
        for(int w=0; w< weiss_field.getNFreq(); w++) {
            fprintf(datap4, "%0.8f", weiss_field.getValue(w) );
            for(int n=0; n< NSpinOrbit_per_atom; n++) {
                for(int m=n; m< NSpinOrbit_per_atom; m++) {
                    fprintf(datap4, "     %0.10f  %0.10f", real( weiss_field.getValue(w, n,m)  ),
                            imag( weiss_field.getValue(w, n,m)));
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


void write_cix_file(int NSector, int maxDim, int Nstate, unsigned long long * dim_for_sector,  std::vector<std::vector<double> > & ImpTotalEnergy, Eigen::MatrixXcd ** f_annMat) {
    FILE * cix = fopen("rutgers_input.cix", "w");
    fprintf(cix, "# Cix file for cluster DMFT with CTQMC\n");
    fprintf(cix, "# cluster_size, number of states, number of baths, maximum matrix size\n");
    fprintf(cix, "1 %d 1 %d\n",NSector, maxDim);
    fprintf(cix, "# baths, dimension, symmetry, global flip\n");
    fprintf(cix, "0  %d\n", Nstate);
    int upper_triangle=0;
    int lower_triangle=0;
    int deltaForRutgersCTQMS[Nstate][Nstate];

    for(int alp=0; alp<Nstate; alp++) {
        for(int bet=0; bet<Nstate; bet++) {
            if (alp==bet) {
                deltaForRutgersCTQMS[alp][bet] = upper_triangle;
                upper_triangle++;
            }
            else if (alp <bet) {
                deltaForRutgersCTQMS[alp][bet] = upper_triangle;
                upper_triangle++;
            }

        }
    }



    for(int alp=0; alp<Nstate; alp++) {
        for(int bet=0; bet<Nstate; bet++) {
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


    fprintf(cix, "# cluster energies for unique baths, eps[k]\n");
    for(int alp=0; alp<Nstate; alp++) {
        fprintf(cix, " 0");
    }
    fprintf(cix, "\n");



    fprintf(cix, "#   N   K   Sz size F^{+,0}, F^{+,1}, F^{+,2}, F^{+,3}, F^{+,4}, F^{+,5}, {Ea,Eb..} ;  {Sa, Sb..}\n");
    for (int Sector=0; Sector < NSector; Sector++) {
        int sectorIndex = Sector+1;
        fprintf(cix, "%d %d 0 0 %d", sectorIndex, Sector,  (int)dim_for_sector[Sector]);

        int nextSector=sectorIndex+1;
        if (nextSector>NSector) nextSector=0;
        for (int alp=0; alp < Nstate; alp++) {
            fprintf(cix, " %d",nextSector);
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
    for (int ketSectorIndex=1; ketSectorIndex <= NSector; ketSectorIndex++) {
        int braSectorIndex = ketSectorIndex+1;

        int ketSector = ketSectorIndex-1;
        int braSector = braSectorIndex-1; //=ketSector+1

        int ketDim = dim_for_sector[ketSectorIndex-1];
        int braDim = ( (braSector>=NSector) ? 0:  dim_for_sector[braSector] );


        braSectorIndex = ( (braSectorIndex>NSector) ? 0:  braSectorIndex );
        for(int alp=0; alp<Nstate; alp++) {
            fprintf(cix, "%d %d %d %d",ketSectorIndex, braSectorIndex, ketDim, braDim);
            Eigen::MatrixXcd  f_dagger;
            if(braSector<NSector)   f_dagger= f_annMat[braSector][alp].adjoint(); // <N+1|f^+ |N>  = <N| f | N+1>.adjoint()
            for (int ketN=0; ketN < ketDim; ketN++) {
                for (int braN=0; braN < braDim; braN++) {
                    fprintf(cix, " %e",real(f_dagger(braN,ketN)) );
                    assert(imag(f_dagger(braN,ketN)) < 1e-5);
                }
                fprintf(cix, "   ");
            }
            fprintf(cix, "\n");
        }
    }
    fprintf(cix, "HB2                # Hubbard-I is used to determine high-frequency\n");

    fprintf(cix, "# UCoulomb : (m1,s1) (m2,s2) (m3,s2) (m4,s1)  Uc[m1,m2,m3,m4]\n");
    int length = (int) Utensor.size();
    int nH=0;
    while(nH < length) {
        if( isOrbitalCorrinHart[Uindex[nH](0) ] and
                isOrbitalCorrinHart[Uindex[nH](1) ] and
                isOrbitalCorrinHart[Uindex[nH](2) ] and
                isOrbitalCorrinHart[Uindex[nH](3) ] ) {
            fprintf(cix, " %d %d %d %d %+0.8f\n",
                    Hart2Corr[Uindex[nH](0)],
                    Hart2Corr[Uindex[nH](1)],
                    Hart2Corr[Uindex[nH](3)],
                    Hart2Corr[Uindex[nH](2)],  real(Utensor[nH]) );
            assert(imag(Utensor[nH]) < 1e-5);
        }
        nH++;
    }


    fprintf(cix, "# number of operators needed\n0");

    fclose(cix);

}








void HamiltonianConstruc(  Eigen::MatrixXcd & Himp,  Eigen::MatrixXi **f_ann, int Nstate , unsigned long long ** SectorToGlobalIndx, int Sector, Eigen::MatrixXcd * f_annMat, Eigen::MatrixXcd Local_Hamiltonian_ED, double muTB, unsigned long long dim, unsigned long long dim_prev   ) {
    unsigned long long  rankf =       choose(Nstate-1,Sector-1);            //# of number state with occupied alp

    Eigen::MatrixXcd numOperat [NSpinOrbit_per_atom];
    for(int alp=0; alp<Nstate; alp++) {
        f_annMat[alp].setZero(dim_prev, dim);
        for (int k=0; k<rankf; k++) {
            f_annMat[alp](f_ann[Sector][alp](k,1),f_ann[Sector][alp](k,0)) =   f_ann[Sector][alp](k,2);
        }
        if(alp < NSpinOrbit_per_atom)        numOperat[alp] = f_annMat[alp].adjoint() * f_annMat[alp];
    }
    Himp.setZero(dim,dim);
    for(int alp=0; alp<NSpinOrbit_per_atom; alp++) {
        for(int bet=0; bet<NSpinOrbit_per_atom; bet++) {
            Himp +=  Local_Hamiltonian_ED(alp,bet) *  (f_annMat[alp].adjoint() * f_annMat[bet]) ; //H_loc e_ab f^+_a *f_b
        }//bet
        Himp -=  muTB *  numOperat[alp] ; //H_loc e_ab f^+_a *f_b
    }//alp

//dd interaction only TODO
    for(int alp=0; alp<NSpinOrbit_per_atom; alp+=2) { //0,2,4...NSpinOrbit_per_atom-2
        Himp += UHubb   * numOperat[alp] *numOperat[alp+1] ;
        for(int bet=0; bet<alp; bet+=2) {
            Himp +=  Uprime         *  numOperat[alp  ]  * numOperat[bet+1]; //H_interaction
            Himp +=  Uprime         *  numOperat[alp+1]  * numOperat[bet];   //H_interaction
            Himp += (Uprime-JHund)  *  numOperat[alp  ]  * numOperat[bet  ]; //H_interaction
            Himp += (Uprime-JHund)  *  numOperat[alp+1]  * numOperat[bet+1]; //H_interaction
        }
    }//alp
}



int F_element(int alpha, unsigned long long state1, unsigned long long state2) {
    if ( state1 == state2-pow(2,alpha) ) {
        int sign,i=0,N=0;
        for(i=0; i<alpha; i++) {
            N+=occupation(i,state2);
        }
        sign = (int)std::pow(-1,N);
        return sign ;
    }
    else return 0;
}

void SectorToGlobalIndxConstruction(unsigned long long ** SectorToGlobalIndx, int Nstate) {
    unsigned long long  i,j, S[Nstate+1];
    unsigned long long  Factorial ;
    for (int N=0; N<Nstate+1; N++) {
        S[N]=0;
    }
    for (i=0; i < pow(2,Nstate); i++) {//sjh
        j=i;
        int N=0;  //N=0,1,,,Nstate = Number of ptl
        while(j>0) {
            if(j%2 ==1) N+=1;
            j=j/2;
        }

        Factorial = choose(Nstate,N);
        if (  N>Nstate || S[N] > Factorial ) {
            ifroot std::cout<<"N:" << N<<"\nSN:" <<S[N] << "\nFactorial:" <<Factorial <<"\nchoose:"<<choose(Nstate,N);
            ifroot std::cout << "Error SOLVER_ED1\n";
            assert(0);
            break;
        }
        /////////////
        SectorToGlobalIndx[N][S[N]] = i;
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


int occupation(int alpha, unsigned long long state) {
    return ( (state/((int)pow(2,alpha))) %2 );
}


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
int F_dagger_sign(int i, unsigned long long alp) {
    int count = 0;
    for(int j=0; j<i; j++)
    {
        count += ((alp >> j) & 1);
    }
    return  std::pow(-1, count);
}
double getTotEng( Eigen::MatrixXcd Local_Hamiltonian_ED, double muTB,  unsigned long long   alp) {

    double TotEng=0;
    for(int i=0; i<NSpinOrbit_per_atom; i++) {
        TotEng +=  ( real(Local_Hamiltonian_ED(i,i))-muTB) *  occupation_Bits(i,alp);
    }//alp
    for(int i=0; i<NSpinOrbit_per_atom; i+=2) { //0,2,4...NSpinOrbit_per_atom-2
        TotEng += UHubb   * occupation_Bits(i,alp) *occupation_Bits(i+1,alp) ;
        for(int j=0; j<i; j+=2) {
            TotEng +=  Uprime         *  occupation_Bits(i  ,alp)  * occupation_Bits(j+1,alp); //H_interaction
            TotEng +=  Uprime         *  occupation_Bits(i+1,alp)  * occupation_Bits(j  ,alp);   //H_interaction
            TotEng += (Uprime-JHund)  *  occupation_Bits(i  ,alp)  * occupation_Bits(j  ,alp); //H_interaction
            TotEng += (Uprime-JHund)  *  occupation_Bits(i+1,alp)  * occupation_Bits(j+1,alp); //H_interaction
        }
    }//alp
    return TotEng;
}
