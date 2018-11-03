#include <time.h>
#include "mpi.h"
#include <fstream>
#include "ImgTimeFtn.h"
#include "tight_common.h"
#include "model.h"
#include "TB.h"
#include <Eigen/Eigenvalues>

//#define  nob(i,j)  ((i)*NumOrbit+(j))


//2015.12.12  jhsim
//2016.03.16  jhsim
//2016.05.29  jhsim
//2016.08.05  jhsim:overlap matrix
//2016.08.15  jhsim:non-orthogonal basis


std::vector<Eigen::MatrixXcd>  DF_CorrBase;

Eigen::VectorXd *KS_eigenEnergy;


///*used in opticalCond*/
//int NumLattx   = 7 ;
//int NumLatty   = 7 ;
//int NumLattz   = 7 ;
//int NumLattxy  = 49 ;
//int NumLattxyz =343  ; //7^3 ;


//low-energy model
std::vector<int> NBAND;
std::vector<std::vector<int> > FromValToKS;
void Time_print();




double **kmesh ;
double *kdist_band;


void TB_allocation(int knum,  std::vector<Eigen::MatrixXcd> &transformMatrix_k) ;
void TB_free();
int TB_allc_downfolding = 0;

void qs_spectra(ImgFreqFtn & SelfE_w,std::vector<Eigen::MatrixXcd>  KS_eigenVectors_orthoBasis,Eigen::VectorXd *KS_eigenEnergy, double  mu,
                std::vector<Eigen::MatrixXcd> H_k_inModelSpace       );

void read_HR(Eigen::MatrixXi &H_Rindex, Eigen::VectorXcd  &H_RMatrix, const std::string &hamiltonian, std::vector<int> & accumulated_Num_SpinOrbital, int NumOrbit);

void  downfolding_ftn
(
    int knum,  int knum_mpiGlobal,
    std::vector<int> &NBAND, std::vector<Eigen::MatrixXcd> &  H_k_inModelSpace, std::vector<Eigen::MatrixXcd> & KS_eigenVectors_orthoBasis, Eigen::VectorXd  * KS_eigenEnergy,
    double muDFT
);



void Find_best_correlated_basis(std::vector<Eigen::MatrixXcd> & H_k_inModelSpace, std::vector<Eigen::MatrixXcd> & SolverBasis, double muDFT) ;
void        NumMat_PWF(  int knum, int knum_mpiGlobal, double muDFT,
                         Eigen::MatrixXcd & NumMatrix, Eigen::VectorXd  * KS_eigenEnergy, std::vector<Eigen::MatrixXcd> & DF_CorrBase  ) ;


int Construct_Hk_Sk(
    int knum, int knum_mpiGlobal,   Eigen::MatrixXi  H_Rindex, Eigen::VectorXcd H_RMatrix, double ** kmesh, std::vector<int> & accumulated_Num_SpinOrbital,
    std::vector<Eigen::MatrixXcd> & H_k_inModelSpace,  std::vector<Eigen::MatrixXcd> & S_overlap,
    std::vector<Eigen::MatrixXcd> & KS_eigenVectors_orthoBasis, Eigen::VectorXd  * KS_eigenEnergy
);



void  ConstructModelHamiltonian
(
    int knum, int knum_mpiGlobal, double ** kmesh, std::vector<int> & accumulated_Num_SpinOrbital,
    std::vector<Eigen::MatrixXcd> & H_k_inModelSpace,  std::vector<Eigen::MatrixXcd> & S_overlap,
    std::vector<Eigen::MatrixXcd> & KS_eigenVectors_orthoBasis, std::vector<Eigen::MatrixXcd> & transformMatrix_k, Eigen::VectorXd  * KS_eigenEnergy,  int overlap_exist,
    double muDFT
);


void SpreadFtn_PWF( int knum,
                    std::vector<Eigen::MatrixXcd> S_overlap, std::vector<Eigen::MatrixXcd> &  transformMatrix_k,std::vector<Eigen::MatrixXcd>  KS_eigenVectors_orthoBasis
                    , std::vector<int> & accumulated_Num_SpinOrbital
                  ) ;



double Nele_non_Inter(
    int knum, int knum_mpiGlobal,
    std::vector<Eigen::MatrixXcd> & H_k_inModelSpace,  std::vector<Eigen::MatrixXcd> & S_overlap
);

double  TightBinding(double mu, const std::string &hamiltonian, ImgFreqFtn & SelfE_w,
                     ImgFreqFtn & weiss_fieldTBweakCorr, ImgFreqFtn & weiss_fieldTBstrongCorr,
                     int mu_adjustTB,  std::vector<Eigen::MatrixXcd> & SolverBasis
                    )
{

//mu_adjustTB == {-2,-1,0,1}
//== -1  Find initial chemical potential setting.
//== -2, just adjusting chemical potential
//==  0, no adjusting chem, Green ftn cal.
//==  1, adjusting chem, Green ftn cal.


    bool k_grid_dos;
    if (  (mode.find(std::string("band")) == std::string::npos ) or SOLVERtype != std::string("TB") or mu_adjustTB==-1 )
        k_grid_dos = true;
    else
        k_grid_dos = false;

    if(k_grid_dos) knum_mpiGlobal = k_pointx*k_pointy*k_pointz;    //dos
    else           knum_mpiGlobal = Nkpath*k_grid+1;               //band
    para_range(0,knum_mpiGlobal-1, mpi_numprocs,mpi_rank, &myksta, &mykend);
    knum= mykend-myksta+1;

    NumSubSpace = 1 ;

//    if (SOCCal == 0 ) NumSubSpace = 1;
//    if (interOrbitHop ==0 ) {
//        if ( SOCCal == 0 ) NumSubSpace = Spin_DegreeOfFreedom;
//        else {
//            std::cout << "InterOrbitalHoppin  option is not optimized yet\n";
//            NumSubSpace = 1;
//        }
//    }
//    if(doublecounting==1) {
//        NumSubSpace = 1;
//        ifroot  std::cout << "doublecounting  option is not optimized yet\n";
//    }

//////////////////////////////////
    std::cout.precision(3);



    std::vector<Eigen::MatrixXcd>  KS_eigenVectors_orthoBasis;
    std::vector<Eigen::MatrixXcd> H_k_inModelSpace;
    std::vector<Eigen::MatrixXcd> S_overlap;

    int    Atom1, Atom2;
    double    HopRe, HopIm;
    cmplx ctemp, ctemp2;
    int ntemp,ltemp,mtemp,a, k0;
    int ntemp2,ltemp2,mtemp2;
    cmplx w;

//read input
    assert(NumAtom>0  );
    assert(NumOrbit>0  );


    std::vector<Eigen::MatrixXcd> transformMatrix_k;
    TB_allocation(knum, transformMatrix_k);
    time_t timeStartTB, timeEndTB, timeStartLocal, timeEndLocal;

    timeStartTB = clock();
////////////////////////////////////////////////////////////
    ifroot std::cout << "******************\n";
    ifroot    Time_print();
    ifroot std::cout << "TB start.... " <<  mu_adjustTB<<"\n";
    ifroot std::cout << "******************\n";


    ifroot std::cout << "Chemical potential : " << mu <<"\n";
////////////////////////////////////////////////////////////
// Construct H(R) and Overlap matrix
////////////////////////////////////////////////////////////
    timeStartLocal = clock();
    /*Read H(R) */
    std::vector<int> accumulated_Num_SpinOrbital(NumAtom+1);
    Eigen::MatrixXi  H_Rindex;
    Eigen::VectorXcd H_RMatrix;
    read_HR(H_Rindex, H_RMatrix, hamiltonian, accumulated_Num_SpinOrbital, NumOrbit);
    MPI_Barrier(MPI_COMM_WORLD);



////////////////////////////////////////////////////////////
// Set k-grid
////////////////////////////////////////////////////////////




    /*FBZ, k-space grid*/
    if (  k_grid_dos ) {//dos
        int kpxy= k_pointx*k_pointy;                             // (kx,ky,kz) = kx + ky*(kpx) + kz*(kpxy)
        double  delx,dely,delz;
        delx       = 2.*pi / (ax * double(k_pointx));
        dely       = 2.*pi / (ay * double(k_pointy));
        delz       = 2.*pi / (az * double(k_pointz));
        for(int k=0; k<knum; k++) {
            kmesh[k][0]= -pi/ax + delx/2. + ((int)  (k+myksta)%k_pointx)               * delx;
            kmesh[k][1]= -pi/ay + dely/2. + ((int)  ((int)((k+myksta)%kpxy))/k_pointx) * dely;
            kmesh[k][2]= -pi/az + delz/2. + ((int)  (k+myksta)/kpxy)                   * delz;
        }
    }
    else   { //band
        for(int k=0; k<knum; k++) {
            int i0=(k+myksta)/k_grid;
            int k_in_path=(k+myksta)%k_grid;
            if(k+myksta < Nkpath*k_grid) {
                kmesh[k][0] = KpathPoint[i0*3+0] *(2.*pi)/ax +  (KpathPoint[(i0+1)*3+0]-KpathPoint[i0*3+0]) * (k_in_path*2.0*pi)/(k_grid*ax);
                kmesh[k][1] = KpathPoint[i0*3+1] *(2.*pi)/ay +  (KpathPoint[(i0+1)*3+1]-KpathPoint[i0*3+1]) * (k_in_path*2.0*pi)/(k_grid*ay);
                kmesh[k][2] = KpathPoint[i0*3+2] *(2.*pi)/az +  (KpathPoint[(i0+1)*3+2]-KpathPoint[i0*3+2]) * (k_in_path*2.0*pi)/(k_grid*az);
            }
            else if (k+myksta == Nkpath*k_grid) {
                kmesh[k][0] = KpathPoint[i0*3+0] *(2.*pi)/ax ;
                kmesh[k][1] = KpathPoint[i0*3+1] *(2.*pi)/ay ;
                kmesh[k][2] = KpathPoint[i0*3+2] *(2.*pi)/az ;
            }
        }
        kdist_band[0] = 0.0;
        for(int k=1; k<knum_mpiGlobal; k++) {
            int k0 = (k-1)/k_grid;
            int k1 = k0+1;
            kdist_band[k] = kdist_band[k-1];
            kdist_band[k]+=(1.0)/(k_grid) *
                           (   (KpathPoint[k1*3+0] *RecipUnitVector_b1 + KpathPoint[k1*3+1] *RecipUnitVector_b2+  KpathPoint[k1*3+2] *RecipUnitVector_b3)
                               -(KpathPoint[k0*3+0] *RecipUnitVector_b1 + KpathPoint[k0*3+1] *RecipUnitVector_b2+  KpathPoint[k0*3+2] *RecipUnitVector_b3) ).norm();
        }
    }
    ifroot printf("kmesh was created ....\n");
////////////////////////////////////////////////////////////
//Construct Hk, overlapmatrix, S(k)*/
////////////////////////////////////////////////////////////
    NBAND.resize(knum);
    H_k_inModelSpace.resize(knum);
    S_overlap.resize(knum);
    DF_CorrBase.resize(knum);
    KS_eigenVectors_orthoBasis.resize(knum);
    TB_allc_downfolding = 1;


    int overlap_exist =  Construct_Hk_Sk(
                             knum, knum_mpiGlobal,   H_Rindex, H_RMatrix,  kmesh,  accumulated_Num_SpinOrbital,
                             H_k_inModelSpace,   S_overlap,
                             KS_eigenVectors_orthoBasis,  KS_eigenEnergy
                         );

    if(mu_adjustTB == -1) {
        mu = Nele_non_Inter(knum, knum_mpiGlobal, H_k_inModelSpace, S_overlap);


        low_energy_subspace_in_KS_basis(knum, knum_mpiGlobal, NBAND, FromValToKS, mu,
                                        KS_eigenVectors_orthoBasis, KS_eigenEnergy);


        ConstructModelHamiltonian (   knum,  knum_mpiGlobal,  kmesh, accumulated_Num_SpinOrbital,
                                      H_k_inModelSpace, S_overlap, KS_eigenVectors_orthoBasis, transformMatrix_k,  KS_eigenEnergy, overlap_exist, mu);


        downfolding_ftn(knum, knum_mpiGlobal, NBAND, H_k_inModelSpace, KS_eigenVectors_orthoBasis, KS_eigenEnergy,  mu);

        Find_best_correlated_basis(H_k_inModelSpace, SolverBasis, mu);
        SpreadFtn_PWF(knum, S_overlap, transformMatrix_k, KS_eigenVectors_orthoBasis, accumulated_Num_SpinOrbital);
        NumMat_PWF(knum, knum_mpiGlobal, mu, NumMatrix, KS_eigenEnergy, DF_CorrBase);

        return mu;
    }


    low_energy_subspace_in_KS_basis(knum, knum_mpiGlobal, NBAND, FromValToKS, muDFT,
                                    KS_eigenVectors_orthoBasis, KS_eigenEnergy);

    ConstructModelHamiltonian (   knum,  knum_mpiGlobal,  kmesh, accumulated_Num_SpinOrbital,
                                  H_k_inModelSpace, S_overlap, KS_eigenVectors_orthoBasis, transformMatrix_k,  KS_eigenEnergy, overlap_exist, muDFT);


    downfolding_ftn(knum, knum_mpiGlobal, NBAND, H_k_inModelSpace, KS_eigenVectors_orthoBasis, KS_eigenEnergy,  muDFT);



////////////////////////////////////////////////////////////
//TB SOLV: solve model hamiltonian, baand, dos, greens'ftn...
////////////////////////////////////////////////////////////

    /*Band & DOS real freq.*/
    if (SOLVERtype==std::string("TB")) {
        qs_spectra( SelfE_w, KS_eigenVectors_orthoBasis,KS_eigenEnergy,    mu, H_k_inModelSpace );
    }
    else {
        std::vector<Eigen::MatrixXcd> densityMatDFT;
        densityMatDFT.resize(knum);
        /*Find chemical potential and Green's ftn*/
        if(mu_adjustTB ==-2 or mu_adjustTB ==1)
            mu = GreenFtn_w_adjust_mu(H_k_inModelSpace,SelfE_w,mu,1,mu_adjustTB);
        if (mu_adjustTB == -2) {
            ifroot  printf("We have found chemical potential mu\n");
            TB_free();
            return mu;
        }
        /*CONSTRUCT diagonal hybridization Ftn delta_w*/
        std::vector<Eigen::MatrixXcd>  Gw;
        GreenFtn_w( NumCluster, NumHartrOrbit_per_cluster,  H_k_inModelSpace,  SelfE_w, Gw, mu, densityMatDFT ) ;
        ifroot  printf("We have found GreenFtn and chemical potential mu\n");



        for(int clust=0; clust < NumCluster; clust++) {
            std::vector<int> weaklyCorr(NumHartrOrbit_per_cluster);
            std::vector<int> strongCorr(NSpinOrbit_per_atom);
            for (int i=0; i<NumHartrOrbit_per_cluster; i++) weaklyCorr[i] = clust* NumHartrOrbit_per_cluster +i;

            Construct_hyb_delta ( NumHartrOrbit_per_cluster, weaklyCorr,  SelfE_w, Gw,  mu, weiss_fieldTBweakCorr, clust, SolverBasis);
            if(NSpinOrbit_per_atom>0) {
                for(int atom=clust*NumAtom_per_cluster; atom < (clust+1)*NumAtom_per_cluster; atom++) {
                    for (int i=0; i<NSpinOrbit_per_atom; i++)  strongCorr[i] = CorrToHartr[atom * NSpinOrbit_per_atom + i ];
                    Construct_hyb_delta (   NSpinOrbit_per_atom, strongCorr, SelfE_w, Gw,  mu, weiss_fieldTBstrongCorr, atom, SolverBasis);
                }
            }
        }



        if(mpi_rank ==0) {
            std::cout << "FILEOUT:Numele.dat\n";
            FILE * OBSERV=fopen("Numele.dat","w");
            /*Numele.dat*/
            for (int orb1=0; orb1<N_peratom_HartrOrbit*NumCorrAtom; orb1++) {
                fprintf(OBSERV, "%d",orb1);
                for (int orb2=0; orb2<N_peratom_HartrOrbit*NumCorrAtom; orb2++) {
                    fprintf(OBSERV, "     %0.3f %0.3f",real(NumMatrix(orb1,orb2)), imag(NumMatrix(orb1,orb2)));
                }
                fprintf(OBSERV, "\n");
            }
            fclose(OBSERV);
        }//mpi_rank
        upfolding_density(densityMatDFT,  KS_eigenVectors_orthoBasis,H_Rindex, mu, S_overlap, transformMatrix_k);
    }

////////////////////////////////////////////////////////////
//Free energy calculation
////////////////////////////////////////////////////////////
//        free_energy_cal( Sw,mu) ;
////////////////////////////////////////////////////////////
    TB_free();
    timeEndTB = clock(); //time
    if (mpi_rank ==0) {
        std::ofstream FINAL("TB_FINISH");
        FINAL << "TB_NORMALLY_DONE";
        FINAL.close() ;

        FILE *tempFile;   //mu.dat
        tempFile = fopen("mu_history.out", "a");
        fprintf(tempFile, "%0.20f\n", mu);
        fclose(tempFile) ;
    }
    return mu;
}






void TB_allocation(int knum,  std::vector<Eigen::MatrixXcd> &transformMatrix_k) {
    /*Memory allocation*/
    kmesh = new double * [knum];
    kdist_band = new double [knum_mpiGlobal];
    KS_eigenEnergy =  new Eigen::VectorXd  [knum];
    transformMatrix_k.resize(knum);
    for(int k=0; k< knum; k++) {
        kmesh[k] = new double [dimension_3d];
        KS_eigenEnergy[k].setZero(NumOrbit);
        transformMatrix_k[k].setIdentity(NumOrbit,NumOrbit);
    }
}


void TB_free() {
    for(int k = 0; k < knum; k++) {
        delete [] kmesh[k];
    }
    delete [] kdist_band;
    delete [] kmesh;
    delete [] KS_eigenEnergy;
    if (TB_allc_downfolding==1) {
        TB_allc_downfolding = 0;
    }
}



void Find_best_correlated_basis(std::vector<Eigen::MatrixXcd> & H_k_inModelSpace, std::vector<Eigen::MatrixXcd> & SolverBasis, double muDFT) {



    for(int ATOM=0 ; ATOM < NumCorrAtom; ATOM++) {
        SolverBasis[ATOM].setIdentity(N_peratom_HartrOrbit, N_peratom_HartrOrbit);
    }



    if(impurityBasisSwitch==1) {


//        Eigen::MatrixXcd  Gloc_w0[NumCorrAtom], Gloc_w0_mpilocal[NumCorrAtom];
        std::vector<Eigen::MatrixXcd>          Gloc_w0(NumCorrAtom);
        std::vector<Eigen::MatrixXcd> Gloc_w0_mpilocal(NumCorrAtom);
        Eigen::MatrixXcd Gkw_w0;
        for(int ATOM=0 ; ATOM < NumCorrAtom; ATOM++) {
            Gloc_w0[ATOM].setZero(NSpinOrbit_per_atom,NSpinOrbit_per_atom);
            Gloc_w0_mpilocal[ATOM].setZero(NSpinOrbit_per_atom,NSpinOrbit_per_atom);
        }

        for(int k=0; k<knum; k++) {
            Gkw_w0 = -H_k_inModelSpace[k];
            for(int i0=0; i0<NBAND[k]; i0++) {
                Gkw_w0(i0,i0) += muDFT;
            }

            Gkw_w0 = Gkw_w0.inverse();
            for(int ATOM=0 ; ATOM < NumCorrAtom; ATOM++) {
                for(int i0=0; i0<NSpinOrbit_per_atom; i0+=1) {
                    for( int m0=0; m0<NSpinOrbit_per_atom; m0+=1) {
                        int  i0F = CorrIndex[ATOM*NSpinOrbit_per_atom+i0];
                        int  m0F = CorrIndex[ATOM*NSpinOrbit_per_atom+m0];
                        Gloc_w0_mpilocal[ATOM](i0,m0) +=Gkw_w0(i0F,m0F);
                    }
                }
            }

        }//k

        for(int ATOM=0 ; ATOM < NumCorrAtom; ATOM++) {
            MPI_Allreduce(Gloc_w0_mpilocal[ATOM].data(), Gloc_w0[ATOM].data(), Gloc_w0[ATOM].size(), MPI_DOUBLE_COMPLEX, MPI_SUM,  MPI_COMM_WORLD);
            Gloc_w0[ATOM] /= knum_mpiGlobal;
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> ces( NSpinOrbit_per_atom );
            ces.compute(Gloc_w0[ATOM]);

            for (int n=0; n<NSpinOrbit_per_atom; n++) {
                for (int m=0; m<NSpinOrbit_per_atom; m++) {
                    int n0 = CorrToHartr[n];
                    int m0 = CorrToHartr[m];

                    SolverBasis.at(ATOM)(n0,m0) =  ces.eigenvectors()(n,m);
                }
            }
        }
        for(int ATOM=0 ; ATOM < NumCorrAtom; ATOM++) {
            ifroot {
                std::cout << "Solver:Gloc_w0:\n" << Gloc_w0[ATOM] <<"\n";
                std::cout <<"\n";
                std::cout << "Solver:Basis:\n" << SolverBasis[ATOM] <<"\n";
            }
        }

    }//impurityBasisSwitch
}

void        NumMat_PWF(  int knum, int knum_mpiGlobal, double muDFT,
                         Eigen::MatrixXcd & NumMatrix, Eigen::VectorXd  * KS_eigenEnergy, std::vector<Eigen::MatrixXcd> & DF_CorrBase  ) {


    Eigen::MatrixXcd NumMatrix_mpilocal;
    NumMatrix_mpilocal.setZero(NumCorrAtom*N_peratom_HartrOrbit, NumCorrAtom*N_peratom_HartrOrbit);
    NumMatrix.setZero(NumCorrAtom*N_peratom_HartrOrbit, NumCorrAtom*N_peratom_HartrOrbit);
    for (int k=0; k<knum; k++) {
        Eigen::VectorXd FermiDirac(NBAND[k]);
        for (int i=0; i<NBAND[k]; i++) FermiDirac(i) =   1./ (1.0 + (std::exp( beta* (KS_eigenEnergy[k][FromValToKS[k][i]] - muDFT  )   )));
        Eigen::MatrixXcd temp;
        temp = DF_CorrBase[k] * FermiDirac.asDiagonal() * DF_CorrBase[k].adjoint();
        NumMatrix_mpilocal += temp.block(0,0,NumCorrAtom*N_peratom_HartrOrbit,NumCorrAtom*N_peratom_HartrOrbit);
    }
    NumMatrix_mpilocal /= knum_mpiGlobal;
    MPI_Allreduce(NumMatrix_mpilocal.data(), NumMatrix.data(), NumMatrix.size(), MPI_DOUBLE_COMPLEX, MPI_SUM,  MPI_COMM_WORLD);

    ifroot  std::cout << "<TB> Occupancy, PWF:\n";
    for(int h1F=0; h1F<N_peratom_HartrOrbit; h1F++) {
        ifroot std::cout << h1F << " " << NumMatrix(h1F,h1F)<<"\n" ;
    }

}
