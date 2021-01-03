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

void qs_spectra(ImgFreqFtn & SelfE_z,std::vector<Eigen::MatrixXcd>  KS_eigenVectors_orthoBasis,Eigen::VectorXd *KS_eigenEnergy, double  mu,
                std::vector<Eigen::MatrixXcd> H_k_inModelSpace, Eigen::MatrixXcd & SolverBasis      );

//void read_HR(Eigen::MatrixXi &H_Rindex, Eigen::VectorXcd  &H_RMatrix, const std::string &hamiltonian, std::vector<int> & accumulated_Num_SpinOrbital, int NumOrbit);

void  downfolding_ftn
(
    int knum,  int knum_mpiGlobal,
    std::vector<int> &NBAND, std::vector<Eigen::MatrixXcd> &  H_k_inModelSpace, std::vector<Eigen::MatrixXcd> & KS_eigenVectors_orthoBasis, Eigen::VectorXd  * KS_eigenEnergy,
    double muDFT
);



void Find_best_correlated_basis (std::vector<Eigen::MatrixXcd> & H_k_inModelSpace, Eigen::MatrixXcd & SolverBasis, double muDFT) ;
void        NumMat_PWF(  int knum, int knum_mpiGlobal, double muDFT,
                         Eigen::MatrixXcd & NumMatrix, Eigen::VectorXd  * KS_eigenEnergy, std::vector<Eigen::MatrixXcd> & DF_CorrBase  ) ;


int Construct_Hk_Sk(  const std::string &hamiltonian,
                      int knum, int knum_mpiGlobal,   Eigen::MatrixXi  H_Rindex, double ** kmesh, std::vector<int> & accumulated_Num_SpinOrbital,
                      std::vector<Eigen::MatrixXcd> & H_k_inModelSpace,  std::vector<Eigen::MatrixXcd> & S_overlap,
                      std::vector<Eigen::MatrixXcd> & KS_eigenVectors_orthoBasis, Eigen::VectorXd  * KS_eigenEnergy
                   );



void  ConstructModelHamiltonian
(
    int knum, int knum_mpiGlobal, double ** kmesh, std::vector<int> & accumulated_Num_SpinOrbital,
    std::vector<Eigen::MatrixXcd> & H_k_inModelSpace,  //std::vector<Eigen::MatrixXcd> & S_overlap,
    std::vector<Eigen::MatrixXcd> & KS_eigenVectors_orthoBasis, std::vector<Eigen::MatrixXcd> & transformMatrix_k, Eigen::VectorXd  * KS_eigenEnergy,  int overlap_exist,
    double muDFT
);


//void SpreadFtn_PWF( int knum,
//                    std::vector<Eigen::MatrixXcd> S_overlap, std::vector<Eigen::MatrixXcd> &  transformMatrix_k,std::vector<Eigen::MatrixXcd>  KS_eigenVectors_orthoBasis
//                    , std::vector<int> & accumulated_Num_SpinOrbital
//                  ) ;
//


double Nele_non_Inter(
    int knum, int knum_mpiGlobal, Eigen::VectorXd  * KS_eigenEnergy
//    std::vector<Eigen::MatrixXcd> & H_k_inModelSpace,  std::vector<Eigen::MatrixXcd> & S_overlap
);

double  TightBinding(double mu, const std::string &hamiltonian, ImgFreqFtn & SelfE_z,
                     ImgFreqFtn & weiss_fieldTBweakCorr, ImgFreqFtn & weiss_fieldTBstrongCorr,
                     int mu_adjustTB, Eigen::MatrixXcd & SolverBasis
                    )
{

//mu_adjustTB == {-2,-1,0,1}
//== -1  Find initial chemical potential (non-interacting system) , and NumMat solverbasis
//== -2, adjusting chem, exit
//==  1, adjusting chem, Green ftn cal.
//==  0, no adjusting chem, Green ftn cal.


    bool k_grid_dos;
    if (  (mode.find(std::string("band")) == std::string::npos ) or SOLVERtype != std::string("TB") or not(mu_adjustTB==0)  )
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



////////////////////////////////////////////////////////////
// Set k-grid
////////////////////////////////////////////////////////////
    /*FBZ, k-space grid*/
// (kx,ky,kz) = kx + ky*(k_pointx) + kz*(kpxy)
// =>
// kx=  k % k_pointx
// ky=  (k % kpxy) / k_pointx
// kx=   k / kpxy

    if (  k_grid_dos ) {//dos
        int kpxy= k_pointx*k_pointy;
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

    // std::vector<Eigen::MatrixXcd> S_overlap;
    // S_overlap.resize(knum);

    DF_CorrBase.resize(knum);
    KS_eigenVectors_orthoBasis.resize(knum);
    TB_allc_downfolding = 1;


    /*Read H(R) S(R) and tromsform to Hk, Sk */
    std::vector<int> accumulated_Num_SpinOrbital(NumAtom+1);
    Eigen::MatrixXi  H_Rindex;
//    Eigen::VectorXcd H_RMatrix;
//    read_HR(H_Rindex, H_RMatrix, hamiltonian, accumulated_Num_SpinOrbital, NumOrbit);
    MPI_Barrier(MPI_COMM_WORLD);
    int overlap_exist =  Construct_Hk_Sk( hamiltonian,
                                          knum, knum_mpiGlobal,   H_Rindex,  kmesh,  accumulated_Num_SpinOrbital,
                                          H_k_inModelSpace,   transformMatrix_k,
                                          KS_eigenVectors_orthoBasis,  KS_eigenEnergy
                                        );         // Here transformMatrix_k = Sk
    // KS_eigenVectors_orthoBasis   =   non-orthogonal AO basis, <AO^i | kn >

    if(mu_adjustTB == -1) {
        muDFT = Nele_non_Inter(knum, knum_mpiGlobal,KS_eigenEnergy);
    }


    ConstructModelHamiltonian (   knum,  knum_mpiGlobal,  kmesh, accumulated_Num_SpinOrbital,
                                  H_k_inModelSpace, KS_eigenVectors_orthoBasis, transformMatrix_k,  KS_eigenEnergy, overlap_exist, muDFT);  //KS_eigenVectors_orthoBasis  written in orthogonal basis

    low_energy_subspace_in_KS_basis(knum, knum_mpiGlobal, NBAND, FromValToKS, muDFT,
                                    KS_eigenVectors_orthoBasis, KS_eigenEnergy);

    downfolding_ftn(knum, knum_mpiGlobal, NBAND, H_k_inModelSpace, KS_eigenVectors_orthoBasis, KS_eigenEnergy,  muDFT);


////////////////////////////////////////////////////////////
//TB SOLV: solve model hamiltonian, baand, dos, greens'ftn...
////////////////////////////////////////////////////////////
    // initialize: SolverBasis
    if(mu_adjustTB == -1) {
//        SpreadFtn_PWF(knum, S_overlap, transformMatrix_k, KS_eigenVectors_orthoBasis, accumulated_Num_SpinOrbital);
        NumMat_PWF(knum, knum_mpiGlobal, muDFT, NumMatrix, KS_eigenEnergy, DF_CorrBase);
        Find_best_correlated_basis(H_k_inModelSpace, SolverBasis, muDFT);

        TB_free();
        return muDFT;
    }
    if(mu_adjustTB ==-2 ) {
        mu = GreenFtn_iw_adjust_mu(H_k_inModelSpace,SelfE_z,mu,1,mu_adjustTB);
        ifroot  printf("We have found chemical potential mu\n");
        TB_free();
        return mu;
    }
    /*Band & DOS real freq.*/
    if (SOLVERtype==std::string("TB")) {
        NumMat_PWF(knum, knum_mpiGlobal, muDFT, NumMatrix, KS_eigenEnergy, DF_CorrBase);
        Find_best_correlated_basis(H_k_inModelSpace, SolverBasis, muDFT);
        qs_spectra( SelfE_z, KS_eigenVectors_orthoBasis,KS_eigenEnergy,    mu, H_k_inModelSpace, SolverBasis );
        TB_free();
        return mu;
    }
    else {
        std::vector<Eigen::MatrixXcd> densityMatDFT;
        densityMatDFT.resize(knum);
        /*Find chemical potential and Green's ftn*/
//        if(mu_adjustTB ==-2 or mu_adjustTB ==1)
        if( mu_adjustTB ==1)
            mu = GreenFtn_iw_adjust_mu(H_k_inModelSpace,SelfE_z,mu,1,mu_adjustTB);
        //if (mu_adjustTB == -2) {
        //    ifroot  printf("We have found chemical potential mu\n");
        //    TB_free();
        //    return mu;
        //}
        /*CONSTRUCT diagonal hybridization Ftn delta_w*/
        std::vector<Eigen::MatrixXcd>  Gw;
        GreenFtn_w( NumCluster, NumHartrOrbit_per_cluster,  H_k_inModelSpace,  SelfE_z, Gw, mu, densityMatDFT ) ;
        ifroot  printf("We have found GreenFtn and chemical potential mu\n");



        for(int clust=0; clust < NumCluster; clust++) {
            ifroot std::cout <<"Construct weiss field for cluster " << clust << std::endl;
            std::vector<int> weaklyCorr(NumHartrOrbit_per_cluster);
            std::vector<int> strongCorr(NSpinOrbit_per_atom);
            for (int i=0; i<NumHartrOrbit_per_cluster; i++) weaklyCorr[i] = clust* NumHartrOrbit_per_cluster +i;

            Construct_hyb_delta ( NumHartrOrbit_per_cluster, weaklyCorr,  SelfE_z, Gw,  mu, weiss_fieldTBweakCorr, clust, SolverBasis, 0);
            if(NSpinOrbit_per_atom>0) {
                for(int atom=clust*NumAtom_per_cluster; atom < (clust+1)*NumAtom_per_cluster; atom++) {
                    ifroot std::cout <<"Construct weiss field for atom " << atom << std::endl;
                    ifroot std::cout <<"Orbital: ";
                    for (int i=0; i<NSpinOrbit_per_atom; i++) {
                        strongCorr[i] = CorrToHartr(atom, i );
                        ifroot std::cout << strongCorr[i] <<" " ;
                    }
                    ifroot std::cout << std::endl;
                    Construct_hyb_delta (   NSpinOrbit_per_atom, strongCorr, SelfE_z, Gw,  mu, weiss_fieldTBstrongCorr, atom, SolverBasis, segmentsolver);
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
//        upfolding_density(densityMatDFT,  KS_eigenVectors_orthoBasis,H_Rindex, mu, S_overlap, transformMatrix_k);
        upfolding_density(densityMatDFT,  KS_eigenVectors_orthoBasis,H_Rindex, mu,  transformMatrix_k);
    }//DMFT

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




void        NumMat_PWF(  int knum, int knum_mpiGlobal, double muDFT,
                         Eigen::MatrixXcd & NumMatrix, Eigen::VectorXd  * KS_eigenEnergy, std::vector<Eigen::MatrixXcd> & DF_CorrBase  ) {


    Eigen::MatrixXcd NumMatrix_mpilocal;
    NumMatrix_mpilocal.setZero(NumCorrAtom*N_peratom_HartrOrbit, NumCorrAtom*N_peratom_HartrOrbit);
    NumMatrix.setZero(NumCorrAtom*N_peratom_HartrOrbit, NumCorrAtom*N_peratom_HartrOrbit);
    for (int k=0; k<knum; k++) {
        Eigen::VectorXd FermiDirac(NBAND[k]);
        for (int i=0; i<NBAND[k]; i++) FermiDirac(i) =   1./ (1.0 + (std::exp( beta_smearing* (KS_eigenEnergy[k][FromValToKS[k][i]] - muDFT  )   )));
        Eigen::MatrixXcd temp;
        temp = DF_CorrBase[k] * FermiDirac.asDiagonal() * DF_CorrBase[k].adjoint();
        NumMatrix_mpilocal += temp.block(0,0,NumCorrAtom*N_peratom_HartrOrbit,NumCorrAtom*N_peratom_HartrOrbit);
    }
    NumMatrix_mpilocal /= knum_mpiGlobal;
    MPI_Allreduce(NumMatrix_mpilocal.data(), NumMatrix.data(), NumMatrix.size(), MPI_DOUBLE_COMPLEX, MPI_SUM,  MPI_COMM_WORLD);

    //ifroot  std::cout << "<TB> Occupancy, PWF:\n";
    //for(int at=0; at<NumCorrAtom; at++) {
    //    double Nat=0;
    //    for(int h1F=at*N_peratom_HartrOrbit; h1F<(at+1)*N_peratom_HartrOrbit; h1F++) {
    //        double Nh1F = real(NumMatrix(h1F,h1F));
    //        Nat += Nh1F;
    //        ifroot std::cout << h1F << " " << Nh1F <<"\n" ;
    //    }
    //    ifroot std::cout << "Atom" << at  << ": " << Nat<<"\n";
    //}
}



void Find_best_correlated_basis(std::vector<Eigen::MatrixXcd> & H_k_inModelSpace, Eigen::MatrixXcd & SolverBasis, double muDFT) {
//
//
    int NumMatBlock ;
    int rotMatBlockSize;
    int BlockDim;
    if(ClusterBasis==2) {
        rotMatBlockSize  = NSpinOrbit_per_atom* NumAtom_per_cluster;
        NumMatBlock = NumCluster;
        BlockDim    = N_peratom_HartrOrbit*NumAtom_per_cluster;
    }
    else
    {
        rotMatBlockSize = NSpinOrbit_per_atom;
        NumMatBlock = NumCorrAtom;
        BlockDim    = N_peratom_HartrOrbit;
    }
//    int spinsize=1;
//    if (magnetism==0 or magnetism==1)  spinsize=2;

    int spinsize=2;
    rotMatBlockSize /= spinsize;

    SolverBasis.setIdentity(NumCorrAtom*N_peratom_HartrOrbit, NumCorrAtom*N_peratom_HartrOrbit);
    std::vector<Eigen::MatrixXcd > SolverBasis_cl(NumMatBlock);
//
    if(impurityBasisSwitch==3 or  impurityBasisSwitch==13) {
/////////////////
//Read SolverBasis
/////////////////
            std::ifstream  input(std::string("./SolverBasis.dat").c_str());
            double reN, imN;
            for(int n=0; n< NumCorrAtom*N_peratom_HartrOrbit;  n++) {
                for(int m=0; m< NumCorrAtom*N_peratom_HartrOrbit; m++) {
                    input >>  reN;
                    input >>  imN;
                    SolverBasis(n,m) = reN + I*imN ;
                }
            }
            ifroot {std::cout << "Solver: basis read from file:\n" ;}
    }
    else if(impurityBasisSwitch != 0 ) {
        for(int cl=0 ; cl < NumMatBlock; cl++) {
            SolverBasis_cl[cl].setIdentity(BlockDim, BlockDim);
        }

        std::vector<Eigen::MatrixXcd>          Gloc_w0(NumMatBlock);
        std::vector<Eigen::MatrixXcd> Gloc_w0_mpilocal(NumMatBlock);
        Eigen::MatrixXcd Gkw_w0;
        for(int cl=0 ; cl < NumMatBlock; cl++) {
            Gloc_w0[cl].setZero(rotMatBlockSize,rotMatBlockSize);
            Gloc_w0_mpilocal[cl].setZero(rotMatBlockSize,rotMatBlockSize);
        }
        if(impurityBasisSwitch==2 or  impurityBasisSwitch==12) {
/////////////////
//Diag Heff
/////////////////
            for(int k=0; k<knum; k++) {
                Gkw_w0 = -H_k_inModelSpace[k]  ;
                for(int i0=0; i0<NBAND[k]; i0++) {
                    Gkw_w0(i0,i0) += muDFT + I*infinitesimal;
                }
                Gkw_w0 = Gkw_w0.inverse();

                for(int cl=0 ; cl < NumMatBlock; cl++) {
                    for(int i0=0; i0<rotMatBlockSize; i0+=1) {
                        for( int m0=0; m0<rotMatBlockSize; m0+=1) {
                            int  i0F = CorrIndex[ (cl*rotMatBlockSize+ i0 )  *  spinsize    ];
                            int  m0F = CorrIndex[ (cl*rotMatBlockSize+ m0 )  *  spinsize    ];
                            Gloc_w0_mpilocal[cl](i0,m0) += Gkw_w0(i0F,m0F);
                            if(spinsize==2) Gloc_w0_mpilocal[cl](i0,m0) += Gkw_w0(i0F+1,m0F+1);
                        }
                    }
                }
            }//k
            for(int cl=0 ; cl < NumMatBlock; cl++) {
                MPI_Allreduce(Gloc_w0_mpilocal[cl].data(), Gloc_w0[cl].data(), Gloc_w0[cl].size(), MPI_DOUBLE_COMPLEX, MPI_SUM,  MPI_COMM_WORLD);
                Gloc_w0[cl] /= knum_mpiGlobal;
                Gloc_w0[cl] *= -1;
            }
            ifroot {
                std::cout << "Solver:Heff diagonal basis:\n" ;
            }
        }
        else if (impurityBasisSwitch==1 or impurityBasisSwitch==11) {

///////////////////
//Or DenMat
///////////////////
            for(int cl=0 ; cl < NumMatBlock; cl++) {
                for(int i0=0; i0<rotMatBlockSize; i0+=1) {
                    for( int m0=0; m0<rotMatBlockSize; m0+=1) {
                        int  i0F =    KS2Hartr[CorrIndex[(cl*rotMatBlockSize+ i0 )  *  spinsize   ]];
                        int  m0F =    KS2Hartr[CorrIndex[(cl*rotMatBlockSize+ m0 )  *  spinsize   ]];
//                    std::cout << i0 <<" " << i0F <<"\n" ;
                        Gloc_w0_mpilocal[cl](i0,m0) += NumMatrix(i0F,m0F);
                        if(spinsize==2) Gloc_w0_mpilocal[cl](i0,m0) += NumMatrix(i0F+1,m0F+1);
                    }
                }
                Gloc_w0[cl] =  -Gloc_w0_mpilocal[cl].inverse();
            }
///////////////////
            ifroot {
                std::cout << "Solver:densMat diagonal basis:\n" ;
            }
        }


        for(int cl=0 ; cl < NumMatBlock; cl++) {


            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> ces( rotMatBlockSize );
            Eigen::MatrixXcd Heff = -(Gloc_w0[cl].inverse() + Gloc_w0[cl].inverse().adjoint())/2.;
//            Eigen::MatrixXcd Heff = NumMatrix;
            ces.compute( Heff);
            if(ClusterBasis==2) {
                for(int at1=cl*NumAtom_per_cluster; at1 < (cl+1)*NumAtom_per_cluster; at1++) {
                    for(int at2=cl*NumAtom_per_cluster; at2 < (cl+1)*NumAtom_per_cluster; at2++) {
                        for (int n=0; n< NSpinOrbit_per_atom; n++) {
                            for (int m=0; m< NSpinOrbit_per_atom; m++) {


                                int n0 =   (CorrToHartr(at1,n) - cl*BlockDim )  /spinsize   ;
                                int m0 =   (CorrToHartr(at2,m) - cl*BlockDim )  /spinsize   ;

                                int n1 =   (n + (at1-cl*NumAtom_per_cluster) * NSpinOrbit_per_atom )/ spinsize  ;
                                int m1 =   (m + (at2-cl*NumAtom_per_cluster) * NSpinOrbit_per_atom )/ spinsize ;


                                SolverBasis_cl.at(cl)(spinsize*n0,spinsize*m0) =  ces.eigenvectors()(n1, m1);
                                if (spinsize==2)           SolverBasis_cl.at(cl)(spinsize*n0+1,spinsize*m0+1) =  ces.eigenvectors()(n1, m1);
                            }
                        }
                    }
                }
            }
            else {
                for (int n=0; n< NSpinOrbit_per_atom; n++) {             //rotMatBlockSize == NSpinOrbit_per_atom / spinsize
                    for (int m=0; m< NSpinOrbit_per_atom; m++) {

                        int n0 =   (CorrToHartr(cl,n) - cl*BlockDim  ) /spinsize  ;
                        int m0 =   (CorrToHartr(cl,m) - cl*BlockDim  ) /spinsize  ;

                        int n1 = n /spinsize  ;
                        int m1 = m /spinsize  ;

                        SolverBasis_cl.at(cl)(spinsize*n0,spinsize*m0) =  ces.eigenvectors()(n1,m1);
                        if (spinsize==2)           SolverBasis_cl.at(cl)(spinsize*n0+1,spinsize*m0+1) =  ces.eigenvectors()(n1,m1);
                    }
                }
            }

            ifroot {
                std::cout << SolverBasis_cl[cl] <<"\n";
            }
            ifroot std::cout <<"orgMat =\n" <<(Heff)  <<"\n";
            ifroot std::cout <<"rotMat =\n" <<( ces.eigenvectors().adjoint()* Heff * ces.eigenvectors() ).diagonal()  <<"\n";
        }//cl
        for(int cl=0 ; cl < NumMatBlock; cl++) {
            SolverBasis.block(cl*BlockDim, cl*BlockDim, BlockDim, BlockDim) = SolverBasis_cl[cl];
        }
    }

        if(mpi_rank ==0) {
            std::cout << "FILEOUT:SolverBasis.dat\n";
            FILE * OBSERV=fopen("SolverBasis.dat","w");
            for (int orb1=0; orb1<N_peratom_HartrOrbit*NumCorrAtom; orb1++) {
                for (int orb2=0; orb2<N_peratom_HartrOrbit*NumCorrAtom; orb2++) {
                    fprintf(OBSERV, "     %0.6f %0.6f",real(SolverBasis(orb1,orb2)), imag(NumMatrix(orb1,orb2)));
                }
                fprintf(OBSERV, "\n");
            }
            fclose(OBSERV);
        }//mpi_rank

    if( impurityBasisSwitch > 10 and  SOLVERtype==std::string("TB")) {
        for(int k=0; k<knum; k++) {
            Eigen::MatrixXcd  temp;
            temp.setIdentity(NBAND[k], NBAND[k]);
            temp.block(0,0, NumHartrOrbit,NumHartrOrbit) = SolverBasis;

            Eigen::MatrixXcd  Hk_sub_d   =  temp.adjoint() * H_k_inModelSpace[k] * temp;
            H_k_inModelSpace[k]= Hk_sub_d;

            Eigen::MatrixXcd  DF_sub_d   =  temp.adjoint() * DF_CorrBase[k];
            DF_CorrBase[k] = DF_sub_d;
        }
        SolverBasis.setIdentity(NumCorrAtom*N_peratom_HartrOrbit, NumCorrAtom*N_peratom_HartrOrbit);
    }
}






double FromHkToNele (double muDFT, Eigen::VectorXd  * KS_eigenEnergy) {

    double TNumEle=0;
    double TNumEle_local=0;
    for(int k=0; k<knum; k++) {
        for(int band=0; band<NumOrbit; band++) {
            TNumEle_local +=  1./(1+std::exp( beta_smearing*(KS_eigenEnergy[k][band]-muDFT) ));
        }
    }
    MPI_Allreduce(&(TNumEle_local), &(TNumEle), 1, MPI_DOUBLE, MPI_SUM,  MPI_COMM_WORLD);
    TNumEle/=knum_mpiGlobal;
    return TNumEle;
}
double Nele_non_Inter(
    int knum, int knum_mpiGlobal,Eigen::VectorXd  * KS_eigenEnergy
//    std::vector<Eigen::MatrixXcd> & H_k_inModelSpace,  std::vector<Eigen::MatrixXcd> & S_overlap
)
{
//find chemical potentiol for non-interacting band

//    std::vector<Eigen::VectorXd> KS_eigenEnergy(knum);
//
//    for(int k=0; k< knum; k++) {
//        Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXcd> ces1(NumOrbit);
//        ces1.compute( H_k_inModelSpace[k], S_overlap[k] );       //  HS \psi = S  \psi E
//        KS_eigenEnergy[k] = ces1.eigenvalues();
//    }

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
