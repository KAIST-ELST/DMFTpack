//How to use : ./a.out [TB hopping parameter]  [SelfEnergy used in imaginary freq calculation, default=0]
#ifndef tight_common_INCLUDE
#define tight_common_INCLUDE
#include <unistd.h>
#include <assert.h>
#include <iomanip>
#include <math.h>
#include "ImgTimeFtn.h"
#include "mpi.h"
#include "dmft_common.h"
#include <Eigen/Core>


//using namespace std;


#ifdef  LAPACK_MKL
extern "C"  void zgetrf_( int*, int* , std::complex<double>* , int*, int* , int* );
extern "C"  void zgetri_( int*, std::complex<double>* , int*, int* , std::complex<double>*, int* , int* );

extern "C"  void zgeev( char* jobvl, char* jobvr, int* n, std::complex<double>* a,
                        int* lda, std::complex<double>* w, std::complex<double>* vl, int* ldvl, std::complex<double>* vr, int* ldvr,
                        std::complex<double>* work, int* lwork, double* rwork, int* info );
#endif




//DFT+DMFT
extern double muTB;
extern double muDFT;
//extern Eigen::MatrixXcd * NumMatDFT;
//extern std::vector<Eigen::MatrixXcd> NumMatDFT;   //KS_eigenVectors[k][i][n]   ith component of eigenstate with energy (E_nk),  direct rep. in non-orthogonal basis,i.e.,  <i|kn>

//downfolding
//extern double  low_energy_window;
extern double  lower_model_window, upper_model_window;
extern double   lower_spectrum_window, upper_spectrum_window;
extern int downfolding;

//
extern int   *Hart2Corr;
extern Eigen::MatrixXi CorrToHartr;
extern int  *HartrIndex_inDFT;


extern int * isOrbitalHartrDFT;
extern int * isOrbitalCorrinHart;
extern int  num_subshell;
extern Eigen::VectorXi  subshell;
extern Eigen::VectorXi  Rydberg_set;
extern Eigen::VectorXi  rot_sym;
extern int  *KS2Hartr;
//extern int  *LongRangeOrder;
//extern int  *LongRangeOrder_DFT;
//extern int  *LongRangeOrder_Hart;
extern int  DFTIt;
extern int  *FromOrbitalToAtom; //0,1,...
extern int  *FromOrbitalToLocalOrbital_DFT; //0,1,...
extern int  *FromOrbitalToAtom_model; //0,1,...
extern int  k_pointx, k_pointy, k_pointz,  k_grid, mu_adjust, SpinP_switch, NumAtom, BraLatt_x,BraLatt_y,BraLatt_z,Measure_nn, NumCorrAtom, NumCluster,NumHartrOrbit;
extern int NumAtom_per_cluster, NumCluster, NumHartrOrbit_per_cluster;
extern int H0_from_OpenMX;
extern double beta, NumberOfElectron, EnergyUnit, doublecounting, Band_renorm,  TotalCorrEle_in_Imp;
extern std::string dctype, localOrbitalType, SOLVERtype, Lowlevel_SOLVERtype;
extern Eigen::Matrix2cd * Zeeman_field_spin;
extern int highFreq, N_freq, N_tau,  magnetism;
extern unsigned long long  Num_MC_steps, THERMALIZATION;
extern std::string mode;
extern double infinitesimal ;
extern int  maxDmftIt, maxDFTIt, maxTime,restart, NSpinOrbit_per_atom, N_peratom_HartrOrbit;
extern double UHubb, Uprime, JHund ,mixing;
extern double nominal_charge;
extern int  mixingFtn;
//extern std::string LongRangeOrder_file;
extern std::string system_name;
extern std::string SOLVERexe, SOLVERdir;
extern cmplx w0, dw;
//extern int SOLVERtype;
extern int  SOCCal,  interOrbitHop;
extern int impurityBasisSwitch;

extern int **HartrRange;
//extern int **HartrRange_DFT;

extern int NsBath;

//Local var
extern Eigen::MatrixXcd impurity_site_Hamiltonian;

extern Eigen::MatrixXcd Sw_doublecounting, NumMatrix;

//extern double  **UMatrix;


extern std::vector<Eigen::VectorXi> Uindex;
extern std::vector<cmplx > Utensor;
extern std::vector<Eigen::VectorXi> Uindex_stronglyCorr;
extern std::vector<cmplx > Utensor_stronglyCorr;

extern cmplx dE, E0;
extern int  Spectral_EnergyGrid ;
extern MPI_Comm localComm;
extern char   processor_name[100];

//Simple Tetragonal Lattice
extern int Nkpath;
extern  Eigen::Vector3d RecipUnitVector_b1;
extern  Eigen::Vector3d RecipUnitVector_b2;
extern  Eigen::Vector3d RecipUnitVector_b3;
extern double KpathPoint[3*30];
//static int const  Nkpath = 5;
//static double const KpathPoint[Nkpath+1][dimension_3d] = {
//    {0.0,  0.0,  0.0},   //G
//    {0.5,  0.0,  0.0},   //X
//    {0.5,  0.5,  0.0},   //M
//    {0.0,  0.0,  0.0},   //G
//    {0.0,  0.0,  0.5},   //Z
//    {0.5,  0.5,  0.5}   //A
//};      //unit = 2pi/a with  reciprocal lattice (b1, b2 basis.)







//void inverse( cmplx *A, int dim);   // eigen library
//void inverse( cmplx **A, int dim );
//void inverse( cmplx **A, cmplx **invA, int dim ) ;
//void inverse( double **A, cmplx **invA, int dim ) ;
//void inverse(cmplx *Greenftn, int n, int NumSubSpace );
//void inverse( cmplx *A, int dim, int * permutations, cmplx * WORK ) ; //lapack



int  diag_jacobi (cmplx **a_hermitian, int dim, double *eigenVal, cmplx **v_comp);
void band(Eigen::VectorXd *KS_eigenEnergy, double muDFT, int knum) ;
void band( std::vector<Eigen::MatrixXcd> & H_k_inModelSpace, double muTB, int knum) ;
void getAsymto_moments (std::vector<Eigen::MatrixXcd> & moments, Eigen::MatrixXcd * delta_w , bool greenFtnAsym = false) ;
void FourierTransform (Eigen::MatrixXcd * delta_w, Eigen::MatrixXcd * delta_t                                       , bool greenFtnAsym   = false, std::string FT_ftn=std::string("NOPRINT"), int print=0);
void FourierTransform (Eigen::MatrixXcd * delta_w, Eigen::MatrixXcd * delta_t , Eigen::MatrixXcd AsymtoV            , bool greenFtnAsym   = false, std::string FT_ftn=std::string("NOPRINT"), int print=0);
void FourierTransform (Eigen::MatrixXcd * delta_w, Eigen::MatrixXcd & delta_t , double tau                          , bool greenFtnAsym   = false, std::string FT_ftn=std::string("NOPRINT"), int print=0);
void FourierTransform (Eigen::MatrixXcd * delta_w, Eigen::MatrixXcd & delta_t , double tau, Eigen::MatrixXcd AsymtoV, bool greenFtnAsym   = false, std::string FT_ftn=std::string("NOPRINT"), int print=0);
void FT_t_to_w (Eigen::MatrixXcd * delta_w, Eigen::MatrixXcd * delta_t, int N_freq );

void dos(Eigen::VectorXd * KS_eigenEnergy,std::vector<Eigen::MatrixXcd> & KS_eigenVectors_orthoBasis , double E_window, double & muDFT ) ;


double read_double(const std::string &Inputfile, const std::string &keyword, bool defalutBool, double dft);
double read_double(const std::string &Inputfile1,const std::string &Inputfile2, const std::string &keyword, double dft) ;
int read_int(const std::string &Inputfile, const std::string &keyword, int dft);
int read_int(const std::string &Inputfile1,const std::string &Inputfile2, const std::string &keyword, int dft) ;
void read_int_array(const std::string &Inputfile, const std::string &keyword,    std::vector<int>  & var,  int length, bool defalutBool, int dft);
void read_double_array(const std::string &Inputfile, const std::string &keyword, std::vector<double>  & var,  int length, bool required, double dft=0);
//std::string read_string(const std::string &Inputfile, const std::string &keyword) ;
std::string read_string(const std::string &Inputfile, const std::string &keyword, bool defalutBool= false, std::string defalutVal = std::string("no def") ) ;
void para_range(int n1, int n2, int nprocs, int myrank, int *mystart, int * myend);
void data_sync(cmplx *A , int startRow, int endRow, int nprocs);
void data_sync(double *A , int startRow, int endRow,  int nprocs) ;
void data_sync(int  *A , int startRow, int endRow,  int nprocs) ;
void data_sync_EigenMat(Eigen::MatrixXcd *A , int startRow, int endRow, int matrix_dim,  int nprocs) ;
void data_sync(cmplx **A , int startRow, int endRow, int lenColumn, int nprocs);
void data_sync(double **A , int startRow, int endRow, int lenColumn,  int nprocs) ;
void data_sync(int **A , int startRow, int endRow, int lenColumn, int nprocs);
void data_sync(cmplx ***A , int startRow, int endRow, int lenColumn, int lenhight,  int nprocs) ;

void read_inputFile( const std::string &hamiltonian) ;
void on_the_fly_control() ;

void cmplxEigenproblem(cmplx **a, int dim, double *eigenVal, cmplx **v_comp);
void cmplxEigenproblem(cmplx *a, int dim, double *eigenVal, cmplx *v_comp);
void cmplxEigenproblem(cmplx **a, int dim, double *eigenVal);
void cmplxEigenproblem(cmplx *a, int dim, double *eigenVal);
double GcmplxEigenproblem(cmplx *a, cmplx *b,  int dim, double *eigenVal);
double GcmplxEigenproblem(cmplx **a,cmplx *b, int dim, double *eigenVal, cmplx **eigenvec);

void optical_conductivity(Eigen::MatrixXcd *H_k_inModelSpace , ImgFreqFtn Sw, double mu, int knum, cmplx ***** H_R, double **kmesh,int alp_direc, int beta_direc) ;
void write_hyb_t ( ImgFreqFtn & weiss_field, Eigen::MatrixXcd * delta_t, int atom) ;
cmplx operator*( int a, cmplx b);
cmplx operator*( cmplx a, int b);
cmplx operator/( cmplx  a, int b);
bool operator!=( cmplx  a, int b);
bool operator==( cmplx  a, int b);
//void ImS_stoch_to_ReS( ImgFreqFtn & SelfEnergy_w);
void rewrite_retardedSw() ;
void free_energy_cal( cmplx **Sw, double mu) ;

//int analytic_selfenergy(int magnetism, const std::string &ImgFtn_data );
//int analytic_stoch(int magnetism, const std::string &ImgFtn_data);



void gen_Uijkl(int n_spinorb, double U, double Uprime, double JH, std::vector<cmplx >  & Utensor, std::vector<Eigen::VectorXi>  & Uindex ) ;
void gen_Uijkl_density_density(int n_spinorb, double U, double Uprime, double JH, std::vector<cmplx > & Utensor, std::vector<Eigen::VectorXi>  & Uindex ) ;
//int read_hdf5 (int index);

void Wait_Run( char fileName[100], int checkTime, int mpi_rank ,int maxTime  ) ;




#endif

