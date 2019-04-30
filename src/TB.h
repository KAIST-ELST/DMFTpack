#ifndef TB_INCLUDE
#define TB_INCLUDE
//#include <complex>
#include <stdlib.h>
//#include <Eigen/Dense>
//using namespace std;
extern int  NumSubSpace, NumOrbit;
//extern int  SOCCal,  interOrbitHop;

//extern std::vector<Eigen::MatrixXcd>  H_k_inModelSpace;
extern std::vector<Eigen::MatrixXcd>  DF_CorrBase;
extern Eigen::MatrixXcd Sw_Hartree;
extern int NumLattx, NumLatty, NumLattz, NumLattxy, NumLattxyz;
//extern std::vector<Eigen::MatrixXcd> densityMatDFT;            
            
extern Eigen::MatrixXcd  weightMatrix;

//extern double  **KS_eigenEnergy;
extern Eigen::VectorXd *KS_eigenEnergy;

int isSameAtom(int , int, int model = 1);
//int nonNeg(int & a);

extern int knum, knum_mpiGlobal, myksta, mykend , kpxy;
//extern int impurityBasisSwitch;


//TB.cpp
extern std::vector<int> NBAND;
//extern int **FromValToKS;
extern std::vector<std::vector<int> > FromValToKS;
//extern Eigen::MatrixXi FromValToKS;
extern double **kmesh;
extern double *kdist_band;
extern Eigen::MatrixXi H_Rindex;




//void Construct_hyb_delta(int impurityDim, std::vector<int> impurityOrbit,
//                         ImgFreqFtn & SelfE_w, std::vector<Eigen::MatrixXcd>   Gw, double mu,
//                         ImgFreqFtn &  weiss_fieldTB, int atom ,     std::vector<Eigen::MatrixXcd> & SolverBasis  ) ;


void Construct_hyb_delta(int impurityDim, std::vector<int> impurityOrbit,
                         ImgFreqFtn & SelfE_w, std::vector<Eigen::MatrixXcd>   Gw, double mu,
                         ImgFreqFtn &  weiss_fieldTB, int atom,     Eigen::MatrixXcd & SolverBasis  ) ;


//void upfolding_k(Eigen::MatrixXcd & Matrix_k, Eigen::MatrixXcd & KS_eigenVectors_k, int k, double mu , std::vector<Eigen::MatrixXcd> & S_overlap);

void upfolding_density(std::vector<Eigen::MatrixXcd> &densityMatDFT, std::vector<Eigen::MatrixXcd> & KS_eigenVectors_orthoBasis ,Eigen::MatrixXi H_Rindex, double mu ,
                       std::vector<Eigen::MatrixXcd> & S_overlap,
                       std::vector<Eigen::MatrixXcd> & transformMatrix_k) ;



double GreenFtn_w_adjust_mu(std::vector<Eigen::MatrixXcd> &  H_k_inModelSpace, ImgFreqFtn & SelfE_w,                 double mu,  double dmu, int mu_adjust) ;
void GreenFtn_w( int NumCorrAtom, int NSpinOrbit_per_atom,   std::vector<Eigen::MatrixXcd> &  H_k_inModelSpace, ImgFreqFtn & SelfE_w, std::vector<Eigen::MatrixXcd>  & Gw, double mu,
                std::vector<Eigen::MatrixXcd> & densityMatDFT
               ) ;
//void retarded_GreenFtn( cmplx *retGkw , std::vector<Eigen::MatrixXcd> &   H_k_inModelSpace,  ImgFreqFtn & SE, double mu, int k, int n, int NumSubSpace) ;
void retarded_GreenFtn2( Eigen::MatrixXcd &retGkw_full , Eigen::MatrixXcd & retGkw,  std::vector<Eigen::MatrixXcd> & H_k_inModelSpace, ImgFreqFtn & SE, double mu, int k, int n) ;


void low_energy_subspace_in_KS_basis(
    int knum,  int knum_mpiGlobal,    
    std::vector<int> & NBAND,  std::vector<std::vector<int> >  & FromValToKS,  double muDFT,
    std::vector<Eigen::MatrixXcd> & KS_eigenVectors_orthoBasis, Eigen::VectorXd  * KS_eigenEnergy
) ;
 

#endif


