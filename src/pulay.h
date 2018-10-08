/*
Pulay Mixing Informations
optimalDensity_{n+1} = \sum_k=0^{s-1}  \beta_k  \inputD_{n-k};
where s=mixing_history and \beta was determined by minimizing
optimalResidue_{n+1} = \sum_k=0^{s-1}  \beta_k  \Residue_{n-k};
Here \beta can be negative or grater than 1. and residual vector = output - input;

Then \inputDensity_{n+1} = optimalDensity_{n+1} + \alpha  * (optimalResidue_{n+1}),
\alpha \in [0,1];
*/
#ifndef pulay_INCLUDE
#define pulay_INCLUDE
#include "ImgTimeFtn.h"
#include "dmft_common.h"
#include <Eigen/Core>


void mixing_checker(double resid, double resid_prev, double & mixing,  double mixing_min,  double mixing_max);

class pulayMixing {
public:
    pulayMixing (int mixing_history_, int start_mixing_,int dim_i, int dim_j, int dim_k, bool parallel_pulay_ = false );
    ~pulayMixing ();
    void mixing(Eigen::MatrixXcd  * inputDensity_n, Eigen::MatrixXcd * outputDensity_n, double  mixingSCGF, int SCGFloop, int mixingStep) ;
//    void my_mixing(Eigen::MatrixXcd  * inputDensity_n, Eigen::MatrixXcd * outputDensity_n, double  mixingSCGF, int SCGFloop, int mixingStep) ;
    void mixing(Eigen::VectorXd   * inputDensity_n, Eigen::VectorXd * outputDensity_n, double  mixingPulay, int SCGFloop, int mixingStep) ;
private:
    int mixing_history;
    int start_mixing;
    int dim_i, dim_j, dim_k;
    Eigen::MatrixXcd ** input_history;
    Eigen::MatrixXcd ** res_history;
    Eigen::MatrixXd OptimalMatrix;
    bool parallel_pulay;
};

#endif
