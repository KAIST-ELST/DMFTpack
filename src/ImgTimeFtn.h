#ifndef ImgTimeFtn_INCLUDE
#define ImgTimeFtn_INCLUDE
#include <stdio.h>
#include <complex>
#include <iostream>
#include "mpi.h"
#include <Eigen/Core>
#include <vector>

extern int mixingType, mixingStart;

class ImgTimeFtn {
public:
    ImgTimeFtn(const std::string &filename, int N_tau,  int NSpinOrbit);
    ImgTimeFtn(double beta, int N_tau,  int NSpinOrbit);
    ~ImgTimeFtn();
    double getValue(int time_indx, int orbit);
    double getValue(int time_indx);
    int getNtau();
    int getNorbit();
    void setValue(int time_indx, int orbit, double value);
    void update(const std::string &filename, double mixing);
    void  dataOut(const std::string &filename, int PrintAtomIdx);

private:
    double  **Ftn;
    double *imgTime;
    int Ntau, Norbit;
};


class ImgFreqFtn {
public:
    int Nfreq_, Norbit, NumSite_, NSpinOrbitPerAtom_;

    Eigen::MatrixXcd Asymto0;
    Eigen::MatrixXcd Asymto1;
    Eigen::MatrixXcd Asymto2;
    Eigen::MatrixXcd Asymto3;
    Eigen::MatrixXcd Asymto4;

    ImgFreqFtn(double beta, int N_freq,  int NSpinOrbit, int NumSite, int mixingType);
    ImgFreqFtn(int a) ;
    ImgFreqFtn(ImgFreqFtn & rhs);
    void Initialize(                                       double beta, int N_freq,  int NSpinOrbitPerAtom, int NumSite,    int mixingType, int startIndx=0, int distMem=0) ;
    void realFreq  (std::complex<double>  E0,std::complex<double> dE,     int N_freq,  int NSpinOrbitPerAtom, int NumSite,    int mixingType, int startIndx=0) ;
    ~ImgFreqFtn();



    int  getNFreq();
    int  getNOrbital();
    double                 getValue(int freq_indx);
    std::complex<double>   getValue(int freq_indx, int orbit1, int orbit2);
    std::complex<double>   getValueSubMat(int freq_indx, int site,  int orbit1, int orbit2);
    Eigen::MatrixXcd       getMatrix(int w) ;
    Eigen::MatrixXcd       getMatrix(int w, int site, int dim) ;
std::vector<Eigen::MatrixXcd> getFtn_data();
    void    setMatrix(int w,            Eigen::MatrixXcd value ) ;
    void    setMatrix(int w, int site,  Eigen::MatrixXcd value ) ;
    void    setValue      (int w,                   int orbit1, int orbit2,  std::complex<double> value);
    void    setValueSubMat(int w, int site, int orbit1, int orbit2,  std::complex<double> value);
    void    setValueSubMat(int w, int site, int dim, int orbit1, int orbit2,  std::complex<double> value);

    void    read_diag  (const std::string &filename );
    void    read_uppertrian  (const std::string &filename, int spindim=0 );
    void    update  (ImgFreqFtn & rhs, double mixing) ;
    void    update  (ImgFreqFtn & FtnOut_c,           double mixing, int updateSite,  int mixingType  ) ;
    void    update  (Eigen::MatrixXcd * FtnOutM,      double mixing,                  int mixingType  ) ;
    void    update_full  (const std::string &filename, double mixing);
void read_full(const std::string &filename,double beta, double beta_prev) ;

//void    generalized_broyden_mixing( std::complex<double> ** FtnOut,               int mixingType=0,   int parameters=0) ;
    void    dataOut (const std::string &filename);
//    void    dataOut (const std::string &filename,  int rank=0);
    void    dataOut_full (const std::string &filename);
    void    dataOut_full_pararell(const std::string &filename) ;
    void    mpiBcast(int root, MPI_Comm comm);
void estimate_asymto(int order);

private:
    std::vector<Eigen::MatrixXcd> Ftn_;
    double *imgFreq;
    int *MEMCHECK;
    int mixingType_;
    int startIndx_;

    /*broyden_mixing options*/
//int ftncount;
//int speed_up_need;
//double speed_now    ;
//double mixingFactor;
//int    speed_up     ;
//int    speed_up_crit;
//double  ** DVm, **DFm;
//double  * Vm, *Fm;
//double  ** B;
};
#endif
