#include <fstream>
#include "ImgTimeFtn.h"
#include <unistd.h>
#include <assert.h>
#include "dmft_common.h"
//#include <vector>
#include <array>



double const  tol = 1e-7;

ImgTimeFtn::ImgTimeFtn(const std::string &filename, int N_tau,  int NSpinOrbit) {
    char next;
    int t,n;
    std::ifstream DATA(filename.c_str());
    Ftn = new double  * [N_tau+1];
    imgTime = new double [N_tau+1];
    for(t=0; t<N_tau+1; t++) {
        Ftn[t] = new double [NSpinOrbit];
    }
    for(t=0; t<N_tau+1; t++) {
        DATA >> imgTime[t];
        for(n=0; n < NSpinOrbit; n++) {
            DATA >> Ftn[t][n];
        }
        while(DATA.get(next)) {
            if (next=='\n') break;
        }
    }
    Ntau=N_tau;
    Norbit=NSpinOrbit;
    DATA.close();
}
ImgTimeFtn::ImgTimeFtn(double beta, int N_tau,  int NSpinOrbit) {
    int t,n;
    Ftn = new double  * [N_tau+1];
    imgTime = new double [N_tau+1];
    for(t=0; t<N_tau+1; t++) {
        Ftn[t] = new double [NSpinOrbit];
    }

    for(t=0; t<N_tau+1; t++) {
        imgTime[t] = (beta/(N_tau*1.0))*t;
        for(n=0; n < NSpinOrbit; n++) {
            Ftn[t][n]=0;
        }
    }
    Ntau=N_tau;
    Norbit=NSpinOrbit;
}
ImgTimeFtn::~ImgTimeFtn() {
    int t;
    for(t=0; t<Ntau+1; t++) {
        delete [] Ftn[t];
    }
    delete [] Ftn;
    delete [] imgTime;
}

double ImgTimeFtn::getValue(int t, int n) {
    return Ftn[t][n];
}
double ImgTimeFtn::getValue(int t) {
    return imgTime[t];
}
int ImgTimeFtn::getNtau() {
    return Ntau;
}
int ImgTimeFtn::getNorbit() {
    return Norbit;
}



void ImgTimeFtn::setValue(int t, int n, double value) {
    Ftn[t][n]= value;
}




void ImgTimeFtn::update(const std::string &filename, double mixing) {
    FILE * Fcheck;
    int fileOpen =0;
    if(mpi_rank==0) {
        for(int i=0; i<20; i++) {
            Fcheck = fopen(filename.c_str(),"r");
            if (Fcheck==NULL and fileOpen<6) {
                fileOpen++;
                sleep(3);
            }
            else if (fileOpen>=6) {
                std::cout << "fail to open" << filename.c_str() <<"\n";
                assert(0);
            }
            else {
                sleep(5);
                fclose(Fcheck);
                break;
            }
        }
        Fcheck = fopen(filename.c_str(),"r");
        assert( Fcheck != NULL );
        fclose(Fcheck);
        sleep(5);

        int readtry=0, check=0;
        do {
            double temp;
            std::ifstream DATA(filename.c_str());
            for(int t=0; t<Ntau+1; t++) {
                DATA>> temp;
                if(fabs(imgTime[t] - temp) > 10e-4) {
                    check=1 ;
                    ifroot    std::cout << filename << " : " <<t <<" " << imgTime[t] << " not eq. to " <<temp <<"\n";
                    exit(1);
                }
                for(int n=0; n < Norbit; n++) {
                    DATA >> temp;
                    Ftn[t][n]=(1-mixing)*Ftn[t][n] + mixing*temp;
                }
            }
            if(check==1) std::cout<< filename << " : Error Img.time Ftn update\n" ;
            DATA.close();
        } while (check==1 and readtry<5)  ;
        assert(check==0);
    }//ifroot
    for(int t=0; t<Ntau+1; t++) {
        MPI_Bcast(&(Ftn[t][0]), Norbit, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
}

void ImgTimeFtn::dataOut(const std::string &filename, int N_ATOMS) {
    if(mpi_rank==0) {
        FILE *datap4 = fopen(filename.c_str(), "w");
        int t,n;
        for(t=0; t<Ntau+1; t++) {
            fprintf(datap4, "%0.5f", imgTime[t]);
            for(n=0; n< Norbit; n++) {
                fprintf(datap4, "     %0.10f", Ftn[t][n]);
            }
            fprintf(datap4,"\n");
        }

        fclose(datap4);
    }
    sleep(3);
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

ImgFreqFtn::ImgFreqFtn(double beta, int N_freq,  int NSpinOrbitPerAtom,  int NumSite, int mixingType) {
    (*this).Initialize(beta,     N_freq,      NSpinOrbitPerAtom,      NumSite,     mixingType);
}


ImgFreqFtn::ImgFreqFtn(int a) {
    imgFreq  = NULL;
    Nfreq_   =  0;
    Norbit  =  0;
    MEMCHECK =  new int;
    *MEMCHECK = magicNumber;
    mixingType_ = 0;
}


void ImgFreqFtn::Initialize(double beta, int N_freq,  int NSpinOrbitPerAtom,int NumSite, int mixingType, int startIndx, int distMem) {
    Nfreq_=N_freq;
    Norbit=NSpinOrbitPerAtom * NumSite;
    NSpinOrbitPerAtom_ = NSpinOrbitPerAtom;
    NumSite_ = NumSite;
    mixingType_ = mixingType;  //should be zero in this version
    startIndx_ = startIndx;

    imgFreq = new double [Nfreq_];
    Ftn_.setZero(Nfreq_+5, (Norbit)* (Norbit));
    for(int w=0; w<Nfreq_; w++) {
        imgFreq[w] = pi*(2.*(w+startIndx_)+1)/beta;
    }


    MEMCHECK =  new int;
    *MEMCHECK = magicNumber;
}

void ImgFreqFtn::realFreq(cmplx E0, cmplx dE, int NumE,  int NSpinOrbitPerAtom, int NumSite, int mixingType,  int startIndx) {
    (*this).Initialize(1., NumE, NSpinOrbitPerAtom, NumSite,mixingType, startIndx);
    for(int w=0; w<NumE; w++) {
        int IDX = w+startIndx;   //IDX = 0,1,...Spectral_EnergyGrid-1
        imgFreq[w] =  real(E0+dE*((double) IDX)) ;
    }
}

ImgFreqFtn::ImgFreqFtn(ImgFreqFtn & rhs) {
    (*this).Initialize(1., rhs.Nfreq_, rhs.Norbit, rhs.mixingType_, rhs.startIndx_);
    for(int w=0; w<rhs.Nfreq_; w++) {
        imgFreq[w] = rhs.imgFreq[w];
        for (int i=0; i<rhs.Norbit; i++) {
            for (int j=0; j<rhs.Norbit; j++) {
                Ftn_(w,i*rhs.Norbit+j) = rhs.Ftn_(w,i*rhs.Norbit+j);
            }
        }
    }
    assert (*(rhs.MEMCHECK) == magicNumber);
    assert (*MEMCHECK == magicNumber);
}

ImgFreqFtn::~ImgFreqFtn() {
    delete [] imgFreq;
    assert (*MEMCHECK == magicNumber);
    delete MEMCHECK;
}


int  ImgFreqFtn::getNFreq() {
    return Nfreq_;
}
int  ImgFreqFtn::getNOrbital() {
    return Norbit;
}
cmplx ImgFreqFtn::getValue(int w, int n, int m) {
    assert(0<= w-startIndx_ );
    assert(w-startIndx_ < Nfreq_+5);
    assert(0<= n );
    assert(0<= m );
    assert(n < Norbit);
    assert(m < Norbit);
    return Ftn_(w-startIndx_,n*Norbit+m);
}
double ImgFreqFtn::getValue( int w) {
    return imgFreq[w-startIndx_];
}
Eigen::MatrixXcd ImgFreqFtn::getMatrix(int w) {
    Eigen::MatrixXcd result (Norbit, Norbit);
    for (int n =0; n<Norbit; n++) {
        for (int m =0; m<Norbit; m++) {
            result(n,m) = Ftn_(w, n*Norbit + m);
        }
    }
    return result;
}

Eigen::MatrixXcd ImgFreqFtn::getMatrix(int w, int site) {
    Eigen::MatrixXcd result (NSpinOrbitPerAtom_, NSpinOrbitPerAtom_);
    for (int n =0; n<NSpinOrbitPerAtom_; n++) {
        for (int m =0; m<NSpinOrbitPerAtom_; m++) {
            int n1 = site *  NSpinOrbitPerAtom_ + n;
            int m1 = site *  NSpinOrbitPerAtom_ + m;
            result(n,m) = Ftn_(w, n1*Norbit + m1);
        }
    }
    return result;
}


void ImgFreqFtn::setMatrix(int w, int site,  Eigen::MatrixXcd value ) {
    if(NSpinOrbitPerAtom_ != value.rows())
        exit(1);

    for (int n =0; n<NSpinOrbitPerAtom_; n++) {
        for (int m =0; m<NSpinOrbitPerAtom_; m++) {
            int n1 = site *  NSpinOrbitPerAtom_ + n;
            int m1 = site *  NSpinOrbitPerAtom_ + m;
            Ftn_(w, n1*Norbit + m1) = value(n,m);
        }
    }
}



void ImgFreqFtn::setMatrix(int w, Eigen::MatrixXcd value ) {
    if(Norbit != value.rows())
        exit(1);

    for (int n =0; n<Norbit; n++) {
        for (int m =0; m<Norbit; m++) {
            Ftn_(w, n*Norbit + m) = value(n,m);
        }
    }
}


void ImgFreqFtn::setValue(int w, int n, int m, cmplx value) {
    Ftn_(w-startIndx_,n*Norbit+m)= value;
}

void  ImgFreqFtn::dataOut(const std::string &filename) {
    if(mpi_rank==0) {
        FILE *datap4 = fopen(filename.c_str(), "w");
        int w,n;
        for(w=0; w<Nfreq_; w++) {
            fprintf(datap4, "%0.8f", imgFreq[w]);
            for(n=0; n< Norbit; n++) {
                fprintf(datap4, "     %0.10f  %0.10f", real(Ftn_(w,n*Norbit+n)), imag(Ftn_(w,n*Norbit+n)));
            }
            fprintf(datap4,"\n");
        }
        //fprintf(datap4,"\n");
        //for(w=Nfreq_; w<Nfreq_+5; w++) {
        //    fprintf(datap4, "%d", w-Nfreq_);
        //    for(n=0; n< Norbit; n++) {
        //        fprintf(datap4, "     %0.10f  %0.10f", real(Ftn_(w,n*Norbit+n)), imag(Ftn_(w,n*Norbit+n)));
        //    }
        //    fprintf(datap4,"\n");
        //}
        fclose(datap4);
        std::cout << "FILEOUT_:" <<  filename <<"\n";
    }
    assert (*MEMCHECK == magicNumber);
}

void  ImgFreqFtn::dataOut_full(const std::string &filename) {
    if(mpi_rank==0) {
        FILE *datap4 = fopen((filename).c_str()  , "w");
        int w,n,m;
        for(w=0; w<Nfreq_; w++) {
            for(m=0; m< Norbit; m++) {
                for(n=0; n< Norbit; n++) {
                    if( std::abs(Ftn_(w,m*Norbit+n)) >tol)  fprintf(datap4, "%d %d %d %+0.8f %+0.8f\n",w,m,n,real(Ftn_(w,m*Norbit+n)), imag(Ftn_(w,m*Norbit+n)));
                }//n
            }//m
        }//w
        for(m=0; m< Norbit; m++) {
            for(n=0; n< Norbit; n++) {
                fprintf(datap4, "-1 %d %d %+0.8f %+0.8f\n",m,n,real(Ftn_(Nfreq_+0,m*Norbit+n)), imag(Ftn_(Nfreq_+0,m*Norbit+n)));
                fprintf(datap4, "-2 %d %d %+0.8f %+0.8f\n",m,n,real(Ftn_(Nfreq_+1,m*Norbit+n)), imag(Ftn_(Nfreq_+1,m*Norbit+n)));
                fprintf(datap4, "-3 %d %d %+0.8f %+0.8f\n",m,n,real(Ftn_(Nfreq_+2,m*Norbit+n)), imag(Ftn_(Nfreq_+2,m*Norbit+n)));
                fprintf(datap4, "-4 %d %d %+0.8f %+0.8f\n",m,n,real(Ftn_(Nfreq_+3,m*Norbit+n)), imag(Ftn_(Nfreq_+3,m*Norbit+n)));
                fprintf(datap4, "-5 %d %d %+0.8f %+0.8f\n",m,n,real(Ftn_(Nfreq_+4,m*Norbit+n)), imag(Ftn_(Nfreq_+4,m*Norbit+n)));
            }//n
        }//m
        fclose(datap4);
        std::cout << "FILEOUT:" << filename.c_str() <<"\n";
    }
    assert (*MEMCHECK == magicNumber);
}
void  ImgFreqFtn::dataOut_full_pararell(const std::string &filename) {
    FILE *datap1;
    ifroot    datap1 = fopen( (filename).c_str()    , "w");
    ifroot    fclose(datap1);

    for(int itsRank=0 ; itsRank<mpi_numprocs; itsRank++) {
        if(mpi_rank==itsRank) {
            datap1 = fopen( (filename).c_str()   , "a");
            for (int w=0; w< Nfreq_; w++) {
                for(int m=0; m< Norbit; m++) {
                    for(int n=0; n< Norbit; n++) {
                        if( std::abs(Ftn_(w,m*Norbit+n)) >tol)  fprintf(datap1, "%d %d %d %+0.8f %+0.8f\n", w+startIndx_, m, n, real(Ftn_(w,m*Norbit+n)), imag(Ftn_(w,m*Norbit+n)));
                    }//n
                }//m
            }//w
            fclose(datap1);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        sleep(1);
    }//itsRank
}

//void ImgFreqFtn::update(const std::string &filename, double mixing, int mixingType) {
//    //file check
//    int check=0;
//    FILE * Fcheck;
//    int Nfileopen=0;
//    while (1) {
//        Nfileopen++;
//        Fcheck = fopen(filename.c_str(),"r");
//        if (Fcheck==NULL) {
//            sleep(1);
//        }
//        else {
//            sleep(1);
//            fclose(Fcheck);
//            break;
//
//        }
//        if (Nfileopen>3) {
//            std::cout << "Cannot open file : " << filename << "\n";
//            return;
//        }
//    }
//
//    Eigen::MatrixXcd  FtnOutM[Nfreq_];
//    for (int w=0; w < Nfreq_+5; w++ ) {
//        FtnOutM[w].setZero(Norbit,Norbit);
//    }
//    double Retemp[Norbit], Imtemp[Norbit], PreRetemp[Norbit], PreImtemp[Norbit], DataFreq, PreDataFreq=0, LinearParameter;
//    cmplx I(0,1);
//    //start to read file
//    char next;
//    int w=0;
//    std::ifstream DATA(filename.c_str());
//    while ( w < Nfreq_ and !DATA.eof() ) {
//        DATA>> DataFreq;
//        for(int n=0; n < Norbit; n++) {
//            DATA >> Retemp[n];
//            DATA >> Imtemp[n];
//        }
//        if(w==0) {
//            for(int n=0; n < Norbit; n++) {
//                PreRetemp[n]=Retemp[n];
//                PreImtemp[n]=Imtemp[n];
//            }
//        }
//
//        if( (imgFreq[w] - DataFreq) > std::min(1., imgFreq[1]-imgFreq[0] )   ) {
//            check =1;
//            if (imgFreq[w] < DataFreq) {
//                while ( imgFreq[w] < DataFreq  and w < Nfreq_ ) {
//                    LinearParameter=(imgFreq[w] - PreDataFreq) / (DataFreq-PreDataFreq);
//                    for(int n=0; n < Norbit; n++) {
//                        FtnOutM[w](n,n)=(  (LinearParameter*Retemp[n]+(1-LinearParameter)*PreRetemp[n])+
//                                           I*(LinearParameter*Imtemp[n]+(1-LinearParameter)*PreImtemp[n]) );
//                    }
//                    w++;
//                }//while
//            }
//            else if (imgFreq[w] > DataFreq) {
//                while ( imgFreq[w] > DataFreq  && !DATA.eof() ) {
//                    DATA>> DataFreq;
//                    for(int n=0; n < Norbit; n++) {
//                        DATA >> Retemp[n];
//                        DATA >> Imtemp[n];
//                        FtnOutM[w](n,n) = ( Retemp[n]+I*Imtemp[n]);
//                    }
//                }
//            }
//        }
//        else {
//            for(int n=0; n < Norbit; n++) {
//                FtnOutM[w](n,n) = ( Retemp[n]+I*Imtemp[n]);
//            }
//            w++;
//        }
//        PreDataFreq = DataFreq;
//        for(int n=0; n < Norbit; n++) {
//            PreRetemp[n]=Retemp[n];
//            PreImtemp[n]=Imtemp[n];
//        }
//        while(DATA.get(next)) {
//            if (next=='\n') break;
//        }
//    }//while
//    if (w < Nfreq_ and DATA.eof()) { //length of the file is too short
//        int wlast = w-1;
//        while (w < Nfreq_ ) {
//            for(int n=0; n < Norbit; n++) {
//                FtnOutM[w](n,n)= FtnOutM[wlast](n,n);
//            }
//            w++;
//        }//while
//    }
//    DATA.close();
//    assert (*MEMCHECK == magicNumber);
//    if(check==1 && mpi_rank==0 ) std::cout<<"#\nWarning2:"<< filename << " : :Img.time Ftn update\n#\n" ;
//    ImgFreqFtn::update(FtnOutM,mixing, mixingType);
//}
//
void ImgFreqFtn::read_uppertrian(const std::string &filename) {
    Eigen::MatrixXcd  * FtnOutM = new Eigen::MatrixXcd [Nfreq_+5];
    for (int w=0; w < Nfreq_+5; w++ ) {
        FtnOutM[w].setZero(Norbit,Norbit);
    }
    double Retemp, Imtemp, DataFreq ;
    cmplx I(0,1);
    //start to read file
    char next;
    int w=0;
    std::ifstream DATA(filename.c_str());
    while ( w < Nfreq_ and !DATA.eof() ) {
        DATA>> DataFreq;

        for(int n=0; n < Norbit; n++) {
            for(int m=n; m < Norbit; m++) {
                DATA >> Retemp;
                DATA >> Imtemp;
                FtnOutM[w](n,m) = ( Retemp+I*Imtemp);
                if(n<m) FtnOutM[w](m,n) = std::conj(FtnOutM[w](n,m));
            }
        }
        w++;
        while(DATA.get(next)) {
            if (next=='\n') break;
        }
    }//while
    DATA.close();
    ImgFreqFtn::update(FtnOutM,1,0);
    delete [] FtnOutM;
}



void ImgFreqFtn::read_diag(const std::string &filename) {
    //file check
    int check=0;
    FILE * Fcheck;
    int Nfileopen=0;

    Eigen::MatrixXcd  * FtnOutM = new Eigen::MatrixXcd [Nfreq_+5];
    for (int w=0; w < Nfreq_+5; w++ ) {
        FtnOutM[w].setZero(Norbit,Norbit);
    }
    double Retemp[Norbit], Imtemp[Norbit],  DataFreq ;
    cmplx I(0,1);
    //start to read file
    char next;
    int w=0;
    std::ifstream DATA(filename.c_str());
    while ( w < Nfreq_ and !DATA.eof() ) {
        DATA>> DataFreq;
        for(int n=0; n < Norbit; n++) {
            DATA >> Retemp[n];
            DATA >> Imtemp[n];
        }

        for(int n=0; n < Norbit; n++) {
            FtnOutM[w](n,n) = ( Retemp[n]+I*Imtemp[n]);
        }
        w++;
        while(DATA.get(next)) {
            if (next=='\n') break;
        }
    }//while
    DATA.close();
    ImgFreqFtn::update(FtnOutM,1,0);
    delete [] FtnOutM;
}



void ImgFreqFtn::update(ImgFreqFtn & FtnOut_c, double mixing, int updateSite, int mixingType) {
    assert(FtnOut_c.Nfreq_ == Nfreq_);
    assert(FtnOut_c.NSpinOrbitPerAtom_ == NSpinOrbitPerAtom_);
    assert(( (FtnOut_c.NumSite_) == NumSite_   and  updateSite==0   ) or
           ( (FtnOut_c.NumSite_) == 1                               ));
    assert ( updateSite < NumSite_);

    int dim = FtnOut_c.Norbit;
    Eigen::MatrixXcd * FtnOutM  = new Eigen::MatrixXcd [Nfreq_+5];
    for (int w=0; w < Nfreq_; w++ ) {
        FtnOutM[w].setZero(dim,dim);
        for(int n=0; n < dim; n++) {
            for(int m=0; m < dim; m++) {
                FtnOutM[w](n,m)=FtnOut_c.Ftn_(w,n*dim+m);
            }
        }
    }//w

    if      (FtnOut_c.NumSite_ == NumSite_)    ImgFreqFtn::update(FtnOutM,mixing, mixingType);
    else if (FtnOut_c.NumSite_ == 1  ) {
        for (int w=0; w < Nfreq_; w++ ) {
            for(int n=0; n < dim; n++) {
                for(int m=0; m < dim; m++) {
                    int n1 = updateSite *  NSpinOrbitPerAtom_ + n;
                    int m1 = updateSite *  NSpinOrbitPerAtom_ + m;
                    Ftn_(w,n1*Norbit+m1) = (1-mixing) * Ftn_(w,n1*Norbit+m1) + mixing * (FtnOutM[w](n,m));
                }
            }
        }
    }
    delete [] FtnOutM;
}

void ImgFreqFtn::update(Eigen::MatrixXcd *  FtnOutM, double mixing, int mixingType) {
    assert(FtnOutM[0].size() == Norbit*Norbit);
    if( mixingType == 0 ) {
        if(mpi_rank==0 and mixing != 1)  std::cout << "Simple Mixing:"<<mixing<<"\n";
        for (int w=0; w < Nfreq_; w++ ) {
            for(int n=0; n < Norbit; n++) {
                for(int m=0; m < Norbit; m++) {
                    Ftn_(w,n*Norbit+m) = (1-mixing) * Ftn_(w,n*Norbit+m) + mixing * (FtnOutM[w](n,m));
                }
            }
        }
    }
}



void ImgFreqFtn::update_full(const std::string &filename, double mixing) {
    cmplx I(0,1);
    /*check file*/
    FILE * Fcheck;
    int Nfileopen=0;
    while (1) {
        Nfileopen++;
        Fcheck = fopen(filename.c_str(),"r");
        if (Fcheck==NULL) {
            sleep(5);
        }
        else {
            fclose(Fcheck);
            break;

        }
        if (Nfileopen>3) {
            std::cout << "Cannot open file : " << filename << "\n";
            assert(0);
        }
    }
    sleep(2);

    /*Start read*/
    int w=0,n,m, check=0, wmax=0;
    double Retemp, Imtemp;
    FILE *dataIN = fopen(filename.c_str(), "r");
    while (  !feof(dataIN) ) {
        fscanf(dataIN, "%d %d %d %lf %lf\n",&w,&m,&n,&Retemp, &Imtemp);
        if(n<0 or m<0 or n>=Norbit or m>=Norbit) exit(1);
//        if(w-startIndx_< Nfreq_ and 0 <= w-startIndx_ ) Ftn_(w-startIndx_, m*Norbit+n) =  -0.0 -I*0.1;
        if(w-startIndx_< Nfreq_ and 0 <= w-startIndx_ ) Ftn_(w-startIndx_, m*Norbit+n) =  Retemp + I* Imtemp;
        if (wmax < w) wmax=w;
    }//while
    fclose(dataIN);

    //length of the file is too short
//    if(wmax+1!=Nfreq_) {
//        std::cout << "WARNING: CHECK FILE SIZE " << filename <<"\n";
//    }
    //for(int w=wmax+1; w<Nfreq_; w++) {
    //    for(n=0; n < Norbit*Norbit; n++) {
    //        Ftn_(w,n)= Ftn_(Nfreq_,n);
    //    }
    //}
    assert (*MEMCHECK == magicNumber);
}


void ImgFreqFtn::read_full(const std::string &filename,double beta, double beta_prev) {
    cmplx I(0,1);
    assert(startIndx_==0);
    /*estimate file size*/
    int w=0,n,m,  wmax=0;
    double Retemp, Imtemp;
    FILE *dataIN = fopen(filename.c_str(), "r");
    while (  !feof(dataIN) ) {
        fscanf(dataIN, "%d %d %d %lf %lf\n",&w,&m,&n,&Retemp, &Imtemp);
        if (wmax < w) wmax=w;
    }
    fclose(dataIN);
    wmax++;
    //
    //Initialize
    //

    Eigen::MatrixXcd * Swread = new Eigen::MatrixXcd [wmax];
    for (int i = 0; i < wmax; i++)
    {
        Swread[i].setZero(Norbit,Norbit);
    }
    std::vector<double> matsubara;
    matsubara.resize(wmax);

    //
    //Start read
    //
    dataIN = fopen(filename.c_str(), "r");
    while (  !feof(dataIN) ) {
        fscanf(dataIN, "%d %d %d %lf %lf\n",&w,&m,&n,&Retemp, &Imtemp);
        if(w>=0) {
            matsubara.at(w) = ((2*w+1)*pi/beta_prev);
            (Swread[w])(m,n) = Retemp+I*Imtemp;
        }
//        else if( w<0 ){
//          (Swread.at(wmax+(std::abs(w)-1)))[m*Norbit+n] =  Retemp + I* Imtemp;
//           Ftn_(Nfreq_+(std::abs(w)-1), m*Norbit+n) =  Retemp + I* Imtemp;
//        }
    }//while
    fclose(dataIN);

    int    wprime=0;
    double this_matsubara;
    for(int w=0; w<Nfreq_; w++) {
        this_matsubara = ((2*w+1)*pi/beta);
        if(this_matsubara <= matsubara.at(0) ) {
            for(int m=0; m < Norbit; m++) {
                for(int n=0; n < Norbit; n++) {
                    Ftn_(w,m*Norbit+n)  = (Swread[1](m,n)-Swread[0](m,n))/(2*pi/beta_prev) *(this_matsubara - matsubara[0])  + Swread[0](m,n);
                }
            }
        }
        else if( matsubara.at(0) < this_matsubara and this_matsubara <= matsubara.at(wmax-1)) {
            while (true) {
                if( matsubara.at(wprime) < this_matsubara and this_matsubara <= matsubara.at(wprime+1)) {
                    for(int m=0; m < Norbit; m++) {
                        for(int n=0; n < Norbit; n++) {
                            Ftn_(w,m*Norbit+n)  = (Swread[wprime+1](m,n)-Swread[wprime](m,n))/(2*pi/beta_prev) *(this_matsubara - matsubara.at(wprime))  + Swread[wprime](m,n);
                        }
                    }
                    break;
                }
                else {
                    wprime++;
                }
            }
        }
        else if(this_matsubara > matsubara.at(wmax-1)) {
            Eigen::MatrixXcd RealSw = (Swread[wmax-1]+Swread[wmax-1].adjoint())/2.0;
            Eigen::MatrixXcd ImagSw = (Swread[wmax-1]-Swread[wmax-1].adjoint())/(2.0*I);
            for(int m=0; m < Norbit; m++) {
                for(int n=0; n < Norbit; n++) {
                    Ftn_(w,m*Norbit+n)= RealSw(m,n) + (matsubara[wmax-1]*ImagSw/(this_matsubara))(m,n);
                }
            }
        }
    }
    delete [] Swread;
}




void ImgFreqFtn::mpiBcast(int root, MPI_Comm comm) {
    MPI_Bcast(imgFreq , Nfreq_, MPI_DOUBLE, root, comm);
    MPI_Bcast(Ftn_.data(), Ftn_.size(), MPI_DOUBLE_COMPLEX, root, comm);
    MPI_Bcast(&startIndx_,1, MPI_INT, root, comm);
}


void ImgFreqFtn::estimate_asymto(int order) {
    Eigen::MatrixXcd tail [3];
    Eigen::MatrixXcd tail_coeff;
    double w0,w1,w2;
    if(order==1) {
        //Asymto, i*S1/w  + ...
        w0 = 1./std::pow((*this).getValue(Nfreq_-3),2);
        w1 = 1./std::pow((*this).getValue(Nfreq_-2),2);
        w2 = 1./std::pow((*this).getValue(Nfreq_-1),2);
        tail[0] = (*this).getMatrix(Nfreq_-3);
        tail[1] = (*this).getMatrix(Nfreq_-2);
        tail[2] = (*this).getMatrix(Nfreq_-1);
        for(int w=0; w<3; w++) {
            tail[w] = (tail[w].adjoint()-tail[w]).eval()/2.;
            tail[w] =  tail[w] / (I* (*this).getValue(Nfreq_-3+w));
        }
        tail_coeff = ( tail[0] + tail[1] + tail[2]) / (w0+w1+w2);
        (*this).setMatrix(Nfreq_+1,  tail_coeff);// sign : 1/i = -i
    }
    else if(order==2) {
        //Asymto S2/(iw^2) + ...
        w0 = 1./std::pow( (*this).getValue(Nfreq_-3)  ,4);
        w1 = 1./std::pow( (*this).getValue(Nfreq_-2)  ,4);
        w2 = 1./std::pow( (*this).getValue(Nfreq_-1)  ,4);

        tail[0] = (*this).getMatrix(Nfreq_-3) ;
        tail[1] = (*this).getMatrix(Nfreq_-2) ;
        tail[2] = (*this).getMatrix(Nfreq_-1) ;
        for(int w=0; w<3; w++) {
            tail[w] =  (tail[w]+tail[w].adjoint()).eval()/2.;
            tail[w] =  -tail[w] / std::pow((*this).getValue(Nfreq_-3+w),2);
        }
        tail_coeff =  (tail[0] + tail[1] + tail[2]) / (w0+w1+w2);
        (*this).setMatrix(Nfreq_+2, tail_coeff);
//        FreqFtn[N_freq+1] = tail_coeff;
    }
    else if(order==3) {
//        //Asymto S2/(iw^2) + ...
//        w0 = 1./std::pow((2*(N_freq-3)+1)*pi/beta,6);
//        w1 = 1./std::pow((2*(N_freq-2)+1)*pi/beta,6);
//        w2 = 1./std::pow((2*(N_freq-1)+1)*pi/beta,6);
//
//        tail[0] = FreqFtn[N_freq-3] -   (FreqFtn[N_freq] / (I*(2*(N_freq-3)+1)*pi/beta)) ;
//        tail[1] = FreqFtn[N_freq-2] -   (FreqFtn[N_freq] / (I*(2*(N_freq-2)+1)*pi/beta)) ;
//        tail[2] = FreqFtn[N_freq-1] -   (FreqFtn[N_freq] / (I*(2*(N_freq-1)+1)*pi/beta)) ;
//        for(int w=0; w<3; w++) {
//            tail[w] =  (tail[w].adjoint()-tail[w]).eval()/2.;
//            tail[w] =  I*tail[w] / std::pow((2*(N_freq-3+w)+1)*pi/beta,3);  //1/i^3=-1/i = i
//        }
//        tail_coeff =  (tail[0] + tail[1] + tail[2]) / (w0+w1+w2);
//        FreqFtn[N_freq+2] = tail_coeff;
    }
}
