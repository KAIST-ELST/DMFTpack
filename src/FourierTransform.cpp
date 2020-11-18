#include <Eigen/Eigenvalues>
#include <unsupported/Eigen/MatrixFunctions>
#include "tight_common.h"
#include "TB.h"
#include "mpi.h"

//FT from iw to \tau
//F(t) = 1/\beta  \sum_w F(iw) e^(-iwt).

/*delta_t shold be zero vectors*/




/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void FourierTransform_tau (Eigen::MatrixXcd * delta_w, Eigen::MatrixXcd  & delta_t, double  tau, std::vector<Eigen::MatrixXcd> moments) {
//
//F(iw_n -> \infty) =  V^+ 1/(iw-H) V = (V^+) (1/iw + H/iw^2 + H^2/iw^3+...)V
//
    int Dim = delta_w[0].rows();

    Eigen::MatrixXcd AsymtoV;
    Eigen::MatrixXcd AsymtoH;
    Eigen::MatrixXcd AsymtoH2;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> ces(moments[1]);
    ces.compute(moments[1]);
    double  a = ces.eigenvalues().minCoeff();
    if(a>0)   {
        AsymtoV = (ces.operatorSqrt()).eval();
        Eigen::MatrixXcd M2prime = AsymtoV.adjoint().inverse() * (2.0 * moments[2]) * AsymtoV.inverse();
        Eigen::MatrixXcd M3prime = AsymtoV.adjoint().inverse() * (2.0 * moments[3]) * AsymtoV.inverse();

        ces.compute( (M3prime/2 - M2prime*M2prime/4) );
        if(  ces.eigenvalues().minCoeff() > 0   ) {
            AsymtoH = M2prime/2.0 + ces.operatorSqrt();
        }
        else   AsymtoH = M2prime/2.0 ;
        AsymtoH2 =  M2prime - AsymtoH;

        double w0local, dwlocal,w ;
        dwlocal = (2*pi)/beta;
        w0local =(pi/beta);
        delta_t.setZero(delta_w[0].rows(), delta_w[0].cols());

        Eigen::MatrixXcd id;
        id.setIdentity( Dim,Dim);

        for(int i=0; i<N_freq; i++) {
            w=w0local+i*dwlocal ;
            Eigen::MatrixXcd tail =  AsymtoV.adjoint() * (( I*w * id  - AsymtoH).inverse()) * AsymtoV
                                     + AsymtoV.adjoint() * (( I*w * id  - AsymtoH2).inverse()) * AsymtoV  ;
            tail /= 2.0;
            delta_t +=   (delta_w[i] - tail ) * std::exp(-I*w*tau);
        }
        delta_t = (delta_t + delta_t.adjoint()).eval(); //alliasing issue and negative freq.sum
        delta_t =  ( delta_t  ) / beta ;


        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> ces2(AsymtoH);
        // G0H = -exp(-\tau H)/(1+exp(-\beta H)
        Eigen::MatrixXcd G0H ;
        G0H.setIdentity(Dim,Dim);
        for(int i=0; i<Dim; i++) {
            G0H(i,i) =  - 1.0 /(std::exp(ces2.eigenvalues()(i) * tau) + std::exp(ces2.eigenvalues()(i) * (tau-beta))   );
        }
        G0H = (ces2.eigenvectors() * G0H * (ces2.eigenvectors().adjoint())).eval();



        Eigen::MatrixXcd G0H2 ;
        G0H2.setIdentity(Dim,Dim);
        ces2.compute(AsymtoH2);
        for(int i=0; i<Dim; i++) {
            G0H2(i,i) =  - 1.0 /(std::exp(ces2.eigenvalues()(i) * tau) + std::exp(ces2.eigenvalues()(i) * (tau-beta))   );
        }
        G0H2 = (ces2.eigenvectors() * G0H2 * (ces2.eigenvectors().adjoint())).eval();


        delta_t += AsymtoV.adjoint() * G0H * AsymtoV /2.0;
        delta_t += AsymtoV.adjoint() * G0H2 * AsymtoV/2.0;


        delta_t = (delta_t + delta_t.adjoint()).eval(); //alliasing issue and negative Hermicity
        delta_t /= 2.0;
    }
    else {
        delta_t  =  (moments[1]*0.5 + moments[2]*(0.5*tau  - beta/4.) )/2.0;
        for(int n=0; n<N_freq; n++) {
            double wn = pi*(2.*n+1.)/beta  ;
            Eigen::MatrixXcd tail = moments[1] / (I*wn)  + moments[2] / (I*wn*I*wn);
            tail /= 2.0;
            delta_t +=    (delta_w[n] - tail) *  std::exp(-I*wn*tau);
        }
        delta_t = (delta_t + delta_t.adjoint()).eval();
        delta_t =  ( delta_t  ) / beta ;
    }
}
Eigen::MatrixXcd  FourierTransform_tau (Eigen::MatrixXcd * delta_w, double  tau, bool greenFtnAsym) {
    Eigen::MatrixXcd   delta_t;
    std::vector<Eigen::MatrixXcd> moments(4);
    getAsymto_moments(moments, delta_w, greenFtnAsym);
    FourierTransform_tau( delta_w, delta_t, tau, moments);
    return delta_t;
}




/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void FourierTransform_(Eigen::MatrixXcd * delta_w, Eigen::MatrixXcd * delta_t, std::vector<Eigen::MatrixXcd> moments)
{
// tree tree


    int Dim = delta_w[0].rows();
    for(int j=0; j<N_tau+1; j++) {
        delta_t[j].setZero(Dim,Dim);
    }

    int  mysta, myend;
    double tau, tau1, tau_prev, mineigval11;
    para_range(0,N_tau, mpi_numprocs, mpi_rank, &mysta, &myend);
    for(int j=mysta; j<=myend; j++) {
        tau=j* beta/N_tau;
        FourierTransform_tau (delta_w, delta_t[j], tau, moments );
    }
    data_sync_EigenMat(delta_t, 0, N_tau, delta_t[0].rows(),  mpi_numprocs);

    int CHECK =0;
    for(int t=0; t<N_tau+1; t++) {
        tau = t*beta/N_tau;

        Eigen::MatrixXcd   test;
        test = delta_t[t];
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> ces00(test); // norm > 0
        ces00.compute(test);
        double mineigval =  ces00.eigenvalues().minCoeff();
        mineigval11 = mineigval;
        if( mineigval >=0 and (t!=0 and t!=N_tau) ) {
            CHECK=1;
            for(int t1=t+1; t1<N_tau+1; t1++) {
                tau1= t1*beta/N_tau;
                tau_prev = (t-1)*beta/N_tau;
                delta_t[t] = (delta_t[t1] - delta_t[t-1]) / (tau1- tau_prev )  * ( tau - tau_prev)  +  delta_t[t-1];

                Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> ces11(delta_t[t]); // norm > 0
                ces11.compute(delta_t[t]);
                mineigval11 =  ces11.eigenvalues().minCoeff();

                if( mineigval11 < 0)     break;
            }//t1
        }
        if(mineigval11 >=0  ) {
            CHECK=1;
            delta_t[t].setIdentity(Dim,Dim);
            delta_t[t] *= -1e-10;
        }
    }
    ifroot   if(CHECK==1)      std::cout <<" Warning : positive Gt\n"   ;


}

//////////////////////////////////////////////////////////////////////////////////////////////////Vy
void FourierTransform (Eigen::MatrixXcd * delta_w, Eigen::MatrixXcd * delta_t, bool greenFtnAsym,  std::string HYB_Ftn, int print) {
    if(mpi_rank==0 and print ==1)
        std::cout << HYB_Ftn <<"\n";

    std::vector<Eigen::MatrixXcd> moments(4);
    getAsymto_moments(moments, delta_w, greenFtnAsym, print);
    FourierTransform_(delta_w, delta_t, moments);
}
void FourierTransform (Eigen::MatrixXcd * delta_w, Eigen::MatrixXcd * delta_t, Eigen::MatrixXcd moments1, bool greenFtnAsym,  std::string HYB_Ftn,   int print) {
    if(mpi_rank==0 and print ==1)
        std::cout << HYB_Ftn <<"\n";

    std::vector<Eigen::MatrixXcd> moments(4);
    getAsymto_moments(moments, delta_w, greenFtnAsym);
    int Dim = delta_w[0].rows();
    moments[0].setZero(Dim,Dim );
    moments[1] = moments1;
    FourierTransform_(delta_w, delta_t, moments);
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void FourierTransform (Eigen::MatrixXcd * delta_w, Eigen::MatrixXcd  & delta_t, double  tau, bool greenFtnAsym,    std::string HYB_Ftn, int print ) {
    if(mpi_rank==0 and print ==1)
        std::cout << HYB_Ftn <<"\n";

    std::vector<Eigen::MatrixXcd> moments(4);
    getAsymto_moments(moments, delta_w, greenFtnAsym);
    FourierTransform_tau(delta_w, delta_t, tau, moments);
}

void FourierTransform (Eigen::MatrixXcd * delta_w, Eigen::MatrixXcd  & delta_t, double  tau, Eigen::MatrixXcd moments1, bool greenFtnAsym,     std::string HYB_Ftn, int print) {
    if(mpi_rank==0 and print ==1)
        std::cout << HYB_Ftn <<"\n";
    std::vector<Eigen::MatrixXcd> moments(4);
    getAsymto_moments(moments, delta_w, greenFtnAsym);
    int Dim = delta_w[0].rows();
    moments[0].setZero(Dim,Dim );
    moments[1] = moments1;
    FourierTransform_tau(delta_w, delta_t, tau, moments );
}







//////////////////////////////////////////////////////////////////////////////////////////////////Vy
void FT_t_to_w (Eigen::MatrixXcd * delta_w, Eigen::MatrixXcd * delta_t, int N_freq) {
    /*Parallelized, Trapezoidal Method for continuous time integral  */
    int N_segment=N_tau;
    std::vector<Eigen::MatrixXcd> at(N_segment);
    std::vector<Eigen::MatrixXcd> bt(N_segment);
    double    dt=beta/N_segment;
    for ( int n=0; n<N_freq; n++) {
        delta_w[n].setZero(delta_t[0].rows(), delta_t[0].cols());
    }

    for(int t=0 ; t<N_segment; t++) {
        double tau0 = t*dt;
        double tau1 = (t+1)*dt  ;
        at[t]=  (delta_t[t+1] - delta_t[t])/ ( tau1 - tau0 );
        bt[t]=  -at[t] * tau0 + delta_t[t];
    }
    int mysta, myend;
    para_range(0,N_freq-1, mpi_numprocs, mpi_rank, &mysta, &myend);
    for ( int n=mysta; n<=myend; n++) {
        double wn = (2*n+1)*pi /beta;
        for(int t=0 ; t< N_segment; t++) {
            double tau0 = t*dt;
            double tau1 = (t+1)*dt  ;
            cmplx  partialI =  ( std::exp(I*tau1*wn) - std::exp(I*tau0*wn) );
            cmplx  partialI2 = ( tau1* std::exp(I*tau1*wn) -  tau0 *std::exp(I*tau0*wn) );
            delta_w[n] += at[t]*(partialI/(std::pow(wn,2)) + partialI2/ (I*wn)) +      bt[t] * partialI/(I*wn);
        }//t



        Eigen::MatrixXcd   test;
        test = (delta_w[n] - delta_w[n].adjoint())/(2*I);
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> ces(test); // norm > 0
        ces.compute(test);
        if(  ces.eigenvalues().maxCoeff() > 0   ) {
            std::cout <<"causality?" << ces.eigenvalues().maxCoeff() <<"\n";
//            assert( ces.eigenvalues().maxCoeff() <= 0  );
        }
    }
    data_sync_EigenMat(delta_w, 0, N_freq-1, delta_w[0].rows(),  mpi_numprocs);
}







/////////////////
void getAsymto_moments (std::vector<Eigen::MatrixXcd> & moments, Eigen::MatrixXcd * delta_w, bool greenFtnAsym, bool rstout ) {


    int Dim = delta_w[0].rows();
    int w_Asymto_start= N_freq/2;



    moments.resize(4);   //   m0/iw + m1/(iw**2) + m2/(iw**3) + m3/(iw**4) =     (-i/w, -1/w^2, i/w^3, 1/w^4)  *  (m0,m1,m2,m3)'
    if (greenFtnAsym == false ) {
//    if ( true  ) {}
        int w_Asymto_end  = N_freq-4;
        int w_Asymto_range = w_Asymto_end - w_Asymto_start +1 ;

        Eigen::VectorXcd ** momentstest = new Eigen::VectorXcd * [Dim*Dim];
        for(int i=0; i<Dim*Dim; i++) {
            momentstest[i] = new Eigen::VectorXcd    [w_Asymto_range];
        }
        Eigen::VectorXd momentstrace[4];
        momentstrace[0].resize(w_Asymto_range);
        momentstrace[1].resize(w_Asymto_range);
        momentstrace[2].resize(w_Asymto_range);
        momentstrace[3].resize(w_Asymto_range);

        for (int w_Asymto = 0;  w_Asymto < w_Asymto_range; w_Asymto++) { //wn=w_Asymto+w_Asymto_start
            int ASlen = N_freq - (w_Asymto + w_Asymto_start) ;   //N_freq-w_Asymto_start ,..., N_freq-w_Asymto_end
            Eigen::MatrixXcd KM;
            KM.setZero(ASlen*2, 4);
            std::vector<Eigen::VectorXcd> GAs(Dim*Dim);
            for(int alp=0; alp<Dim; alp++) {
                for(int bet=0; bet<Dim; bet++) {
                    GAs[alp*Dim+bet].setZero(ASlen*2);
                }
            }
            momentstrace[0][w_Asymto]=0.;
            momentstrace[1][w_Asymto]=0.;
            momentstrace[2][w_Asymto]=0.;
            momentstrace[3][w_Asymto]=0.;


            for (int jj =0; jj < ASlen; jj++) {
                int wn= w_Asymto+w_Asymto_start+jj   ;// = w_Asymto + w_Asymto_start , ..., N_freq-1 = w_Asymto_start+jj,..., w_Asymto_end+jj ;
                double z = (2*wn+1)*pi/beta;
                KM(2*jj, 0) =  1.0;                     //Re
                KM(2*jj+1, 1) = -1./std::pow(z,1);      //Im
                KM(2*jj, 2) = -1./std::pow(z,2);        //Re
                KM(2*jj+1, 3) =  1./std::pow(z,3);      //Im
                for(int alp=0; alp<Dim; alp++) {
                    for(int bet=0; bet<Dim; bet++) {
                        GAs[alp*Dim+bet](2*jj+0)= (delta_w[wn](alp,bet) + std::conj(delta_w[wn](bet,alp)))/2;         //Re
                        GAs[alp*Dim+bet](2*jj+1)= (delta_w[wn](alp,bet) - std::conj(delta_w[wn](bet,alp)))/(2*I);     //Im
                    }
                }
            }
            //G=K * M =>   K^+ G = K^+ K M   =>   (K'K)^-1 K' G = M
            Eigen::MatrixXcd KM_solve  = (KM.adjoint() * KM).inverse() * KM.adjoint();
            for(int alp=0; alp<Dim; alp++) {
                for(int bet=0; bet<Dim; bet++) {
                    momentstest[alp*Dim+bet][w_Asymto] = KM_solve * GAs[alp*Dim+bet];
                }
                momentstrace[0][w_Asymto]+= real(momentstest[alp*Dim +alp][w_Asymto][0]);
                momentstrace[1][w_Asymto]+= real(momentstest[alp*Dim +alp][w_Asymto][1]);
                momentstrace[2][w_Asymto]+= real(momentstest[alp*Dim +alp][w_Asymto][2]);
                momentstrace[3][w_Asymto]+= real(momentstest[alp*Dim +alp][w_Asymto][3]);
            }
        }//w_Asymto

        int Nv = N_freq/20;
        if(Nv==0) Nv=1;


        std::vector<Eigen::MatrixXcd > moments0_avg;
        std::vector<Eigen::MatrixXcd > moments1_avg;
        std::vector<Eigen::MatrixXcd > moments2_avg;
        std::vector<Eigen::MatrixXcd > moments3_avg;
        for(int w_Asymto_center=Nv; w_Asymto_center < w_Asymto_range-Nv; w_Asymto_center++) {
            Eigen::MatrixXcd a, b,c,d;
            a.setZero(Dim,Dim);
            b.setZero(Dim,Dim);
            c.setZero(Dim,Dim);
            d.setZero(Dim,Dim);
            moments0_avg.push_back(a);
            moments1_avg.push_back(b);
            moments2_avg.push_back(c);
            moments3_avg.push_back(d);
            for (int iw=w_Asymto_center-Nv; iw<=w_Asymto_center+Nv; iw++) {
                for(int alp=0; alp<Dim; alp++) {
                    for(int bet=0; bet<Dim; bet++) {
                        moments0_avg[w_Asymto_center-Nv](alp,bet)+= momentstest[alp*Dim+bet][iw][0];
                        moments1_avg[w_Asymto_center-Nv](alp,bet)+= momentstest[alp*Dim+bet][iw][1];
                        moments2_avg[w_Asymto_center-Nv](alp,bet)+= momentstest[alp*Dim+bet][iw][2];
                        moments3_avg[w_Asymto_center-Nv](alp,bet)+= momentstest[alp*Dim+bet][iw][3];
                    }
                }
            }
            moments0_avg[w_Asymto_center-Nv] /= (2*Nv+1);
            moments1_avg[w_Asymto_center-Nv] /= (2*Nv+1);
            moments2_avg[w_Asymto_center-Nv] /= (2*Nv+1);
            moments3_avg[w_Asymto_center-Nv] /= (2*Nv+1);
        }


//        moments.resize(4);   //  m0 +  m1/iw + m2/(iw**2) + m3/(iw**3) =     (-i/w, -1/w^2, i/w^3, 1/w^4)  *  (m0,m1,m2,m3)'
        for (int mom=0; mom<4; mom++) {
            std::vector<double> var_j;
            for(int w_Asymto_center=Nv; w_Asymto_center < w_Asymto_range-Nv; w_Asymto_center++) {
                Eigen::ArrayXd vec = momentstrace[mom].segment(w_Asymto_center-Nv, 2*Nv+1).array();
                var_j.push_back(   std::sqrt( (vec - vec.mean()).square().sum()/(vec.size()) )   );
                //w_Asymto_center = Nv, Nv+1, ... w_Asymto_range-Nv-1   =>   (w_Asymto_range - 2*Nv)  loops
            }

            int minPos = var_j.size()-1 ;    // =w_Asymto_range-2*Nv-1
            for (int i = 0; i < var_j.size(); i++)
            {
                Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> ces(moments1_avg[i]); // norm > 0
                ces.compute(moments1_avg[i]);

                if (ces.eigenvalues().minCoeff() >0 ) {
                    Eigen::MatrixXcd V_inv = (ces.operatorSqrt()).inverse();



                    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> ces2(moments1_avg[i]);
                    ces2.compute(     V_inv*moments3_avg[i]*V_inv
                                      - (V_inv*moments2_avg[i]*V_inv) * (V_inv*moments2_avg[i]*V_inv) );   //   <w^2> -  <w>^2  > 0

                    double  a = ces.eigenvalues().minCoeff();
                    double  b = ces2.eigenvalues().minCoeff();
                    if (var_j[i] < var_j[minPos] and a>0 and b>0)         // Found a smaller min
                        minPos = i;                       //w_Asymto_center = minPos+Nv
                }
            }
            int minPos_w_Asymto_center = minPos +Nv;   // <= w_Asymto_range-Nv-1


            moments[mom].setZero(Dim,Dim);
            for(int w_Asymto=minPos_w_Asymto_center-Nv; w_Asymto <= minPos_w_Asymto_center+Nv; w_Asymto++) {
                assert(w_Asymto < w_Asymto_range);
                for(int alp=0; alp<Dim; alp++) {
                }
                for(int alp=0; alp<Dim; alp++) {
                    for(int bet=0; bet<Dim; bet++) {
                        moments[mom](alp,bet) += (momentstest[alp*Dim+bet][w_Asymto][mom] +  std::conj(momentstest[bet*Dim+alp][w_Asymto][mom]))/2.0;
                    }
                }
            }
            moments[mom] /= (2*Nv+1);
        }
        for(int i=0; i<Dim*Dim; i++) {
            delete [] momentstest[i]  ;
        }
        delete  [] momentstest;
    }
    else if (greenFtnAsym == true ) {
        int w_Asymto_end  = N_freq-2;
        int w_Asymto_range = w_Asymto_end - w_Asymto_start +1 ;

        Eigen::VectorXcd ** momentstest = new Eigen::VectorXcd * [Dim*Dim];
        for(int i=0; i<Dim*Dim; i++) {
            momentstest[i] = new Eigen::VectorXcd    [w_Asymto_range];
        }
        Eigen::VectorXd momentstrace[2];
        momentstrace[0].resize(w_Asymto_range);
        momentstrace[1].resize(w_Asymto_range);

        for (int w_Asymto = 0;  w_Asymto < w_Asymto_range; w_Asymto++) { //wn=w_Asymto+w_Asymto_start
            int ASlen = N_freq - (w_Asymto + w_Asymto_start) ;   //N_freq-w_Asymto_start ,..., N_freq-w_Asymto_end
            Eigen::MatrixXcd KM;
            KM.setZero(ASlen*2, 2);
            std::vector<Eigen::VectorXcd> GAs(Dim*Dim);
            for(int alp=0; alp<Dim; alp++) {
                for(int bet=0; bet<Dim; bet++) {
                    GAs[alp*Dim+bet].setZero(ASlen*2);
                }
            }
            momentstrace[0][w_Asymto]=0.;
            momentstrace[1][w_Asymto]=0.;


            for (int jj =0; jj < ASlen; jj++) {
                int wn= w_Asymto+w_Asymto_start+jj   ;// = w_Asymto + w_Asymto_start , ..., N_freq-1 = w_Asymto_start+jj,..., w_Asymto_end+jj ;
                double z = (2*wn+1)*pi/beta;
                KM(2*jj,   0) = -1./std::pow(z,2);        //Re
                KM(2*jj+1, 1) =  1./std::pow(z,3);      //Im
                for(int alp=0; alp<Dim; alp++) {
                    for(int bet=0; bet<Dim; bet++) {
                        GAs[alp*Dim+bet](2*jj+0)= (delta_w[wn](alp,bet) + std::conj(delta_w[wn](bet,alp)))/2;         //Re
                        GAs[alp*Dim+bet](2*jj+1)= (delta_w[wn](alp,bet) - std::conj(delta_w[wn](bet,alp)))/(2*I) ;     //ImG  = - 1 /wn  + ? /(wn)^3
                    }
                    GAs[alp*Dim+alp](2*jj+1) += 1.0 / z;     //ImG   +  1 /wn  =  ? /(wn)^3
                }
            }
            //G=K * M =>   K^+ G = K^+ K M   =>   (K'K)^-1 K' G = M
            Eigen::MatrixXcd KM_solve  = (KM.adjoint() * KM).inverse() * KM.adjoint();
            for(int alp=0; alp<Dim; alp++) {
                for(int bet=0; bet<Dim; bet++) {
                    momentstest[alp*Dim+bet][w_Asymto] = KM_solve * GAs[alp*Dim+bet];
                }
                momentstrace[0][w_Asymto]+= real(momentstest[alp*Dim +alp][w_Asymto][0]);
                momentstrace[1][w_Asymto]+= real(momentstest[alp*Dim +alp][w_Asymto][1]);
            }
        }//w_Asymto

        int Nv = N_freq/20;
        if(Nv==0) Nv=1;


        std::vector<Eigen::MatrixXcd > moments0_avg;
        std::vector<Eigen::MatrixXcd > moments1_avg;
        for(int w_Asymto_center=Nv; w_Asymto_center < w_Asymto_range-Nv; w_Asymto_center++) {
            Eigen::MatrixXcd a, b,c,d;
            a.setZero(Dim,Dim);
            b.setZero(Dim,Dim);
            moments0_avg.push_back(a);
            moments1_avg.push_back(b);
            for (int iw=w_Asymto_center-Nv; iw<=w_Asymto_center+Nv; iw++) {
                for(int alp=0; alp<Dim; alp++) {
                    for(int bet=0; bet<Dim; bet++) {
                        moments0_avg[w_Asymto_center-Nv](alp,bet)+= momentstest[alp*Dim+bet][iw][0];
                        moments1_avg[w_Asymto_center-Nv](alp,bet)+= momentstest[alp*Dim+bet][iw][1];
                    }
                }
            }
            moments0_avg[w_Asymto_center-Nv] /= (2*Nv+1);
            moments1_avg[w_Asymto_center-Nv] /= (2*Nv+1);
        }


        for (int mom=0; mom<2; mom++) {
            std::vector<double> var_j;
            for(int w_Asymto_center=Nv; w_Asymto_center < w_Asymto_range-Nv; w_Asymto_center++) {
                Eigen::ArrayXd vec = momentstrace[mom].segment(w_Asymto_center-Nv, 2*Nv+1).array();
                var_j.push_back(   std::sqrt( (vec - vec.mean()).square().sum()/(vec.size()) )   );
                //w_Asymto_center = Nv, Nv+1, ... w_Asymto_range-Nv-1   =>   (w_Asymto_range - 2*Nv)  loops
            }

            int minPos = var_j.size()-1 ;    // =w_Asymto_range-2*Nv-1
            for (int i = 0; i < var_j.size(); i++)
            {
                Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> ces(moments1_avg[i]);
                ces.compute(moments1_avg[i]  - moments0_avg[i]*moments0_avg[i]);

//                Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> ces2(moments1_avg[i]);
//                ces2.compute( moments3_avg[i] - moments1_avg[i]*moments1_avg[i]);

                double  a = ces.eigenvalues().minCoeff();
//                double  b = ces2.eigenvalues().minCoeff();
//                if (var_j[i] < var_j[minPos] and a>0 and b>0)         // Found a smaller min
                if (var_j[i] < var_j[minPos] and a>1  )         // Found a smaller min
                    minPos = i;                       //w_Asymto_center = minPos+Nv
            }
            int minPos_w_Asymto_center = minPos +Nv;   // <= w_Asymto_range-Nv-1


            moments[mom+2].setZero(Dim,Dim);
            for(int w_Asymto=minPos_w_Asymto_center-Nv; w_Asymto <= minPos_w_Asymto_center+Nv; w_Asymto++) {
                assert(w_Asymto < w_Asymto_range);
                for(int alp=0; alp<Dim; alp++) {
                }
                for(int alp=0; alp<Dim; alp++) {
                    for(int bet=0; bet<Dim; bet++) {
                        moments[mom+2](alp,bet) += (momentstest[alp*Dim+bet][w_Asymto][mom] +  std::conj(momentstest[bet*Dim+alp][w_Asymto][mom]))/2.0;
                    }
                }
            }
            moments[mom+2] /= (2*Nv+1);
        }//for_mom
        moments[0].setZero(Dim,Dim);
        moments[1].setIdentity(Dim,Dim);
        for(int i=0; i<Dim*Dim; i++) {
            delete [] momentstest[i]  ;
        }
        delete  [] momentstest;
    }//greenFtn
    moments[0] = (moments[0]+moments[0].adjoint())/2.0;
    moments[1] = (moments[1]+moments[1].adjoint())/2.0;
    moments[2] = (moments[2]+moments[2].adjoint())/2.0;
    moments[3] = (moments[3]+moments[3].adjoint())/2.0;
    if (rstout==true    and mpi_rank==0) {
        std::cout << "moments0\n" << moments[0]  << "\n";
        std::cout << "moments1\n" << moments[1]  << "\n";
        std::cout << "moments2\n" << moments[2]  << "\n";
        std::cout << "moments3\n" << moments[3]  << "\n";
    }




}
