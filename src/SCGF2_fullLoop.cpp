//arXiv::1511.03911
//Implemented by J.H. Sim
//2017/10/16
//
//
#include <iostream>
#include "ImgTimeFtn.h"
#include "pulay.h"
#include "tight_common.h"
#include <Eigen/LU>
#include <Eigen/QR>
#include <Eigen/Eigenvalues>

#define   history                  2

int stdoutLevel=  2 ;

void getoccMat(Eigen::MatrixXcd * Gwimp, Eigen::MatrixXcd & occMat, int solverDim ) ;
void getGwimp(Eigen::MatrixXcd * Gwimp, Eigen::MatrixXcd projimpurity_site_Hamiltonian, Eigen::MatrixXcd Fock,  Eigen::MatrixXcd * Swimp_secondOrder,
              Eigen::MatrixXcd * delta_w,double muTB, int solverDim) ;
void getFockOperator( Eigen::MatrixXcd & Fock, Eigen::MatrixXcd & occMat, int solverDim , std::vector<Eigen::VectorXi> projUindex, std::vector<cmplx > projUtensor        ) ;
void getStimp(Eigen::MatrixXcd * Stimp, Eigen::MatrixXcd * Gtimp, bool skeleton, int solverDim,std::vector<Eigen::VectorXi> projUindex, std::vector<cmplx > projUtensor       ) ;




void writeResults( Eigen::MatrixXcd & projNumMatrix, Eigen::MatrixXcd * Swimp_secondOrder, Eigen::MatrixXcd Fock, Eigen::MatrixXcd * Gwimp,  Eigen::MatrixXcd occMat,
                   ImgFreqFtn & SE_out, ImgFreqFtn & Gwimp_out, int solverDim) ;

void mixing_check_outer(  double normGwimpDiffOuter, double & mixingSCGF);
void mixing_check_inner(double  normOccDiffInner, double & mixingHF);
//void paramagnetic_solution(Eigen::MatrixXcd & occMat ) ;
void estimate_asymto(Eigen::MatrixXcd * FreqFtn, int order);
double MatsubaraFtnnorm(Eigen::MatrixXcd * Swimp_secondOrder, int N_freq);
double MatsubaraFtnnorm(Eigen::MatrixXcd * Ftn1, Eigen::MatrixXcd * Ftn2, int N_freq);


void  HFDC(   Eigen::MatrixXcd *  Gwimp, std::vector<Eigen::VectorXi> Uindex_stronglyCorr, std::vector<cmplx > Utensor_stronglyCorr ,
              std::vector<int> strongCorr,
              Eigen::MatrixXcd & weakFock ) {

    int strongDim = strongCorr.size();

    Eigen::MatrixXcd  occMat(strongDim, strongDim);



    Eigen::MatrixXcd Gw0[N_freq+4];
    Eigen::MatrixXcd Fock(strongDim, strongDim);


    for (int n=0; n<N_freq; n++) {
        Gw0[n].setZero(strongDim, strongDim);
        for (int i=0; i<strongDim; i++) {
            for (int j=0; j<strongDim; j++) {
                Gw0[n](i,j) = Gwimp[n](strongCorr[i],  strongCorr[j]);
            }
        }

    }




    Gw0[N_freq].setIdentity(strongDim, strongDim);
    estimate_asymto(Gw0,2);
    estimate_asymto(Gw0,3);
    Gw0[N_freq+3].setZero(strongDim,strongDim);




    getoccMat(Gw0, occMat,strongDim);
    getFockOperator( Fock , occMat, strongDim, Uindex_stronglyCorr, Utensor_stronglyCorr) ;


    for (int i=0; i<strongDim; i++) {
        for (int j=0; j<strongDim; j++) {

            weakFock(strongCorr[i], strongCorr[j]) -= Fock(i,j);
        }
    }
}


void  SecOrdDC(   Eigen::MatrixXcd *  Gwimp, std::vector<Eigen::VectorXi> Uindex_stronglyCorr, std::vector<cmplx > Utensor_stronglyCorr ,
                  std::vector<int> strongCorr,
                  Eigen::MatrixXcd & weakFock, Eigen::MatrixXcd  * weakSwimp_secondOrder ) {

    int strongDim = strongCorr.size();

    Eigen::MatrixXcd  occMat(strongDim, strongDim);

    Eigen::MatrixXcd Gt0[N_tau+1];
    Eigen::MatrixXcd Stimp[N_tau+1];


    Eigen::MatrixXcd Gw0[N_freq+4];
    Eigen::MatrixXcd Swimp_secondOrder[N_freq+3];
    Eigen::MatrixXcd Fock(strongDim, strongDim);


    for (int n=0; n<N_freq; n++) {
        Gw0[n].setZero(strongDim, strongDim);
        for (int i=0; i<strongDim; i++) {
            for (int j=0; j<strongDim; j++) {
                Gw0[n](i,j) = Gwimp[n](strongCorr[i],  strongCorr[j]);
            }
        }

        Swimp_secondOrder[n].setZero(strongDim, strongDim);
    }
    Swimp_secondOrder[N_freq+0].setZero(strongDim, strongDim);
    Swimp_secondOrder[N_freq+1].setZero(strongDim, strongDim);
    Swimp_secondOrder[N_freq+2].setZero(strongDim, strongDim);


    for(int t =0 ; t<N_tau+1; t++) {
        Gt0[t].setZero(strongDim, strongDim);
        Stimp[t].setZero(strongDim,strongDim);
    }


    Gw0[N_freq].setIdentity(strongDim, strongDim);
    estimate_asymto(Gw0,2);
    estimate_asymto(Gw0,3);
    Gw0[N_freq+3].setZero(strongDim,strongDim);




    getoccMat(Gw0, occMat,strongDim);
    getFockOperator( Fock , occMat, strongDim, Uindex_stronglyCorr, Utensor_stronglyCorr) ;


    FourierTransform( Gw0, Gt0, Gw0[N_freq] );
    getStimp ( Stimp, Gt0, false, strongDim, Uindex_stronglyCorr, Utensor_stronglyCorr);
    FT_t_to_w(Swimp_secondOrder, Stimp, N_freq);



    for (int i=0; i<strongDim; i++) {
        for (int j=0; j<strongDim; j++) {
            for (int n=0; n<N_freq; n++) {
                weakSwimp_secondOrder[n](strongCorr[i], strongCorr[j]) -= Swimp_secondOrder[n](i,j);
            }

            weakFock(strongCorr[i], strongCorr[j]) -= Fock(i,j);
        }
    }
}









//void  SecOrdDC( ImgFreqFtn & Gwimp, Eigen::MatrixXcd & projNumMatrix, std::vector<Eigen::VectorXi> projUindex, std::vector<cmplx > projUtensor, int solverDim,
//                ImgFreqFtn & SE_out ) {
//
//    Eigen::MatrixXcd  occMat(solverDim, solverDim);
//
//    Eigen::MatrixXcd Gt0[N_tau+1];
//    Eigen::MatrixXcd Stimp[N_tau+1];
//
//
//    Eigen::MatrixXcd Gw0[N_freq+4];
//    Eigen::MatrixXcd Swimp_secondOrder[N_freq+3];
//    Eigen::MatrixXcd Fock(solverDim, solverDim);
//
//
//    for (int n=0; n<N_freq; n++) {
//        Gw0[n] = Gwimp.getMatrix(n);
//        Swimp_secondOrder[n].setZero(solverDim, solverDim);
//    }
//    Swimp_secondOrder[N_freq+0].setZero(solverDim, solverDim);
//    Swimp_secondOrder[N_freq+1].setZero(solverDim, solverDim);
//    Swimp_secondOrder[N_freq+2].setZero(solverDim, solverDim);
//
//
//    for(int t =0 ; t<N_tau+1; t++) {
//        Gt0[t].setZero(solverDim, solverDim);
//        Stimp[t].setZero(solverDim,solverDim);
//    }
//
//
//    Gw0[N_freq].setIdentity(solverDim, solverDim);
//    estimate_asymto(Gw0,2);
//    estimate_asymto(Gw0,3);
//    Gw0[N_freq+3].setZero(solverDim,solverDim);
//
//
//
//
//    getoccMat(Gw0, occMat,solverDim);
//    getFockOperator( Fock , occMat, solverDim, projUindex, projUtensor) ;
//
//
//    FourierTransform( Gw0, Gt0, Gw0[N_freq] );
//    getStimp ( Stimp, Gt0, false, solverDim, projUindex, projUtensor);
//    FT_t_to_w(Swimp_secondOrder, Stimp, N_freq);
//
//    writeResults(projNumMatrix,  Swimp_secondOrder, Fock, Gw0,   projNumMatrix,  SE_out, Gwimp, solverDim);
//}





void non_ingeracting_Gw ( Eigen::MatrixXcd * Gw0, int solverDim,  Eigen::MatrixXcd projimpurity_site_Hamiltonian, double muAdju, double NumTotInitial, Eigen::MatrixXcd * delta_w)  {
    Eigen::MatrixXcd Swimp_secondOrder[N_freq+3];
    for (int n=0; n<N_freq; n++) {
        Swimp_secondOrder[n].setZero(solverDim, solverDim);
    }
    Swimp_secondOrder[N_freq].setZero(solverDim, solverDim);
    Swimp_secondOrder[N_freq+1].setZero(solverDim, solverDim);
    Swimp_secondOrder[N_freq+2].setZero(solverDim, solverDim);


    Eigen::MatrixXcd  occMatHF(solverDim, solverDim);
    Eigen::MatrixXcd Fock(solverDim, solverDim);
    Fock.setZero(solverDim, solverDim);

    double Nele[history];
    for (int i=0; i<history; i++)  Nele[i]=-1;
    double muRD;
    double mu_next;
    double muUB=0, muLB=0;
    int nearAns=0;
    double dmu=1;
    muAdju-=dmu;
    int step=0;

    while (  std::abs(Nele[0] - NumTotInitial) > std::min(1e-5, 1e-5/((Fock.norm())+1e-5)  )  ) { //   and  std::abs(dmu)>1e-5 ) {}
        muAdju += dmu;

        getGwimp( Gw0,projimpurity_site_Hamiltonian,Fock, Swimp_secondOrder, delta_w, muAdju, solverDim) ;
        getoccMat(Gw0, occMatHF,solverDim);

        for(int i=history-1; i>0 ; i--) Nele[i]=Nele[i-1];
        Nele[0] = real(occMatHF.trace());

        if (nearAns ==0) {
            if ( Nele[0] < NumTotInitial) {
                muLB = muAdju;
                if (Nele[1]>NumTotInitial && Nele[1] >0 ) nearAns = 1;
                else dmu = fabs(dmu);
            }
            if ( Nele[0] > NumTotInitial) {
                muUB = muAdju;
                if (Nele[1] < NumTotInitial && Nele[1] > 0)  nearAns =1;
                else dmu = -1*fabs(dmu);
            }
            if ( Nele[0] == NumTotInitial) break;
        }
        else if (nearAns ==1) {
            if (Nele[0] > NumTotInitial) muUB = muAdju;
            else if(Nele[0] < NumTotInitial) muLB =muAdju;
            mu_next = (muUB+muLB)/2.;
            dmu = mu_next - muAdju;
        }
    }//while
}



void SecondOrderPerturbation  ( int solverDim,Eigen::MatrixXcd projimpurity_site_Hamiltonian,  Eigen::MatrixXcd & projNumMatrix,
                                ImgFreqFtn & SE_out,      ImgFreqFtn & Gwimp_out,
                                std::vector<Eigen::VectorXi> projUindex, std::vector<cmplx > projUtensor,  bool high_level_solver ,
                                std::vector<Eigen::VectorXi> Uindex_stronglyCorr, std::vector<cmplx > Utensor_stronglyCorr, std::vector<int> strongCorr       ) {
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> ces(solverDim);
    ifroot std::cout <<"Himp (diag):";
    ces.compute(projimpurity_site_Hamiltonian);
    for(int n=0; n<solverDim; n+=1) {
        ifroot std::cout <<std::fixed << std::setprecision(3)<< std::fixed   <<  ( ces.eigenvalues()[n]  ) <<" " ;
    }
    ifroot std::cout <<"\n";

    /*Initial Setting*/
    Eigen::MatrixXcd Fock(solverDim, solverDim);
    Eigen::MatrixXcd  occMat0(solverDim, solverDim);
    Eigen::MatrixXcd  occMatHF(solverDim, solverDim);
    Eigen::MatrixXcd *  Swimp_secondOrder =  new Eigen::MatrixXcd [N_freq+3];

    Eigen::MatrixXcd * Gwimp = new Eigen::MatrixXcd [N_freq+4];
    Eigen::MatrixXcd * Gt0   = new Eigen::MatrixXcd [N_tau+1];
    Eigen::MatrixXcd * Stimp = new Eigen::MatrixXcd [N_tau+1];




    Fock.setZero(solverDim,solverDim);
    for(int t =0 ; t<N_tau+1; t++) {
        Gt0[t].setZero(solverDim, solverDim);
        Stimp[t].setZero(solverDim,solverDim);
    }
    for (int n=0; n<N_freq; n++) {
        Gwimp[n] = Gwimp_out.getMatrix(n);
        Swimp_secondOrder[n].setZero(solverDim, solverDim);
    }
    Swimp_secondOrder[N_freq+0].setZero(solverDim, solverDim);
    Swimp_secondOrder[N_freq+1].setZero(solverDim, solverDim);
    Swimp_secondOrder[N_freq+2].setZero(solverDim, solverDim);

    Gwimp[N_freq+0].setIdentity(solverDim, solverDim);
    estimate_asymto(Gwimp,2);
    estimate_asymto(Gwimp,3);
    Gwimp[N_freq+3].setZero(solverDim, solverDim);





    getoccMat(Gwimp, occMat0,solverDim);
    getFockOperator( Fock , occMat0, solverDim, projUindex, projUtensor) ;



    FourierTransform( Gwimp, Gt0, Gwimp[N_freq] );
    getStimp ( Stimp, Gt0, false, solverDim, projUindex, projUtensor);
    FT_t_to_w(Swimp_secondOrder, Stimp, N_freq);

    if(high_level_solver==false) {
        SecOrdDC(  Gwimp, Uindex_stronglyCorr,  Utensor_stronglyCorr ,
                   strongCorr,
                   Fock,  Swimp_secondOrder ) ;
    }
    /*write result to NumMatrix and Sw*/
    writeResults(projNumMatrix,  Swimp_secondOrder, Fock, Gwimp,   occMat0, SE_out, Gwimp_out, solverDim);
    MPI_Barrier(MPI_COMM_WORLD);
    delete []  Swimp_secondOrder ;

    delete [] Gwimp ;
    delete [] Gt0   ;
    delete [] Stimp ;
}





void SCHF (int solverDim, Eigen::MatrixXcd projimpurity_site_Hamiltonian, Eigen::MatrixXcd & projNumMatrix,
           ImgFreqFtn & SE_out,      ImgFreqFtn & Gwimp_out,  ImgFreqFtn & weiss_field, double muTB,
           std::vector<Eigen::VectorXi> projUindex, std::vector<cmplx > projUtensor, bool high_level_solver,
           std::vector<Eigen::VectorXi> Uindex_stronglyCorr, std::vector<cmplx > Utensor_stronglyCorr, std::vector<int> strongCorr, ImgFreqFtn & SE_strong       ) {

    /*Initial Setting*/
    Eigen::MatrixXcd Fock(solverDim, solverDim);
    Eigen::MatrixXcd  occMat(solverDim, solverDim);
    Eigen::MatrixXcd * delta_w           = new Eigen::MatrixXcd [N_freq+2];
    Eigen::MatrixXcd * Swimp_secondOrder = new Eigen::MatrixXcd [N_freq+3];

    Eigen::MatrixXcd * Gwimp_in =  new Eigen::MatrixXcd [N_freq+4];
    Eigen::MatrixXcd * Gwimp    =  new Eigen::MatrixXcd [N_freq+4];


    Fock.setZero(solverDim,solverDim);
    for (int n=0; n<N_freq; n++) {
        Gwimp[n].setZero(solverDim, solverDim);
        Swimp_secondOrder[n].setZero(solverDim, solverDim);
        delta_w[n] = weiss_field.getMatrix(n);
    }
    Swimp_secondOrder[N_freq].setZero(solverDim, solverDim);
    Swimp_secondOrder[N_freq+1].setZero(solverDim, solverDim);
    Swimp_secondOrder[N_freq+2].setZero(solverDim, solverDim);
    delta_w[N_freq]   = weiss_field.getMatrix(N_freq+1 );
    delta_w[N_freq+1] = weiss_field.getMatrix(N_freq+2 );
    Gwimp[N_freq].setZero(solverDim, solverDim);
    Gwimp[N_freq+1].setZero(solverDim, solverDim);
    Gwimp[N_freq+2].setZero(solverDim, solverDim);
    Gwimp[N_freq+3].setZero(solverDim, solverDim);

    for(int i = 0 ; i<  solverDim ; i++) {
        for(int j = 0 ; j<  solverDim ; j++) {
            occMat(i,j)= projNumMatrix(i,j);
        }
    }
    if(high_level_solver==false) {
        int strongDim = strongCorr.size();
        for (int i=0; i<strongDim; i++) {
            for (int j=0; j<strongDim; j++) {
                for (int n=0; n<N_freq; n++) {
                    Swimp_secondOrder[n](strongCorr[i], strongCorr[j]) =  (SE_strong.getValue(n,i,j) - SE_strong.getValue(N_freq,i,j));
                }
            }
        }
    }

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> ces(solverDim);
    ifroot std::cout <<"Himp (diag):";
    ces.compute(projimpurity_site_Hamiltonian);
    for(int n=0; n<solverDim; n+=1) {
        ifroot std::cout <<std::fixed << std::setprecision(3)<< std::fixed   <<  ( ces.eigenvalues()[n]  ) <<" " ;
    }
    ifroot std::cout <<"\n";




    //non_interacting_GreenFtn
    double NumTotInitial = real(occMat.trace());
    double muAdju = real(projimpurity_site_Hamiltonian.trace())/solverDim;
    non_ingeracting_Gw(Gwimp, solverDim,projimpurity_site_Hamiltonian, muAdju, NumTotInitial, delta_w);
    getoccMat(Gwimp, occMat, solverDim);



    /*SCGF2-loop*/
    int SCGFloop=0;
    double normGwimpDiffOuter=0;
    double mixingSCGF = mixing/10.;
    double mixingHF = mixing/10.;
    int  RDOccdec =0;
    int  RDOccinc =0;
    int  RDSwdec  =0;
    int  RDSwinc  =0;
    double normUnit =  MatsubaraFtnnorm(Gwimp,N_freq);
    double normOccMat=0;
    pulayMixing Mixing_SCGF(3, 100, N_freq , solverDim, solverDim );
    while( (SCGFloop<10 or(normGwimpDiffOuter > 1e-6) or normOccMat>1e-6) and SCGFloop <2000 ) {

        //transfer previous results
        for (int n=0; n<N_freq+4; n++) {
            Gwimp_in[n] = Gwimp[n];
        }
        /*get occupation and Gt*/
        Gwimp[N_freq].setIdentity(solverDim, solverDim);
        getoccMat(Gwimp, occMat, solverDim);
        Eigen::MatrixXcd occMat_in = occMat;


        /*and HF*/
        getFockOperator( Fock , occMat, solverDim, projUindex, projUtensor) ;





        /*obtain results*/
        getGwimp( Gwimp,projimpurity_site_Hamiltonian,Fock, Swimp_secondOrder, delta_w, muTB,solverDim) ;
        getoccMat(Gwimp, occMat, solverDim);

        /*cal. RD*/
        normOccMat = (occMat-occMat_in).norm();
        normGwimpDiffOuter= MatsubaraFtnnorm(Gwimp_in , Gwimp, N_freq)  /  normUnit;

        /*mixing, Gwimp is used in the next iteration*/
        if(SCGFloop%500==0) mixingSCGF /= 1.5;
        if(SCGFloop < 1500) Mixing_SCGF.mixing( Gwimp_in, Gwimp, mixingSCGF, SCGFloop, 1);
        else {
            for(int n=0; n<N_freq; n++) {
                Gwimp[n] =  (1-mixingSCGF) * Gwimp_in[n] + mixingSCGF * Gwimp[n];
            }
        }
        Gwimp[N_freq].setIdentity(solverDim, solverDim);
        estimate_asymto(Gwimp,2);
        estimate_asymto(Gwimp,3);

        if(stdoutLevel>=2 and mpi_rank==0) {


            std::cout <<"oc:";
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> ces(solverDim);
            ces.compute(occMat);
            for(int n=0; n<solverDim; n+=1) {
                std::cout <<std::fixed << std::setprecision(5)<< std::fixed   <<  ( ces.eigenvalues()[n]  ) <<" " ;
            }
            std::cout <<": " << occMat.trace() <<"\n";

        }


        /*mixing check*/
        if (SCGFloop%10 == 0) mixing_check_outer(  normGwimpDiffOuter,  mixingSCGF);

        SCGFloop++;
    }//outer loop

    /*Write resutls*/
    ifroot std::cout << "SCGF loop:" << SCGFloop;
    ifroot printf      ("  GwRD:%e" , normGwimpDiffOuter);
    ifroot printf      ("  mixing:%e\n:" , mixingSCGF);




    /*write result to NumMatrix and Sw*/

    if(high_level_solver==false) {
        int strongDim = strongCorr.size();
        for (int i=0; i<strongDim; i++) {
            for (int j=0; j<strongDim; j++) {
                for (int n=0; n<N_freq; n++) {
                    Swimp_secondOrder[n](strongCorr[i], strongCorr[j]) =  0.0 ;
                }
                Fock(strongCorr[i], strongCorr[j]) -=  SE_strong.getValue(N_freq, i,j);
            }
        }
    }
    writeResults(projNumMatrix, Swimp_secondOrder, Fock, Gwimp,   occMat, SE_out, Gwimp_out, solverDim);

    MPI_Barrier(MPI_COMM_WORLD);

    delete [] delta_w;
    delete [] Swimp_secondOrder;
    delete [] Gwimp_in;
    delete  [] Gwimp;
}//SCHF




void SCGF2 (int solverDim, Eigen::MatrixXcd projimpurity_site_Hamiltonian, Eigen::MatrixXcd & projNumMatrix,
            ImgFreqFtn & SE_out,      ImgFreqFtn & Gwimp_out,  ImgFreqFtn & weiss_field, double muTB,
            std::vector<Eigen::VectorXi> projUindex, std::vector<cmplx > projUtensor, bool high_level_solver,
            std::vector<Eigen::VectorXi> Uindex_stronglyCorr, std::vector<cmplx > Utensor_stronglyCorr, std::vector<int> strongCorr, ImgFreqFtn & SE_strong       ) {

    /*Initial Setting*/
    Eigen::MatrixXcd Fock(solverDim, solverDim);
    Eigen::MatrixXcd  occMat(solverDim, solverDim);
    Eigen::MatrixXcd Stimp[N_tau+1];
    Eigen::MatrixXcd delta_w[N_freq+2];
    Eigen::MatrixXcd Swimp_secondOrder[N_freq+3];

    Eigen::MatrixXcd Gtimp[N_tau+1];
    Eigen::MatrixXcd Gwimp_in[N_freq+4];
    Eigen::MatrixXcd Gwimp[N_freq+4];


    Fock.setZero(solverDim,solverDim);
    for(int t =0 ; t<N_tau+1; t++) {
        Gtimp[t].setZero(solverDim, solverDim);
        Stimp[t].setZero(solverDim,solverDim);
    }
    for (int n=0; n<N_freq; n++) {
        Gwimp[n].setZero(solverDim, solverDim);
        Swimp_secondOrder[n].setZero(solverDim, solverDim);
        delta_w[n] = weiss_field.getMatrix(n);
    }
    Swimp_secondOrder[N_freq].setZero(solverDim, solverDim);
    Swimp_secondOrder[N_freq+1].setZero(solverDim, solverDim);
    Swimp_secondOrder[N_freq+2].setZero(solverDim, solverDim);
    delta_w[N_freq]   = weiss_field.getMatrix(N_freq+1 );
    delta_w[N_freq+1] = weiss_field.getMatrix(N_freq+2 );
    Gwimp[N_freq].setZero(solverDim, solverDim);
    Gwimp[N_freq+1].setZero(solverDim, solverDim);
    Gwimp[N_freq+2].setZero(solverDim, solverDim);
    Gwimp[N_freq+3].setZero(solverDim, solverDim);

    for(int i = 0 ; i<  solverDim ; i++) {
        for(int j = 0 ; j<  solverDim ; j++) {
            occMat(i,j)= projNumMatrix(i,j);
        }
    }

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> ces(solverDim);
    ifroot std::cout <<"Himp (diag):";
    ces.compute(projimpurity_site_Hamiltonian);
    for(int n=0; n<solverDim; n+=1) {
        ifroot std::cout <<std::fixed << std::setprecision(3)<< std::fixed   <<  ( ces.eigenvalues()[n]  ) <<" " ;
    }
    ifroot std::cout <<"\n";


    //non_interacting_GreenFtn
    double NumTotInitial = real(occMat.trace());
    double muAdju = real(projimpurity_site_Hamiltonian.trace())/solverDim;
    non_ingeracting_Gw(Gwimp, solverDim,projimpurity_site_Hamiltonian, muAdju, NumTotInitial, delta_w);
    getoccMat(Gwimp, occMat, solverDim);


    /*SCGF2-loop*/
    int SCGFloop=0;
    double normGwimpDiffOuter=0;
    double mixingSCGF = mixing/10.;
    double mixingHF = mixing/10.;
    int  RDOccdec =0;
    int  RDOccinc =0;
    int  RDSwdec  =0;
    int  RDSwinc  =0;
    double normUnit =  MatsubaraFtnnorm(Gwimp,N_freq);
    double normOccMat=0;
    pulayMixing Mixing_SCGF(3, 100, N_freq , solverDim, solverDim );
    while( (SCGFloop<10 or(normGwimpDiffOuter > 1e-6) or normOccMat>1e-6) and SCGFloop <2000 ) {

        //transfer previous results
        for (int n=0; n<N_freq+4; n++) {
            Gwimp_in[n] = Gwimp[n];
        }
        getoccMat(Gwimp, occMat, solverDim);
        Eigen::MatrixXcd occMat_in = occMat;

        /*get occupation and Gt*/
//        getGwimp( Gwimp,projimpurity_site_Hamiltonian,Fock, Swimp_secondOrder, delta_w, muTB,solverDim) ;
        Gwimp[N_freq].setIdentity(solverDim, solverDim);
        FourierTransform( Gwimp, Gtimp, Gwimp[N_freq] );


        /*and HF*/
        getFockOperator( Fock , occMat, solverDim, projUindex, projUtensor) ;
        /*get  Sigma_w,  */
        getStimp ( Stimp, Gtimp, true, solverDim, projUindex, projUtensor);
        FT_t_to_w(Swimp_secondOrder, Stimp, N_freq);


//embedding strong Sw
        if(high_level_solver==false) {
            SecOrdDC(  Gwimp, Uindex_stronglyCorr,  Utensor_stronglyCorr ,
                       strongCorr,
                       Fock,  Swimp_secondOrder ) ;
            int strongDim = strongCorr.size();
            for (int i=0; i<strongDim; i++) {
                for (int j=0; j<strongDim; j++) {
                    for (int n=0; n<N_freq; n++) {
                        Swimp_secondOrder[n](strongCorr[i], strongCorr[j]) +=  (SE_strong.getValue(n,i,j) - SE_strong.getValue(N_freq,i,j));
                    }

                    Fock(strongCorr[i], strongCorr[j]) += SE_strong.getValue(N_freq,i,j);
                }
            }
        }



        /*obtain results*/
        getoccMat(Gwimp, occMat, solverDim);
        getGwimp( Gwimp,projimpurity_site_Hamiltonian,Fock, Swimp_secondOrder, delta_w, muTB,solverDim) ;

        /*cal. RD*/
        normOccMat = (occMat-occMat_in).norm();
        normGwimpDiffOuter= MatsubaraFtnnorm(Gwimp_in , Gwimp, N_freq)  /  normUnit;

        /*mixing, Gwimp is used in the next iteration*/
        if(SCGFloop%500==0) mixingSCGF /= 1.5;
        if(SCGFloop < 1500) Mixing_SCGF.mixing( Gwimp_in, Gwimp, mixingSCGF, SCGFloop, 1);
        else {
            for(int n=0; n<N_freq; n++) {
                Gwimp[n] =  (1-mixingSCGF) * Gwimp_in[n] + mixingSCGF * Gwimp[n];
            }
        }
        Gwimp[N_freq].setIdentity(solverDim, solverDim);
        estimate_asymto(Gwimp,2);
        estimate_asymto(Gwimp,3);

        if(stdoutLevel>=2 and mpi_rank==0) {


            std::cout <<"gt:";
            ces.compute(-Gtimp[N_tau]);
            for(int n=0; n<solverDim; n+=1) {
                std::cout <<std::fixed << std::setprecision(5)<< std::fixed   <<  ( ces.eigenvalues()[n]  ) <<" " ;
            }
            std::cout <<": " << -(Gtimp[N_tau].trace()) <<"\n";

            std::cout <<"oc:";
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> ces(solverDim);
            ces.compute(occMat);
            for(int n=0; n<solverDim; n+=1) {
                std::cout <<std::fixed << std::setprecision(5)<< std::fixed   <<  ( ces.eigenvalues()[n]  ) <<" " ;
            }
            std::cout <<": " << occMat.trace() <<"\n";

        }


        /*mixing check*/
        if (SCGFloop%10 == 0) mixing_check_outer(  normGwimpDiffOuter,  mixingSCGF);

        SCGFloop++;
    }//outer loop

    /*Write resutls*/
    ifroot std::cout << "SCGF loop:" << SCGFloop;
    ifroot printf      ("  SwRD:%e" , normGwimpDiffOuter);
    ifroot printf      ("  mixing:%e\n:" , mixingSCGF);




    /*write result to NumMatrix and Sw*/

    if(high_level_solver==false) {
        SecOrdDC(  Gwimp, Uindex_stronglyCorr,  Utensor_stronglyCorr ,
                   strongCorr,
                   Fock,  Swimp_secondOrder ) ;
    }

    writeResults(projNumMatrix, Swimp_secondOrder, Fock, Gwimp,   occMat, SE_out, Gwimp_out, solverDim);

    MPI_Barrier(MPI_COMM_WORLD);
}//SCGF2





void getoccMat(Eigen::MatrixXcd * Gwimp, Eigen::MatrixXcd & occMat, int solverDim) {


//occMat is not really occupation, but density matrix.
//i.e. occMat(a,b) = G(a,b; 0-) = -G(a,b; beta-) =   <c^\dagger_b,  c_a> =! <c^\dagger_a,  c_b>
//

    Eigen::MatrixXcd Id;
    Id.setIdentity(solverDim,solverDim);
    FourierTransform(Gwimp, occMat, beta, Id);
    occMat *= -1;



//    Eigen::MatrixXcd occMat_local(solverDim,solverDim);
//    occMat.setZero(solverDim,solverDim);
//    int  mysta, myend;
//    double tail1=0, tail2=0;
//    para_range(0,N_freq-1, mpi_numprocs, mpi_rank, &mysta, &myend);
//    for(int n=mysta; n<=myend; n++) {
//        double wn=  (2*n+1)*pi/beta;
//        occMat_local +=  (-1.*( Gwimp[n]    )   ); //where sign -= cos(wn*\beta), \sum_iwn G(iwn)*exp(iwn\beta) = \sum n>0 G(iwn)* (exp(iwn\beta) + c.c.)   // Gwimp[N_freq]/(i*wn) terms give no contribution. since G_Nf = Hermition
//        tail1 +=    1./ std::pow( wn ,2);          //sign + = (i)^2 * cos(wn*\beta)
//        tail2 +=   -1./std::pow(wn,4);             //sign - = (i)^4 * cos(wn*\beta)
//    }
//    occMat_local -= tail1 * Gwimp[N_freq+1];
//    occMat_local -= tail2 * Gwimp[N_freq+3];
//    MPI_Allreduce(occMat_local.data(), occMat.data(), solverDim*solverDim, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
//
//    occMat = (occMat + occMat.adjoint()).eval();
//    occMat /=  beta ;
//    occMat +=  ((Gwimp[N_freq])*(-0.5)  );                   // \sum 1/iw = -0.5
//    occMat +=  ((Gwimp[N_freq+1])*(beta /4.));               // 1/1+1/9+1/25+ ... = pi^2/8 ==> \sum 1/(w_n)^2  = beta^2/4
//    occMat += -((Gwimp[N_freq+3])*(std::pow(beta,3)/48.));     // \sum 1/(2n+1)^4  for n=0,1,..  = (pi^4)/96  , where - sign = i^4*cos(wn*\beta)
//    occMat *= (-1.0);
}

void getGwimp(Eigen::MatrixXcd * Gwimp, Eigen::MatrixXcd projimpurity_site_Hamiltonian,  Eigen::MatrixXcd Fock, Eigen::MatrixXcd * Swimp_secondOrder, Eigen::MatrixXcd * delta_w,  double muTB, int solverDim) {
    Eigen::MatrixXcd id;
    id.setIdentity(solverDim,solverDim);
    int mysta, myend;
    para_range(0,N_freq-1, mpi_numprocs, mpi_rank, &mysta, &myend);
    for (int n=mysta; n<=myend; n++) {
        double wn = (pi*(2.*n+1.))/beta;
        Gwimp[n] =  (I*wn + muTB)*id - projimpurity_site_Hamiltonian - delta_w[n] - Fock- Swimp_secondOrder[n] ;
        Gwimp[n] = Gwimp[n].inverse();
    }
    data_sync_EigenMat(Gwimp, 0, N_freq-1, solverDim, mpi_numprocs);
    Gwimp[N_freq].setIdentity(solverDim, solverDim);
    estimate_asymto(Gwimp,2);
    estimate_asymto(Gwimp,3);
}




void getFockOperator( Eigen::MatrixXcd & Fock, Eigen::MatrixXcd & occMat , int solverDim,
                      std::vector<Eigen::VectorXi> projUindex, std::vector<cmplx > projUtensor
                    ) {
    Eigen::MatrixXcd newFock( solverDim, solverDim);
    for(int idx=0;  idx <  projUtensor.size(); idx++) {
        // Utenser = U_{abcd} = <ab|U|cd>    //a<-c, b<-d
        int  a = projUindex[idx](0);
        int  b = projUindex[idx](1);
        int  c = projUindex[idx](2);
        int  d = projUindex[idx](3);
        newFock(b,d) +=   projUtensor[idx] * occMat(c,a); //H
        newFock(b,c) -=   projUtensor[idx] * occMat(d,a); //F
    }
    Fock = newFock;
}


//void getFockOperator( Eigen::MatrixXcd & Fock, Eigen::MatrixXcd & occMat ) {
//    Eigen::MatrixXcd newFock(NSpinOrbit_per_atom,NSpinOrbit_per_atom);
//    for(int idx=0;  idx <  Utensor.size(); idx++) {
//        if (    isOrbitalCorrinHart[Uindex[idx](0)] and
//                isOrbitalCorrinHart[Uindex[idx](1)] and
//                isOrbitalCorrinHart[Uindex[idx](2)] and
//                isOrbitalCorrinHart[Uindex[idx](3)] ) {
//            int  a = Hart2Corr[Uindex[idx](0)];
//            int  b = Hart2Corr[Uindex[idx](1)];
//            int  c = Hart2Corr[Uindex[idx](2)];
//            int  d = Hart2Corr[Uindex[idx](3)];
//            // Utenser = U_{abcd} = <ab|U|cd>    //a<-c, b<-d
//            newFock(b,d) +=   Utensor[idx] * occMat(c,a); //H
//            newFock(b,c) -=   Utensor[idx] * occMat(d,a); //F
//        }
//    }
//    Fock = newFock;
//}


//void getStimp(Eigen::MatrixXcd * Stimp, Eigen::MatrixXcd * Gtimp, bool skeleton, int solverDim  ) {
//    Eigen::MatrixXcd  newStimp;
//    int  mysta, myend;
//    para_range(0,N_tau, mpi_numprocs, mpi_rank, &mysta, &myend);
//    for(int t =mysta ; t<=myend; t++) {
//        newStimp.setZero(solverDim, solverDim);
//        for(int idx1=0;  idx1 <  Utensor.size(); idx1++) {
//            if (    isOrbitalCorrinHart[Uindex(idx1,0)] and
//                    isOrbitalCorrinHart[Uindex(idx1,1)] and
//                    isOrbitalCorrinHart[Uindex(idx1,2)] and
//                    isOrbitalCorrinHart[Uindex(idx1,3)] ) {
//                int  i = Hart2Corr[Uindex(idx1,0)];
//                int  q = Hart2Corr[Uindex(idx1,1)];
//                int  m = Hart2Corr[Uindex(idx1,2)];
//                int  k = Hart2Corr[Uindex(idx1,3)];
//
//                for(int idx2=0;  idx2 <  Utensor.size(); idx2++) {
//                    if (    isOrbitalCorrinHart[Uindex(idx2,0)] and
//                            isOrbitalCorrinHart[Uindex(idx2,1)] and
//                            isOrbitalCorrinHart[Uindex(idx2,2)] and
//                            isOrbitalCorrinHart[Uindex(idx2,3)] ) {
//                        int  l = Hart2Corr[Uindex(idx2,0)];
//                        int  n = Hart2Corr[Uindex(idx2,1)];
//                        int  p = Hart2Corr[Uindex(idx2,2)];
//                        int  j = Hart2Corr[Uindex(idx2,3)];
//                        //Utensor[idx1] = U_iqmk  = <iq|U|mk>  :    i<-m, q<-k,   Here m, k is at time t=0^-  and i,q are time t=0;
//                        //Utensor[idx2] = U_lnpj  = <ln|U|pj>  :    l<-p, n<-j
//
//                        newStimp(i,j) +=    (-1) * Gtimp[t](k,l) * Gtimp[t](m,n) *(-Gtimp[N_tau-t](p,q)) * Utensor[idx1] * Utensor[idx2];  //Bubble, sign = one-buble*Gt(-t)
//                        newStimp(i,j) +=          Gtimp[t](k,n) * Gtimp[t](m,l) * (-Gtimp[N_tau-t](p,q)) * Utensor[idx1] * Utensor[idx2];  //2ndOrder exchange, sign = Gt(-t)
//
//                    }//if
//                }//idx2
//            }//if
//        }//idx1
//        Stimp[t] = newStimp;
//    }//t
//    data_sync_EigenMat(Stimp, 0, N_tau, solverDim, mpi_numprocs);
//}

void getStimp(Eigen::MatrixXcd * Stimp, Eigen::MatrixXcd * Gtimp, bool skeleton, int solverDim,
              std::vector<Eigen::VectorXi> projUindex, std::vector<cmplx > projUtensor        ) {
    Eigen::MatrixXcd  newStimp;
    int  mysta, myend;
    para_range(0,N_tau, mpi_numprocs, mpi_rank, &mysta, &myend);
    for(int t =mysta ; t<=myend; t++) {
        newStimp.setZero(solverDim, solverDim);
        for(int idx1=0;  idx1 <  projUtensor.size(); idx1++) {
            int  i = projUindex[idx1](0);
            int  q = projUindex[idx1](1);
            int  m = projUindex[idx1](2);
            int  k = projUindex[idx1](3);

            for(int idx2=0;  idx2 <  projUtensor.size(); idx2++) {
                int  l = projUindex[idx2](0);
                int  n = projUindex[idx2](1);
                int  p = projUindex[idx2](2);
                int  j = projUindex[idx2](3);
                //Utensor[idx1] = U_iqmk  = <iq|U|mk>  :    i<-m, q<-k,   Here m, k is at time t=0^-  and i,q are time t=0;
                //Utensor[idx2] = U_lnpj  = <ln|U|pj>  :    l<-p, n<-j
                newStimp(i,j) +=    (-1) * Gtimp[t](k,l) * Gtimp[t](m,n) *(-Gtimp[N_tau-t](p,q)) * projUtensor[idx1] * projUtensor[idx2];  //Bubble, sign = one-buble*Gt(-t)
                newStimp(i,j) +=          Gtimp[t](k,n) * Gtimp[t](m,l) * (-Gtimp[N_tau-t](p,q)) * projUtensor[idx1] * projUtensor[idx2];  //2ndOrder exchange, sign = Gt(-t)

            }//idx2
        }//idx1
        Stimp[t] = newStimp;
    }//t
    data_sync_EigenMat(Stimp, 0, N_tau, solverDim, mpi_numprocs);
}

void writeResults( Eigen::MatrixXcd & projNumMatrix, Eigen::MatrixXcd * Swimp_secondOrder, Eigen::MatrixXcd Fock, Eigen::MatrixXcd * Gwimp,  Eigen::MatrixXcd occMat,
                   ImgFreqFtn & SE_out, ImgFreqFtn & Gwimp_out, int solverDim) {


    estimate_asymto(Swimp_secondOrder  ,1);
    estimate_asymto(Swimp_secondOrder  ,2);

    for(int i=0; i<solverDim; i++) {
        for(int j=0; j<solverDim; j++) {
            for(int n=0; n<N_freq; n++) {
                SE_out.setValue(n, i, j, Swimp_secondOrder[n](i,j) + Fock(i,j) );
//
                Gwimp_out.setValue(n, i, j, Gwimp[n](i,j));
            }
            SE_out.setValue(N_freq,  i,j, Fock(i,j));
            SE_out.setValue(N_freq+1,i,j, Swimp_secondOrder[N_freq+1](i,j));
            SE_out.setValue(N_freq+2,i,j, Swimp_secondOrder[N_freq+2](i,j));


            Gwimp_out.setValue(N_freq+1, i, j, Gwimp[N_freq](i,j));
            Gwimp_out.setValue(N_freq+2, i, j, Gwimp[N_freq+1](i,j));
            Gwimp_out.setValue(N_freq+3, i, j, Gwimp[N_freq+2](i,j));
            Gwimp_out.setValue(N_freq+4, i, j, Gwimp[N_freq+3](i,j));
        }
    }//n

//    ifroot std::cout << "Num ele (SOLVER, 2ndOrder) = " ;
//    for(int alp=0; alp<solverDim; alp++) {
//        ifroot std::cout << " " << real(occMat(alp,alp));
//        for(int bet=0; bet<solverDim; bet++) {
//            projNumMatrix(alp, bet) = occMat(alp,bet);
//        }
//    }
//    ifroot std::cout <<"\n";
    sleep(5);
}

void mixing_check_outer( double normGwimpDiffOuter, double & mixingSCGF) {
    static double normGwimpDiffOuter_prev = 0;
    if       ( normGwimpDiffOuter_prev > normGwimpDiffOuter ) mixingSCGF = std::min (1.07 *mixingSCGF ,   0.5);
    else if  ( normGwimpDiffOuter_prev < normGwimpDiffOuter ) mixingSCGF = std::max (1./1.09 *mixingSCGF , 1e-3);
    normGwimpDiffOuter_prev = normGwimpDiffOuter;
}

void mixing_check_inner(double  normOccDiffInner, double & mixingHF) {
    static double normOccDiffInner_prev  = 0;
    if      (normOccDiffInner < normOccDiffInner_prev)  mixingHF = std::min(1.075*mixingHF , 0.5);
    else if (normOccDiffInner > normOccDiffInner_prev)  mixingHF = std::max(mixingHF /1.09, 1e-4);
    normOccDiffInner_prev = normOccDiffInner;
}

//void paramagnetic_solution(Eigen::MatrixXcd & occMat ) {
//    if(magnetism==0) {
//        for(int n =0; n<occMat.rows() ; n+=2) {
//            for(int m =0; m<occMat.cols(); m+=2) {
//                cmplx  temp = (occMat(n,m) + occMat(n+1,m+1))/2.;
//                occMat(n  ,m  ) = temp;
//                occMat(n+1,m+1) = temp;
//                occMat(n  ,m+1) = 0;
//                occMat(n+1,m  ) = 0;
//            }
//        }
//    }
//
//}
void estimate_asymto(Eigen::MatrixXcd * FreqFtn, int order) {
    Eigen::MatrixXcd tail [3];
    Eigen::MatrixXcd tail_coeff;
    double w0,w1,w2;
    if(order==1) {
        //Asymto, i*S1/w  + ...
        w0 = 1./std::pow((2*(N_freq-3)+1)*pi/beta,2);
        w1 = 1./std::pow((2*(N_freq-2)+1)*pi/beta,2);
        w2 = 1./std::pow((2*(N_freq-1)+1)*pi/beta,2);
        tail[0] = FreqFtn[N_freq-3];
        tail[1] = FreqFtn[N_freq-2];
        tail[2] = FreqFtn[N_freq-1];
        for(int w=0; w<3; w++) {
            tail[w] = (tail[w].adjoint()-tail[w]).eval()/2.;
            tail[w] =  tail[w] / (I*(2*(N_freq-3+w)+1)*pi/beta );
        }
        tail_coeff = ( tail[0] + tail[1] + tail[2]) / (w0+w1+w2);
        FreqFtn[N_freq+0] = tail_coeff;// sign : 1/i = -i
    }
    else if(order==2) {
        //Asymto S2/(iw^2) + ...
        w0 = 1./std::pow((2*(N_freq-3)+1)*pi/beta,4);
        w1 = 1./std::pow((2*(N_freq-2)+1)*pi/beta,4);
        w2 = 1./std::pow((2*(N_freq-1)+1)*pi/beta,4);

        tail[0] = FreqFtn[N_freq-3] ;
        tail[1] = FreqFtn[N_freq-2] ;
        tail[2] = FreqFtn[N_freq-1] ;
        for(int w=0; w<3; w++) {
            tail[w] =  (tail[w]+tail[w].adjoint()).eval()/2.;
            tail[w] =  -tail[w] / std::pow((2*(N_freq-3+w)+1)*pi/beta,2);
        }
        tail_coeff =  (tail[0] + tail[1] + tail[2]) / (w0+w1+w2);
        FreqFtn[N_freq+1] = tail_coeff;
    }
    else if(order==3) {
        //Asymto S2/(iw^3) + ...
        w0 = 1./std::pow((2*(N_freq-3)+1)*pi/beta,6);
        w1 = 1./std::pow((2*(N_freq-2)+1)*pi/beta,6);
        w2 = 1./std::pow((2*(N_freq-1)+1)*pi/beta,6);

        tail[0] = FreqFtn[N_freq-3] -   (FreqFtn[N_freq] / (I*(2*(N_freq-3)+1)*pi/beta)) ;
        tail[1] = FreqFtn[N_freq-2] -   (FreqFtn[N_freq] / (I*(2*(N_freq-2)+1)*pi/beta)) ;
        tail[2] = FreqFtn[N_freq-1] -   (FreqFtn[N_freq] / (I*(2*(N_freq-1)+1)*pi/beta)) ;
        for(int w=0; w<3; w++) {
            tail[w] =  (tail[w].adjoint()-tail[w]).eval()/2.;
            tail[w] =  I*tail[w] / std::pow((2*(N_freq-3+w)+1)*pi/beta,3);  //1/i^3=-1/i = i
        }
        tail_coeff =  (tail[0] + tail[1] + tail[2]) / (w0+w1+w2);
        FreqFtn[N_freq+2] = tail_coeff;
    }
}



double MatsubaraFtnnorm(Eigen::MatrixXcd * Swimp_secondOrder, int N_freq) {
    double normValue = 0;
    for (int n=0; n<N_freq; n++) {
        double wn=  (2*n+1)*pi/beta;
        normValue += std::pow(Swimp_secondOrder[n].norm(),2)/wn;
    }
    return std::sqrt(normValue);
}
double MatsubaraFtnnorm(Eigen::MatrixXcd * Ftn1, Eigen::MatrixXcd * Ftn2, int N_freq) {
    double normValue = 0;
    for (int n=0; n<N_freq; n++) {
        double wn=  (2*n+1)*pi/beta;
        normValue += std::pow( (Ftn1[n]-Ftn2[n]).norm(),2)/wn;
    }
    return std::sqrt(normValue);
}














///*

void IPT( int solverDim,Eigen::MatrixXcd projimpurity_site_Hamiltonian,  ImgFreqFtn & weiss_field, Eigen::MatrixXcd & projNumMatrix,
          ImgFreqFtn & SE_out,      ImgFreqFtn & Gwimp_out, double muTB,
          std::vector<Eigen::VectorXi> projUindex, std::vector<cmplx > projUtensor        ) {

//See M. Ptthoff, T. Wegner, and W. Nolting, PRB, 16132 (1997)
//For more information,
//Georges and Kotliar PRB 45, 6479 (1992); (HF+2PT)
//H. Kajueter and G. Kotliar, PRL 77, 131 (1996); modification
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> ces(solverDim);
    ifroot std::cout <<"Himp (diag):";
    ces.compute(projimpurity_site_Hamiltonian);
    for(int n=0; n<solverDim; n+=1) {
        ifroot std::cout <<std::fixed << std::setprecision(3)<< std::fixed   <<  ( ces.eigenvalues()[n]  ) <<" " ;
    }
    ifroot std::cout <<"\n";

    //initial Setting
    Eigen::MatrixXcd Fock_in(solverDim, solverDim);
    Eigen::MatrixXcd Fock(solverDim, solverDim);
    Eigen::MatrixXcd occMat(solverDim, solverDim);
    Eigen::MatrixXcd occMatHF(solverDim, solverDim);
    Eigen::MatrixXcd delta_w[N_freq+2];
    Eigen::MatrixXcd Swimp_secondOrder[N_freq+3];

    Eigen::MatrixXcd Gwimp[N_freq+4];
    Eigen::MatrixXcd Gw0[N_freq+4];
    Eigen::MatrixXcd Gt0[N_tau+1];

    Eigen::MatrixXcd Stimp[N_tau+1];


    Fock_in.setZero(solverDim,solverDim);
    Fock.setZero(solverDim,solverDim);
    for(int t =0 ; t<N_tau+1; t++) {
        Gt0[t].setZero(solverDim, solverDim);
        Stimp[t].setZero(solverDim,solverDim);
    }
    for (int n=0; n<N_freq; n++) {
        Gwimp[n] = Gwimp_out.getMatrix(n);
        Gw0[n] =Gwimp_out.getMatrix(n);
//        Gw0[n].setZero(solverDim, solverDim);
        Swimp_secondOrder[n].setZero(solverDim, solverDim);
        delta_w[n] = weiss_field.getMatrix(n);
    }
    Swimp_secondOrder[N_freq].setZero(solverDim, solverDim);
    Swimp_secondOrder[N_freq+1].setZero(solverDim, solverDim);
    Swimp_secondOrder[N_freq+2].setZero(solverDim, solverDim);
    delta_w[N_freq]   = weiss_field.getMatrix(N_freq+1);
    delta_w[N_freq+1] = weiss_field.getMatrix(N_freq+2);

    Gwimp[N_freq].setIdentity(solverDim, solverDim);
    Gwimp[N_freq+1].setZero(solverDim, solverDim);
    Gwimp[N_freq+2].setZero(solverDim, solverDim);
    Gwimp[N_freq+3].setZero(solverDim, solverDim);

    Gw0[N_freq].setIdentity(solverDim, solverDim);
    Gw0[N_freq+1].setZero(solverDim, solverDim);
    Gw0[N_freq+2].setZero(solverDim, solverDim);
    Gw0[N_freq+3].setZero(solverDim, solverDim);

    for(int i = 0 ; i<  solverDim ; i++) {
        for(int j = 0 ; j<  solverDim ; j++) {
            occMat(i,j)= projNumMatrix(i,j);
        }
    }
    occMatHF = occMat;
    //
    //Initial condition: non-interacting (or HF with n>=2 loop) initial  Greens Ftn
//    double NumTotInitial = real(occMat.trace());
//    double muHF = real(projimpurity_site_Hamiltonian.trace())/ solverDim;
//    non_ingeracting_Gw(Gw0, solverDim,projimpurity_site_Hamiltonian, muHF, NumTotInitial, delta_w);

    double muHF = muTB;
    getFockOperator( Fock , occMatHF, solverDim, projUindex, projUtensor) ;
    getGwimp( Gw0, projimpurity_site_Hamiltonian, Fock, Swimp_secondOrder, delta_w, muHF, solverDim) ;

// HF

    pulayMixing Mixing_HF(3, 10, 1, solverDim, solverDim );
    double HFloop=0;
    double HFmixing=0.05;
    Fock_in=Fock;
    do {
        HFloop++;
        //mixing, Gwimp is used in the next iteration
        if(HFloop < 1500) Mixing_HF.mixing( &Fock_in, &Fock, HFmixing, HFloop, 1);
        else {
            Fock =  (1-HFmixing) * Fock_in + HFmixing * Fock;
        }
        Fock_in = Fock;
        getGwimp( Gw0, projimpurity_site_Hamiltonian,Fock, Swimp_secondOrder, delta_w, muTB, solverDim) ; //non-interacting (or HF with n>=2 loop)  Greens Ftn
        getoccMat(Gw0, occMatHF, solverDim);
        if(magnetism==0) {
            for (int n=0; n<solverDim; n+=2) {
                int up = n, dn =n+1;
                cmplx avgN = (occMatHF(up,up)+occMatHF(dn,dn))/2.0;
                occMatHF(up,up)=avgN;
                occMatHF(dn,dn)=avgN;
                occMatHF(dn,up)=0.0;
                occMatHF(up,dn)=0.0;

            }
        }

        /*
        ///////////////////////////////////////////////////////////////////////////////////////////////
        ///Adjust Chem.Pot
                double Nele[history];
                for (int i=0; i<history; i++)  Nele[i]=-1;
                Nele[0] = real(occMatHF.trace());
                double muRD;
                double mu_next;
                double muUB=0, muLB=0;
                int nearAns=0;
                double dmu=1;
                muHF-=dmu;
                int step=0;
                while (  std::abs(Nele[0] - NumTotInitial) > std::min(1e-5, 1e-5/((Fock.norm())+1e-5)  )  ) { //   and  std::abs(dmu)>1e-5 ) {}
                    muHF += dmu;

                    getGwimp( Gw0,projimpurity_site_Hamiltonian,Fock, Swimp_secondOrder, delta_w,  muHF, solverDim) ;
                    estimate_asymto(Gw0,2);
                    estimate_asymto(Gw0,3);
                    getoccMat(Gw0, occMatHF,solverDim);

                    for(int i=history-1; i>0 ; i--) Nele[i]=Nele[i-1];
                    Nele[0] = real(occMatHF.trace());

                    if (nearAns ==0) {
                        if ( Nele[0] < NumTotInitial) {
                            muLB = muHF;
                            if (Nele[1]>NumTotInitial && Nele[1] >0 ) nearAns = 1;
                            else dmu = fabs(dmu);
                        }
                        if ( Nele[0] > NumTotInitial) {
                            muUB = muHF;
                            if (Nele[1] < NumTotInitial && Nele[1] > 0)  nearAns =1;
                            else dmu = -1*fabs(dmu);
                        }
                        if ( Nele[0] == NumTotInitial) break;
        //                assert( std::abs(muHF) <  100)   ;
                    }
                    else if (nearAns ==1) {
                        if (Nele[0] > NumTotInitial) muUB = muHF;
                        else if(Nele[0] < NumTotInitial) muLB =muHF;
                        mu_next = (muUB+muLB)/2.;
                        dmu = mu_next - muHF;
                    }
                }//while
        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        //// */
        getFockOperator( Fock , occMatHF, solverDim, projUindex, projUtensor) ;
        ifroot std::cout << "HF_iter: "<<HFloop<<"\n occMatHF:\n " << occMatHF << "\n Fock_RD:\n" << (Fock_in-Fock).norm() << "\n\n";
    } while ( (Fock_in-Fock).norm() > 1e-4   ) ;


    getFockOperator( Fock , occMatHF, solverDim, projUindex, projUtensor) ;
    FourierTransform( Gw0, Gt0, Gw0[N_freq]);
    getStimp ( Stimp, Gt0, false, solverDim, projUindex, projUtensor);
    FT_t_to_w(Swimp_secondOrder, Stimp, N_freq);


    double a[2];
    double b[2];
    double c[2];
    a[0] = real((UHubb*UHubb*occMat(1,1) *(1.0-occMat(1,1))) );
    a[1] = real((UHubb*UHubb*occMat(0,0) *(1.0-occMat(0,0))) );
    b[0] = real( (muHF- muTB +UHubb*(1.0-2.0*occMat(1,1))))  ;
    b[1] = real( (muHF- muTB +UHubb*(1.0-2.0*occMat(0,0))))  ;
    c[0] = real((UHubb*UHubb*occMatHF(1,1) *(1.0-occMatHF(1,1))) );
    c[1] = real((UHubb*UHubb*occMatHF(0,0) *(1.0-occMatHF(0,0))) );
    ifroot std::cout <<"b:" << b[0] <<" " << b[1] <<"\n";
    for (int n=0; n<N_freq; n++) {
        Swimp_secondOrder[n](0,0) = a[0]*Swimp_secondOrder[n](0,0) / (c[0]- b[0]*Swimp_secondOrder[n](0,0) );
        Swimp_secondOrder[n](1,1) = a[1]*Swimp_secondOrder[n](1,1) / (c[1]- b[1]*Swimp_secondOrder[n](1,1) );
    }

    getGwimp( Gwimp, projimpurity_site_Hamiltonian , Fock, Swimp_secondOrder, delta_w, muTB, solverDim) ;
    estimate_asymto(Gwimp,2);
    estimate_asymto(Gwimp,3);


    //write result to NumMatrix and Sw
    writeResults(projNumMatrix,  Swimp_secondOrder, Fock, Gwimp,   occMat, SE_out, Gwimp_out, solverDim);
    MPI_Barrier(MPI_COMM_WORLD);
}
// */
