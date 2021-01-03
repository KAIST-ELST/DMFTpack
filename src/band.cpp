#include "tight_common.h"
#include "TB.h"
#include <Eigen/Eigenvalues>
//Non-interaction Hamiltonian band canculation
//160128 Jae-Hoon Sim, KAIST
void band(Eigen::VectorXd  *KS_eigenEnergy, double muDFT, int knum) {

    ifroot printf("********************************\n");
    ifroot printf("             Band  \n");
    ifroot printf("********************************\n");

    double * temp  = new double           [knum*NumOrbit];
    double * band  = new double [knum_mpiGlobal*NumOrbit];
    for (int k=0; k< knum; k++) {
        for(int j=0; j<NumOrbit; j++) {
            temp[k*NumOrbit+j] = KS_eigenEnergy[k][j]-muDFT;
        }
    }

    int ircnt[mpi_numprocs], idisp[mpi_numprocs];
    for (int itsRank=0; itsRank<mpi_numprocs; itsRank++) {
        int itskend, itsksta;
        para_range(0,knum_mpiGlobal-1, mpi_numprocs,itsRank, &itsksta, &itskend);
        ircnt[itsRank] =  (itskend-itsksta+1)*NumOrbit;
        idisp[itsRank]  = itsksta*NumOrbit;
    }
    MPI_Gatherv (temp, knum*NumOrbit, MPI_DOUBLE, band, ircnt, idisp, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    FILE *datap1;
    ifroot datap1 = fopen("band.dat", "w");
    /*band*/
    if(mpi_rank==0) {
        for(int j=0; j<NumOrbit; j++) {
            for (int k=0; k< knum_mpiGlobal; k++) {
                fprintf(datap1, "%d  %0.8f    %0.8f\n",k, kdist_band[k],  band[k*NumOrbit+j] );
//                fprintf(datap1, "%0.8f    %0.8f\n" , k*0.1,  band[k][j] );
            }
            fprintf(datap1, "\n");
            fflush(datap1);
        }
    }/*band*/
    ifroot fclose(datap1);
    delete [] temp;
    delete [] band;
}
void band(    std::vector<Eigen::MatrixXcd> & H_k_inModelSpace, double mu, int knum) {


    ifroot printf("********************************\n");
    ifroot printf("             Band  \n");
    ifroot printf("********************************\n");

    int minNBAND=999;
    for(int k = 0;  k < knum; k++) {
        if( minNBAND > NBAND[k]) minNBAND = NBAND[k];
    }

    double Model_eigenEnergy[knum][minNBAND];
    for(int k = 0;  k < knum; k++) {
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> ces(NBAND[k]);
        ces.compute(H_k_inModelSpace[k],0);
        for(int i =0; i<minNBAND; i++) {
            Model_eigenEnergy[k][i] = (ces.eigenvalues()[i]) ;  //(V[0][0], V[1][0], V[2][0],..V[n][0]) = 0th eigenvector
        }
        if (mpi_rank ==0 and k==0  )          std::cout << "Effective Hamiltonian in model space = H_DFT - doublecounting + Sw_Hartree\n";
    }


    //double temp[knum][minNBAND];
    //double band[knum_mpiGlobal][minNBAND];
    double * temp  = new double           [knum*minNBAND];
    double * band  = new double [knum_mpiGlobal*minNBAND];
    for (int k=0; k< knum; k++) {
        for(int j=0; j<minNBAND; j++) {
            temp[k*minNBAND+j] = Model_eigenEnergy[k][j]-mu;
        }
    }

    int ircnt[mpi_numprocs], idisp[mpi_numprocs];
    for (int itsRank=0; itsRank<mpi_numprocs; itsRank++) {
        int itskend, itsksta;
        para_range(0,knum_mpiGlobal-1, mpi_numprocs,itsRank, &itsksta, &itskend);
        ircnt[itsRank] =  (itskend-itsksta+1)*minNBAND;
        idisp[itsRank]  = itsksta*minNBAND;
    }
    MPI_Gatherv (temp, knum*minNBAND, MPI_DOUBLE, band, ircnt, idisp, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    FILE *datap1;
    ifroot datap1 = fopen("band.dat", "w");
    /*band*/
    for(int j=0; j<minNBAND; j++) {
        if(mpi_rank==0) {
            for (int k=0; k< knum_mpiGlobal; k++) {
                fprintf(datap1, "%d  %0.8f    %0.8f\n", k, kdist_band[k],  band[k*minNBAND+j] );
            }
            fprintf(datap1, "\n");
            fflush(datap1);
        }
    }/*band*/
    ifroot fclose(datap1);
    delete [] temp;
    delete [] band;
}
