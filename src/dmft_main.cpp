#if SLEPc_ENABLED_dmft
#include <slepceps.h>
#endif

#include <fstream>
#include <cstring>
#include <stdlib.h>
#include <time.h>
#include "ImgTimeFtn.h"
#include "tight_common.h"
#include "pulay.h"


#include <Eigen/Eigenvalues>


//#define debug_jhs 1

bool dmft_scf_check ( Eigen::MatrixXcd NumMatrixLatt, Eigen::MatrixXcd NumMatrixImp, time_t timeStartIt, time_t timeEndIt, int currentIt ) ;

double  TightBinding(double mu, const std::string &hamiltonian, ImgFreqFtn & SelfE_w,
                     ImgFreqFtn & weiss_fieldTB, ImgFreqFtn & weiss_fieldTBCorr,
                     int mu_adjustTB,  Eigen::MatrixXcd & SolverBasis
                    );



//void analysis_example(std::string scfout_file);


void weak_solver(
    std::string SOLVERtype,
    int solverDim,
    ImgFreqFtn & SE_out,  ImgFreqFtn & Gwimp_in_out, ImgFreqFtn & GwHF,
    std::vector<Eigen::VectorXi> projUindex, std::vector<cmplx > projUtensor
);


void SOLVER(
    std::string SOLVERtype,
    int solverDim, int solver_block,  std::vector<int> impurityOrbit,
    Eigen::MatrixXcd  impurity_site_Hamiltonian,  ImgFreqFtn & weiss_field,
    ImgFreqFtn & SE_out,  ImgFreqFtn & Gwimp_in_out,
    std::vector<Eigen::VectorXi> projUindex, std::vector<cmplx > projUtensor,
    Eigen::MatrixXcd  Sw_doublecounting,            std::vector<Eigen::MatrixXcd >  dc_weakCorr,
    double muTB
);

void Time_print();
//MPI_Comm localComm;
void write_gnuplot(int writingFile) ;
void write_results(int DFTIt, int currentIt, std::string system_name, int  NumCluster);

double muTB;
double muDFT;

void dc_for_dmft( Eigen::MatrixXcd  & Sw_doublecounting, std::vector<Eigen::VectorXi> & Uindex, std::vector<cmplx >  & Utensor, Eigen::MatrixXcd & NelectronDFT,
                  int NumCorrAtom, int initial_dc, ImgFreqFtn & SelfEnergy_w) ;
//void set_HFself_energy( Eigen::MatrixXcd & NumMatrix, Eigen::MatrixXi & Uindex, Eigen::VectorXcd & Utensor, ImgFreqFtn & SelfEnergy_w   );
//void set_Hartree_self_energy( Eigen::MatrixXcd & NumMatrix, std::vector<Eigen::VectorXi >  Uindex, std::vector<cmplx >  Utensor,ImgFreqFtn & SelfEnergy_w   );




int main(int argc, char *argv[]) {
/////////////////////////////////////////////
//set mpi_env
/////////////////////////////////////////////



    char chr[100];
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Get_processor_name(processor_name, &mpi_namelen);
    printf("Process %d on %s out of %d\n", mpi_rank, processor_name, mpi_numprocs);

//    int key, color=0, s=0;
//    while(processor_name[s]) {
//        if(
//            processor_name[s] == '0' or
//            processor_name[s] == '1' or
//            processor_name[s] == '2' or
//            processor_name[s] == '3' or
//            processor_name[s] == '4' or
//            processor_name[s] == '5' or
//            processor_name[s] == '6' or
//            processor_name[s] == '7' or
//            processor_name[s] == '8' or
//            processor_name[s] == '9'
//        )
//            color += ((int)(processor_name[s]-'0'))+color*10;
//        s++;
//    }
//    key = mpi_rank;
//    MPI_Comm_split(MPI_COMM_WORLD, color, key,  &localComm);
//    MPI_Comm_size(localComm, &node_local_numprocs);
//    MPI_Comm_rank(localComm, &node_local_rank);
//    printf("Local node Process %d on %s (%d) out of %d\n", node_local_rank, processor_name,color, node_local_numprocs);


/////////////////////////////////////////////
//read_inputFile & Initializing
/////////////////////////////////////////////

//std::string hamiltonian = std::string(argv[1]);
    read_inputFile( std::string("Hk.HWR") );
//    read_inputFile( std::string(argv[1]) );
//    read_inputFile( std::string(argv[1]) );
//    if(H0_from_OpenMX!=0  ) {
//        ifroot analysis_example(std::string(argv[1]));
//        MPI_Barrier(MPI_COMM_WORLD);
//        sleep(5);
//    }

    DFTIt =0;
    if (DFTIt >0 and restart ==0 ) {
        restart =1;
        ifroot std::cout <<"Restart set to 1\n";
    }


    ifroot std::cout << "This is the end of the reading input parameters...\n";
    std::cout << std::fixed << std::setprecision(4);
    std::cout << std::showpos ;


    Eigen::MatrixXcd SolverBasis;
    SolverBasis.setIdentity(NumCorrAtom*N_peratom_HartrOrbit, NumCorrAtom*N_peratom_HartrOrbit);
    ImgFreqFtn weiss_field_weakCorr(0);
    ImgFreqFtn weiss_field_strongCorr(0);
    ImgFreqFtn weiss_fieldTB_weakCorr(0);
    ImgFreqFtn weiss_fieldTB_strongCorr(0);
    ImgFreqFtn SelfEnergy_w(0);
    ImgFreqFtn SelfEnergy_E(0);

    int mystaE, myendE;
    para_range(0,Spectral_EnergyGrid-1, mpi_numprocs, mpi_rank, &mystaE, &myendE);

    weiss_fieldTB_weakCorr.Initialize(  beta, N_freq, NumHartrOrbit_per_cluster,NumCluster, mixingType);
    weiss_field_weakCorr.Initialize(    beta, N_freq, NumHartrOrbit_per_cluster,NumCluster, mixingType);
    weiss_fieldTB_strongCorr.Initialize(beta, N_freq, NSpinOrbit_per_atom,NumCorrAtom, mixingType);
    weiss_field_strongCorr.Initialize(  beta, N_freq, NSpinOrbit_per_atom,NumCorrAtom, mixingType);




    Eigen::MatrixXcd  Sw_doublecounting;
    ifroot std::cout << "This is the end of the reading input parameters...\n";

/////////////////////////////////////////////
//////File arrange
/////////////////////////////////////////////
//    /* Arrange restart files.. */
    if(mpi_rank==0) {
        std::string cp_comm;
        std::ostringstream mkdir_comm;
        if(restart == 0) {
            system("mkdir Restart");
            system("rm mu_history.out DM_DFT_DMFT.HWR");
            system("rm dual_DM_direct.HWR");
            std::cout << "\nMaking work directory\n";
        }
        else {
            sprintf(chr,"%d",restart);
            system_name = system_name + "_restart" +chr;
            system("cp ./Restart/*  ./");
            std::cout << "\nRestart...";
        }


        mkdir_comm <<"mkdir "<<  system_name;
        system(mkdir_comm.str().c_str());

        cp_comm = std::string("mkdir Inputfiles");
        system(cp_comm.c_str());

        cp_comm = std::string("cp   ")+  "Hk.HWR input.parm   input.solver  OverlapMatrix*.HWR  PARAMS ./Inputfiles"   ;
        system(cp_comm.c_str());

        FILE * betaWrite = fopen("./Restart/beta.dat","w");
        fprintf(betaWrite, "%e\n", beta);
        fclose(betaWrite);
        sleep(1);

        cp_comm = std::string("cp -r ")+ "./Restart"    + " ./"+system_name+ "/";
        system(cp_comm.c_str());

        if(SOLVERtype!=std::string("TB")) {
            ifroot write_gnuplot(0);
        }
    } //mpi_rank==0

    ifroot std::cout <<"Files are arranged..\n";




/////////////////////////////////////////////////////
//  occupation, muDFT
////////////////////////////////////////////////////

    /*DFT results and double counting */
    muDFT = 0.0;
    muDFT    =  TightBinding (muDFT, std::string("Hk.HWR"), SelfEnergy_w,   weiss_fieldTB_weakCorr, weiss_fieldTB_strongCorr,  -1, SolverBasis);
    muTB = muDFT;
    ifroot std::cout << "First run, Initial Chemical potential, we have muDFT: " << muDFT <<"\n";
    FILE *tempFile;   //mu.dat
    tempFile = fopen("muDFT.out", "w");
    fprintf(tempFile, "%0.20f\n", muDFT);
    fclose(tempFile) ;

/////////////////////////////////////////////////////
//  Initial Self energy, double counting,  and Sw_inf and muTB
////////////////////////////////////////////////////
    /*retarded self-energy, S_E, for quasi-ptl calculation */
    if(SOLVERtype !=std::string("TB")) {
        SelfEnergy_w.Initialize(         beta, N_freq,  NumHartrOrbit_per_cluster, NumCluster, mixingType);
        if(restart!=0 and N_peratom_HartrOrbit >0) {
            /*Matsubara self-energy, S_w*/
            ifroot std::cout <<"Reading Self-energy...\n";
            double beta_prev;
            std::ifstream inputbeta("./beta.dat");
            inputbeta >> beta_prev;
            SelfEnergy_w.read_full(std::string("./Restart/Sw_SOLVER.full.dat"),beta, beta_prev);
        }
        SelfEnergy_w.dataOut(std::string("InitialSw.dat"));
    }
    else if(SOLVERtype== std::string("TB"))
        SelfEnergy_E.realFreq(       E0,  real(dE), (myendE-mystaE+1),  NumHartrOrbit_per_cluster, NumCluster, 0, mystaE    );

    if (SOLVERtype == std::string("TB") and mode!=std::string("band") and mode!=std::string("dos")  ) {
        if (restart<=0 ) {
            std::cout << "Please, use restart option\n";
            exit(1);
        }
        if(restart!=0) {
            /*occupation matrix*/
            std::ifstream  input(std::string("./Restart/Numele.dat").c_str());
            double reN, imN;
            for(int n=0; n< NumCluster*NumHartrOrbit_per_cluster; n++) {
                int i;
                input >> i;
                assert(n==i);
                for(int m=0; m< NumCluster*NumHartrOrbit_per_cluster; m++) {
                    input >>  reN;
                    input >>  imN;
                    NumMatrix(n,m) = reN + I*imN ;
                }
            }
            FILE * Chem = fopen("./Restart/mu_history.out","r");
            while(!feof(Chem)) {
                fscanf(Chem, "%lf\n",&muTB);
            }
            fclose(Chem);
        }//restart
        rewrite_retardedSw();

        MPI_Barrier(MPI_COMM_WORLD);
        if( N_peratom_HartrOrbit > 0) SelfEnergy_E.update_full(std::string("realFreq_Sw.full.dat"),1);
        if( N_peratom_HartrOrbit > 0) SelfEnergy_E.dataOut_full_pararell(std::string("realFreq_Sw.full.dat_"));
    }

    ifroot std::cout << "Initial Self-energy was constructed," << mpi_rank <<"\n";

/////////////////////////////////////////////
//single TB run  (band, dos...)
/////////////////////////////////////////////
    if (SOLVERtype==std::string("TB")) {
        ifroot std::cout << "Here we start TB.\n";
        muTB = TightBinding (muTB, std::string("Hk.HWR"), SelfEnergy_E, weiss_fieldTB_weakCorr, weiss_fieldTB_strongCorr, 0, SolverBasis);
        ifroot std::cout << "Solver=0 loop finished\n";
        MPI_Finalize();
        return 0;
    }

    /*******************************************
    /////////////////////////////////////////////
    //Charge density loop//
    ////////////////////////////////////////////////////////////////////////////////////////////////
    *******************************************/

    DFTIt++;
    ifroot     std::cout <<"\n#############################################\n";
    ifroot     std::cout <<"         DFT, Charge iteration (mu) : "<<DFTIt <<"\n";
    ifroot     std::cout <<  "#############################################\n";

    /*******************************************************************************************************
    OPENMX HK construction HERE
    //        Sw_Hartree_new.setZero(NumHartrOrbit_per_cluster*NumCluster, NumHartrOrbit_per_cluster*NumCluster);
    //        Sw_doublecounting.setZero(N_peratom_HartrOrbit*NumCorrAtom, N_peratom_HartrOrbit*NumCorrAtom);
    //        ifroot std::cout <<"Find double counting from Hk with sigma=0\n";
    //        muTB = TightBinding (muTB, argv[1], SelfEnergy_zero,  Sw_doublecounting, Sw_Hartree_new   , weiss_fieldTB,   1);
    //        NelectronDFT = NumMatrix;
    //        if(doublecounting==1) dc_for_dmft( Sw_doublecounting,  Uindex,  Utensor,  NelectronDFT);
    *******************************************************************************************************/


    muTB = TightBinding (muTB, std::string("Hk.HWR"), SelfEnergy_w,  weiss_fieldTB_weakCorr, weiss_fieldTB_strongCorr, 1, SolverBasis);
    dc_for_dmft( Sw_doublecounting,  Uindex,  Utensor, NumMatrix, NumCorrAtom, 0, SelfEnergy_w);

    ifroot {
        for(int at=0; at<NumCorrAtom; at++) {
            std::cout << "Num ele (decomp, Initial) atom"<<at<<" = ";
            for(int n=0; n<N_peratom_HartrOrbit; n+=1) {
                int n0 = at* N_peratom_HartrOrbit + n;
                std::cout <<std::fixed << std::setprecision(4)<< std::fixed <<    real(NumMatrix(n0,n0)) <<" " ;
            }
            std::cout <<"\n";
        }
//        for(int at=0; at<NumCorrAtom; at++) {
//            Eigen::MatrixXcd NumMat_atom1(N_peratom_HartrOrbit, N_peratom_HartrOrbit);
//            for(int n=0; n<N_peratom_HartrOrbit; n+=1) {
//                for(int m=0; m<N_peratom_HartrOrbit; m+=1) {
//                    NumMat_atom1(n,m) = NumMatrix(at*N_peratom_HartrOrbit+n, at*N_peratom_HartrOrbit+m);
//                }
//            }
//        }
        for(int at=0; at<NumCorrAtom; at++) {
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> ces(N_peratom_HartrOrbit);
            Eigen::MatrixXcd NumMat_atom1(N_peratom_HartrOrbit, N_peratom_HartrOrbit);
            for(int n=0; n<N_peratom_HartrOrbit; n+=1) {
                for(int m=0; m<N_peratom_HartrOrbit; m+=1) {
                    NumMat_atom1(n,m) = NumMatrix(at*N_peratom_HartrOrbit+n, at*N_peratom_HartrOrbit+m);
                }
            }
            ces.compute(NumMat_atom1);
            std::cout << "Num ele (diagon, Initial) atom"<<at<<" = ";
            for(int n=0; n<N_peratom_HartrOrbit; n+=1) {
                std::cout <<std::fixed << std::setprecision(4)<< std::fixed   <<  ( ces.eigenvalues()[n]  ) <<" " ;
            }
            std::cout <<": " << NumMat_atom1.trace() <<"\n";
        }
    }

    /*************************************************
     DMFT SC loop
    ***************************************************/
    int currentIt =1;

    pulayMixing mixing_self_energy(2, 5, N_freq, NumHartrOrbit_per_cluster*NumCluster, NumHartrOrbit_per_cluster*NumCluster, true );
    while(true) {
        Eigen::MatrixXcd NumMatrixImp;
        Eigen::MatrixXcd NumMatrixLatt;
        time_t timeStartIt, timeEndIt;
        timeStartIt = clock();
        ifroot     std::cout <<"\n#########################################\n";
        ifroot     std::cout <<"    DFT : "<<DFTIt <<"    DMFT : SCF  "<< currentIt <<"\n";
        ifroot     std::cout <<"#########################################\n";

        on_the_fly_control();
        if(  (currentIt >=maxDmftIt/2 and DFTIt != 1) or (currentIt >=(maxDmftIt) and DFTIt == 1)   ) {
            ifroot std::cout <<"DMFT: The maximum number of iterations has been reached." << mpi_rank <<"\n";
            break;
        }

        if(mixingFtn==0 and currentIt != 1) {
            weiss_field_weakCorr.update(weiss_fieldTB_weakCorr,      mixing, 0, mixingType);   //alps, diag selfenergy update
            weiss_field_strongCorr.update(weiss_fieldTB_strongCorr,  mixing, 0, mixingType);   //alps, diag selfenergy update
        }
        else {
            weiss_field_weakCorr.update(weiss_fieldTB_weakCorr,      1,0, 0);
            weiss_field_strongCorr.update(weiss_fieldTB_strongCorr,  1,0, 0);
        }


//        weiss_field_strongCorr.dataOut(std::string("delta_w.dat"));
//        weiss_field_strongCorr.dataOut_full(std::string("delta_w.full.dat"));
//        ifroot std::cout << "FILEOUT:delta_w.dat\n" ;


        /***********************************************************
        run impurity
        in  = delta_w, delta_t, mu, onsite_Hamiltonian, SelfEnergy_w(initial) ;
        out = SelfEnergy_w, GwImp_DMFT, NumMat
        note : Self-energy = \Sigma_Hartree + \Simga_imp -\Simga_dc
        *************************************************************/

        ImgFreqFtn SelfEnergy_w_out(0);
        SelfEnergy_w_out.Initialize(         beta, N_freq, NumHartrOrbit_per_cluster, NumCluster, mixingType);
        for(int cl=0; cl < NumCluster;  cl++) {
            ImgFreqFtn SE_out(0);
            ImgFreqFtn SE_lowlevel(0);
            ImgFreqFtn SE_strong(0);

            SE_out.Initialize(         beta, N_freq,  NumHartrOrbit_per_cluster, 1, mixingType);
            SE_lowlevel.Initialize(beta, N_freq, NumHartrOrbit_per_cluster,1,0);
            SE_strong.Initialize(beta,   N_freq, NSpinOrbit_per_atom,      1,0);

            std::vector<int> strongCorr(NSpinOrbit_per_atom);

            ImgFreqFtn Gw_weak(beta, N_freq, NumHartrOrbit_per_cluster, 1,0);
            Gw_weak.update_full(std::string("Gw_loc.full.dat") + intToString(cl),1);


            ImgFreqFtn GwHF_weak(beta, N_freq, NumHartrOrbit_per_cluster, 1,0);
            GwHF_weak.update_full(std::string("GwHF_loc.full.dat") + intToString(cl),1);


            ImgFreqFtn Gw_strong(0);
            Gw_strong.Initialize(beta, N_freq, NSpinOrbit_per_atom,1,0);

            ImgFreqFtn GwHF_strong(0);
            GwHF_strong.Initialize(beta, N_freq, NSpinOrbit_per_atom,1,0);
            for(int at=cl*NumAtom_per_cluster; at < (cl+1)*NumAtom_per_cluster; at++) {
                for (int n=0; n<N_freq+3; n++) {
                    for (int i=0; i<NSpinOrbit_per_atom; i++) {
                        for (int j=0; j<NSpinOrbit_per_atom; j++) {
                            int iF = CorrToHartr(at,i)  - cl * NumHartrOrbit_per_cluster ;
                            int jF = CorrToHartr(at,j)  - cl * NumHartrOrbit_per_cluster ;
                            Gw_strong.setValue(n,i,j, Gw_weak.getValue(n,iF,jF));
                            GwHF_strong.setValue(n,i,j, GwHF_weak.getValue(n,iF,jF));
                        }
                    }
                }
            }

            /*Weakly correlated subspace in the impurity site*/
            if(NumHartrOrbit_per_cluster>NSpinOrbit_per_atom) {
                ifroot std::cout <<"Run low level solver\n";
                weak_solver(Lowlevel_SOLVERtype,
                            NumHartrOrbit_per_cluster,
                            SE_lowlevel, Gw_weak, GwHF_weak,
                            Uindex,Utensor);

                SE_lowlevel.dataOut(std::string("Bare_Sw_lowlevel.dat") + intToString(cl));

                /*remove double counting contribution*/
                for(int at=cl*NumAtom_per_cluster; at < (cl+1)*NumAtom_per_cluster; at++) {
                    ifroot std::cout <<"Run low level solver for active space\n";
                    ImgFreqFtn SE_lowlevel_local(0);
                    SE_lowlevel_local.Initialize(beta, N_freq, NumHartrOrbit_per_cluster,1,0);
                    weak_solver(Lowlevel_SOLVERtype,
                                NSpinOrbit_per_atom,
                                SE_lowlevel_local, Gw_strong, GwHF_strong,
                                Uindex_stronglyCorr,Utensor_stronglyCorr);

                    for (int n=0; n<N_freq+1; n++) {
                        for (int i=0; i<NSpinOrbit_per_atom; i++) {
                            for (int j=0; j<NSpinOrbit_per_atom; j++) {
                                int iF = CorrToHartr(at,i) - cl*NumHartrOrbit_per_cluster;
                                int jF = CorrToHartr(at,j) - cl*NumHartrOrbit_per_cluster;
                                SE_lowlevel.setValue(n,iF,jF,
                                                     SE_lowlevel.getValue(n,iF,jF)-  SE_lowlevel_local.getValue(n,i,j)   );
                            }
                        }
                    }
                }//at
                for (int n=0; n<N_freq+2; n++) {
                    SE_out.setMatrix(n,SE_lowlevel.getMatrix(n));
                }

            }//if, weak solver

            //strong solver
            for(int at=cl*NumAtom_per_cluster; at < (cl+1)*NumAtom_per_cluster; at++) {

                //double counting
                std::vector<Eigen::MatrixXcd >  dc_weakCorr(N_freq+1);
                for (int n=0; n<N_freq+1; n++) {
                    dc_weakCorr[n].setZero(NSpinOrbit_per_atom,NSpinOrbit_per_atom);
                    for (int i=0; i<NSpinOrbit_per_atom; i++) {
                        for (int j=0; j<NSpinOrbit_per_atom; j++) {
                            int iF = CorrToHartr(at, i)  - cl * NumHartrOrbit_per_cluster;
                            int jF = CorrToHartr(at, j)  - cl * NumHartrOrbit_per_cluster;
                            dc_weakCorr[n](i,j) = SE_lowlevel.getValue(n,iF,jF)   ;
                        }
                    }
                }
                for (int n=0; n<N_freq; n++) {
                    dc_weakCorr[n] -= dc_weakCorr[N_freq];
                }

                //strongly correlated subspace
                for (int i=0; i<NSpinOrbit_per_atom; i++)  strongCorr[i] = CorrToHartr(at, i ) ;
                ifroot std::cout << "High-level solver for atom " << at <<"\n";

                SOLVER(SOLVERtype,
                       NSpinOrbit_per_atom, at, strongCorr,
                       impurity_site_Hamiltonian, weiss_field_strongCorr,
                       SE_strong, Gw_strong,
                       Uindex_stronglyCorr, Utensor_stronglyCorr,
                       Sw_doublecounting, dc_weakCorr,
                       muTB);



                for (int n=0; n<N_freq+3; n++) {
                    for (int i=0; i<NSpinOrbit_per_atom; i++) {
                        for (int j=0; j<NSpinOrbit_per_atom; j++) {
                            int iF = CorrToHartr(at,i) - cl*NumHartrOrbit_per_cluster;
                            int jF = CorrToHartr(at,j) - cl*NumHartrOrbit_per_cluster;
                            Gw_weak.setValue(n,iF,jF, Gw_strong.getValue(n,i,j));   // G^-1 = G_weak^-1 + G_strong^-1
                        }
                    }
                }

                /*combine all result (Sw_weak, Sw_strong, Sw_DFTdc)  to SE_out*/
                for (int n=0; n<N_freq+2; n++) {
                    for (int i=0; i<NSpinOrbit_per_atom; i++) {
                        for (int j=0; j<NSpinOrbit_per_atom; j++) {
                            int iF = CorrToHartr(at,i) - cl *NumHartrOrbit_per_cluster  ;
                            int jF = CorrToHartr(at,j) - cl *NumHartrOrbit_per_cluster  ;
                            SE_out.setValue(n,iF,jF,        SE_lowlevel.getValue(n,iF,jF)  + SE_strong.getValue(n,i,j));
                        }
                    }
                }
                SE_strong.dataOut(std::string("Bare_Sw_strong.dat")+intToString(at));
            }//at

            Gw_weak.dataOut_full((std::string("Gw_imp.full.dat")+ intToString(cl)  ) );
            Gw_weak.dataOut(std::string("Gw_imp.dat"+ intToString(cl)));

            SE_out.dataOut(std::string("Bare_Sw_SOLVER.dat") + intToString(cl));
            ifroot std::cout << "FILEOUT:Bare_Sw_SOLVER.dat for cluster " << cl <<"\n" ;


            if(doublecounting==1) {
                for (int  w=0; w<N_freq+1; w++) {
                    SE_out.setMatrix(w,
                                     SE_out.getMatrix(w) - Sw_doublecounting.block(cl*NumHartrOrbit_per_cluster, cl*NumHartrOrbit_per_cluster, NumHartrOrbit_per_cluster, NumHartrOrbit_per_cluster  ));
                }
            }
            SelfEnergy_w_out.update(     SE_out,   1,    cl, 0 );
        }//cl, solver

        //self-energy mixing.
        if(mixingFtn==1)   {
            mixing_self_energy.mixing ( SelfEnergy_w, SelfEnergy_w_out, mixing, currentIt, 1);
            MPI_Barrier(MPI_COMM_WORLD);
            SelfEnergy_w.update(SelfEnergy_w_out,1.0);
        }


        SelfEnergy_w.dataOut_full(std::string("Sw_SOLVER.full.dat"));
        SelfEnergy_w.dataOut(std::string("Sw_SOLVER.dat"));
        NumMatrixImp = NumMatrix;



        /***********************************************************
        run tight binding
        in= SelfEnergy_w;
        out = delta_w, delta_t, mu, onsite_Hamiltonian, GwLoc_DMFT
        Sw = Sw_LDA - Sw_DC  + Sw_HF + Sw_imp  =  Sw_static + Sw_imp  . Here  we belive Sw_LDA \approx Sw_DC
        Then, H + Sw =  (H+Sw_static) + Sw_imp == HR + Sw_imp   (see TB.cpp),
        where HR = H_DFT-Sw_DC+Sw_HF
        *************************************************************/
        ifroot    std::cout << "*********************************\n";
        ifroot    std::cout << "****TB: Constructing G_local*****\n";
        ifroot    std::cout << "*********************************\n";
        double muIN = muTB;
        /*Solve TB hamiltonian to get, Gloc, mu, ...*/
        muTB = TightBinding (muTB, std::string("Hk.HWR"), SelfEnergy_w,
                             weiss_fieldTB_weakCorr, weiss_fieldTB_strongCorr,1, SolverBasis);
        dc_for_dmft( Sw_doublecounting,  Uindex,  Utensor, NumMatrix, NumCorrAtom, 1, SelfEnergy_w);
        NumMatrixLatt = NumMatrix;

        /***********************************************************
        //print and convergence check
        ***********************************************************/
        bool converg = dmft_scf_check( NumMatrixLatt, NumMatrixImp, timeStartIt, timeEndIt, currentIt );
        if (converg) {
            std::cout << "DMFT: The DMFT calculation has reached convergence." << mpi_rank <<"\n" ;
            break;
        }
        MPI_Barrier(MPI_COMM_WORLD);
        currentIt++;
    }// while, dmft iterations


//    }/*charge density iterations*/
    for(int i=0; i<NumCorrAtom; i++) {
        delete [] HartrRange[i];
    }
    delete [] HartrRange;
    delete [] KS2Hartr;

    delete [] FromOrbitalToAtom;
    delete [] FromOrbitalToLocalOrbital_DFT;
    delete [] HartrIndex_inDFT;

#if SLEPc_ENABLED_dmft
    PetscFinalize();
#endif
    MPI_Finalize();
    return 0;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////
void dc_for_dmft( Eigen::MatrixXcd  & Sw_doublecounting, std::vector<Eigen::VectorXi> & Uindex, std::vector<cmplx >  & Utensor, Eigen::MatrixXcd & NelectronDFT,
                  int NumCorrAtom, int initial_dc, ImgFreqFtn & SelfEnergy_w) {
    ifroot std::cout << "*************************\n";
    ifroot std::cout << "Calculate double counting\n";
    ifroot std::cout << "*************************\n";


    std::vector<Eigen::MatrixXcd> Sw_doublecounting_atom(NumCorrAtom);
    for (int  at=0; at< NumCorrAtom; at++) {
        Sw_doublecounting_atom[at].setZero(N_peratom_HartrOrbit, N_peratom_HartrOrbit);
    }
    double averU =( N_peratom_HartrOrbit * UHubb + N_peratom_HartrOrbit*(N_peratom_HartrOrbit-1)*Uprime) /
                  (N_peratom_HartrOrbit*N_peratom_HartrOrbit);
    double averJ = JHund;

    if (dctype == "sigma0" and initial_dc==1) {
        Eigen::MatrixXcd Id;
        Id.setIdentity(N_peratom_HartrOrbit,N_peratom_HartrOrbit);
        for(int at=0; at<NumCorrAtom; at++) {
            double temp = real((SelfEnergy_w.getMatrix(0,at, N_peratom_HartrOrbit)).trace() / N_peratom_HartrOrbit);
            Sw_doublecounting_atom[at] =
                Sw_doublecounting.block(at*N_peratom_HartrOrbit,at*N_peratom_HartrOrbit, N_peratom_HartrOrbit, N_peratom_HartrOrbit)
                + temp *Id;
        }
    }
    else {
        if(dctype.find(std::string("fll")) != std::string::npos or dctype.find(std::string("sigma0")) != std::string::npos )  {
            /*FLL DC*/
            if( dctype.find(std::string("Uprime")) != std::string::npos)  averU = Uprime;

            double sumSpinUp=0, sumSpinDown=0;

            for(int at=0; at<NumCorrAtom; at++) {
                sumSpinUp=0;
                sumSpinDown=0;
                for(int n=0; n<N_peratom_HartrOrbit; n+=2) {
                    int nH =    at*N_peratom_HartrOrbit + n;
                    sumSpinUp   += real(NelectronDFT(nH,nH));
                    sumSpinDown += real(NelectronDFT(nH+1,nH+1));
                }
                for(int i=0; i<N_peratom_HartrOrbit; i+=2) {
                    Sw_doublecounting_atom[at](i,i  )=     averU*(sumSpinDown+sumSpinUp-0.5)-averJ*(sumSpinUp-0.5);
                    Sw_doublecounting_atom[at](i+1,i+1)=     averU*(sumSpinDown+sumSpinUp-0.5)-averJ*(sumSpinDown-0.5);
                }
                ifroot std::cout << "Spin up:" << sumSpinUp << " Spin down:" << sumSpinDown <<"\n";
            }//at
        }
        else if(dctype =="amf") {
            /*AMF DC*/
            double sumSpinUp=0, sumSpinDown=0;

            for(int at=0; at<NumCorrAtom; at++) {
                sumSpinUp=0;
                sumSpinDown=0;
                for(int n=0; n<N_peratom_HartrOrbit; n+=2) {
                    int nH =    at*N_peratom_HartrOrbit + n;
                    sumSpinUp   += real(NelectronDFT(nH,nH));
                    sumSpinDown += real(NelectronDFT(nH+1,nH+1));
                }
                double n0 = (sumSpinDown+sumSpinDown)/N_peratom_HartrOrbit;

                for(int i=0; i<N_peratom_HartrOrbit; i++) {
                    Sw_doublecounting_atom[at](i,i  )=     averU*n0*(N_peratom_HartrOrbit-1) - averJ*n0*(N_peratom_HartrOrbit/2 - 1);
                }
            }//at
        }
//        if(dctype.find(std::string("fll")) != std::string::npos or dctype.find(std::string("sigma0")) != std::string::npos )  {
        else if ( dctype.find(std::string("nominal")) != std::string::npos      ) {
            /*nominal DC*/
            double sumSpinUp=0, sumSpinDown=0;
            if( dctype.find(std::string("Uprime")) != std::string::npos)  averU = Uprime;

            for(int at=0; at<NumCorrAtom; at++) {
                sumSpinUp=0;
                sumSpinDown=0;
                for(int i=0; i<N_peratom_HartrOrbit; i+=2) {
                    Sw_doublecounting_atom[at](i,i  )=     averU*(nominal_charge-0.5) - averJ*(nominal_charge/2.0 - 0.5);
                    Sw_doublecounting_atom[at](i+1,i+1)=     averU*(nominal_charge-0.5) - averJ*(nominal_charge/2.0 - 0.5);
                }
                ifroot std::cout << "Spin up:" << sumSpinUp << " Spin down:" << sumSpinDown <<"\n";
            }//at
        }
    }
    Sw_doublecounting.setZero(NumCorrAtom*N_peratom_HartrOrbit, NumCorrAtom*N_peratom_HartrOrbit);
    for (int at=0; at<NumCorrAtom; at++) {
        Sw_doublecounting.block(at*N_peratom_HartrOrbit,at*N_peratom_HartrOrbit, N_peratom_HartrOrbit, N_peratom_HartrOrbit) = Sw_doublecounting_atom[at];
    }


    //print
    ifroot {
        std::cout << "double counting Self Energy:\n";
        for(int at=0; at<NumCorrAtom; at++) {
            for(int i=0; i<N_peratom_HartrOrbit; i+=2) {
                std::cout <<  Sw_doublecounting_atom[at](i,i) <<"  ";
                std::cout <<  Sw_doublecounting_atom[at](i+1,i+1) <<"\n";
            }
        }
        FILE * DC = fopen("DC.dat","w");
        for (int at=0; at<NumCorrAtom; at++) {
            for (int orb1=0; orb1<N_peratom_HartrOrbit; orb1++) {
                for (int orb2=0; orb2<N_peratom_HartrOrbit; orb2++) {
                    fprintf(DC, "     %0.5f %0.5f",real(Sw_doublecounting_atom[at](orb1,orb2)),
                            imag(Sw_doublecounting_atom[at](orb1,orb2)));
                }
                fprintf(DC, "\n");
            }
        }
        fclose(DC);
        sleep(1);
    }



}
///////////////////////////////////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////////////////////////////////
void Time_print() {
    struct tm *t;
    time_t timer;
    char s[20];
    timer = time(NULL);
    t = localtime(&timer);
    sprintf(s,"%04d-%02d-%02d %02d:%02d:%02d",
            t->tm_year + 1900, t->tm_mon + 1, t->tm_mday,
            t->tm_hour, t->tm_min, t->tm_sec);
    printf("%s\n",s);
}
///////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////
void write_gnuplot(int writingFile) {
    //and write .gnuplot files
    FILE * GNU = fopen("gw_loc_Im.gnuplot", "w");
    fprintf(GNU,"set term x11 dashed\n p");
    if(SOLVERtype==std::string("ALPS_CTSEG")) { /*CT-HYB, alps*/
        fprintf(GNU,"\\\n\"./Gw_loc.dat0\" u 1:%d   w l  lw 2      lc rgb  \"black\"   lt 1      title \"%d\"", 2*0+3, 0);
        for ( int i=1; i< N_peratom_HartrOrbit; i++) {
            fprintf(GNU,",\\\n\"./Gw_loc.dat0\" u 1:%d   w l  lw 2      lc rgb  \"black\"   lt 1      title \"%d\"", 2*i+3, i);
        }
        for ( int i=0; i< N_peratom_HartrOrbit; i++) {
            fprintf(GNU,",\\\n\"./Gw.dat\" u 1:%d   w l  lw 2      lc rgb  \"red\"   lt 1      title \"%d\"", 2*i+3, i);
        }
    }
    else { /*ED*/
        for ( int i=0; i< N_peratom_HartrOrbit; i++) {
            fprintf(GNU,"\\\n\"./Gw_loc.dat0\" u 1:%d   w l  lw 2      lc rgb  \"black\"   lt 1      title \"%d\" ,", 2*i+3, i);
        }
        for ( int i=0; i< N_peratom_HartrOrbit; i++) {
            fprintf(GNU,"\\\n\"./Gw_imp.dat0\" u 1:%d   w l  lw 2      lc rgb  \"red\"   lt 1      title \"%d\" ,", 2*i+3, i);
        }
    }
    fprintf(GNU,"\npause -1");
    fclose(GNU);

    //write Sw_Im.gnuplot
    GNU = fopen("sw_Im.gnuplot", "w");
    fprintf(GNU,"set term x11 dashed\n p");
    if(SOLVERtype==std::string("ALPS_CTSEG")) { /*CT-HYB, alps*/

        fprintf(GNU,"\\\n\"./Sw.dat\" u 1:($%d)   w l  lw 2      lc rgb  \"red\"   lt 1      title \"%d\"", 2*0+3, 0);
        for ( int i=1; i< N_peratom_HartrOrbit; i++) {
            fprintf(GNU,",\\\n\"./Sw.dat\" u 1:($%d)   w l  lw 2      lc rgb  \"red\"   lt 1      title \"%d\"", 2*i+3, i);
        }
        for ( int i=0; i< N_peratom_HartrOrbit; i++) {
            fprintf(GNU,",\\\n\"./Swl.dat\" u 1:($%d)   w l  lw 2      lc rgb  \"black\"   lt 1      title \"%d\"", 2*i+3, i);
        }
        for ( int i=0; i< N_peratom_HartrOrbit; i++) {
            fprintf(GNU,",\\\n\"./Sw_SOLVER.dat\" u 1:($%d)   w l  lw 2      lc rgb  \"gray\"   lt 1      title \"%d\"", 2*i+3, i);
        }
    }
    else { /*ED*/
        for ( int i=0; i< N_peratom_HartrOrbit; i++) {
            fprintf(GNU,"\\\n\"./Sw_SOLVER.dat\" u 1:($%d)   w l  lw 2      lc rgb  \"red\"   lt 1      title \"%d\" ,", 2*i+3, i);
        }
    }
    fprintf(GNU,"\npause -1");
    fclose(GNU);
}
///////////////////////////////////////////////////////////////////////////////////////////////////////







///////////////////////////////////////////////////////////////////////////////////////////////////////
//void set_Hartree_self_energy( Eigen::MatrixXcd & NumMatrix, std::vector<Eigen::VectorXi > Uindex, std::vector<cmplx >  Utensor, ImgFreqFtn & SelfEnergy_w   ) {
//
//    for(int b=0; b<NumCorrAtom*N_peratom_HartrOrbit; b++) {
//        for(int d=0; d<NumCorrAtom*N_peratom_HartrOrbit; d++) {
//            SelfEnergy_w.setValue(N_freq, b,d, 0);
//        }
//    }
//
//    for(int at=0; at<NumCorrAtom; at++) {
//        for(int idx=0;  idx <  Utensor.size(); idx++) {
//            int aH = at* N_peratom_HartrOrbit + Uindex[idx](0);
//            int bH = at* N_peratom_HartrOrbit + Uindex[idx](1);
//            int cH = at* N_peratom_HartrOrbit + Uindex[idx](2);
//            int dH = at* N_peratom_HartrOrbit + Uindex[idx](3);
//            if (not( isOrbitalCorrinHart[aH] and
//                     isOrbitalCorrinHart[bH] and
//                     isOrbitalCorrinHart[cH] and
//                     isOrbitalCorrinHart[dH] ) ) {
//                SelfEnergy_w.setValue(N_freq, bH,dH, SelfEnergy_w.getValue(N_freq,bH,dH) + Utensor[idx] * NumMatrix(cH,aH)   );
//                SelfEnergy_w.setValue(N_freq, bH,cH, SelfEnergy_w.getValue(N_freq,bH,cH) - Utensor[idx] * NumMatrix(dH,aH)   );
//            }
//        }
//    }
//    for(int n=0; n<N_freq; n++) {
//        for(int b=0; b<NumCorrAtom*N_peratom_HartrOrbit; b++) {
//            for(int d=0; d<NumCorrAtom*N_peratom_HartrOrbit; d++) {
//                SelfEnergy_w.setValue(n, b,d, SelfEnergy_w.getValue(N_freq,b,d));
//            }
//        }
//
//    }
//    FILE * DC = fopen("SwHartree.dat","w");
//    for (int orb1=0; orb1<N_peratom_HartrOrbit*NumCorrAtom; orb1++) {
//        for (int orb2=0; orb2<N_peratom_HartrOrbit*NumCorrAtom; orb2++) {
//            fprintf(DC, "     %0.5f %0.5f",real(SelfEnergy_w.getValue(N_freq,orb1,orb2)), imag(SelfEnergy_w.getValue(N_freq,orb1,orb2)));
//        }
//        fprintf(DC, "\n");
//    }
//    fclose(DC);
//    sleep(1);
//}
///////////////////////////////////////////////////////////////////////////////////////////////////////




///////////////////////////////////////////////////////////////////////////////////////////////////////
void write_results(int DFTIt, int currentIt, std::string system_name, int NumCorrAtom) {
    ifroot {
        FILE * DC;
        std::string cp_comm;
        char chr[100];
        system("cp   Sw_SOLVER.full.dat Sw_SOLVER_weak.full.dat  Numele.dat  mu_history.out DM_DFT_DMFT.HWR dual_DM_direct.HWR muDFT.out   ./Restart");



        sprintf(chr,"%d_%d",DFTIt,currentIt);


        for(int at=0; at<NumCorrAtom; at++) {
            cp_comm = std::string("cp Gw_loc.dat") + intToString(at) +  std::string(" ./") +system_name+ "/" + chr + "Gw_loc.dat" + intToString(at);
            system(cp_comm.c_str());
        }


        cp_comm = std::string("cp std.out  ./") +system_name+ "/" +  "std.out";
        system(cp_comm.c_str());
        cp_comm = std::string("cp delta_t.dat  ./") +system_name+ "/" + chr + "delta_t.dat";
        system(cp_comm.c_str());
        cp_comm = std::string("cp ctqmc.log  ./") +system_name+ "/" + chr + "ctqmc.log";
        system(cp_comm.c_str());

//        cp_comm = std::string("cp delta_w.dat  ./") +system_name+ "/" + chr + "delta_w.dat";
//        system(cp_comm.c_str());
        cp_comm = std::string("cp Numele.dat  ./") +system_name+ "/" + chr + "Numele.dat";
        system(cp_comm.c_str());
        cp_comm = std::string("cp ctqmc.log  ./") +system_name+ "/" + chr + "ctqmc.log";
        system(cp_comm.c_str());
        cp_comm = std::string("cp Sw_SOLVER.dat  ./") +system_name+ "/" + chr + "Sw_SOLVER.dat";
        system(cp_comm.c_str());
        cp_comm = std::string("cp mu_vector.alps.dat  ./") +system_name+ "/" + chr + "mu_vector.alps.dat";
        system(cp_comm.c_str());
        cp_comm = std::string("rm -r ")+ "  ./"+system_name+ "/Restart";
        system(cp_comm.c_str());
        cp_comm = std::string("cp -r ")+ "./Restart"    + "  ./"+system_name+ "/Restart";
        system(cp_comm.c_str());
        if(SOLVERtype==std::string("ALPS_CTSEG")) {
            cp_comm = std::string("cp Gw.dat  ./") +system_name+ "/" + chr + "Gw.dat";
            system(cp_comm.c_str());
            cp_comm = std::string("cp orders.dat  ./") +system_name+ "/" + chr + "orders.dat";
            system(cp_comm.c_str());
            cp_comm = std::string("cp Gwl.dat  ./") +system_name+ "/" + chr + "Gwl.dat";
            system(cp_comm.c_str());
            cp_comm = std::string("cp Ft.dat  ./") +system_name+ "/" + chr + "Ft.dat";
            system(cp_comm.c_str());
            cp_comm = std::string("cp Gt.dat  ./") +system_name+ "/" + chr + "Gt.dat";
            system(cp_comm.c_str());
            cp_comm = std::string("cp Gtl.dat  ./") +system_name+ "/" + chr + "Gtl.dat";
            system(cp_comm.c_str());
            cp_comm = std::string("cp Swl.dat  ./") +system_name+ "/" + chr +"Swl.dat";
            system(cp_comm.c_str());
            cp_comm = std::string("cp Sw.dat  ./") +system_name+ "/" + chr +"Sw.dat";
            system(cp_comm.c_str());
            cp_comm = std::string("cp observables.dat  ./") +system_name+ "/" + chr +"observables.dat";
            system(cp_comm.c_str());
            cp_comm = std::string("cp simulation.dat  ./") +system_name+ "/" + chr +"simulation.dat";
            system(cp_comm.c_str());
            cp_comm = std::string("cp sector_statistics.dat ./") +system_name+ "/" + chr +"sector_statistics.dat";
            system(cp_comm.c_str());
        }
    }/*write results, ifroot*/
}
///////////////////////////////////////////////////////////////////////////////////////////////////////





bool dmft_scf_check( Eigen::MatrixXcd NumMatrixLatt, Eigen::MatrixXcd NumMatrixImp, time_t timeStartIt, time_t timeEndIt, int currentIt ) {


    ifroot {

        std::cout<<"\nElectron Number matrix:\n";
        for(int at=0; at<NumCorrAtom; at++) {
            for(int  n=at*N_peratom_HartrOrbit; n<at*N_peratom_HartrOrbit+N_peratom_HartrOrbit; n++) {
                std::cout << "         ";
                for(int m=at*N_peratom_HartrOrbit; m<at*N_peratom_HartrOrbit+N_peratom_HartrOrbit; m++) {
                    std::cout <<NumMatrixLatt(n,m) <<"   ";
                }
                std::cout << "\n";
            }
        }

        for(int at=0; at<NumCorrAtom; at++) {
            std::cout << "\nNum ele (decomp,solver)  = ";
            double sum=0;
            for(int n=0; n<N_peratom_HartrOrbit; n+=1) {
                std::cout <<std::fixed << std::setprecision(4)<< std::fixed   <<  real(NumMatrixImp(at*N_peratom_HartrOrbit+n,at*N_peratom_HartrOrbit+n)) <<" " ;
                sum+= real(NumMatrixImp(at*N_peratom_HartrOrbit+n,at*N_peratom_HartrOrbit+n));
            }
            std::cout << std::fixed << std::setprecision(4)<< std::fixed   << ";  (total)  = "<<  sum <<"\n" ;
        }

        for(int at=0; at<NumCorrAtom; at++) {
            std::cout << "\nNum ele (decomp,lattice)  = ";
            double sum=0;
            for(int n=0; n<N_peratom_HartrOrbit; n+=1) {
                std::cout <<std::fixed << std::setprecision(4)<< std::fixed   <<  real(NumMatrixLatt(at*N_peratom_HartrOrbit+n,at*N_peratom_HartrOrbit+n)) <<" " ;
                sum+= real(NumMatrixLatt(at*N_peratom_HartrOrbit+n,at*N_peratom_HartrOrbit+n));
            }
            std::cout << std::fixed << std::setprecision(4)<< std::fixed   << ";  (total)  = "<<  sum <<"\n" ;

            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> ces(N_peratom_HartrOrbit);
            Eigen::MatrixXcd NumMat_atom1(N_peratom_HartrOrbit, N_peratom_HartrOrbit);
            for(int n=0; n<N_peratom_HartrOrbit; n+=1) {
                for(int m=0; m<N_peratom_HartrOrbit; m+=1) {
                    NumMat_atom1(n,m) = NumMatrixLatt(at*N_peratom_HartrOrbit+n, at*N_peratom_HartrOrbit+m);
                }
            }
            ces.compute(NumMat_atom1);
            std::cout << "Num ele (diagon, lattice) atom"<<at<<" = ";
            for(int n=0; n<N_peratom_HartrOrbit; n+=1) {
                std::cout <<std::fixed << std::setprecision(4)<< std::fixed   <<  ( ces.eigenvalues()[n]  ) <<" " ;
            }
            std::cout <<"\n";
        }

        // Write results
        write_results(DFTIt, currentIt, system_name, NumCorrAtom);
        timeEndIt = clock();
        std::cout <<"Computaional time for single SC-loop:\n"
                  << ((timeEndIt-timeStartIt)/(CLOCKS_PER_SEC))/3600
                  <<":"<< (((timeEndIt-timeStartIt)/(CLOCKS_PER_SEC))%3600)/60
                  <<":"<<(((timeEndIt-timeStartIt)/(CLOCKS_PER_SEC))%3600)%60 <<  "\n";
        std::cout <<"using "  <<mpi_numprocs<< " processors\n";
    }//ifroot


    /*SCF check; */
    ImgFreqFtn GwDMFT_Loc(beta, N_freq, N_peratom_HartrOrbit, NumCorrAtom,0);
    ImgFreqFtn GwDMFT_Imp(beta, N_freq, N_peratom_HartrOrbit, NumCorrAtom,0);
    if(N_peratom_HartrOrbit>0) GwDMFT_Loc.update_full(std::string("Gw_loc.full.dat0"),1);
    if(N_peratom_HartrOrbit>0) GwDMFT_Imp.update_full(std::string("Gw_imp.full.dat0"),1);
    double        GwNorm=1e-10;
    double        GwlowNorm=1e-10;
    double        GwRD=0;
    double        GwRD_low=0;
    double        GwRD_max=0;
    double GwRDMAXnorm=1;
    double NumMatRD = 0;
    for (int n=0; n<N_freq; n++) {
        double GwRD_n=0, GwNorm_n=1e-10;
        for (int i=0; i<N_peratom_HartrOrbit; i++) {
            for (int j=0; j<N_peratom_HartrOrbit; j++) {
                if  (GwDMFT_Imp.getValue(n) < std::min(UHubb,3.) ) {
                    GwRD_low  +=   (pow(std::abs(GwDMFT_Imp.getValue(n,i,j) - GwDMFT_Loc.getValue(n,i,j)),2));
                    GwlowNorm +=   (pow(std::abs(GwDMFT_Loc.getValue(n,i,j)                             ),2));
                }
                GwRD_n   += (pow(std::abs(GwDMFT_Imp.getValue(n,i,j) - GwDMFT_Loc.getValue(n,i,j)),2));
                GwNorm_n += (pow(std::abs(GwDMFT_Loc.getValue(n,i,j)                             ),2));
                GwRD   += (pow(std::abs(GwDMFT_Imp.getValue(n,i,j) - GwDMFT_Loc.getValue(n,i,j)),2));
                GwNorm += (pow(std::abs(GwDMFT_Loc.getValue(n,i,j)                             ),2));
            }//forj
        }
        if ( GwRD_n/ GwNorm_n> GwRD_max ) {
            GwRD_max    = GwRD_n/GwNorm_n;
        }
    }
    GwRD = std::sqrt(GwRD/GwNorm);
    GwRD_low = std::sqrt(GwRD_low/GwlowNorm);
    GwRD_max = std::sqrt(GwRD_max);
    NumMatRD = ((NumMatrixLatt-NumMatrixImp).norm()/NumMatrixLatt.norm());


    ifroot printf(  "GwRD       =%e\n", GwRD);
    ifroot printf(  "GwRD_low   =%e\n", GwRD_low);
    ifroot printf(  "GwRD_max   =%e\n", GwRD_max);
//    ifroot printf(  "muRD_fluct =%e\n", muRDFluct);
    ifroot printf(  "NumMatRD   =%e\n", NumMatRD);

    if(  GwRD        < 1e-5
            and GwRD_low    < 1e-5
            and GwRD_max    < 1e-5
            and NumMatRD   <  1e-3
      ) {
//            and muRDFluct   < 0.05
        return true;
    }
    else return false;
}
