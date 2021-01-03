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

double  TightBinding(double mu, const std::string &hamiltonian, ImgFreqFtn & SelfE_z,
                     ImgFreqFtn & weiss_fieldTB, ImgFreqFtn & weiss_fieldTBCorr,
                     int mu_adjustTB,  Eigen::MatrixXcd & SolverBasis
                    );

void rot_Uijkl(
    std::vector<cmplx > & Utensor, std::vector<Eigen::VectorXi>  & Uindex,
    std::vector<cmplx > & rotUtensor, std::vector<Eigen::VectorXi>  & rotUindex,
    Eigen::MatrixXcd & SolverBasis, int n_spinorb
) ;
void proj_Uijkl(
    std::vector<cmplx > & Utensor, std::vector<Eigen::VectorXi>  & Uindex,
    std::vector<cmplx > & projUtensor, std::vector<Eigen::VectorXi>  & projUindex,
    std::vector<int> & sub2full, int n_spinorb
) ;


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
void write_restart(int DFTIt, int currentIt, std::string system_name, int  NumCluster);

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
    ImgFreqFtn weiss_field_weakCorr(0);
    ImgFreqFtn weiss_field_strongCorr(0);
    ImgFreqFtn weiss_fieldTB_weakCorr(0);
    ImgFreqFtn weiss_fieldTB_strongCorr(0);
    ImgFreqFtn SelfEnergy_w(0);
    ImgFreqFtn SelfEnergy_E(0);

//    int mystaE, myendE;
//    para_range(0,Spectral_EnergyGrid-1, mpi_numprocs, mpi_rank, &mystaE, &myendE);

    weiss_fieldTB_weakCorr.Initialize(  beta, N_freq, NumHartrOrbit_per_cluster,NumCluster, mixingType);
    weiss_field_weakCorr.Initialize(    beta, N_freq, NumHartrOrbit_per_cluster,NumCluster, mixingType);
    weiss_fieldTB_strongCorr.Initialize(beta, N_freq, NSpinOrbit_per_atom,NumCorrAtom, mixingType);
    weiss_field_strongCorr.Initialize(  beta, N_freq, NSpinOrbit_per_atom,NumCorrAtom, mixingType);



//    Eigen::MatrixXcd  Sw_doublecounting;
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
//  Initial Self energy,   and Sw_inf and muTB
////////////////////////////////////////////////////
    /*retarded self-energy, S_E, for quasi-ptl calculation */
    SelfEnergy_w.Initialize(         beta, N_freq,  NumHartrOrbit_per_cluster, NumCluster, mixingType);
    if(SOLVERtype !=std::string("TB")) {


        for(int n=0; n<N_freq; n++) {
            for(int at=0; at<NumCorrAtom; at++) {
                for(int  i=at*N_peratom_HartrOrbit; i<at*N_peratom_HartrOrbit+N_peratom_HartrOrbit; i+=2) {
                    for(int  j=at*N_peratom_HartrOrbit; j<at*N_peratom_HartrOrbit+N_peratom_HartrOrbit; j+=2) {
                        SelfEnergy_w.setValue(n, i+0,j+0, Zeeman_field_spin_corratom[at](0,0));
                        SelfEnergy_w.setValue(n, i+0,j+1, Zeeman_field_spin_corratom[at](0,1));
                        SelfEnergy_w.setValue(n, i+1,j+0, Zeeman_field_spin_corratom[at](1,0));
                        SelfEnergy_w.setValue(n, i+1,j+1, Zeeman_field_spin_corratom[at](1,1));
                    }
                }
            }

        }

        if(restart!=0 and N_peratom_HartrOrbit >0) {
            /*Matsubara self-energy, S_w*/
            double beta_prev;
            std::ifstream inputbeta("./beta.dat");
            inputbeta >> beta_prev;
            ifroot std::cout <<"Reading Self-energy...  " << beta_prev << "\n";
            SelfEnergy_w.read_full(std::string("./Restart/Sw_SOLVER.full.dat"),beta, beta_prev);
        }
        SelfEnergy_w.dataOut(std::string("InitialSw.dat"));
        SelfEnergy_w.dataOut_full(std::string("InitialSw.full.dat"));
    }
    else if(SOLVERtype == std::string("TB")) {
        SelfEnergy_E.realFreq(       E0,  real(dE),  Spectral_EnergyGrid,  NumHartrOrbit_per_cluster, NumCluster, 0, 0    );
        if ((mode==std::string("qsband") or mode==std::string("qsdos"))  ) {
            if (restart<=0 ) {
                std::cout << "Please, use restart option\n";
                exit(1);
            }
            FILE * Chem = fopen("./Restart/mu_history.out","r");
            while(!feof(Chem)) {
                fscanf(Chem, "%lf\n",&muTB);
            }
            fclose(Chem);
//            Chem = fopen("./Restart/muDFT.out","r");
//            while(!feof(Chem)) {
//                fscanf(Chem, "%lf\n",&muDFT);
//            }
            fclose(Chem);
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


            /*Matsubara self-energy, S_w used to adjust chemical potential*/
            if( 0 == access("Sw_SOLVER.full_fromRetardedSw.dat", F_OK) ) {
                double beta_prev;
                std::ifstream inputbeta("./beta.dat");
                inputbeta >> beta_prev;
                ifroot std::cout <<"Reading Self-energy... to adjust chem. pot.  from beta_prev" << beta_prev << "\n";
                SelfEnergy_w.read_full(std::string("Sw_SOLVER.full_fromRetardedSw.dat"),beta, beta_prev);
                SelfEnergy_w.dataOut(std::string("Sw_SOLVER_fromRetardedSw.dat"));
                muTB = TightBinding (muTB, std::string("Hk.HWR"), SelfEnergy_w, weiss_fieldTB_weakCorr, weiss_fieldTB_strongCorr, -2, SolverBasis);
            }
            else {
                double beta_prev;
                std::ifstream inputbeta("./beta.dat");
                inputbeta >> beta_prev;
                SelfEnergy_w.read_full(std::string("./Restart/Sw_SOLVER.full.dat"),beta, beta_prev);
            }
            rewrite_retardedSw( SelfEnergy_w );
            if( N_peratom_HartrOrbit > 0) SelfEnergy_E.update_full(std::string("realFreq_Sw.full.dat"),1);
            SelfEnergy_E.dataOut(std::string("realFreq_Sw.dat"));
        }
    }
    ifroot std::cout << "Initial Self-energy was constructed," << mpi_rank <<"\n";

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
    *******************************************************************************************************/
    /////////////////////////////////////////////////////
    //  occupation, muDFT
    ////////////////////////////////////////////////////
    /*DFT results and double counting */
    muDFT = 0.0;
    muDFT    =  TightBinding (muDFT, std::string("Hk.HWR"), SelfEnergy_w,   weiss_fieldTB_weakCorr, weiss_fieldTB_strongCorr,  -1, SolverBasis);
    ifroot std::cout << "First run, Initial Chemical potential, we have muDFT: " << muDFT <<"\n";
    FILE *tempFile;   //mu.dat
//    tempFile = fopen("muDFT.out", "w");
//    fprintf(tempFile, "%0.20f\n", muDFT);
//    fclose(tempFile) ;
    /////////////////////////////////////////////
    //single TB run  (band, dos...)
    /////////////////////////////////////////////
    if (SOLVERtype==std::string("TB")) {
        ifroot std::cout << "Here we start TB.\n";
        muTB = TightBinding (muTB, std::string("Hk.HWR"), SelfEnergy_E, weiss_fieldTB_weakCorr, weiss_fieldTB_strongCorr, 0, SolverBasis);   //both muTB (qsdos, qsband)  and muDFT (dos, band) are used .
        ifroot std::cout << "Solver=0 loop finished\n";
        MPI_Finalize();
        return 0;
    }



    /*************************************************
     DMFT SC loop :   Σ -TB->  Δ, Σdc, mu   -solver->  Σ
    ***************************************************/
    int currentIt =1;
    dc_for_dmft( Sw_doublecounting,  Uindex,  Utensor, NumMatrix, NumCorrAtom, 0, SelfEnergy_w);
    muTB = muDFT;
    pulayMixing mixing_self_energy(pulay_hist, pulay_start, N_freq, NumHartrOrbit_per_cluster*NumCluster, NumHartrOrbit_per_cluster*NumCluster, true );
    while(true) {
        Eigen::MatrixXcd NumMatrixImp;
        Eigen::MatrixXcd NumMatrixLatt;
        time_t timeStartIt, timeEndIt;
        timeStartIt = clock();
        ifroot     std::cout <<"\n#########################################\n";
        ifroot     std::cout <<"    DFT : "<<DFTIt <<"    DMFT : SCF  "<< currentIt <<"\n";
        ifroot     std::cout <<"#########################################\n";

        on_the_fly_control();


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
//        double muIN = muTB;
        /*Solve TB hamiltonian to get, Gloc, mu, ...*/
        muTB = TightBinding (muTB, std::string("Hk.HWR"), SelfEnergy_w, weiss_fieldTB_weakCorr, weiss_fieldTB_strongCorr,1, SolverBasis);
        NumMatrixLatt = NumMatrix;




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


        weiss_field_weakCorr.dataOut(std::string("delta_weak_w.dat"));
        weiss_field_strongCorr.dataOut(std::string("delta_w.dat"));
        weiss_field_strongCorr.dataOut_full(std::string("delta_w.full.dat"));
        ifroot std::cout << "FILEOUT:delta_w.dat\n" ;
        /*double counting*/
        if( dctype.find(std::string("dmft"))  != std::string::npos )
            dc_for_dmft( Sw_doublecounting,  Uindex,  Utensor, NumMatrix, NumCorrAtom, 0, SelfEnergy_w);
        /////////////////
        //Delta[Sigma], Sdc[Sigma]

        /***********************************************************
        run impurity
        in  = delta_w, delta_t, mu, onsite_Hamiltonian, SelfEnergy_w(initial) ;
        out = SelfEnergy_w, GwImp_DMFT, NumMat
        note : Self-energy = \Sigma_Hartree + \Simga_imp -\Simga_dc
        *************************************************************/

        ImgFreqFtn SelfEnergy_w_out(0);
        SelfEnergy_w_out.Initialize(         beta, N_freq, NumHartrOrbit_per_cluster, NumCluster, mixingType);

        NumMatrix =                                     (SolverBasis.adjoint() * NumMatrix                 * SolverBasis).eval();
        Eigen::MatrixXcd  rotimpurity_site_Hamiltonian = SolverBasis.adjoint() * impurity_site_Hamiltonian * SolverBasis;
        Eigen::MatrixXcd  rotSw_doublecounting =         SolverBasis.adjoint() * Sw_doublecounting         * SolverBasis;


        Epot_weak = 0.0;
        Epot = 0.0;
        for(int cl=0; cl < NumCluster;  cl++) {

            Eigen::MatrixXcd  SolverBasis_cl;
            SolverBasis_cl =  SolverBasis.block(cl*NumHartrOrbit_per_cluster, cl*NumHartrOrbit_per_cluster, NumHartrOrbit_per_cluster, NumHartrOrbit_per_cluster) ;

            std::vector<cmplx > rotUtensor;
            std::vector<Eigen::VectorXi> rotUindex;
            rot_Uijkl(Utensor, Uindex, rotUtensor, rotUindex, SolverBasis_cl, NumHartrOrbit_per_cluster);


            ImgFreqFtn Gw_weak(beta, N_freq, NumHartrOrbit_per_cluster, 1,0);
            ImgFreqFtn GwHF_weak(beta, N_freq, NumHartrOrbit_per_cluster, 1,0);
            ImgFreqFtn rotGw_weak(beta, N_freq, NumHartrOrbit_per_cluster, 1,0);
            ImgFreqFtn rotGwHF_weak(beta, N_freq, NumHartrOrbit_per_cluster, 1,0);

            Gw_weak.update_full(std::string("Gw_loc.full.dat") + intToString(cl),1);
            GwHF_weak.update_full(std::string("GwHF_loc.full.dat") + intToString(cl),1);
            for(int n =0; n<N_freq; n++) {
                rotGw_weak.setMatrix(n,    (SolverBasis_cl).adjoint()* Gw_weak.getMatrix(n)   *SolverBasis_cl);
                rotGwHF_weak.setMatrix(n,  (SolverBasis_cl).adjoint()* GwHF_weak.getMatrix(n) *SolverBasis_cl);
            }


            ImgFreqFtn SE_out(0);
            ImgFreqFtn rotSE_lowlevel(0);
            ImgFreqFtn SE_strong(0);

            SE_out.Initialize(         beta, N_freq,  NumHartrOrbit_per_cluster, 1, mixingType);
            rotSE_lowlevel.Initialize(beta, N_freq, NumHartrOrbit_per_cluster,1,0);
            SE_strong.Initialize(beta,   N_freq, NSpinOrbit_per_atom,      1,0);

            std::vector<int> strongCorr(NSpinOrbit_per_atom);
            std::vector<int> strongCorrinCl(NSpinOrbit_per_atom);

            if(NumHartrOrbit_per_cluster>NSpinOrbit_per_atom) {
                weaksolver_run=1;
            }

            ///////////////////////////////////////////////////
            /*Weakly correlated subspace in the impurity site*/
            ///////////////////////////////////////////////////
            if(weaksolver_run == 1 ) {
                ifroot std::cout <<"Run low level solver\n";
                weak_solver(Lowlevel_SOLVERtype,
                            NumHartrOrbit_per_cluster,
                            rotSE_lowlevel, Gw_weak, GwHF_weak,
                            Uindex,Utensor);



                Eigen::MatrixXcd * Sw_Eig = new Eigen::MatrixXcd [N_freq];
                Eigen::MatrixXcd * Gw_Eig = new Eigen::MatrixXcd [N_freq];
                std::vector<Eigen::MatrixXcd >  Swmoments;
                for(int n=0; n < N_freq; n++) {
                    Sw_Eig[n] = rotSE_lowlevel.getMatrix(n);
                    Gw_Eig[n] = Gw_weak.getMatrix(n);
                }
                getAsymto_moments(Swmoments, Sw_Eig);

                Epot_weak += real(  (Swmoments[0] *  (-1.*FourierTransform_tau(Gw_Eig, beta, true ) )).trace())     *0.5 ;
                for(int n =0; n<N_freq; n++) {
                    Epot_weak +=  2.0* real(((Sw_Eig[n]-Swmoments[0]) * Gw_Eig[n] ).trace())/beta    *0.5;
                    rotSE_lowlevel.setMatrix(n,  (SolverBasis_cl).adjoint()* rotSE_lowlevel.getMatrix(n) *SolverBasis_cl);
                }
                rotSE_lowlevel.setMatrix(N_freq,  (SolverBasis_cl).adjoint()* rotSE_lowlevel.getMatrix(N_freq) *SolverBasis_cl);
                delete [] Gw_Eig;
                delete [] Sw_Eig;

            }//if, weak solver

            ////////////////////////////////
            /*strongly correlated subspace*/
            ///////////////////////////////
            ImgFreqFtn rotGw_strong(0);
            ImgFreqFtn rotGwHF_strong(0);
            rotGw_strong.Initialize(beta, N_freq, NSpinOrbit_per_atom,1,0);
            rotGwHF_strong.Initialize(beta, N_freq, NSpinOrbit_per_atom,1,0);
            for(int at=cl*NumAtom_per_cluster; at < (cl+1)*NumAtom_per_cluster; at++) {
                for (int n=0; n<N_freq+3; n++) {
                    for (int i=0; i<NSpinOrbit_per_atom; i++) {
                        int iF = CorrToHartr(at,i)  - cl * NumHartrOrbit_per_cluster ;
                        if(segmentsolver==1) {
                            rotGw_strong.setValue(n,i,i, rotGw_weak.getValue(n,iF,iF));
                            rotGwHF_strong.setValue(n,i,i, rotGwHF_weak.getValue(n,iF,iF));
                        }
                        else {
                            for (int j=0; j<NSpinOrbit_per_atom; j++) {
                                int jF = CorrToHartr(at,j)  - cl * NumHartrOrbit_per_cluster ;
                                rotGw_strong.setValue(n,i,j,   rotGw_weak.getValue(n,iF,jF));
                                rotGwHF_strong.setValue(n,i,j, rotGwHF_weak.getValue(n,iF,jF));
                            }
                        }
                    }
                }
//            {}
//            for(int at=cl*NumAtom_per_cluster; at < (cl+1)*NumAtom_per_cluster; at++) {}
                for (int i=0; i<NSpinOrbit_per_atom; i++)  strongCorr[i] = CorrToHartr(at, i ) ;
                for (int i=0; i<NSpinOrbit_per_atom; i++)  strongCorrinCl[i] = CorrToHartr(at,i) - cl * NumHartrOrbit_per_cluster ;
                std::vector<cmplx >          projUtensor;
                std::vector<Eigen::VectorXi> projUindex;

                //proj U into subspace  (dd term only).
                proj_Uijkl( rotUtensor,  rotUindex,
                            projUtensor, projUindex, strongCorrinCl, NSpinOrbit_per_atom);

                //double counting
                if(weaksolver_run == 1) {
                    /*remove double counting contribution*/
                    ifroot std::cout <<"Run low level solver for active space\n";
                    ImgFreqFtn SE_lowlevel_local(0);
                    SE_lowlevel_local.Initialize(beta, N_freq, NSpinOrbit_per_atom,1,0);
                    weak_solver(Lowlevel_SOLVERtype,
                                NSpinOrbit_per_atom,
                                SE_lowlevel_local, rotGw_strong, rotGwHF_strong,
                                projUindex,projUtensor);



                    Eigen::MatrixXcd * Sw_Eig = new Eigen::MatrixXcd [N_freq];
                    Eigen::MatrixXcd * Gw_Eig = new Eigen::MatrixXcd [N_freq];
                    std::vector<Eigen::MatrixXcd >  Swmoments;
                    for(int n=0; n < N_freq; n++) {
                        Sw_Eig[n] =SE_lowlevel_local.getMatrix(n);
                        Gw_Eig[n] =rotGw_strong.getMatrix(n);
                    }
                    getAsymto_moments(Swmoments, Sw_Eig);
                    Epot_weak -= real(  (Swmoments[0] *  (-1.*FourierTransform_tau(Gw_Eig, beta, true) )).trace())   *0.5;
                    for(int n =0; n<N_freq; n++) {
                        Epot_weak -=  2.0* real(((Sw_Eig[n]-Swmoments[0]) * Gw_Eig[n] ).trace()) /beta   *0.5;
                    }
                    delete [] Gw_Eig;
                    delete [] Sw_Eig;




                    for (int n=0; n<N_freq+1; n++) {
                        for (int i=0; i<NSpinOrbit_per_atom; i++) {
                            for (int j=0; j<NSpinOrbit_per_atom; j++) {
                                int iF = CorrToHartr(at,i) - cl*NumHartrOrbit_per_cluster;
                                int jF = CorrToHartr(at,j) - cl*NumHartrOrbit_per_cluster;
                                rotSE_lowlevel.setValue(n,iF,jF,
                                                        rotSE_lowlevel.getValue(n,iF,jF)-  SE_lowlevel_local.getValue(n,i,j)   );
                            }
                        }
                    }
                }//if, weak solver
                std::vector<Eigen::MatrixXcd >  dc_weakCorr(N_freq+1);
                for (int n=0; n<N_freq+1; n++) {
                    dc_weakCorr[n].setZero(NSpinOrbit_per_atom,NSpinOrbit_per_atom);
                    for (int i=0; i<NSpinOrbit_per_atom; i++) {
                        for (int j=0; j<NSpinOrbit_per_atom; j++) {
                            int iF = CorrToHartr(at, i)  - cl * NumHartrOrbit_per_cluster;
                            int jF = CorrToHartr(at, j)  - cl * NumHartrOrbit_per_cluster;
                            dc_weakCorr[n](i,j) = rotSE_lowlevel.getValue(n,iF,jF)   ;
                        }
                    }
                }
                for (int n=0; n<N_freq; n++) {
                    dc_weakCorr[n] -= dc_weakCorr[N_freq];
                }

                ifroot std::cout << "High-level solver for atom " << at <<"\n";
                MPI_Barrier(MPI_COMM_WORLD);
                SOLVER(SOLVERtype,
                       NSpinOrbit_per_atom, at, strongCorr,
                       rotimpurity_site_Hamiltonian, weiss_field_strongCorr,
                       SE_strong, rotGw_strong,
                       projUindex, projUtensor,
                       rotSw_doublecounting, dc_weakCorr,
                       muTB);





                if (SOLVERtype==std::string("RUTGERS_CTSEG")) {


                    std::string s, line;
                    std::ifstream fs;

                    fs.open( "ctqmc.log" );

                    if(fs.is_open())
                    {
                        while( getline(fs, line, ' ') ) {
//                            std::cout <<"line " << line <<"\n";
                            std::stringstream ss(line);
                            while(getline(ss, s, '=')) {
//                                std::cout << "s " << s << std::endl;
                                if(s==std::string("TrSigmaG")) {
                                    getline(ss, s, '=');
                                    //                             std::cout << "rst: " << s << std::endl;
                                    Epot += std::stod(s);
                                }
                            }
                        }
                    }
                }
                else {

                    Eigen::MatrixXcd * Sw_Eig = new Eigen::MatrixXcd [N_freq];
                    Eigen::MatrixXcd * Gw_Eig = new Eigen::MatrixXcd [N_freq];
                    std::vector<Eigen::MatrixXcd >  Swmoments;
                    for(int n=0; n < N_freq; n++) {
                        Sw_Eig[n] =SE_strong.getMatrix(n);
                        Gw_Eig[n] =rotGw_strong.getMatrix(n);
                    }
                    getAsymto_moments(Swmoments, Sw_Eig);
                    Epot   += real(  (Swmoments[0] *  (-1.*FourierTransform_tau(Gw_Eig, beta, true) )).trace())   *0.5;
                    for(int n =0; n<N_freq; n++) {
                        Epot +=  2.0* real(((Sw_Eig[n]-Swmoments[0]) * Gw_Eig[n] ).trace()) /beta    *0.5;
                    }
                    delete [] Gw_Eig;
                    delete [] Sw_Eig;
                }





                for (int n=0; n<N_freq+3; n++) {
                    for (int i=0; i<NSpinOrbit_per_atom; i++) {
                        if(segmentsolver==1) {
                            int iF = CorrToHartr(at,i) - cl*NumHartrOrbit_per_cluster;
                            rotGw_weak.setValue(n,iF,iF, rotGw_strong.getValue(n,i,i));   // G^-1 = G_weak^-1 + G_strong^-1
                        }
                        else {
                            for (int j=0; j<NSpinOrbit_per_atom; j++) {
                                int iF = CorrToHartr(at,i) - cl*NumHartrOrbit_per_cluster;
                                int jF = CorrToHartr(at,j) - cl*NumHartrOrbit_per_cluster;
                                rotGw_weak.setValue(n,iF,jF, rotGw_strong.getValue(n,i,j));   // G^-1 = G_weak^-1 + G_strong^-1
                            }
                        }
                    }
                }
                /*combine all result (Sw_weak, Sw_strong, Sw_DFTdc)  to SE_out*/
                for (int n=0; n<N_freq+2; n++) {
                    for (int i=0; i<NSpinOrbit_per_atom; i++) {
                        for (int j=0; j<NSpinOrbit_per_atom; j++) {
                            int iF = CorrToHartr(at,i) - cl *NumHartrOrbit_per_cluster  ;
                            int jF = CorrToHartr(at,j) - cl *NumHartrOrbit_per_cluster  ;
                            SE_out.setValue(n,iF,jF,        rotSE_lowlevel.getValue(n,iF,jF)  + SE_strong.getValue(n,i,j));
                        }
                    }
                }
                SE_strong.dataOut(std::string("Bare_Sw_strong.dat")+intToString(at));
                if(mpi_rank==0) {
                    std::cout << "MO occ = ";
                    double sum=0;
                    for(int n=0; n<N_peratom_HartrOrbit; n+=1) {
                        std::cout <<std::fixed << std::setprecision(4)<< std::fixed   <<  real(NumMatrix(at*N_peratom_HartrOrbit+n,at*N_peratom_HartrOrbit+n)) <<" " ;
                        sum+= real(NumMatrix(at*N_peratom_HartrOrbit+n,at*N_peratom_HartrOrbit+n));
                    }
                    std::cout << std::fixed << std::setprecision(4)<< std::fixed   << ";  (total)  = "<<  sum <<"\n" ;
                }
                MPI_Barrier(MPI_COMM_WORLD);
            }//at


            rotSE_lowlevel.dataOut(std::string("Bare_Sw_lowlevel.dat") + intToString(cl));
            SE_out.dataOut(std::string("Bare_Sw_SOLVER.dat") + intToString(cl));
//            if(doublecounting==1) {
//                for (int  w=0; w<N_freq+1; w++) {
//                    SE_out.setMatrix(w,
//                                     SE_out.getMatrix(w) - rotSw_doublecounting.block(cl*NumHartrOrbit_per_cluster, cl*NumHartrOrbit_per_cluster, NumHartrOrbit_per_cluster, NumHartrOrbit_per_cluster  ));
//                }
//            }
            SelfEnergy_w_out.update(     SE_out,   1,    cl, 0 );

            for(int n =0; n<N_freq; n++) {
                Gw_weak.setMatrix(n,  (SolverBasis_cl)* rotGw_weak.getMatrix(n) *SolverBasis_cl.adjoint());
            }
            Gw_weak.dataOut_full((std::string("Gw_imp.full.dat")+ intToString(cl)  ) );
            Gw_weak.dataOut(std::string("Gw_imp.dat"+ intToString(cl)));
        }//cl, solver


        NumMatrix = (SolverBasis * NumMatrix * SolverBasis.adjoint()).eval();
        for(int n =0; n<N_freq+1; n++) {
            SelfEnergy_w_out.setMatrix(n,  SolverBasis* SelfEnergy_w_out.getMatrix(n) * SolverBasis.adjoint() );
        }

        /*double counting*/
        if(doublecounting==1) {
//       if( dctype.find(std::string("dmft"))  != std::string::npos )
//            dc_for_dmft( Sw_doublecounting,  Uindex,  Utensor, NumMatrix, NumCorrAtom, 0, SelfEnergy_w_out);
            for (int  w=0; w<N_freq+1; w++) {
                SelfEnergy_w_out.setMatrix(w,
                                           SelfEnergy_w_out.getMatrix(w) - Sw_doublecounting            );
            }
        }



        //self-energy mixing.
        if(mixingFtn==1)   {
            mixing_self_energy.mixing ( SelfEnergy_w, SelfEnergy_w_out, mixing, currentIt, 1);
            SelfEnergy_w.update(SelfEnergy_w_out,1.0);
        }


        SelfEnergy_w.dataOut_full(std::string("Sw_SOLVER.full.dat"));
        SelfEnergy_w.dataOut(std::string("Sw_SOLVER.dat"));
        NumMatrixImp = NumMatrix;


        ifroot  std::cout << "Kinetic Energy(EkinDFT) : " << EkinDFT  << "\n";
        ifroot  std::cout << "Kinetic Energy(Ekin)    : " << Ekin  << "\n";
        ifroot  std::cout << "Coulomb Energy(Epot)    : " << Epot + Epot_weak  << "\n";
        ifroot  std::cout << "Coulomb Energy(Epot_strong)   : " << Epot       << "\n";
        ifroot  std::cout << "Coulomb Energy(Edc)      : "  << Edc  << "\n";
        EngDMFT=Ekin-EkinDFT+Epot-Edc ;
        ifroot  std::cout << "Energy Correction(E_DMFT): " << EngDMFT << "\n";




        /***********************************************************
        //print and convergence check
        ***********************************************************/
        bool converg = dmft_scf_check( NumMatrixLatt, NumMatrixImp, timeStartIt, timeEndIt, currentIt );
        if (converg and std::abs(EngDMFT-EngDMFT_prev)<1e-3) {
            std::cout << "DMFT: The DMFT calculation has reached convergence." << mpi_rank <<"\n" ;
            break;
        }
        EngDMFT_prev = EngDMFT;
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
    Edc = 0;
    int Norbitals = N_peratom_HartrOrbit/2;



    double averU  = ( Norbitals * UHubb + Norbitals*(Norbitals-1)*Uprime) / ( Norbitals * Norbitals);
    double averJ  = -(Uprime-JHund - averU);             // = J_{avg} = JHund + 2JHund/N  ,          See pavarini, The LDA+DMFT approach to strongly correlated materials, Chap6

    ifroot std::cout << "Uavg:" << averU << " Javg:" << averJ <<"\n";
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
        //if(dctype.find(std::string("fll2")) != std::string::npos) {

        //    double sumSpinUp=0, sumSpinDown=0, sumTot;
        //    double Ubar =  (UHubb + (Norbitals-1)*Uprime + (Norbitals-1) * (Uprime-JHund)) /  (2*Norbitals-1);



        //    for(int at=0; at<NumCorrAtom; at++) {
        //        sumSpinUp=0;
        //        sumSpinDown=0;
        //        for(int n=0; n<N_peratom_HartrOrbit; n+=2) {
        //            int nH =    at*N_peratom_HartrOrbit + n;
        //            sumSpinUp   += real(NelectronDFT(nH,nH));
        //            sumSpinDown += real(NelectronDFT(nH+1,nH+1));
        //        }
        //        sumTot = sumSpinUp+sumSpinDown;
        //        for(int i=0; i<N_peratom_HartrOrbit; i+=2) {
        //            Sw_doublecounting_atom[at](i,i  )=     Ubar*(sumTot-0.5) ;
        //            Sw_doublecounting_atom[at](i+1,i+1)=   Ubar*(sumTot-0.5) ;
        //        }
        //        ifroot std::cout << "Spin up:" << sumSpinUp << " Spin down:" << sumSpinDown <<"\n";
        //        Edc +=    (
        //                      +  0.5* Ubar * sumTot    * (sumTot    -1)   );
        //    }//at
        //}
        if(dctype.find(std::string("fllsk")) != std::string::npos) {

            double sumSpinUp=0, sumSpinDown=0, sumTot;

            for(int at=0; at<NumCorrAtom; at++) {
                sumSpinUp=0;
                sumSpinDown=0;
                for(int n=0; n<N_peratom_HartrOrbit; n+=2) {
                    int nH =    at*N_peratom_HartrOrbit + n;
                    sumSpinUp   += real(NelectronDFT(nH,nH));
                    sumSpinDown += real(NelectronDFT(nH+1,nH+1));
                }
                sumTot = sumSpinUp+sumSpinDown;
                for(int i=0; i<N_peratom_HartrOrbit; i+=2) {
                    Sw_doublecounting_atom[at](i,i  )=     UHubb*0.5 + Uprime*(sumTot-1.0) - JHund*(sumTot/2. - 0.5)   ;
                    Sw_doublecounting_atom[at](i+1,i+1)=   UHubb*0.5 + Uprime*(sumTot-1.0) - JHund*(sumTot/2. - 0.5)  ;
                }
                ifroot std::cout << "Spin up:" << sumSpinUp << " Spin down:" << sumSpinDown <<"\n";
                Edc +=    (             UHubb * sumTot/2.
                                        +  0.5*Uprime * sumTot    * (sumTot    -2)
                                        - 0.5*JHund  * sumTot/2. * (sumTot/2. -1)
                                        - 0.5*JHund  * sumTot/2. * (sumTot/2. -1)                           );
            }//at
        }
        else if(dctype.find(std::string("fll")) != std::string::npos or dctype.find(std::string("sigma0")) != std::string::npos )  {
            /*FLL DC*/
            double sumSpinUp=0, sumSpinDown=0, sumTot;

            for(int at=0; at<NumCorrAtom; at++) {
                sumSpinUp=0;
                sumSpinDown=0;
                for(int n=0; n<N_peratom_HartrOrbit; n+=2) {
                    int nH =    at*N_peratom_HartrOrbit + n;
                    sumSpinUp   += real(NelectronDFT(nH,nH));
                    sumSpinDown += real(NelectronDFT(nH+1,nH+1));
                }
                sumTot=sumSpinUp+sumSpinDown;
                sumSpinUp = sumTot/2.;
                sumSpinDown = sumTot/2.;
                for(int i=0; i<N_peratom_HartrOrbit; i+=2) {
                    Sw_doublecounting_atom[at](i,i  )=       averU*(sumSpinDown+sumSpinUp-0.5)-averJ*(sumSpinUp-0.5);
                    Sw_doublecounting_atom[at](i+1,i+1)=     averU*(sumSpinDown+sumSpinUp-0.5)-averJ*(sumSpinDown-0.5);
                }
                ifroot std::cout << "Spin up:" << sumSpinUp << " Spin down:" << sumSpinDown <<"\n";
                Edc +=  0.5*  (    averU * (sumSpinDown+sumSpinUp) * (sumSpinDown+sumSpinUp-1)
                                   - averJ * (sumSpinDown) * (sumSpinDown-1)
                                   - averJ * (sumSpinUp)   * (sumSpinUp-1)                           );
            }//at
        }
        else if(dctype.find(std::string("nominal_amf")) != std::string::npos) {
            /*AMF DC   [ Ylvisaker and Pickett, PHYSICAL REVIEW B 79, 035103  2009 ] */
            double sumSpinUp=0, sumSpinDown=0, sumTot;

            for(int at=0; at<NumCorrAtom; at++) {
                sumSpinUp=nominal_charge/2.0;
                sumSpinDown=nominal_charge/2.0;
//                for(int n=0; n<N_peratom_HartrOrbit; n+=2) {
//                    int nH =    at*N_peratom_HartrOrbit + n;
//                    sumSpinUp   += real(NelectronDFT(nH,nH));
//                    sumSpinDown += real(NelectronDFT(nH+1,nH+1));
//                }
                double n0 = (sumSpinUp+sumSpinDown)/N_peratom_HartrOrbit;

                sumTot=sumSpinUp+sumSpinDown;
                sumSpinUp = sumTot/2.;
                sumSpinDown = sumTot/2.;
                if(nominal_charge<0) {
                    nominal_charge =  sumTot;
                }

                for(int i=0; i<N_peratom_HartrOrbit; i++) {
                    Sw_doublecounting_atom[at](i,i  )=     averU*n0*(N_peratom_HartrOrbit-1) - averJ*n0*(N_peratom_HartrOrbit/2 - 1);
                }
                ifroot std::cout << "Spin up:" << sumSpinUp << " Spin down:" << sumSpinDown <<"\n";
                Edc  += 0.5* averU                                    * (sumSpinDown+sumSpinUp  )     * (sumSpinDown+sumSpinUp  )
                        -0.5* (averU + (Norbitals-1)* averJ)/Norbitals * (sumSpinUp*sumSpinUp    )
                        -0.5* (averU + (Norbitals-1)* averJ)/Norbitals * (sumSpinDown*sumSpinDown) ;
            }//at
        }
        else if(dctype.find(std::string("amf")) != std::string::npos) {
            /*AMF DC   [ Ylvisaker and Pickett, PHYSICAL REVIEW B 79, 035103  2009 ] */
            double sumSpinUp=0, sumSpinDown=0, sumTot;

            for(int at=0; at<NumCorrAtom; at++) {
                sumSpinUp=0;
                sumSpinDown=0;
                for(int n=0; n<N_peratom_HartrOrbit; n+=2) {
                    int nH =    at*N_peratom_HartrOrbit + n;
                    sumSpinUp   += real(NelectronDFT(nH,nH));
                    sumSpinDown += real(NelectronDFT(nH+1,nH+1));
                }
                double n0 = (sumSpinUp+sumSpinDown)/N_peratom_HartrOrbit;

                sumTot=sumSpinUp+sumSpinDown;
                sumSpinUp = sumTot/2.;
                sumSpinDown = sumTot/2.;

                for(int i=0; i<N_peratom_HartrOrbit; i++) {
                    Sw_doublecounting_atom[at](i,i  )=     averU*n0*(N_peratom_HartrOrbit-1) - averJ*n0*(N_peratom_HartrOrbit/2 - 1);
                }
                ifroot std::cout << "Spin up:" << sumSpinUp << " Spin down:" << sumSpinDown <<"\n";
                Edc  += 0.5* averU                                    * (sumSpinDown+sumSpinUp  )     * (sumSpinDown+sumSpinUp  )
                        -0.5* (averU + (Norbitals-1)* averJ)/Norbitals * (sumSpinUp*sumSpinUp    )
                        -0.5* (averU + (Norbitals-1)* averJ)/Norbitals * (sumSpinDown*sumSpinDown) ;
            }//at
        }
        else if ( dctype.find(std::string("nominal")) != std::string::npos      ) {
            /*nominal DC*/
            double sumSpinUp=0, sumSpinDown=0, sumTot;
            for(int at=0; at<NumCorrAtom; at++) {
                sumSpinUp=0;
                sumSpinDown=0;
                for(int n=0; n<N_peratom_HartrOrbit; n+=2) {
                    int nH =    at*N_peratom_HartrOrbit + n;
                    sumSpinUp   += real(NelectronDFT(nH,nH));
                    sumSpinDown += real(NelectronDFT(nH+1,nH+1));
                }
                sumTot=sumSpinUp+sumSpinDown;
                if(nominal_charge<0) {
                    nominal_charge =  sumTot;
                }
//                sumSpinUp = sumTot/2.;
//                sumSpinDown = sumTot/2.;

                for(int i=0; i<N_peratom_HartrOrbit; i+=2) {
                    Sw_doublecounting_atom[at](i,i  )=     averU*(nominal_charge-0.5) - averJ*(nominal_charge/2.0 - 0.5);
                    Sw_doublecounting_atom[at](i+1,i+1)=     averU*(nominal_charge-0.5) - averJ*(nominal_charge/2.0 - 0.5);
                }
                ifroot std::cout << "Spin up:" << sumSpinUp << " Spin down:" << sumSpinDown <<"\n";
                Edc +=  0.5*averU * ( sumTot*(nominal_charge-1) + nominal_charge*(sumTot-1) - nominal_charge*(nominal_charge-1)                         )
                        -averJ * ( ( sumTot/2)*(nominal_charge/2-1) +( nominal_charge/2)*(sumTot/2-1)-( nominal_charge/2)*(nominal_charge/2-1)       );
            }//at
        }
    }
    if( dctype.find(std::string("unif")) != std::string::npos) {

        cmplx Sw_up=0, Sw_dn=0;
        for(int at=0; at<NumCorrAtom; at++) {
            Sw_up += Sw_doublecounting_atom[at](0,0) ;
            Sw_dn += Sw_doublecounting_atom[at](1,1) ;
        }
        Sw_up /= NumCorrAtom;
        Sw_dn /= NumCorrAtom;
        for(int at=0; at<NumCorrAtom; at++) {
            for(int i=0; i<N_peratom_HartrOrbit; i+=2) {
                Sw_doublecounting_atom[at](i,i  )=      Sw_up;
                Sw_doublecounting_atom[at](i+1,i+1)=    Sw_dn;
            }
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
void set_Hartree_self_energy( Eigen::MatrixXcd & NumMatrix, std::vector<Eigen::VectorXi > Uindex, std::vector<cmplx >  Utensor, ImgFreqFtn & SelfEnergy_w   ) {

    for(int b=0; b<NumCorrAtom*N_peratom_HartrOrbit; b++) {
        for(int d=0; d<NumCorrAtom*N_peratom_HartrOrbit; d++) {
            SelfEnergy_w.setValue(N_freq, b,d, 0);
        }
    }

    for(int at=0; at<NumCorrAtom; at++) {
        for(int idx=0;  idx <  Utensor.size(); idx++) {
            int aH = at* N_peratom_HartrOrbit + Uindex[idx](0);
            int bH = at* N_peratom_HartrOrbit + Uindex[idx](1);
            int cH = at* N_peratom_HartrOrbit + Uindex[idx](2);
            int dH = at* N_peratom_HartrOrbit + Uindex[idx](3);
            if (not( isOrbitalCorrinHart[aH] and
                     isOrbitalCorrinHart[bH] and
                     isOrbitalCorrinHart[cH] and
                     isOrbitalCorrinHart[dH] ) ) {
                SelfEnergy_w.setValue(N_freq, bH,dH, SelfEnergy_w.getValue(N_freq,bH,dH) + Utensor[idx] * NumMatrix(cH,aH)   );
                SelfEnergy_w.setValue(N_freq, bH,cH, SelfEnergy_w.getValue(N_freq,bH,cH) - Utensor[idx] * NumMatrix(dH,aH)   );
            }
        }
    }
    for(int n=0; n<N_freq; n++) {
        for(int b=0; b<NumCorrAtom*N_peratom_HartrOrbit; b++) {
            for(int d=0; d<NumCorrAtom*N_peratom_HartrOrbit; d++) {
                SelfEnergy_w.setValue(n, b,d, SelfEnergy_w.getValue(N_freq,b,d));
            }
        }

    }
    FILE * DC = fopen("SwHartree.dat","w");
    for (int orb1=0; orb1<N_peratom_HartrOrbit*NumCorrAtom; orb1++) {
        for (int orb2=0; orb2<N_peratom_HartrOrbit*NumCorrAtom; orb2++) {
            fprintf(DC, "     %0.5f %0.5f",real(SelfEnergy_w.getValue(N_freq,orb1,orb2)), imag(SelfEnergy_w.getValue(N_freq,orb1,orb2)));
        }
        fprintf(DC, "\n");
    }
    fclose(DC);
    sleep(1);
}
///////////////////////////////////////////////////////////////////////////////////////////////////////




///////////////////////////////////////////////////////////////////////////////////////////////////////
void write_restart (int DFTIt, int currentIt, std::string system_name, int NumCorrAtom) {
    ifroot {
        FILE * DC;
        std::string cp_comm;
        char chr[100];
//        system("cp   Sw_SOLVER.full.dat Sw_SOLVER_weak.full.dat  Numele.dat  mu_history.out DM_DFT_DMFT.HWR dual_DM_direct.HWR muDFT.out   ./Restart");
        system("cp   Sw_SOLVER.full.dat Sw_SOLVER_weak.full.dat  Numele.dat  mu_history.out DM_DFT_DMFT.HWR   ./Restart");
    }
}
void write_results (int DFTIt, int currentIt, std::string system_name, int NumCorrAtom) {
    ifroot {
        FILE * DC;
        std::string cp_comm;
        char chr[100];
//        system("cp   Sw_SOLVER.full.dat Sw_SOLVER_weak.full.dat  Numele.dat  mu_history.out DM_DFT_DMFT.HWR dual_DM_direct.HWR muDFT.out   ./Restart");
        system("cp   Sw_SOLVER.full.dat Sw_SOLVER_weak.full.dat  Numele.dat  mu_history.out DM_DFT_DMFT.HWR  ./Restart");



        sprintf(chr,"%d_%d",DFTIt,currentIt);
        cp_comm = std::string("mkdir  ./") +system_name+ "/" + chr ;
        system(cp_comm.c_str());


        for(int at=0; at<NumCorrAtom; at++) {
            cp_comm = std::string("cp Gw_loc.dat") + intToString(at) +  std::string(" ./") +system_name+ "/" + chr + "/" ;
            system(cp_comm.c_str());

            cp_comm = std::string("cp ctqmc.log") + intToString(at) +  std::string(" ./") +system_name+ "/" + chr + "/" ;
            system(cp_comm.c_str());

            cp_comm = std::string("cp Probability.dat") + intToString(at) +  std::string(" ./") +system_name+ "/" + chr + "/" ;
            system(cp_comm.c_str());
        }



        cp_comm = std::string("cp delta_t.dat  ./") +system_name+ "/" + chr + "/";
        system(cp_comm.c_str());

        cp_comm = std::string("cp delta_w.dat  ./") +system_name+ "/" + chr + "/";
        system(cp_comm.c_str());
        cp_comm = std::string("cp delta_weak_w.dat  ./") +system_name+ "/" + chr + "/";
        system(cp_comm.c_str());
        cp_comm = std::string("cp projdelta_w.dat* ./") +system_name+ "/" + chr + "/";
        system(cp_comm.c_str());
        cp_comm = std::string("cp Numele.dat  ./") +system_name+ "/" + chr + "/";
        system(cp_comm.c_str());

        cp_comm = std::string("cp Sw_SOLVER.dat  ./") +system_name+ "/" + chr + "/";
        system(cp_comm.c_str());

        cp_comm = std::string("cp Bare_Sw_lowlevel.dat*  ./") +system_name+ "/" + chr + "/";
        system(cp_comm.c_str());
        cp_comm = std::string("rm -r ")+ "  ./"+system_name+ "/Restart";
        system(cp_comm.c_str());
        cp_comm = std::string("cp -r ")+ "./Restart"    + "  ./"+system_name+ "/Restart";
        system(cp_comm.c_str());

        if(SOLVERtype==std::string("ALPS_CTSEG")) {
            cp_comm = std::string("cp Gw.dat  ./") +system_name+ "/" + chr + "/";
            system(cp_comm.c_str());
            cp_comm = std::string("cp orders.dat  ./") +system_name+ "/" + chr + "/";
            system(cp_comm.c_str());
            cp_comm = std::string("cp Gwl.dat  ./") +system_name+ "/" + chr + "/";
            system(cp_comm.c_str());
            cp_comm = std::string("cp Gt.dat  ./") +system_name+ "/" + chr + "/";
            system(cp_comm.c_str());
            cp_comm = std::string("cp Gtl.dat  ./") +system_name+ "/" + chr + "/";
            system(cp_comm.c_str());

            cp_comm = std::string("cp observables.dat  ./") +system_name+ "/" + chr +"/";
            system(cp_comm.c_str());
            cp_comm = std::string("cp simulation.dat  ./") +system_name+ "/" + chr +"/";
            system(cp_comm.c_str());
            cp_comm = std::string("cp sector_statistics.dat ./") +system_name+ "/" + chr +"/";
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

        for(int at=0; at<NumCorrAtom; at++) {
            std::cout << "Num ele (decomp,solver)  = ";
            double sum=0;
            for(int n=0; n<N_peratom_HartrOrbit; n+=1) {
                std::cout <<std::fixed << std::setprecision(4)<< std::fixed   <<  real(NumMatrixImp(at*N_peratom_HartrOrbit+n,at*N_peratom_HartrOrbit+n)) <<" " ;
                sum+= real(NumMatrixImp(at*N_peratom_HartrOrbit+n,at*N_peratom_HartrOrbit+n));
            }
            std::cout << std::fixed << std::setprecision(4)<< std::fixed   << ";  (total)  = "<<  sum <<"\n" ;
        }


        // Write results
//        write_results(DFTIt, currentIt, system_name, NumCorrAtom);
        write_restart(DFTIt, currentIt, system_name, NumCorrAtom);
        timeEndIt = clock();
        std::cout <<"Computaional time for single SC-loop:\n"
                  << ((timeEndIt-timeStartIt)/(CLOCKS_PER_SEC))/3600
                  <<":"<< (((timeEndIt-timeStartIt)/(CLOCKS_PER_SEC))%3600)/60
                  <<":"<<(((timeEndIt-timeStartIt)/(CLOCKS_PER_SEC))%3600)%60 <<  "\n";
        std::cout <<"using "  <<mpi_numprocs<< " processors\n";
    }//ifroot
    MPI_Barrier(MPI_COMM_WORLD);


    /*SCF check; */
    ImgFreqFtn GwDMFT_Loc(beta, N_freq, N_peratom_HartrOrbit, NumCorrAtom,0);
    ImgFreqFtn GwDMFT_Imp(beta, N_freq, N_peratom_HartrOrbit, NumCorrAtom,0);
//        std::cout <<"using "  <<mpi_numprocs<< " processors\n";
    if(N_peratom_HartrOrbit>0) GwDMFT_Loc.update_full(std::string("Gw_loc.full.dat0"),1);
    if(N_peratom_HartrOrbit>0) GwDMFT_Imp.update_full(std::string("Gw_imp.full.dat0"),1);
    double        GwNorm=1e-10;
    double        GwlowNorm=1e-10;
    double        GwRD=0;
    double        GwRD_low=0;
    double        GwRD_max=0;
    double GwRDMAXnorm=1;
    double NumMatRD = 0;
//        std::cout <<"using "  <<mpi_numprocs<< " processors\n";
    for (int n=0; n<N_freq; n++) {
        double GwRD_n=0, GwNorm_n=1e-10;
        for (int i=0; i<N_peratom_HartrOrbit; i++) {
            for (int j=0; j<N_peratom_HartrOrbit; j++) {
                if  (GwDMFT_Imp.getValue(n) < std::min( beta*UHubb/3.0, N_freq/2.0) ) {
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
//        std::cout <<"using "  <<mpi_numprocs<< " processors\n";
    GwRD = std::sqrt(GwRD/GwNorm);
    GwRD_low = std::sqrt(GwRD_low/GwlowNorm);
    GwRD_max = std::sqrt(GwRD_max);
//    NumMatRD = ((NumMatrixLatt-NumMatrixImp).norm()/NumMatrixLatt.norm());
    NumMatRD = std::abs(real((NumMatrixLatt-NumMatrixImp).trace()/NumCorrAtom));


    ifroot printf(  "GwRD       =%e\n", GwRD);
    ifroot printf(  "GwRD_low   =%e\n", GwRD_low);
    ifroot printf(  "GwRD_max   =%e\n", GwRD_max);
//    ifroot printf(  "muRD_fluct =%e\n", muRDFluct);
    ifroot printf(  "NumMatRD   =%e\n", NumMatRD);

    if(      GwRD        < 1e-3
             and GwRD_low    < 1e-3
             and GwRD_max    < 1e-3
             and NumMatRD    <  1e-3
      ) {
        return true;
    }
    else return false;
}
