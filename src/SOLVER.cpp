//#include <time.h>
#include "ImgTimeFtn.h"
#include "tight_common.h"
#include <cstring>
#include <stdlib.h>
#include <Eigen/LU>
#include <Eigen/Eigenvalues>
#ifdef READ_HDF5_dmft
#include <H5Cpp.h>
#endif


void rot_Uijkl(
    std::vector<cmplx > & Utensor, std::vector<Eigen::VectorXi>  & Uindex,
    std::vector<cmplx > & rotUtensor, std::vector<Eigen::VectorXi>  & rotUindex,
    Eigen::MatrixXcd & SolverBasis, int n_spinorb
) ;
void rot_Uijkl_dd(
    std::vector<cmplx > & Utensor, std::vector<Eigen::VectorXi>  & Uindex,
    std::vector<cmplx > & rotUtensor, std::vector<Eigen::VectorXi>  & rotUindex,
    Eigen::MatrixXcd & SolverBasis, int n_spinorb
) ;

void IPT( int solverDim,Eigen::MatrixXcd projimpurity_site_Hamiltonian,  ImgFreqFtn & weiss_field, Eigen::MatrixXcd & projNumMatrix,
          ImgFreqFtn & SE_out,      ImgFreqFtn & Gwimp_out, double muTB,
          std::vector<Eigen::VectorXi> projUindex, std::vector<cmplx > projUtensor        ) ;
//

void on_shot_HF (int solverDim,
                 ImgFreqFtn & SE_out,      ImgFreqFtn & Gwimp_out,
                 std::vector<Eigen::VectorXi> projUindex, std::vector<cmplx > projUtensor
                ) ;

void SC2PT_weak  ( int solverDim,
                   ImgFreqFtn & SE_out,      ImgFreqFtn & Gwimp_out, ImgFreqFtn & GwHF, bool SC,
                   std::vector<Eigen::VectorXi> projUindex, std::vector<cmplx > projUtensor);



void SecondOrderPerturbation  ( int solverDim,   Eigen::MatrixXcd projimpurity_site_Hamiltonian,  ImgFreqFtn & weiss_field, Eigen::MatrixXcd & projNumMatrix,
                                ImgFreqFtn & SE_out,      ImgFreqFtn & Gwimp_out, double muTB,
                                std::vector<Eigen::VectorXi> projUindex, std::vector<cmplx > projUtensor
                              ) ;


void SCHF (int solverDim, Eigen::MatrixXcd projimpurity_site_Hamiltonian, Eigen::MatrixXcd & projNumMatrix,
           ImgFreqFtn & SE_out,      ImgFreqFtn & Gwimp_out,  ImgFreqFtn & weiss_field, double muTB,
           std::vector<Eigen::VectorXi> projUindex, std::vector<cmplx > projUtensor
          );

void SCGF2 (int solverDim, Eigen::MatrixXcd projimpurity_site_Hamiltonian, Eigen::MatrixXcd & projNumMatrix,
            ImgFreqFtn & SE_out,      ImgFreqFtn & Gwimp_out,  ImgFreqFtn & weiss_field, double muTB,
            std::vector<Eigen::VectorXi> projUindex, std::vector<cmplx > projUtensor);


void ctqmc_rutgers(  Eigen::MatrixXcd Local_Hamiltonian_ED, double muTB,  ImgFreqFtn & weiss_field, std::vector<cmplx > Utensor, std::vector<Eigen::VectorXi> Uindex, int solverDim  ) ;
void ctqmc_rutgers_seg(  Eigen::MatrixXcd Local_Hamiltonian_ED, double muTB,  ImgFreqFtn & weiss_field, std::vector<cmplx > Utensor, std::vector<Eigen::VectorXi> Uindex, int solverDim) ;






void proj_to_site( int solverDim, int solver_block, std::vector<int> impurityOrbit,
                   Eigen::MatrixXcd  impurity_site_Hamiltonian,  Eigen::MatrixXcd NumMatrix,  ImgFreqFtn & weiss_field,
                   Eigen::MatrixXcd & projimpurity_site_Hamiltonian,Eigen::MatrixXcd & projSolverBasis,Eigen::MatrixXcd & projNumMatrix, ImgFreqFtn & projweiss_field,
                   Eigen::MatrixXcd  Sw_doublecounting, std::vector<Eigen::MatrixXcd >  dc_weakCorr
                 ) ;



void weak_solver(
    std::string SOLVERtype,
    int solverDim,
    ImgFreqFtn & SE_out,  ImgFreqFtn & Gwimp_in_out, ImgFreqFtn & GwHF,
    std::vector<Eigen::VectorXi> projUindex, std::vector<cmplx > projUtensor
)
{
    if(SOLVERtype == std::string("HF")) {

        on_shot_HF( solverDim,
                    SE_out, Gwimp_in_out,
                    projUindex, projUtensor);

    }//SCHF
    else if(SOLVERtype == std::string("2PT")) {
        ifroot std::cout << "\n2PT solver\n";
        SC2PT_weak( solverDim,
                    SE_out, Gwimp_in_out, GwHF, false,
                    projUindex, projUtensor);
    }//2PT
    else if(SOLVERtype == std::string("SC2PT")   ) {
        ifroot std::cout << "\nSC2PT solver\n";
        SC2PT_weak( solverDim,
                    SE_out, Gwimp_in_out, GwHF, true,
                    projUindex, projUtensor);
    }//SC2PT
    else {
        ifroot std::cout << "Please set  Lowlevel_SOLVERTYPE = HF or SC2PT\n";
    }

}


void SOLVER(
    std::string SOLVERtype,
    int solverDim, int solver_block,  std::vector<int> impurityOrbit,
    Eigen::MatrixXcd  impurity_site_Hamiltonian,  ImgFreqFtn & weiss_field,
    ImgFreqFtn & SE_out,  ImgFreqFtn & Gwimp_in_out,
    std::vector<Eigen::VectorXi> projUindex, std::vector<cmplx > projUtensor,
    Eigen::MatrixXcd  Sw_doublecounting,            std::vector<Eigen::MatrixXcd >  dc_weakCorr,
    double muTB
)
{


    ifroot std::cout << "<SOLVER>\n";


    Eigen::MatrixXcd   projimpurity_site_Hamiltonian, projSolverBasis, projNumMatrix ;
    ImgFreqFtn projweiss_field(0);
    projweiss_field.Initialize(beta, N_freq, solverDim,1,0);

    projimpurity_site_Hamiltonian.setZero(solverDim, solverDim);
    projNumMatrix.setZero(solverDim, solverDim);
    projSolverBasis.setIdentity(solverDim, solverDim);

    proj_to_site ( solverDim, solver_block, impurityOrbit,
                   impurity_site_Hamiltonian, NumMatrix, weiss_field,
                   projimpurity_site_Hamiltonian, projSolverBasis, projNumMatrix, projweiss_field,
                   Sw_doublecounting, dc_weakCorr);







    if ( SOLVERtype.find(std::string("CT")) != std::string::npos ) {
        /*diagonalize onsite-Hamiltonian*/
        //SOLVER basis transformation
        Eigen::MatrixXcd  projimpurity_site_Hamiltonian_diag ;
        if ( SOLVERtype.find(std::string("SEG")) != std::string::npos ) {
            projimpurity_site_Hamiltonian_diag.setZero(solverDim,solverDim);
            for(int h1=0; h1<solverDim; h1++) {
                projimpurity_site_Hamiltonian_diag(h1,h1) = ( projimpurity_site_Hamiltonian )(h1,h1);
            }
        }


        if (SOLVERtype==std::string("ALPS_CTSEG")) {//ALPS_CTSEG
            /*hyb_t write*/
            Eigen::MatrixXcd * delta_t  = new Eigen::MatrixXcd  [N_tau+1];
            write_hyb_t (projweiss_field, delta_t, 0);



            /*write mu_alps for ALPS SOLVER.*/
            double mu_alps[solverDim];
            for(int h1=0; h1<solverDim; h1++) {
                mu_alps[h1] = muTB -  real( projimpurity_site_Hamiltonian_diag(h1,h1) ) ;
            }

            if (mpi_rank==0) {
                FILE *tempFile2;   //mu.dat
                tempFile2 = fopen("mu_vector.alps.dat", "w");
                for(int h1=0; h1<solverDim; h1++) {
                    fprintf(tempFile2, " %0.8f ", mu_alps[h1]);
                }
                fclose(tempFile2);

                FILE *datap4 = fopen("delta_t.dat", "w");
                Eigen::MatrixXcd    temp;
                for(int tau=0; tau<=N_tau; tau++) {
                    fprintf(datap4, "%0.5f",tau*beta/N_tau);
                    for(int h1=0; h1<solverDim; h1++) {
                        fprintf(datap4, "    %0.10f", real(  delta_t[tau](h1,h1)));
                    }
                    fprintf(datap4, "\n");
                }
                fclose(datap4);
            }
            sleep(3);
            MPI_Barrier(MPI_COMM_WORLD);

            /*set solver command*/
            ifroot std::cout << "ALPS, CT-SEG solver\n";
            MPI_Comm communication;

            std::ostringstream impSolver_hyb_comm;
            char *hyb_In [] = {"alps_input.h5", NULL};
//            impSolver_hyb_comm << (solver_bin_dir + std::string("hybridization")).c_str();   //   alps_input.h5 >> std.hyb.out"; C2
            impSolver_hyb_comm << ( SOLVERdir + SOLVERexe).c_str();   //   alps_input.h5 >> std.hyb.out"; C2
            char * command_hyb = new char[impSolver_hyb_comm.str().size()+1];
            std::strcpy(command_hyb,impSolver_hyb_comm.str().c_str());

            /*neglect possitive data*/
            int t, CHECK, t1;
            ImgTimeFtn DeltDiag(beta,N_tau, solverDim);
            DeltDiag.update(std::string ("delta_t.dat"), 1);

            for(t=0; t<N_tau+1; t++) {
                for(int n=0; n<DeltDiag.getNorbit(); n++) {
                    if(DeltDiag.getValue(t,n)>=0 and (t!=0 and t!=N_tau) ) {
                        CHECK=1;
                        for( t1=t+1; t1<N_tau+1; t1++) {
                            DeltDiag.setValue(t,n, ( ((DeltDiag.getValue(t1,n) - DeltDiag.getValue(t-1,n))/( DeltDiag.getValue(t1) - DeltDiag.getValue(t-1)) )
                                                     *(DeltDiag.getValue(t)-DeltDiag.getValue(t-1)) + DeltDiag.getValue(t-1,n))  );
                            if(DeltDiag.getValue(t,n) < 0)     break;
                        }//t1
                    }
                    if(DeltDiag.getValue(t,n)>=0  ) {
                        CHECK=1;
                        DeltDiag.setValue(t,n, -1e-10);
                    }
                }//n
            }
            ifroot   if(CHECK==1)      std::cout <<" Warning : positive DeltDiag\n"   ;

            /*run alps solver*/
            ifroot   DeltDiag.dataOut("delta_t.dat",NumAtom);
            ifroot  std::cout << "****Solver type=ALPS_CTSEG....***\n";
            ifroot  std::cout <<  impSolver_hyb_comm.str().c_str() <<"\n" ;

            int errors[mpi_numprocs], checkTime=10;
            MPI_Comm_spawn(command_hyb, hyb_In, mpi_numprocs,MPI_INFO_NULL, 0, MPI_COMM_WORLD, &communication, errors);
            Wait_Run("Sw.dat", checkTime, mpi_rank, maxTime);


            /*read output*/
            ImgTimeFtn  Gt_imp(beta, N_tau, solverDim);

            Gt_imp.update(std::string("Gt.dat"), 1);
            SE_out.read_diag(std::string("Swl.dat"));
            Gwimp_in_out.read_diag(std::string("Gw.dat"));

            MPI_Barrier(MPI_COMM_WORLD);
            ifroot  std::cout << "Reading output data..\n";

            for (int l=0; l<solverDim; l++) {
                projNumMatrix(l,l) = - Gt_imp.getValue(N_tau,l) ;
            }

            ifroot  std::cout << "Writing output data..\n";
            delete[] command_hyb;
            delete [] delta_t;
        }//ALPS_CTSEG
        else if (SOLVERtype==std::string("RUTGERS_CTSEG")) {//RUTGERS_CTSEG
            ifroot std::cout << "Rutgers(K. Haule), CT-HYB(SEG) solver\n";

            Eigen::MatrixXcd * delta_t   = new Eigen::MatrixXcd  [N_tau+1];
            write_hyb_t (projweiss_field, delta_t, 0);
            ifroot std::cout <<"FILE OUT: delta_t.dat\n";
            ifroot{
                FILE *datap4 = fopen("delta_t.dat", "w");
                for(int tau=0; tau<=N_tau; tau++) {
                    fprintf(datap4, "%0.5f",tau*beta/N_tau);
                    for(int h1=0; h1<solverDim; h1++) {
                        fprintf(datap4, "    %0.10f", real(  delta_t[tau](h1,h1)));
                    }
                    fprintf(datap4, "\n");
                }
                fclose(datap4);
            }
            sleep(2);


            std::vector<cmplx > rotUtensor;
            std::vector<Eigen::VectorXi> rotUindex;
            rot_Uijkl_dd(projUtensor, projUindex, rotUtensor, rotUindex, projSolverBasis, solverDim);

            ctqmc_rutgers_seg(  projimpurity_site_Hamiltonian_diag,  muTB, projweiss_field, rotUtensor, rotUindex, solverDim   );

            /*set solver command*/
            MPI_Comm communication;

//            sleep(5);
            std::ostringstream impSolver_hyb_comm;
            char *hyb_In [] = {"PARAMS", NULL};
            impSolver_hyb_comm << ( SOLVERdir + SOLVERexe).c_str();   //   alps_input.h5 >> std.hyb.out"; C2
            char * command_hyb = new char[impSolver_hyb_comm.str().size()+1];
            std::strcpy(command_hyb,impSolver_hyb_comm.str().c_str());

            /*run alps solver*/
            ifroot  std::cout << "****Solver type=Rutgers CT-HYB(SEG)  ***\n";
            ifroot  std::cout <<  impSolver_hyb_comm.str().c_str() <<"\n" ;

            int errors[mpi_numprocs], checkTime=10;
            MPI_Comm_spawn(command_hyb, hyb_In, mpi_numprocs,MPI_INFO_NULL, 0, MPI_COMM_WORLD, &communication, errors);
            Wait_Run("Gcoeff.dat", checkTime, mpi_rank, maxTime);


            /*read output*/
            SE_out.read_diag(std::string("Sw.dat"));
            Gwimp_in_out.read_diag(std::string("Gw.dat"));
            MPI_Barrier(MPI_COMM_WORLD);
            ifroot  std::cout << "Reading output data..\n";


            Eigen::MatrixXcd Gwimp_EIG [N_freq];
            for(int n=0; n<N_freq; n++) {
                Gwimp_EIG[n].setZero(NSpinOrbit_per_atom,NSpinOrbit_per_atom);
                for(int h1=0; h1<NSpinOrbit_per_atom; h1++) {
                    for(int h2=0; h2<NSpinOrbit_per_atom; h2++) {
                        Gwimp_EIG[n](h1,h2) = Gwimp_in_out.getValue(n,h1,h2);
                    }
                }//h1
            }
            Eigen::MatrixXcd Id;
            Id.setIdentity(solverDim,solverDim);
            std::string print_msg = std::string("Gwimp_EIG in SOLVER.cpp");
            FourierTransform (Gwimp_EIG, projNumMatrix, beta, Id) ;
            projNumMatrix = -projNumMatrix;

            ifroot  std::cout << "Writing output data..\n";
            std::stringstream ss;
            ss << solver_block;
            ifroot system(  (std::string("cp ")+"rutgers_input.cix  rutgers_input.cix" +ss.str()).c_str());
            delete[] command_hyb;
            delete [] delta_t;
        }//RUTGERS_CTSEG
        else if (SOLVERtype==std::string("RUTGERS_CTHYB")) {//Rutgers
            /*diagonalize onsite-Hamiltonian*/
            //SOLVER basis transformation

            std::vector<cmplx > rotUtensor;
            std::vector<Eigen::VectorXi> rotUindex;
            rot_Uijkl_dd(projUtensor, projUindex, rotUtensor, rotUindex, projSolverBasis, solverDim);


            ifroot "write rutgers_input file\n";
            ctqmc_rutgers(  projimpurity_site_Hamiltonian,  muTB, projweiss_field, rotUtensor, rotUindex, solverDim);

            /*set solver command*/
            ifroot std::cout << "Rutgers(K. Haule), CT-HYB solver\n";

            sleep(3);
            MPI_Comm communication;
            std::ostringstream impSolver_hyb_comm;
            char *hyb_In [] = {"PARAMS", NULL};
            impSolver_hyb_comm << ( SOLVERdir + SOLVERexe).c_str();   //   alps_input.h5 >> std.hyb.out"; C2
            char * command_hyb = new char[impSolver_hyb_comm.str().size()+1];
            std::strcpy(command_hyb,impSolver_hyb_comm.str().c_str());

            /*run alps solver*/
            ifroot  std::cout << "****Solver type=Rutgers CT-HYB  ***\n";
            ifroot  std::cout <<  impSolver_hyb_comm.str().c_str() <<"\n" ;

            int errors[mpi_numprocs], checkTime=10;
            MPI_Comm_spawn(command_hyb, hyb_In, mpi_numprocs,MPI_INFO_NULL, 0, MPI_COMM_WORLD, &communication, errors);
            Wait_Run("Gcoeff.dat", checkTime, mpi_rank, maxTime);


            /*read output*/

            int spindim=1;
            if(magnetism==0 or magnetism==1) spindim=2;
            SE_out.read_uppertrian(std::string("Sw.dat") , spindim   );
            Gwimp_in_out.read_uppertrian(std::string("Gw.dat") , spindim );


            MPI_Barrier(MPI_COMM_WORLD);
            ifroot  std::cout << "Reading output data..\n";



            Eigen::MatrixXcd Gwimp_EIG [N_freq];
            for(int n=0; n<N_freq; n++) {
                Gwimp_EIG[n].setZero(NSpinOrbit_per_atom,NSpinOrbit_per_atom);
                for(int h1=0; h1<NSpinOrbit_per_atom; h1++) {
                    for(int h2=0; h2<NSpinOrbit_per_atom; h2++) {
                        Gwimp_EIG[n](h1,h2) = Gwimp_in_out.getValue(n,h1,h2);
                    }
                }//h1
            }
            Eigen::MatrixXcd Id;
            Id.setIdentity(solverDim,solverDim);
            std::string print_msg = std::string("Gwimp_EIG in SOLVER.cpp");
            FourierTransform (Gwimp_EIG, projNumMatrix, beta, Id) ;
            projNumMatrix *= -1;



            ifroot  std::cout << "Writing output data..\n";
            std::stringstream ss;
            ss << solver_block;
            system(  (std::string("cp ")+"rutgers_input.cix  rutgers_input.cix" +ss.str()).c_str());
            delete[] command_hyb;

        }//Rutgers

    }
    else if(SOLVERtype == std::string("SC2PT")) {

        SCGF2( solverDim,projimpurity_site_Hamiltonian, projNumMatrix,
               SE_out, Gwimp_in_out, projweiss_field, muTB,
               projUindex, projUtensor);

    }//SC2PT
    else if(SOLVERtype == std::string("SCHF")) {

        SCHF( solverDim,projimpurity_site_Hamiltonian, projNumMatrix,
              SE_out, Gwimp_in_out, projweiss_field, muTB,
              projUindex, projUtensor);

    }//SCHF
    else if(SOLVERtype == std::string("2PT")) {
        ifroot std::cout << "\n2PT solver\n";
        SecondOrderPerturbation( solverDim, projimpurity_site_Hamiltonian, projweiss_field, projNumMatrix,SE_out, Gwimp_in_out, muTB, projUindex, projUtensor);
    }//2PT
    else if(SOLVERtype == std::string("IPT")) {
        ifroot std::cout << "\nIPT solver\n";
        IPT( solverDim, projimpurity_site_Hamiltonian, projweiss_field, projNumMatrix,SE_out, Gwimp_in_out, muTB, projUindex, projUtensor);
    }//IPT
    else if(SOLVERtype == std::string("HF")) {

        on_shot_HF( solverDim,
                    SE_out, Gwimp_in_out,
                    projUindex, projUtensor);
    }//HF
    else {
        std::cout << "Please check <SOLVER_TYPE> \n";
        exit(1);
    }

    for(int n =0; n<N_freq; n++) {
        SE_out.setMatrix(n,  projSolverBasis* SE_out.getMatrix(n) *(projSolverBasis).adjoint());
    }

    for(int n =0; n<N_freq; n++) {
        Gwimp_in_out.setMatrix(n,  projSolverBasis* Gwimp_in_out.getMatrix(n) *(projSolverBasis).adjoint());
    }

    projNumMatrix = projSolverBasis * projNumMatrix * (projSolverBasis.adjoint());



///end of solver main part and ...

    if ( magnetism == 0) {
        for (int n=0; n<solverDim; n+=2) {
            int upH = n, dnH=n+1;
            cmplx avgN = (projNumMatrix(upH,upH)+ projNumMatrix(dnH,dnH))/2.;
            projNumMatrix(upH,upH) = avgN;
            projNumMatrix(dnH,dnH) = avgN;
            projNumMatrix(upH,dnH) = 0;
            projNumMatrix(dnH,upH) = 0;
            for ( int w=0; w<N_freq+5; w++) {
                cmplx avgS = ( SE_out.getValue(w,n,n) + SE_out.getValue(w,n+1,n+1) ) /2.; //alps CT-SEG,
                SE_out.setValue(w,n,n,avgS);
                SE_out.setValue(w,n+1,n+1,avgS);
                SE_out.setValue(w,n,n+1,0);
                SE_out.setValue(w,n+1,n,0);
            }//w
        }//n
    }//magnetism
    for (int n=0; n<solverDim; n++) {
        for (int m=0; m<solverDim; m++) {
            int n0 = impurityOrbit[n] ;
            int m0 = impurityOrbit[m] ;
            NumMatrix(n0,m0) = projNumMatrix(n,m);  // alps CT-SEG num matrix.
        }
    }





//
//    //Asymto, S0 + S1/iw + ...
//    //
//    //S0 part
//    Eigen::MatrixXcd tail [3];
//    Eigen::MatrixXcd tail_coeff;
//    Eigen::MatrixXcd temp_eig [3];
//    double w0,w1,w2;
//    tail[0] = SE_out.getMatrix(N_freq-3);
//    tail[1] = SE_out.getMatrix(N_freq-2);
//    tail[2] = SE_out.getMatrix(N_freq-1);
//    for(int w=0; w<3; w++) {
//        tail[w] = (tail[w]+tail[w].adjoint()).eval()/2.;
//    }
//    tail_coeff = ( tail[0] + tail[1] + tail[2]) /3.;
//    SE_out.setMatrix(N_freq, tail_coeff);
//
//
//    //S1,
//    w0 = 1./std::pow(SE_out.getValue(N_freq-3),2);
//    w1 = 1./std::pow(SE_out.getValue(N_freq-2),2);
//    w2 = 1./std::pow(SE_out.getValue(N_freq-1),2);
//    tail[0] = SE_out.getMatrix(N_freq-3);
//    tail[1] = SE_out.getMatrix(N_freq-2);
//    tail[2] = SE_out.getMatrix(N_freq-1);
//    for(int w=0; w<3; w++) {
//        tail[w] = (tail[w]-tail[w].adjoint()).eval()/2.;
//        tail[w] = tail[w] / (I*SE_out.getValue(N_freq-3+w));
//    }
//    tail_coeff = -( tail[0] + tail[1] + tail[2]) / (w0+w1+w2);
//    SE_out.setMatrix(N_freq+1, tail_coeff);


}


void proj_to_site( int solverDim, int solver_block, std::vector<int> impurityOrbit,
                   Eigen::MatrixXcd  impurity_site_Hamiltonian,  Eigen::MatrixXcd NumMatrix,  ImgFreqFtn & weiss_field,
                   Eigen::MatrixXcd & projimpurity_site_Hamiltonian,Eigen::MatrixXcd & projSolverBasis,Eigen::MatrixXcd & projNumMatrix, ImgFreqFtn & projweiss_field,
                   Eigen::MatrixXcd  Sw_doublecounting, std::vector<Eigen::MatrixXcd >  dc_weakCorr
                 ) {



    ifroot std::cout << "proj_to_site\n";
    for(int h1=0; h1< solverDim; h1++) {
        for(int h2=0; h2< solverDim; h2++) {
            int h1F = impurityOrbit.at(h1) ;
            int h2F = impurityOrbit.at(h2) ;
            projimpurity_site_Hamiltonian(h1,h2) = impurity_site_Hamiltonian(h1F,h2F) - Sw_doublecounting(h1F,h2F) +dc_weakCorr.at(N_freq)(h1,h2) ;
            projNumMatrix(h1,h2) = NumMatrix(h1F,h2F);
        }
    }

    ifroot std::cout << "impurity site info\n";
    for (int n=0; n<N_freq; n++) {
        projweiss_field.setMatrix(n,  weiss_field.getMatrix(n,solver_block, solverDim ) +dc_weakCorr.at(n)    );
    }

    ifroot std::cout << "Imp H0:\n"  << std::fixed << std::setprecision(6)<< projimpurity_site_Hamiltonian <<"\n";
    ifroot std::cout << "Re[D([w=0)]:\n" << std::fixed << std::setprecision(6)<<( (projweiss_field.getMatrix(0)) + (projweiss_field.getMatrix(0)).adjoint() ) /2 <<"\n";
    std::cout << std::fixed << std::setprecision(4);


    Eigen::MatrixXcd projimp_off  =   projimpurity_site_Hamiltonian.diagonal().asDiagonal();
    Eigen::MatrixXcd weiss_off  =  (projweiss_field.getMatrix(0)).diagonal().asDiagonal();
    double projimp_off_norm =  (projimpurity_site_Hamiltonian - projimp_off).norm();
    double weiss_off_norm   =  (projweiss_field.getMatrix(0) - weiss_off).norm();




    projSolverBasis.setIdentity(solverDim, solverDim);
    if(impurityBasisSwitch) {
        std::cout << std::fixed << std::setprecision(4);

        if(  ((  projimpurity_site_Hamiltonian  ).imag()).norm()  <  1e-5 and (magnetism==0 or magnetism==1) ) {
            Eigen::MatrixXcd weiss0_re(solverDim/2, solverDim/2);
            for (int i=0; i<solverDim; i+=2) {
                for (int j=0; j<solverDim; j+=2) {
                    cmplx temp = ( (projweiss_field.getMatrix(0))(i,j)  + (projweiss_field.getMatrix(0))(i+1,j+1) ) /2;
                    temp += (projimpurity_site_Hamiltonian(i,j) + projimpurity_site_Hamiltonian(i+1, j+1)) /2;
                    weiss0_re( i/2, j/2   ) = temp;
                }
            }//n
            weiss0_re  = (weiss0_re+weiss0_re.adjoint()).eval();
            weiss0_re /= 2.0;

            Eigen::MatrixXd temp =  ((weiss0_re).real());
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ces( solverDim/2 );
            ces.compute(    temp   );

            for (int i=0; i<solverDim; i+=2) {
                for (int j=0; j<solverDim; j+=2) {
                    projSolverBasis(i ,j ) = (ces.eigenvectors())(i/2,j/2);
                    projSolverBasis(i+1,j+1) = (ces.eigenvectors())(i/2,j/2);
                }
            }//i
        }
        else {
            Eigen::MatrixXcd weiss0_reSOC =  projweiss_field.getMatrix(0) ;
            weiss0_reSOC = (weiss0_reSOC + weiss0_reSOC.adjoint()).eval();
            weiss0_reSOC /= 2.0;

            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> ces( solverDim );
            ces.compute(    projimpurity_site_Hamiltonian   );
            projSolverBasis = ces.eigenvectors();



            Eigen::MatrixXcd temp =  projSolverBasis.adjoint() * (projimpurity_site_Hamiltonian + weiss0_reSOC) * projSolverBasis;
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ces2( solverDim );
            ces2.compute(    temp.real()  );

            projSolverBasis = (  projSolverBasis *  ces2.eigenvectors()  ).eval();

        }
    }



    ifroot std::cout << "Solver:Basis:\n" << projSolverBasis<<"\n";
    projimpurity_site_Hamiltonian = (projSolverBasis.adjoint() * projimpurity_site_Hamiltonian* projSolverBasis);
    for(int n =0; n<N_freq; n++) {
        projweiss_field.setMatrix(n,  (projSolverBasis).adjoint()* projweiss_field.getMatrix(n) *(projSolverBasis));
    }

    ifroot std::cout << "Imp H0:\n"  << std::fixed << std::setprecision(6)<< projimpurity_site_Hamiltonian <<"\n";
    ifroot std::cout << "Re[D([w=0)]:\n" << std::fixed << std::setprecision(6)<<( (projweiss_field.getMatrix(0)) + (projweiss_field.getMatrix(0)).adjoint() ) /2 <<"\n";
    projimp_off  =   projimpurity_site_Hamiltonian.diagonal().asDiagonal();
    weiss_off  =  (projweiss_field.getMatrix(0)).diagonal().asDiagonal();
    double projimp_off_norm2 =  (projimpurity_site_Hamiltonian - projimp_off).norm();
    double weiss_off_norm2   =  (projweiss_field.getMatrix(0) - weiss_off).norm();
    ifroot std::cout<< "offdiagonal part of the impurity hamilonian ( delta_0) is reduced from " << projimp_off_norm << " (" <<weiss_off_norm<<") to " << projimp_off_norm2 << " (" <<weiss_off_norm2<<")\n";

    projNumMatrix = projSolverBasis.adjoint() * projNumMatrix * projSolverBasis ;

}
