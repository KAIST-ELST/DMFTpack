#include "mpi.h"
#include <fstream>
#include "ImgTimeFtn.h"
#include "tight_common.h"
#include "model.h"
#include "TB.h"
#include <Eigen/Eigenvalues>

int read_OverlapMat(Eigen::MatrixXi &S_overlap_Rindex, Eigen::VectorXcd  &S_overlap_RMatrix, const std::string &hamiltonian, std::vector<int> & accumulated_Num_SpinOrbital);

void SpreadFtn( int knum,
                std::vector<Eigen::MatrixXcd> S_overlap, std::vector<Eigen::MatrixXcd> &  transformMatrix_k,
                std::vector<int> & accumulated_Num_SpinOrbital    ) {



    Eigen::MatrixXi  x_overlap_Rindex,y_overlap_Rindex, z_overlap_Rindex   ;
    Eigen::VectorXcd x_overlap_RMatrix,y_overlap_RMatrix, z_overlap_RMatrix    ;
    int temp;
    temp =          read_OverlapMat(x_overlap_Rindex, x_overlap_RMatrix, std::string("OverlapMatrix_x.HWR"), accumulated_Num_SpinOrbital);
    temp =          read_OverlapMat(y_overlap_Rindex, y_overlap_RMatrix, std::string("OverlapMatrix_y.HWR"), accumulated_Num_SpinOrbital);
    temp =          read_OverlapMat(z_overlap_Rindex, z_overlap_RMatrix, std::string("OverlapMatrix_z.HWR"), accumulated_Num_SpinOrbital);





    Eigen::MatrixXcd x_overlap[knum];
    Eigen::MatrixXcd y_overlap[knum];
    Eigen::MatrixXcd z_overlap[knum];

    /*DFT_KS hamiltonian in k-space*/
    for(int k = 0;  k < knum; k++) {

        /*S(k),   H*S |\psi> = E |\psi>      */
        x_overlap[k].setZero(NumOrbit,NumOrbit);
        y_overlap[k].setZero(NumOrbit,NumOrbit);
        z_overlap[k].setZero(NumOrbit,NumOrbit);
        for(int indx=0; indx<x_overlap_RMatrix.size(); indx++) {
            x_overlap[k](x_overlap_Rindex(indx,3)*2, x_overlap_Rindex(indx,4)*2) += x_overlap_RMatrix(indx)
                    * exp ( -I*( (kmesh[k][0]*ax*x_overlap_Rindex(indx,0))+(kmesh[k][1]*ay*x_overlap_Rindex(indx,1))+(kmesh[k][2]*az*x_overlap_Rindex(indx,2))) )  ;
        }
        for(int indx=0; indx<y_overlap_RMatrix.size(); indx++) {
            y_overlap[k](y_overlap_Rindex(indx,3)*2, y_overlap_Rindex(indx,4)*2) += y_overlap_RMatrix(indx)
                    * exp ( -I*( (kmesh[k][0]*ax*y_overlap_Rindex(indx,0))+(kmesh[k][1]*ay*y_overlap_Rindex(indx,1))+(kmesh[k][2]*az*y_overlap_Rindex(indx,2))) )  ;
        }
        for(int indx=0; indx<z_overlap_RMatrix.size(); indx++) {
            z_overlap[k](z_overlap_Rindex(indx,3)*2, z_overlap_Rindex(indx,4)*2) += z_overlap_RMatrix(indx)
                    * exp ( -I*( (kmesh[k][0]*ax*z_overlap_Rindex(indx,0))+(kmesh[k][1]*ay*z_overlap_Rindex(indx,1))+(kmesh[k][2]*az*z_overlap_Rindex(indx,2))) )  ;
        }
        int Norbital_spatial = NumOrbit/Spin_DegreeOfFreedom;
        for (int i0=0; i0<Norbital_spatial; i0++) {
            for (int m0=0; m0<Norbital_spatial; m0++) {
                x_overlap[k]((2*i0+1), (2*m0+1)) = x_overlap[k](2*i0,2*m0);
                y_overlap[k]((2*i0+1), (2*m0+1)) = y_overlap[k](2*i0,2*m0);
                z_overlap[k]((2*i0+1), (2*m0+1)) = z_overlap[k](2*i0,2*m0);
            }
        }
    }//k
















    Eigen::MatrixXcd    x_overlap2_R0_mpiLocal;
    Eigen::MatrixXcd    y_overlap2_R0_mpiLocal;
    Eigen::MatrixXcd    z_overlap2_R0_mpiLocal;
    Eigen::MatrixXcd    x_overlap2_R0;
    Eigen::MatrixXcd    y_overlap2_R0;
    Eigen::MatrixXcd    z_overlap2_R0;

    Eigen::MatrixXcd    x_overlap_R0_mpiLocal;
    Eigen::MatrixXcd    y_overlap_R0_mpiLocal;
    Eigen::MatrixXcd    z_overlap_R0_mpiLocal;
    Eigen::MatrixXcd    x_overlap_R0;
    Eigen::MatrixXcd    y_overlap_R0;
    Eigen::MatrixXcd    z_overlap_R0;



    x_overlap2_R0_mpiLocal.setZero(NumOrbit, NumOrbit);
    y_overlap2_R0_mpiLocal.setZero(NumOrbit, NumOrbit);
    z_overlap2_R0_mpiLocal.setZero(NumOrbit, NumOrbit);
    x_overlap2_R0.setZero(NumOrbit,NumOrbit);
    y_overlap2_R0.setZero(NumOrbit, NumOrbit);
    z_overlap2_R0.setZero(NumOrbit, NumOrbit);

    x_overlap_R0_mpiLocal.setZero(NumOrbit, NumOrbit);
    y_overlap_R0_mpiLocal.setZero(NumOrbit, NumOrbit);
    z_overlap_R0_mpiLocal.setZero(NumOrbit, NumOrbit);
    x_overlap_R0.setZero(NumOrbit,NumOrbit);
    y_overlap_R0.setZero(NumOrbit, NumOrbit);
    z_overlap_R0.setZero(NumOrbit, NumOrbit);




    if (localOrbitalType == std::string("direct")  ) {
        for(int k = 0;  k < knum; k++) {
            x_overlap[k] = transformMatrix_k[k].adjoint() * x_overlap[k] * transformMatrix_k[k];
            y_overlap[k] = transformMatrix_k[k].adjoint() * y_overlap[k] * transformMatrix_k[k];
            z_overlap[k] = transformMatrix_k[k].adjoint() * z_overlap[k] * transformMatrix_k[k];

            x_overlap2_R0_mpiLocal += x_overlap[k].adjoint() * (S_overlap[k]).inverse()* x_overlap[k];
            y_overlap2_R0_mpiLocal += y_overlap[k].adjoint() * (S_overlap[k]).inverse()* y_overlap[k];
            z_overlap2_R0_mpiLocal += z_overlap[k].adjoint() * (S_overlap[k]).inverse()* z_overlap[k] ;
            x_overlap_R0_mpiLocal += x_overlap[k] ;
            y_overlap_R0_mpiLocal += y_overlap[k] ;
            z_overlap_R0_mpiLocal += z_overlap[k] ;
        }
    }
    else {
        for(int k = 0;  k < knum; k++) {
            x_overlap[k] = transformMatrix_k[k].adjoint() * x_overlap[k] * transformMatrix_k[k];
            y_overlap[k] = transformMatrix_k[k].adjoint() * y_overlap[k] * transformMatrix_k[k];
            z_overlap[k] = transformMatrix_k[k].adjoint() * z_overlap[k] * transformMatrix_k[k];
            x_overlap2_R0_mpiLocal += x_overlap[k].adjoint() * x_overlap[k];
            y_overlap2_R0_mpiLocal += y_overlap[k].adjoint() * y_overlap[k];
            z_overlap2_R0_mpiLocal += z_overlap[k].adjoint() * z_overlap[k];
            x_overlap_R0_mpiLocal += x_overlap[k] ;
            y_overlap_R0_mpiLocal += y_overlap[k] ;
            z_overlap_R0_mpiLocal += z_overlap[k] ;
        }
    }



    MPI_Allreduce(x_overlap2_R0_mpiLocal.data(), x_overlap2_R0.data(), x_overlap2_R0.size(), MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(y_overlap2_R0_mpiLocal.data(), y_overlap2_R0.data(), y_overlap2_R0.size(), MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(z_overlap2_R0_mpiLocal.data(), z_overlap2_R0.data(), z_overlap2_R0.size(), MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);

    MPI_Allreduce(x_overlap_R0_mpiLocal.data(), x_overlap_R0.data(), x_overlap_R0.size(), MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(y_overlap_R0_mpiLocal.data(), y_overlap_R0.data(), y_overlap_R0.size(), MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(z_overlap_R0_mpiLocal.data(), z_overlap_R0.data(), z_overlap_R0.size(), MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);

    x_overlap2_R0/=knum_mpiGlobal;
    y_overlap2_R0/=knum_mpiGlobal;
    z_overlap2_R0/=knum_mpiGlobal;

    x_overlap_R0/=knum_mpiGlobal;
    y_overlap_R0/=knum_mpiGlobal;
    z_overlap_R0/=knum_mpiGlobal;


    if(mpi_rank==0 )   std::cout << "<TB> Spread function:\n";



    for(int at=0; at<NumAtom; at++) {
        for(int h1F=accumulated_Num_SpinOrbital[at]; h1F<accumulated_Num_SpinOrbital[at+1]; h1F++) {
            if(isOrbitalHartrDFT[h1F]) {
                double Dx = real(x_overlap2_R0(h1F,h1F) - x_overlap_R0(h1F,h1F) * x_overlap_R0(h1F,h1F) );
                double Dy = real(y_overlap2_R0(h1F,h1F) - y_overlap_R0(h1F,h1F) * y_overlap_R0(h1F,h1F) );
                double Dz = real(z_overlap2_R0(h1F,h1F) - z_overlap_R0(h1F,h1F) * z_overlap_R0(h1F,h1F) );
                ifroot std::cout << at+1<< " " << h1F-accumulated_Num_SpinOrbital[at]+1 <<" " << std::sqrt(Dx+Dy+Dz)<<"\n" ;
            }
        }
    }


}




void SpreadFtn_PWF( int knum,
                    std::vector<Eigen::MatrixXcd> S_overlap,
                    std::vector<Eigen::MatrixXcd> &  transformMatrix_k,
                    std::vector<Eigen::MatrixXcd>  KS_eigenVectors_orthoBasis
                    ,  std::vector<int> & accumulated_Num_SpinOrbital) {



    ifroot std::cout << "Calculate Spread function..\n";
    Eigen::MatrixXi  x_overlap_Rindex,y_overlap_Rindex, z_overlap_Rindex   ;
    Eigen::VectorXcd x_overlap_RMatrix,y_overlap_RMatrix, z_overlap_RMatrix    ;
    int temp;
    temp =          read_OverlapMat(x_overlap_Rindex, x_overlap_RMatrix, std::string("OverlapMatrix_x.HWR"), accumulated_Num_SpinOrbital);
    temp =          read_OverlapMat(y_overlap_Rindex, y_overlap_RMatrix, std::string("OverlapMatrix_y.HWR"), accumulated_Num_SpinOrbital);
    temp =          read_OverlapMat(z_overlap_Rindex, z_overlap_RMatrix, std::string("OverlapMatrix_z.HWR"), accumulated_Num_SpinOrbital);
    int NumCorrOrbit = NumCorrAtom * N_peratom_HartrOrbit;





    Eigen::MatrixXcd x_overlap_k;
    Eigen::MatrixXcd y_overlap_k;
    Eigen::MatrixXcd z_overlap_k;






    Eigen::MatrixXcd    x_overlap2_R0_mpiLocal;
    Eigen::MatrixXcd    y_overlap2_R0_mpiLocal;
    Eigen::MatrixXcd    z_overlap2_R0_mpiLocal;
    Eigen::MatrixXcd    x_overlap2_R0;
    Eigen::MatrixXcd    y_overlap2_R0;
    Eigen::MatrixXcd    z_overlap2_R0;

    Eigen::MatrixXcd    x_overlap_R0_mpiLocal;
    Eigen::MatrixXcd    y_overlap_R0_mpiLocal;
    Eigen::MatrixXcd    z_overlap_R0_mpiLocal;
    Eigen::MatrixXcd    x_overlap_R0;
    Eigen::MatrixXcd    y_overlap_R0;
    Eigen::MatrixXcd    z_overlap_R0;



    x_overlap2_R0_mpiLocal.setZero(NumCorrOrbit, NumCorrOrbit);
    y_overlap2_R0_mpiLocal.setZero(NumCorrOrbit, NumCorrOrbit);
    z_overlap2_R0_mpiLocal.setZero(NumCorrOrbit, NumCorrOrbit);
    x_overlap2_R0.setZero(NumCorrOrbit, NumCorrOrbit);
    y_overlap2_R0.setZero(NumCorrOrbit, NumCorrOrbit);
    z_overlap2_R0.setZero(NumCorrOrbit, NumCorrOrbit);

    x_overlap_R0_mpiLocal.setZero(NumCorrOrbit, NumCorrOrbit);
    y_overlap_R0_mpiLocal.setZero(NumCorrOrbit, NumCorrOrbit);
    z_overlap_R0_mpiLocal.setZero(NumCorrOrbit, NumCorrOrbit);
    x_overlap_R0.setZero(NumCorrOrbit, NumCorrOrbit);
    y_overlap_R0.setZero(NumCorrOrbit, NumCorrOrbit);
    z_overlap_R0.setZero(NumCorrOrbit, NumCorrOrbit);


    Eigen::MatrixXcd  eye;
    eye.setIdentity(NumOrbit, NumOrbit);
    Eigen::MatrixXcd Bra_ProjOrbitals_Ket_KSenergy;

    for(int k = 0;  k < knum; k++) {
        /*S(k),   H*S |\psi> = E |\psi>      */
        x_overlap_k.setZero(NumOrbit,NumOrbit);
        y_overlap_k.setZero(NumOrbit,NumOrbit);
        z_overlap_k.setZero(NumOrbit,NumOrbit);
        for(int indx=0; indx<x_overlap_RMatrix.size(); indx++) {
            x_overlap_k(x_overlap_Rindex(indx,3)*2, x_overlap_Rindex(indx,4)*2) += x_overlap_RMatrix(indx)
                    * exp ( -I*( (kmesh[k][0]*ax*x_overlap_Rindex(indx,0))+(kmesh[k][1]*ay*x_overlap_Rindex(indx,1))+(kmesh[k][2]*az*x_overlap_Rindex(indx,2))) )  ;
        }
        for(int indx=0; indx<y_overlap_RMatrix.size(); indx++) {
            y_overlap_k(y_overlap_Rindex(indx,3)*2, y_overlap_Rindex(indx,4)*2) += y_overlap_RMatrix(indx)
                    * exp ( -I*( (kmesh[k][0]*ax*y_overlap_Rindex(indx,0))+(kmesh[k][1]*ay*y_overlap_Rindex(indx,1))+(kmesh[k][2]*az*y_overlap_Rindex(indx,2))) )  ;
        }
        for(int indx=0; indx<z_overlap_RMatrix.size(); indx++) {
            z_overlap_k(z_overlap_Rindex(indx,3)*2, z_overlap_Rindex(indx,4)*2) += z_overlap_RMatrix(indx)
                    * exp ( -I*( (kmesh[k][0]*ax*z_overlap_Rindex(indx,0))+(kmesh[k][1]*ay*z_overlap_Rindex(indx,1))+(kmesh[k][2]*az*z_overlap_Rindex(indx,2))) )  ;
        }
        int Norbital_spatial = NumOrbit/Spin_DegreeOfFreedom;
        for (int i0=0; i0<Norbital_spatial; i0++) {
            for (int m0=0; m0<Norbital_spatial; m0++) {
                x_overlap_k((2*i0+1), (2*m0+1)) = x_overlap_k(2*i0,2*m0);
                y_overlap_k((2*i0+1), (2*m0+1)) = y_overlap_k(2*i0,2*m0);
                z_overlap_k((2*i0+1), (2*m0+1)) = z_overlap_k(2*i0,2*m0);
            }
        }

        Bra_ProjOrbitals_Ket_KSenergy = eye;
        for (int i=0; i<NBAND[k]; i++) {
            for (int j=0; j<NBAND[k]; j++) {
                Bra_ProjOrbitals_Ket_KSenergy(FromValToKS[k][i], FromValToKS[k][j]) = DF_CorrBase[k](i,j);

            }
        }
        Eigen::MatrixXcd KS_eigenVectors_nonOrtho_DFT_dual_Basis
            =  S_overlap[k].inverse() *  (   (transformMatrix_k[k].adjoint().inverse())* KS_eigenVectors_orthoBasis[k]); // <\phi^\alpha |  kn>

//        Eigen::MatrixXcd KS_eigenVectors_nonOrtho_DFT_dual_Basis
//            =  transformMatrix_k[k] * KS_eigenVectors_orthoBasis[k];


        // <nk|r|nk>
        x_overlap_k = (KS_eigenVectors_nonOrtho_DFT_dual_Basis.adjoint()*x_overlap_k*KS_eigenVectors_nonOrtho_DFT_dual_Basis).eval();
        y_overlap_k = (KS_eigenVectors_nonOrtho_DFT_dual_Basis.adjoint()*y_overlap_k*KS_eigenVectors_nonOrtho_DFT_dual_Basis).eval();
        z_overlap_k = (KS_eigenVectors_nonOrtho_DFT_dual_Basis.adjoint()*z_overlap_k*KS_eigenVectors_nonOrtho_DFT_dual_Basis).eval();

        //  <PWF|r|PWF>
        x_overlap_k = (Bra_ProjOrbitals_Ket_KSenergy*x_overlap_k*Bra_ProjOrbitals_Ket_KSenergy.adjoint()).eval();
        y_overlap_k = (Bra_ProjOrbitals_Ket_KSenergy*y_overlap_k*Bra_ProjOrbitals_Ket_KSenergy.adjoint()).eval();
        z_overlap_k = (Bra_ProjOrbitals_Ket_KSenergy*z_overlap_k*Bra_ProjOrbitals_Ket_KSenergy.adjoint()).eval();

        x_overlap2_R0_mpiLocal += (x_overlap_k.adjoint() * x_overlap_k ).block(FromValToKS[k][0], FromValToKS[k][0], NumCorrOrbit, NumCorrOrbit);
        y_overlap2_R0_mpiLocal += (y_overlap_k.adjoint() * y_overlap_k ).block(FromValToKS[k][0], FromValToKS[k][0], NumCorrOrbit, NumCorrOrbit);
        z_overlap2_R0_mpiLocal += (z_overlap_k.adjoint() * z_overlap_k ).block(FromValToKS[k][0], FromValToKS[k][0], NumCorrOrbit, NumCorrOrbit);
        x_overlap_R0_mpiLocal  += (x_overlap_k                         ).block(FromValToKS[k][0], FromValToKS[k][0], NumCorrOrbit, NumCorrOrbit);
        y_overlap_R0_mpiLocal  += (y_overlap_k                         ).block(FromValToKS[k][0], FromValToKS[k][0], NumCorrOrbit, NumCorrOrbit);
        z_overlap_R0_mpiLocal  += (z_overlap_k                         ).block(FromValToKS[k][0], FromValToKS[k][0], NumCorrOrbit, NumCorrOrbit);
    }//k



    MPI_Allreduce(x_overlap2_R0_mpiLocal.data(), x_overlap2_R0.data(), x_overlap2_R0.size(), MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(y_overlap2_R0_mpiLocal.data(), y_overlap2_R0.data(), y_overlap2_R0.size(), MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(z_overlap2_R0_mpiLocal.data(), z_overlap2_R0.data(), z_overlap2_R0.size(), MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);

    MPI_Allreduce(x_overlap_R0_mpiLocal.data(), x_overlap_R0.data(), x_overlap_R0.size(), MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(y_overlap_R0_mpiLocal.data(), y_overlap_R0.data(), y_overlap_R0.size(), MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(z_overlap_R0_mpiLocal.data(), z_overlap_R0.data(), z_overlap_R0.size(), MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);

    x_overlap2_R0/=knum_mpiGlobal;
    y_overlap2_R0/=knum_mpiGlobal;
    z_overlap2_R0/=knum_mpiGlobal;

    x_overlap_R0/=knum_mpiGlobal;
    y_overlap_R0/=knum_mpiGlobal;
    z_overlap_R0/=knum_mpiGlobal;

    if(mpi_rank==0 )   std::cout << "<TB> center (Angs), PWF:\n";
    for(int h1F=0; h1F<NumCorrOrbit; h1F++) {
        double Dx = real( x_overlap_R0(h1F,h1F) );
        double Dy = real( y_overlap_R0(h1F,h1F) );
        double Dz = real( z_overlap_R0(h1F,h1F) );
        ifroot std::cout << h1F << ": " <<(Dx) <<" " <<   (Dy) <<" " << (Dz) <<"\n";
    }

    if(mpi_rank==0 )   std::cout << "<TB> Spread function, PWF:\n";
    for(int h1F=0; h1F<NumCorrOrbit; h1F++) {
        double Dx = real(x_overlap2_R0(h1F,h1F) - x_overlap_R0(h1F,h1F) * x_overlap_R0(h1F,h1F) );
        double Dy = real(y_overlap2_R0(h1F,h1F) - y_overlap_R0(h1F,h1F) * y_overlap_R0(h1F,h1F) );
        double Dz = real(z_overlap2_R0(h1F,h1F) - z_overlap_R0(h1F,h1F) * z_overlap_R0(h1F,h1F) );
        ifroot std::cout << h1F << ": " <<std::sqrt(Dx) <<" " <<   std::sqrt(Dy) <<" " << std::sqrt(Dz) <<" : "   << std::sqrt(Dx+Dy+Dz)<<"\n" ;
    }
}
