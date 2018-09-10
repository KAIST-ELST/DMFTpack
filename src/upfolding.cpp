#include "mpi.h"
#include "ImgTimeFtn.h"
#include "tight_common.h"
#include "model.h"
#include "TB.h"
#include <Eigen/Eigenvalues>


//from model basis to DFT basis.
void upfolding_k(Eigen::MatrixXcd & densityMatDFT_k, Eigen::MatrixXcd & KS_eigenVectors_nonOrtho_DFTBasis, int k,double mu)
{
    /*densityMatDFT_k from Model to KS, |nk>, basis*/
    densityMatDFT_k =  DF_CorrBase[k].adjoint() * densityMatDFT_k * DF_CorrBase[k];

    Eigen::VectorXcd densityMatDFT_k_diagonal_in_KSSpace(NumOrbit);
    for(int l0=0; l0<NumOrbit; l0++) {
        if     (  (beta*(KS_eigenEnergy[k][l0] - mu)) >  100.0 ) densityMatDFT_k_diagonal_in_KSSpace[l0] = 0.0;
        else if(  (beta*(KS_eigenEnergy[k][l0] - mu)) < -100.0 ) densityMatDFT_k_diagonal_in_KSSpace[l0] = 1.0;
        else     densityMatDFT_k_diagonal_in_KSSpace[l0] = 1./( 1.+std::exp(  beta*(KS_eigenEnergy[k][l0]-mu) ));
    }


    /* set Matrix in Full KS base*/
    Eigen::MatrixXcd Matrix_KSBase;
    Matrix_KSBase.setZero(NumOrbit,NumOrbit);
    for(int i=0 ; i< NumOrbit; i++)  Matrix_KSBase(i,i) = densityMatDFT_k_diagonal_in_KSSpace[i];

    for(int i=0 ; i< NBAND[k]; i++) {
        for(int j=0 ; j< NBAND[k]; j++) {
            Matrix_KSBase(FromValToKS[k][i],FromValToKS[k][j]) = densityMatDFT_k(i,j);
        }
    }
    /*We het Matrix in DFT base*/
    densityMatDFT_k = KS_eigenVectors_nonOrtho_DFTBasis * Matrix_KSBase * KS_eigenVectors_nonOrtho_DFTBasis.adjoint();
}


void upfolding_density(std::vector<Eigen::MatrixXcd> &densityMatDFT, std::vector<Eigen::MatrixXcd> & KS_eigenVectors_orthoBasis ,Eigen::MatrixXi H_Rindex, double mu ,
                       std::vector<Eigen::MatrixXcd> & S_overlap,
                       std::vector<Eigen::MatrixXcd> & transformMatrix_k) {

    Eigen::MatrixXcd dual_denMat_dual[knum];
    Eigen::MatrixXcd dual_denMat_direct[knum];
    for(int k=0; k<knum; k++) {
        if(downfolding == 1 ) {
            Eigen::MatrixXcd KS_eigenVectors_nonOrtho_DFT_direct_Basis = ((transformMatrix_k[k].adjoint().inverse())* KS_eigenVectors_orthoBasis[k]);
            upfolding_k(densityMatDFT[k], KS_eigenVectors_nonOrtho_DFT_direct_Basis, k, mu);

            dual_denMat_dual[k] =     S_overlap[k].inverse() *densityMatDFT[k] * S_overlap[k].inverse();
            dual_denMat_direct[k] =  (S_overlap[k].inverse() *densityMatDFT[k]  ).eval();
        }
    }//k

    Eigen::VectorXcd dual_DMR_dual_loc, dual_DMR_dual;
    dual_DMR_dual_loc.setZero(H_Rindex.rows());
    dual_DMR_dual.setZero(H_Rindex.rows());

    Eigen::VectorXcd dual_DM_direct_DFTR_loc, dual_DM_direct_DFTR;
    dual_DM_direct_DFTR_loc.setZero(H_Rindex.rows());
    dual_DM_direct_DFTR.setZero(H_Rindex.rows());

    Eigen::VectorXcd direct_DMR_direct_loc, direct_DMR_direct;
    direct_DMR_direct_loc.setZero(H_Rindex.rows());
    direct_DMR_direct.setZero(H_Rindex.rows());


    for(int index=0; index<H_Rindex.rows(); index++) {
        int n=  H_Rindex(index,0) ;
        int l=  H_Rindex(index,1) ;
        int m=  H_Rindex(index,2) ;
        int i0= H_Rindex(index,3) ;
        int m0= H_Rindex(index,4) ;
        for(int k=0; k<knum; k++) {
            dual_DMR_dual_loc(index)      +=  dual_denMat_dual[k](i0,m0) * exp ( I*( (kmesh[k][0]*ax*n)+(kmesh[k][1]*ay*l)+(kmesh[k][2]*az*m)) )  ;
            dual_DM_direct_DFTR_loc(index) += dual_denMat_direct[k](i0,m0) * exp ( I*( (kmesh[k][0]*ax*n)+(kmesh[k][1]*ay*l)+(kmesh[k][2]*az*m)) )  ;
            direct_DMR_direct_loc(index)   += densityMatDFT[k](i0,m0) * exp ( I*( (kmesh[k][0]*ax*n)+(kmesh[k][1]*ay*l)+(kmesh[k][2]*az*m)) )  ;
        }

    }
    MPI_Allreduce(dual_DMR_dual_loc.data(), dual_DMR_dual.data(), dual_DMR_dual.size(), MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(dual_DM_direct_DFTR_loc.data(), dual_DM_direct_DFTR.data(), dual_DM_direct_DFTR.size(), MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(direct_DMR_direct_loc.data(), direct_DMR_direct.data(), direct_DMR_direct.size(), MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);

    dual_DMR_dual /= knum_mpiGlobal;
    dual_DM_direct_DFTR/= knum_mpiGlobal;
    direct_DMR_direct  /= knum_mpiGlobal;


    ifroot{
        FILE *fp_jh3;
        fp_jh3 = fopen("DM_DFT_DMFT.HWR","w");
        for(int index=0; index<H_Rindex.rows(); index++) {
            int n=  H_Rindex(index,0) ;
            int l=  H_Rindex(index,1) ;
            int m=  H_Rindex(index,2) ;
            int i0= H_Rindex(index,3) ;
            int m0= H_Rindex(index,4) ;
            fprintf(fp_jh3, " %+d %+d %+d %d %d %d %d      %e    %e  \n",
            n,l,m   ,
            FromOrbitalToAtom[i0]+1,             FromOrbitalToAtom[m0]+1,
            FromOrbitalToLocalOrbital_DFT[i0]+1, FromOrbitalToLocalOrbital_DFT[m0]+1,
            real(dual_DMR_dual(index)), imag(dual_DMR_dual(index)));  //Collinear.
        }
        std::cout << "FILEOUT:DM_DFT_DMFT.HWR\n" ;
//        std::cout << "DFT+DMFT: Collinear calculation in current version!!!\n";
        fclose(fp_jh3);

        fp_jh3 = fopen("dual_DM_direct.HWR","w");
        for(int index=0; index<H_Rindex.rows(); index++) {
            int n=  H_Rindex(index,0) ;
            int l=  H_Rindex(index,1) ;
            int m=  H_Rindex(index,2) ;
            int i0= H_Rindex(index,3) ;
            int m0= H_Rindex(index,4) ;
            fprintf(fp_jh3, " %+d %+d %+d %d %d %d %d      %e    %e  \n",
            n,l,m   ,
            FromOrbitalToAtom[i0]+1,             FromOrbitalToAtom[m0]+1,
            FromOrbitalToLocalOrbital_DFT[i0]+1, FromOrbitalToLocalOrbital_DFT[m0]+1,
            real(dual_DM_direct_DFTR(index)), imag(dual_DM_direct_DFTR(index)));  //Collinear.
        }
        fclose(fp_jh3);

        fp_jh3 = fopen("direct_DM_direct.HWR","w");
        for(int index=0; index<H_Rindex.rows(); index++) {
            int n=  H_Rindex(index,0) ;
            int l=  H_Rindex(index,1) ;
            int m=  H_Rindex(index,2) ;
            int i0= H_Rindex(index,3) ;
            int m0= H_Rindex(index,4) ;
            fprintf(fp_jh3, " %+d %+d %+d %d %d %d %d      %e    %e  \n",
            n,l,m   ,
            FromOrbitalToAtom[i0]+1,             FromOrbitalToAtom[m0]+1,
            FromOrbitalToLocalOrbital_DFT[i0]+1, FromOrbitalToLocalOrbital_DFT[m0]+1,
            real(direct_DMR_direct(index)), imag(direct_DMR_direct(index)));  //Collinear.
        }
        std::cout << "FILEOUT:DM_DFTDMFT_dual_direct.HWR\n" ;
        fclose(fp_jh3);
    }
}
