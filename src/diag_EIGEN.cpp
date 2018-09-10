#include <Eigen/Eigenvalues>
#include "tight_common.h"
void cmplxEigenproblem(cmplx **a, int dim, double *eigenVal, cmplx **v_comp)
{
    //Compute all eigenvalues and eigenvectors of a general cmplx  matrix A[1..n][1..n].
    //On output, elements of a above the diagonal are destroyed. d[1...n] returns the eigenvalues of a.
    //V[1..n][1..n] is a matrix whose coulumns contain, on output, the normalized eigenvectors of a. i.e.  (V[0][0], V[0][1], V[0][2],..V[0][n]) = 0th eigenvector
    //http://eigen.tuxfamily.org/dox/GettingStarted.html
    Eigen::MatrixXcd aMat(dim, dim);  //where X=size(dynamic), cd=cmplx double
    int i,j,k, check=1;

    for(i=0; i<dim; i++) {
        for(j=0; j<dim; j++) {
            aMat.row(i)[j] = a[i][j];
        }
    }
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> ces(dim);
    ces.compute(aMat);
    // ces.computeDirect(aMat);  faster  but it might also be less accurate.
    for(i =0; i<dim; i++) {
        eigenVal[i] = ces.eigenvalues()[i] ;
        if(check==1 and  eigenVal[i] != 0 ) check=0;  //zero-matrix?
        for(j =0; j<dim; j++) {
            v_comp[i][j] = ces.eigenvectors().col(i)[j];
        }
    }
    if(check==1) {
        for(i =0; i<dim; i++) {
            for(j =0; j<dim; j++) {
                v_comp[i][j] =0;
            }
            v_comp[i][i] =1;
        }
    }
}
void cmplxEigenproblem(cmplx *a, int dim, double *eigenVal, cmplx *v_comp)
{
    //Compute all eigenvalues and eigenvectors of a general cmplx  matrix A[1..n][1..n].
    //On output, elements of a above the diagonal are destroyed. d[1...n] returns the eigenvalues of a.
    //V[1..n][1..n] is a matrix whose coulumns contain, on output, the normalized eigenvectors of a. i.e.  (V[0][0], V[0][1], V[0][2],..V[0][n]) = 0th eigenvector
    //http://eigen.tuxfamily.org/dox/GettingStarted.html
    Eigen::MatrixXcd aMat(dim, dim);  //where X=size(dynamic), cd=cmplx double
    int i,j,k, check=0;

    for(i=0; i<dim; i++) {
        for(j=0; j<dim; j++) {
            aMat.row(i)[j] = a[i*dim+j];
        }
    }
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> ces(dim);
    ces.compute(aMat);
    for(i =0; i<dim; i++) {
        eigenVal[i] = ces.eigenvalues()[i] ;
        if(check==1 and  eigenVal[i] != 0 ) check=0;  //zero-matrix?
        for(j =0; j<dim; j++) {
            v_comp[i*dim+j] = ces.eigenvectors().col(i)[j];
        }
    }
}



void cmplxEigenproblem(cmplx **a, int dim, double *eigenVal)
{
    if(dim>0) {
        Eigen::MatrixXcd aMat(dim, dim);  //where X=size(dynamic), cd=cmplx double
        int i,j,k, check=1;

        for(i=0; i<dim; i++) {
            for(j=0; j<dim; j++) {
                aMat.row(i)[j] = a[i][j];
            }
        }
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> ces(dim);
        ces.compute(aMat,0);
        for(i =0; i<dim; i++) {
            eigenVal[i] = ces.eigenvalues()[i] ;
        }
    }
}

void cmplxEigenproblem(cmplx *a, int dim, double *eigenVal)
{
    Eigen::MatrixXcd aMat(dim, dim);  //where X=size(dynamic), cd=cmplx double
    int i,j,k, check=1;

    for(i=0; i<dim; i++) {
        for(j=0; j<dim; j++) {
            aMat.row(i)[j] = a[i*dim+j];
        }
    }


    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> ces(dim);

    ces.compute(aMat,0);

    for(i =0; i<dim; i++) {
        eigenVal[i] = ces.eigenvalues()[i] ;
    }
}




int cmpfunc (const void * a, const void * b)
{
    if (( *((double*)a) - *((double*)b) )>0) return 1;
    else return -1;
}
double GcmplxEigenproblem(cmplx *a, cmplx *b, int dim, double *eigenVal)
{
    Eigen::MatrixXcd aMat(dim, dim);  //where X=size(dynamic), cd=cmplx double
    Eigen::MatrixXcd bMat(dim, dim);  //where X=size(dynamic), cd=cmplx double
    int i,j,k, check=1;
    double maxImg=0;

    for(i=0; i<dim; i++) {
        for(j=0; j<dim; j++) {
            aMat.row(i)[j] = a[i*dim+j];
            bMat.row(i)[j] = b[i*dim+j];
        }
    }
    //Eigen::ComplexEigenSolver<Eigen::MatrixXcd> ces(dim);
    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXcd> ces(dim);

    ces.compute(aMat,bMat);
    for(i =0; i<dim; i++) {
        eigenVal[i] = (ces.eigenvalues()[i]) ;
//        if(maxImg<imag(ces.eigenvalues()[i])) maxImg = imag(ces.eigenvalues()[i]);
    }
    qsort(eigenVal, dim, sizeof(double), cmpfunc);
    return maxImg;
}

double GcmplxEigenproblem(cmplx **a, cmplx *b, int dim, double *eigenVal, cmplx **eigenVec)
{
//Sol  eigen problem,   av= \lambda b v
    Eigen::MatrixXcd aMat(dim, dim);  //where X=size(dynamic), cd=cmplx double
    Eigen::MatrixXcd bMat(dim, dim);  //where X=size(dynamic), cd=cmplx double
    double maxImg=0;

    for(int i=0; i<dim; i++) {
        for(int j=0; j<dim; j++) {
            aMat.row(i)[j] = a[i][j];
            bMat.row(i)[j] = b[i*dim+j];
        }
    }
    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXcd> ces(dim);
    ces.compute(aMat,bMat);
    for(int i =0; i<dim; i++) {
        eigenVal[i] = (ces.eigenvalues()[i]) ;
        for(int j =0; j<dim; j++) {
            eigenVec[i][j] = ces.eigenvectors().col(i)[j];   //(V[0][0], V[0][1], V[0][2],..V[0][n]) = 0th eigenvector
        }
    }
    return maxImg;
}
