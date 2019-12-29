#include "tight_common.h"
#include "mpi.h"


int min(int a,  int b) {
    if (a>=b) return b;
    else  return a;
}

void para_range(int n1, int n2, int nprocs, int myrank, int *mystart, int * myend) {
//    MPI_Barrier(MPI_COMM_WORLD);
    int work1, work2;
//assert( (n2-n1+1) >= nprocs);

    work1 = (n2-n1+1)/nprocs;
    work2 = (n2-n1+1)% nprocs;

    *mystart = myrank * work1 + n1 + min(myrank,work2);
    *myend = *mystart + work1 -1;
    if ( myrank < work2) *myend=*myend +1;
// Example
// rank = 0,1,...,work2-1 ==> (work1+1) job
// rank = work2 ... nprocs ==> work1 job
//
// if nprocs =3,  n1=1, n2=11 (for i=1,2,...11)  ==> work1=3, work2=2
// rank0 = mystart = 1, myend = 4
// rank1 = mystart = 5, myend = 8
// rank2 = mystart = 9, myend = 11
// Now, (*myend - *mystart) +1  = work job
//
// if nprocs = 10, n1=1, n2=5  ==> work1 = 0,  work2 = 5
// rank0 = mystart = 1, myend = 1
// rank1 = mystart = 2, myend = 2
// rank2 = mystart = 3, myend = 3
// rank3 = mystart = 4, myend = 4
// rank4 = mystart = 5, myend = 5
// rank5 = mystart = 5, myend = 4
// rank6 = mystart = 5, myend = 4
// ...
// rank9 = mystart = 5, myend = 4
//
// if nprocs = 10, n1=1, n2=1  ==> work1 = 0,  work2 = 1
// rank0 = mystart = 1, myend = 1
// rank1 = mystart = 2, myend = 1
// rank2 = mystart = 2, myend = 1
// rank3 = mystart = 2, myend = 1
// rank4 = mystart = 2, myend = 1
// rank5 = mystart = 2, myend = 1
// rank6 = mystart = 2, myend = 1
// ...
// rank9 = mystart = 2, myend = 1
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
void data_sync(cmplx *A, int startRow, int endRow,  int nprocs) {
    int itsRank,istart,iend, myrank, row;

    int * start = new int [nprocs];
    int * end = new int [nprocs];
    for(itsRank=0 ; itsRank<nprocs; itsRank++) {
        para_range(startRow,endRow,nprocs,itsRank, &istart, &iend);
        start[itsRank]=istart;
        end[itsRank]=iend;
    }
    for(myrank=0; myrank<nprocs; myrank++) {
        MPI_Bcast(&A[0], end[myrank]-start[myrank]+1, MPI_DOUBLE_COMPLEX, myrank, MPI_COMM_WORLD);
        //for(row=start[myrank]; row <=end[myrank] ; row++ ) {
        //    MPI_Bcast(&A[row], 1, MPI_DOUBLE_COMPLEX, myrank, MPI_COMM_WORLD);
        //}
    }

    delete [] start;
    delete [] end;
}
void data_sync(double *A, int startRow, int endRow,  int nprocs) {
    int itsRank,istart,iend, myrank, row;

    int * start = new int [nprocs];
    int * end =   new int [nprocs];
    for(itsRank=0 ; itsRank<nprocs; itsRank++) {
        para_range(startRow,endRow,nprocs,itsRank, &istart, &iend);
        start[itsRank]=istart;
        end[itsRank]=iend;
    }
    for(myrank=0; myrank<nprocs; myrank++) {
        MPI_Bcast(&A[start[myrank]], end[myrank]-start[myrank]+1, MPI_DOUBLE, myrank, MPI_COMM_WORLD);
    }

    delete [] start;
    delete [] end;
}
void data_sync(int  *A, int startRow, int endRow,  int nprocs) {
    int itsRank,istart,iend, myrank, row;

    int * start = new int [nprocs];
    int * end = new int [nprocs];
    for(itsRank=0 ; itsRank<nprocs; itsRank++) {
        para_range(startRow,endRow,nprocs,itsRank, &istart, &iend);
        start[itsRank]=istart;
        end[itsRank]=iend;
    }
    for(myrank=0; myrank<nprocs; myrank++) {
        MPI_Bcast(&(A[start[myrank]]), end[myrank]-start[myrank]+1, MPI_INT, myrank, MPI_COMM_WORLD);
    }

    delete [] start;
    delete [] end;
}
void data_sync_EigenMat(Eigen::MatrixXcd *A, int startRow, int endRow, int matrix_dim,  int nprocs) {
    int itsRank,istart,iend, myrank, row;

    int * start = new int [nprocs];
    int * end = new int [nprocs];
    for(itsRank=0 ; itsRank<nprocs; itsRank++) {
        para_range(startRow,endRow,nprocs,itsRank, &istart, &iend);
        start[itsRank]=istart;
        end[itsRank]=iend;
    }
    int lenRow = endRow-startRow+1;
    int mat_size= matrix_dim*matrix_dim;
    cmplx * AA = new cmplx [lenRow * mat_size ];
    for(row=0; row <lenRow ; row++ ) {
        for(int i=0; i<matrix_dim; i++) {
            for(int j=0; j<matrix_dim; j++) {
                AA[row*mat_size + i*matrix_dim + j] = A[row+startRow](i,j);
            }
        }
    }

    for(myrank=0; myrank<nprocs; myrank++) {
        MPI_Bcast(&(AA[(start[myrank]-startRow)*mat_size]), (end[myrank]-start[myrank]+1)*(matrix_dim*matrix_dim), MPI_DOUBLE_COMPLEX, myrank, MPI_COMM_WORLD);
        //for(row=start[myrank]; row <=end[myrank] ; row++ ) {
        //    MPI_Bcast(A[row].data(), A[row].size(), MPI_DOUBLE_COMPLEX, myrank, MPI_COMM_WORLD);
        //}
    }


    for(row=0; row <lenRow ; row++ ) {
        for(int i=0; i<matrix_dim; i++) {
            for(int j=0; j<matrix_dim; j++) {
                A[row+startRow](i,j) = AA[row*mat_size+i*matrix_dim + j] ;
            }
        }
    }

    delete [] AA;
    delete [] start;
    delete [] end;
}


void data_sync(cmplx **A, int startRow, int endRow, int lenColumn,  int nprocs) {
    int itsRank,istart,iend, myrank, row;

    int * start = new int [nprocs];
    int * end = new int [nprocs];
    for(itsRank=0 ; itsRank<nprocs; itsRank++) {
        para_range(startRow,endRow,nprocs,itsRank, &istart, &iend);
        start[itsRank]=istart;
        end[itsRank]=iend;
    }
    int lenRow = endRow-startRow+1;
    cmplx * AA = new cmplx [lenRow * lenColumn];
    for(int i=0; i<lenRow; i++) {
        for(int j=0; j<lenColumn; j++) {
            AA[i*lenColumn + j] = A[i+startRow][j];
        }
    }

    for(myrank=0; myrank<nprocs; myrank++) {
        MPI_Bcast(&(AA[ (start[myrank]-startRow)*lenColumn]), (end[myrank]-start[myrank]+1)*lenColumn, MPI_DOUBLE_COMPLEX, myrank, MPI_COMM_WORLD);
    }
    for(int i=0; i<lenRow; i++) {
        for(int j=0; j<lenColumn; j++) {
            A[i+startRow][j]=AA[i*lenColumn+j];
        }
    }
    delete [] AA;
    delete [] start;
    delete [] end;
}

void data_sync(double **A, int startRow, int endRow, int lenColumn,  int nprocs) {
    int itsRank,istart,iend, myrank, row;
    int * start = new int [nprocs];
    int * end = new int [nprocs];

    for(itsRank=0 ; itsRank<nprocs; itsRank++) {
        para_range(startRow,endRow,nprocs,itsRank, &istart, &iend);
        start[itsRank]=istart;
        end[itsRank]=iend;
    }
    int lenRow = endRow-startRow+1;
    double * AA = new double [lenRow * lenColumn];
    for(int i=0; i<lenRow; i++) {
        for(int j=0; j<lenColumn; j++) {
            AA[i*lenColumn+j] = A[i+startRow][j];
        }
    }
    for(myrank=0; myrank<nprocs; myrank++) {
        MPI_Bcast(&(AA[(start[myrank]-startRow) * lenColumn]), (end[myrank]-start[myrank]+1)*lenColumn, MPI_DOUBLE, myrank, MPI_COMM_WORLD);
    }
    for(int i=0; i<lenRow; i++) {
        for(int j=0; j<lenColumn; j++) {
            A[i+startRow][j]=AA[i*lenColumn+j];
        }
    }

    delete [] AA;
    delete [] start;
    delete [] end;
}

void data_sync(int **A, int startRow, int endRow, int lenColumn,  int nprocs) {
    int itsRank,istart,iend, myrank, row;

    int * start = new int [nprocs];
    int * end = new int [nprocs];
    for(itsRank=0 ; itsRank<nprocs; itsRank++) {
        para_range(startRow,endRow,nprocs,itsRank, &istart, &iend);
        start[itsRank]=istart;
        end[itsRank]=iend;
    }
    int lenRow = endRow-startRow+1;
    int * AA = new int [lenRow * lenColumn];
    for(int i=0; i<lenRow; i++) {
        for(int j=0; j<lenColumn; j++) {
            AA[i*lenColumn +j] = A[i+startRow][j];
        }
    }
    for(myrank=0; myrank<nprocs; myrank++) {
        MPI_Bcast(&(AA[(start[myrank]-startRow)*lenColumn]), (end[myrank]-start[myrank]+1)*lenColumn, MPI_INT, myrank, MPI_COMM_WORLD);
    }
    for(int i=0; i<lenRow; i++) {
        for(int j=0; j<lenColumn; j++) {
            A[i+startRow][j]=AA[i*lenColumn +j];
        }
    }
    delete [] AA;
    delete [] start;
    delete [] end;
}

void data_sync(cmplx ***A, int startRow, int endRow, int lenColumn, int lenhight,  int nprocs) {
    int itsRank,istart,iend, myrank, row;

    int * start = new int [nprocs];
    int * end = new int [nprocs];
    for(itsRank=0 ; itsRank<nprocs; itsRank++) {
        para_range(startRow,endRow,nprocs,itsRank, &istart, &iend);
        start[itsRank]=istart;
        end[itsRank]=iend;
    }
    int lenRow = endRow-startRow+1;
    int mat_size= lenColumn*lenhight;
    cmplx * AA = new cmplx [lenRow * mat_size ];
    for(int i=0; i<lenRow; i++) {
        for(int j=0; j<lenColumn; j++) {
            for(int k=0; k<lenhight; k++) {
                AA[i*mat_size+ j*lenhight+k] = A[i+startRow][j][k];
            }
        }
    }
    for(myrank=0; myrank<nprocs; myrank++) {
        MPI_Bcast(&(AA[(start[myrank]-startRow)*mat_size]), (end[myrank]-start[myrank]+1)*(lenColumn*lenhight), MPI_DOUBLE_COMPLEX, myrank, MPI_COMM_WORLD);
//        for(row=start[myrank]; row <=end[myrank] ; row++ ) {
//            MPI_Bcast(&A[row][0][0], lenColumn*lenhight, MPI_DOUBLE_COMPLEX, myrank, MPI_COMM_WORLD);
//        }
    }
    for(int i=0; i<lenRow; i++) {
        for(int j=0; j<lenColumn; j++) {
            for(int k=0; k<lenhight; k++) {
                A[i+startRow][j][k] = AA[i*mat_size +j*lenhight+k] ;
            }
        }
    }
    delete [] AA;
    delete [] start;
    delete [] end;
}
