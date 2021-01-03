#include <fstream>
#include <Eigen/Core>
#include "tight_common.h"
//#include "model.h"
//#include "TB.h"
void read_HR(Eigen::MatrixXi &H_Rindex, Eigen::VectorXcd  &H_RMatrix, const std::string &hamiltonian, std::vector<int> & accumulated_Num_SpinOrbital, int NumOrbit) {

    int number_nonzero_element=0;
    FILE *dataIN;
    ifroot{
        dataIN = fopen(hamiltonian.c_str(), "r");
        if (dataIN == NULL) {
            printf("Cannot open Hopping data file...\n");
            exit(1);
        }
        accumulated_Num_SpinOrbital[0]=0;
        for(int i=1; i<NumAtom+1; i++) {
            int temp1;
            fscanf(dataIN,"Norbital of Atom = %d\n", &temp1);
            accumulated_Num_SpinOrbital[i]=temp1+accumulated_Num_SpinOrbital[i-1];
        }
        std::string line;
        std::ifstream hoppingfile(hamiltonian.c_str());
        while (std::getline(hoppingfile, line)) {
            ++number_nonzero_element;
        }
        number_nonzero_element -= NumAtom;
    }

    MPI_Bcast(accumulated_Num_SpinOrbital.data(), accumulated_Num_SpinOrbital.size(),  MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast( &number_nonzero_element, 1, MPI_INT, 0, MPI_COMM_WORLD   );
    H_Rindex.setZero(number_nonzero_element,5);
    H_RMatrix.setZero(number_nonzero_element);

    ifroot{

        int index=0;
        while(!feof(dataIN)) {
            int n,l,m, Atom1, Atom2,  i,j;
            double HopRe, HopIm;
            fscanf(dataIN,"%d %d %d   %d %d    %d %d     %lf  %lf\n",&n,&l,&m,  &Atom1, &Atom2,    &i,&j,   &HopRe,  &HopIm);
            int i0=accumulated_Num_SpinOrbital[Atom1-1]+i-1;
            int m0=accumulated_Num_SpinOrbital[Atom2-1]+j-1;
            H_Rindex(index,0) = n;
            H_Rindex(index,1) = l;
            H_Rindex(index,2) = m;

            H_Rindex(index,3) = i0;
            H_Rindex(index,4) = m0;

            H_RMatrix(index) = HopRe + I*HopIm;



            if(i0 > NumOrbit or m0>NumOrbit) {
                printf("out of range, Hk.HWR: %d %d %d   %d %d    %d %d     %lf  %lf\n",n,l,m,  Atom1, Atom2,    i,j,   HopRe,  HopIm);
                fflush(stdout);
            }

            index++;
        }
        fclose(dataIN);
    }

//    MPI_Bcast( &number_nonzero_element, 1, MPI_INT, 0 , MPI_COMM_WORLD   );
//    H_Rindex.setZero(number_nonzero_element,5);
//    H_RMatrix.setZero(number_nonzero_element);

    MPI_Bcast(H_Rindex.data(), H_Rindex.size(), MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(H_RMatrix.data(), H_RMatrix.size(), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

    ifroot std::cout << "Hk.HWR file was normally read.\n";
    sleep(10);
}

int read_OverlapMat(Eigen::MatrixXi &S_overlap_Rindex, Eigen::VectorXcd  &S_overlap_RMatrix, const std::string &readfile, std::vector<int> & accumulated_Num_SpinOrbital) {
    double S_max_ofdiagonal=0.;
    int overlap_exist;
    int number_nonzero_element=0;

    ifroot{
        if (  access(readfile.c_str(), F_OK)!= 0) {
            ifroot  std::cout << "Cannot found Overlap matrix\n"  ;
            overlap_exist =0;
        }
        else {
            overlap_exist=1;
            std::string line;
            std::ifstream overlapfile(readfile.c_str());
            while (std::getline(overlapfile, line)) {
                ++number_nonzero_element;
            }
            S_overlap_Rindex.setZero(number_nonzero_element,5);
            S_overlap_RMatrix.setZero(number_nonzero_element);
            int    index=0;


            FILE *Overlap= fopen(readfile.c_str(),"r");

            int n_max=0,l_max=0,m_max=0, Atom1_max=0, Atom2_max=0, i_max=0,j_max=0;
            while(!feof(Overlap)) {
                int n,l,m, Atom1, Atom2, i,j;
                double HopRe, HopIm;
                fscanf(Overlap,"%d %d %d   %d %d    %d %d     %lf  %lf\n",&n,&l,&m,  &Atom1, &Atom2,    &i,&j,   &HopRe, &HopIm);
                int i0=(accumulated_Num_SpinOrbital[Atom1-1]/Spin_DegreeOfFreedom) +(i-1);
                int m0=(accumulated_Num_SpinOrbital[Atom2-1]/Spin_DegreeOfFreedom) +(j-1);
                S_overlap_Rindex(index,0)=n;
                S_overlap_Rindex(index,1)=l;
                S_overlap_Rindex(index,2)=m;
                S_overlap_Rindex(index,3)=i0;
                S_overlap_Rindex(index,4)=m0;
                S_overlap_RMatrix(index) = HopRe + I * HopIm;
                if((Atom1!=Atom2 or n!=0 or l!=0 or m!=0) and std::abs(HopRe+I*HopIm) > S_max_ofdiagonal  ) {
                    n_max = n;
                    l_max = l;
                    m_max = m;
                    Atom1_max = Atom1;
                    Atom2_max = Atom2;
                    i_max = i;
                    j_max = j;
                    S_max_ofdiagonal=  std::abs(HopRe+I*HopIm);
                }
                index++;
            }

            if (mpi_rank ==0)  std::cout << "We have Overlap matrix"  <<"\n";
            if (mpi_rank ==0)  std::cout <<"OverlapMatrix_max_offdiagonal: "
                                             <<n_max<< ","
                                             <<l_max<<","
                                             <<m_max<< "; "
                                             << Atom1_max <<","
                                             << Atom2_max <<"; "
                                             <<i_max <<","
                                             <<j_max <<"= "
                                             << S_max_ofdiagonal<<"\n";
        }
    }//ifroot

    MPI_Bcast(&overlap_exist, 1, MPI_INT,0, MPI_COMM_WORLD);
    if(overlap_exist) {
        MPI_Bcast(&number_nonzero_element, 1, MPI_INT,0, MPI_COMM_WORLD);
        if(mpi_rank!=0) {
            S_overlap_Rindex.setZero(number_nonzero_element,5);
            S_overlap_RMatrix.setZero(number_nonzero_element);
        }
        MPI_Bcast(S_overlap_Rindex.data(),  S_overlap_Rindex.size(), MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(S_overlap_RMatrix.data(), S_overlap_RMatrix.size(), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    }


    return overlap_exist;
}





int  read_DensityMat(Eigen::MatrixXi &NumMat_Rindex, Eigen::VectorXcd  &NumMat_RMatrix, const std::string &hamiltonian, std::vector<int> & accumulated_Num_SpinOrbital) {
    FILE *dataIN = fopen(hamiltonian.c_str(), "r");

    int overlap_exist;
    if (  access(hamiltonian.c_str(), F_OK)!= 0) {
        overlap_exist =0;
        ifroot printf("Cannot open Occupation Mat file...\n");
    }
    else {
        overlap_exist=1;

        int number_nonzero_element=0;
        std::string line;
        std::ifstream hoppingfile(hamiltonian.c_str());
        while (std::getline(hoppingfile, line)) {
            ++number_nonzero_element;
        }

        NumMat_Rindex.setZero(number_nonzero_element,5);
        NumMat_RMatrix.setZero(number_nonzero_element);



        int index=0;
        while(!feof(dataIN)) {
            int n,l,m, Atom1, Atom2,  i,j;
            double HopRe, HopIm;
            fscanf(dataIN,"%d %d %d   %d %d    %d %d     %lf  %lf\n",&n,&l,&m,  &Atom1, &Atom2,    &i,&j,   &HopRe,  &HopIm);
            int i0=accumulated_Num_SpinOrbital[Atom1-1]+i-1;
            int m0=accumulated_Num_SpinOrbital[Atom2-1]+j-1;
            NumMat_Rindex(index,0) = n;
            NumMat_Rindex(index,1) = l;
            NumMat_Rindex(index,2) = m;
            NumMat_Rindex(index,3) = i0;
            NumMat_Rindex(index,4) = m0;

            NumMat_RMatrix(index) = HopRe + I*HopIm;

            index++;
        }
        fclose(dataIN);
        if (mpi_rank ==0)  std::cout << "We have Occupation matrix with " <<number_nonzero_element <<"lines..."  <<"\n";
    }
    return overlap_exist;
}
