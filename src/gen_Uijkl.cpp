#include <Eigen/Dense>
#include <fstream>
#include "tight_common.h"
// H_inter = 1/2  \sum_{ijkl}  U_{ijkl} c^+_i c^+_j c_l c_k ,
//  and U_{ijkl} = \int dr1 dr2   \i*(r1) \j*(r2) V(r1,r2) \k(r1) \l(r2)

void write_Uijkl(int size, std::vector<cmplx >  & Utensor,std::vector<Eigen::VectorXi>  & Uindex, int indexC) ;

void gen_Uijkl(int n_spinorb, double U, double Uprime, double JH, std::vector<cmplx >  & Utensor, std::vector<Eigen::VectorXi>  & Uindex ) {
//    for(int estimate_dim=0; estimate_dim<2; estimate_dim++) {
    int index=0, indexC=0;
    for(int i=0; i<n_spinorb; i++) {
        for(int j=0; j<n_spinorb; j++) {
            for(int k=0; k<n_spinorb; k++) {
                for(int l=0; l<n_spinorb; l++) {
                    int iorb =  i/2;
                    int ispin = i%2;
                    int jorb =  j/2;
                    int jspin = j%2;
                    int korb =  k/2;
                    int kspin = k%2;
                    int lorb =  l/2;
                    int lspin = l%2;
                    double coeff=0;
//                    if(i!=j and k!=l) { //c^+_i c^+_j  = c_k c_l = 0.
                    //density-density terms
                    if(i==k and j==l) {
//                            if (iorb==jorb and ispin!=jspin) coeff = U;             //intra-orbital
                        if (iorb==jorb) coeff = U;             //intra-orbital
                        if (iorb!=jorb and ispin!=jspin) coeff = Uprime;        //inter-orbital
                        if (iorb!=jorb and ispin==jspin) coeff = Uprime - JH;   //inter-orbital
                    }
//                        //spin-flip
//                        else if (iorb  == lorb  and jorb  == korb  and iorb  !=  jorb and
//                                 ispin == kspin and jspin == lspin and ispin !=  jspin    )   coeff = JH;
//                        //pair-hopping
//                        else if (iorb  == jorb  and korb  == lorb  and iorb  !=  lorb and
//                                 ispin == kspin and jspin == lspin and ispin !=  jspin    )   coeff = JH;
                    //spin-flip
                    else if (korb  == iorb  and lorb  == jorb  and korb  !=  lorb and
                             kspin != ispin and lspin != jspin and kspin !=  lspin    )   coeff = JH;
                    //pair-hopping
                    else if (korb  == lorb  and korb  != iorb  and iorb  ==  jorb and
                             kspin != lspin and kspin == ispin and lspin ==  jspin    )   coeff = JH;
//                   }
                    if(coeff!=0) {
                        Eigen::VectorXi temp(4);

                        temp(0) = i;
                        temp(1) = j;
                        temp(2) = k;
                        temp(3) = l;
                        Uindex.push_back(temp);
                        Utensor.push_back(coeff);

                        index++;
                        if ( isOrbitalCorrinHart[i] and  isOrbitalCorrinHart[j] and  isOrbitalCorrinHart[k] and  isOrbitalCorrinHart[l] ) indexC++;
                    }
                }
            }
        }
    }
    write_Uijkl(index, Utensor, Uindex, indexC);
}

void gen_Uijkl_density_density(int n_spinorb, double U, double Uprime, double JH, std::vector<cmplx > & Utensor, std::vector<Eigen::VectorXi>  & Uindex ) {

//    for(int estimate_dim=0; estimate_dim<2; estimate_dim++) {
    int index=0, indexC=0;
    for(int  n=0;   n<n_spinorb; n++) {
        for(int m=0; m<n_spinorb; m++) {
            if(n!=m) {
                Eigen::VectorXi temp(4);

                if(  n/2==m/2 and n%2!=m%2)     Utensor.push_back( U          )  ;              //intra-orbital
                if(  n/2!=m/2 and n%2==m%2)     Utensor.push_back( Uprime - JH);      //inter-orbital,  F
                if(  n/2!=m/2 and n%2!=m%2)     Utensor.push_back( Uprime     );           //inter-orbital, AF
                temp(0) = n;
                temp(1) = m;
                temp(2) = n;
                temp(3) = m;
                Uindex.push_back(temp);
                index++;
                if ( isOrbitalCorrinHart[n] and  isOrbitalCorrinHart[m]  ) indexC++;
            }
        }
    }
    write_Uijkl(n_spinorb*n_spinorb, Utensor, Uindex, indexC);
}


//void write_Uijkl(int size,  Eigen::VectorXcd & Utensor, Eigen::MatrixXi & Uindex, int indexC) {
void write_Uijkl(int size, std::vector<cmplx >  & Utensor,std::vector<Eigen::VectorXi>  & Uindex, int indexC) {
    int length = (int) Utensor.size();
    FILE * Uijk_file = fopen("Uijkl.dat", "w");
    FILE * Uijk_full = fopen("Uijkl_Hart.dat", "w");
    fprintf(Uijk_file, "%d\n", indexC);
    fprintf(Uijk_full, "%d\n", length);
    int nH=0, nC=0;
    while(nH < length) {
        fprintf(Uijk_full, "%d %d %d %d %d %+0.8f %+0.8f\n",
                nH, Uindex[nH](0),
                Uindex[nH](1),
                Uindex[nH](2),
                Uindex[nH](3),  real(Utensor.at(nH)), imag(Utensor.at(nH)) );
        if( isOrbitalCorrinHart[Uindex[nH](0) ] and
                isOrbitalCorrinHart[Uindex[nH](1) ] and
                isOrbitalCorrinHart[Uindex[nH](2) ] and
                isOrbitalCorrinHart[Uindex[nH](3) ] ) {
            fprintf(Uijk_file, "%d %d %d %d %d %+0.8f %+0.8f\n",
                    nC, Hart2Corr[Uindex[nH](0)],
                    Hart2Corr[Uindex[nH](1)],
                    Hart2Corr[Uindex[nH](2)],
                    Hart2Corr[Uindex[nH](3)],  real(Utensor.at(nH)), imag(Utensor.at(nH)) );
            nC++;
        }
        nH++;
    }
    for(int n=0; n<length; n++) {
    }
    fclose(Uijk_file);
    fclose(Uijk_full);

}

void read_Uijkl( std::vector<cmplx > & Utensor, std::vector<Eigen::VectorXi>  & Uindex) {
    int index;
    double realPart, imagPart;
    int i,j,k,l;
    FILE * Uijk_full = fopen("Uijkl_Hart.dat", "r");
    fscanf(Uijk_full, "%d\n", &index);
    int nH=0;
    while(nH < index) {
        fscanf(Uijk_full, "%d %d %d %d %d %lf %lf\n",
               &nH,    &i, &j, &k, &l,  &realPart, &imagPart );
        Eigen::VectorXi temp(4);
        temp(0)=i;
        temp(1)=j;
        temp(2)=k;
        temp(3)=l;
        Uindex.push_back(temp);
        Utensor.push_back(realPart +I*imagPart);
        nH++;
    }
    fclose(Uijk_full);
}

void rot_Uijkl(
    std::vector<cmplx > & Utensor, std::vector<Eigen::VectorXi>  & Uindex,
    std::vector<cmplx > & rotUtensor, std::vector<Eigen::VectorXi>  & rotUindex,
    Eigen::MatrixXcd & SolverBasis, int n_spinorb
) {
    ifroot std::cout << "original U matrix\n";
    Eigen::MatrixXcd SolverBasis_adj = SolverBasis.adjoint();
    int length = (int) Utensor.size();
    int nH=0;
    while(nH < length) {
        ifroot std::cout << Uindex[nH](0) <<" " <<Uindex[nH](1)<<" " <<Uindex[nH](2)<<" " <<Uindex[nH](3) <<" : " << Utensor[nH]<<"\n";
        nH++;
    }
    ifroot std::cout << "\n";


    ifroot std::cout << "rot U matrix\n";
    for(int i=0; i<n_spinorb; i++) {
        for(int j=0; j<n_spinorb; j++) {
            for(int k=0; k<n_spinorb; k++) {
                for(int l=0; l<n_spinorb; l++) {
                    Eigen::VectorXi temp(4);
                    temp(0)=i;
                    temp(1)=j;
                    temp(2)=k;
                    temp(3)=l;
                    int nH=0;
                    cmplx  interaction=0;
                    while(nH < length) {
                        if(
                            Uindex[nH](0) >= n_spinorb or
                            Uindex[nH](1) >= n_spinorb or
                            Uindex[nH](2) >= n_spinorb or
                            Uindex[nH](3) >= n_spinorb)
                        {
                            std::cout <<  Uindex[nH](0)   << "\n";
                            std::cout <<  Uindex[nH](1)   << "\n";
                            std::cout <<  Uindex[nH](2)   << "\n";
                            std::cout <<  Uindex[nH](3)   << "\n";
                        }
                        interaction +=  SolverBasis_adj(i, Uindex[nH](0)) * SolverBasis_adj(j, Uindex[nH](1)) *
                                        Utensor[nH] *
                                        SolverBasis(Uindex[nH](2), k) *  SolverBasis(Uindex[nH](3),l);
                        nH++;
                    }
                    if(std::abs(interaction) > 1e-5) {
                        ifroot std::cout << i <<" " <<j<<" " <<k<<" " <<l <<" : " << interaction<<"\n";
                        rotUindex.push_back(temp);
                        rotUtensor.push_back(interaction);
                    }
                }
            }
        }
    }
//    write_Uijkl(index, Utensor, Uindex, indexC);

}
