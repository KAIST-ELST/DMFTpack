#include <Eigen/Dense>
#include <fstream>
#include "tight_common.h"
// H_inter = 1/2  \sum_{ijkl}  U_{ijkl} c^+_is1 c^+_js2 c_ls2 c_ks1
//  and U_{ijkl} = \int dr1 dr2   \i*(r1) \j*(r2) V(r1,r2) \k(r1) \l(r2) = < ij | V | kl >

void write_Uijkl(int size, std::vector<cmplx >  & Utensor,std::vector<Eigen::VectorXi>  & Uindex, int indexC) ;

void gen_Uijkl(int n_spinorb, double U, double Uprime, double JH, double JX, double JP, std::vector<cmplx >  & Utensor, std::vector<Eigen::VectorXi>  & Uindex ) {
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
                    if(i!=j and k!=l) { //c^+_i c^+_j  = c_k c_l = 0.
                        if ( ispin== kspin and jspin==lspin) {
                            //density-density terms
                            if(i==k and j==l) {
                                // c^+_is1 c^+_js2 c_js2 c_is1 =  n_is1 n_js2 ,
                                if (iorb==jorb) coeff = U;                  //intra-orbital
                                else if (iorb!=jorb) coeff = Uprime;        //inter-orbital
                                if(ispin == jspin) coeff -= JH;             //  c^+_is1 c^+_js1 c_is1 c_js1 =  - n_is1 n_js1 ,
                            }
                            //spin-flip
                            else if (korb  == jorb  and lorb  == iorb  and korb  !=  lorb and
                                     ispin !=  jspin    )   coeff = JX;             //  J c^+_ls c^+_ks' c_ls' c_ks ,
                            //pair-hopping
                            else if (korb  == lorb  and iorb  ==  jorb and korb  != iorb  and
                                     ispin != jspin     )   coeff = JP;             //  J c^+_is1 c^+_is2 c_ks2 c_ks1
                        }
                    }
                    if(coeff >= 1e-5) {
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

    int index=0, indexC=0;
    for(int  n=0;   n<n_spinorb; n++) {
        for(int m=0; m<n_spinorb; m++) {
            Eigen::VectorXi temp(4);
            double interaction =0 ;

            if(  n/2==m/2      )        interaction = U;            //intra-orbital
            else if(  n/2!=m/2 )   interaction = Uprime;           //inter-orbital,  F

            if (n%2 == m%2) interaction -= JH;
            Utensor.push_back(interaction);

            temp(0) = n;
            temp(1) = m;
            temp(2) = n;
            temp(3) = m;
            Uindex.push_back(temp);
            index++;
            if ( isOrbitalCorrinHart[n] and  isOrbitalCorrinHart[m]  ) indexC++;
        }
    }
    write_Uijkl(n_spinorb*n_spinorb, Utensor, Uindex, indexC);
}

void gen_Uijkl_sphd(int n_spinorb, double U, double J, std::vector<cmplx >  & Utensor, std::vector<Eigen::VectorXi>  & Uindex ) {
//ordered basis = {dz2, dx2, dxy, dxz, dyz}
    int index=0, indexC=0;
    double U0= U+8.0/7.0*J;
    double J1= 0.77167*J;
    double J2= 0.88645*J;
    double J3= 0.42735*J;
    double J4= 0.54212*J;
    Eigen::MatrixXd  UMat(5,5);
    Eigen::MatrixXd  JMat(5,5);


    UMat << U0, U0, U0, U0, U0,
         U0, U0, U0, U0, U0,
         U0, U0, U0, U0, U0,
         U0, U0, U0, U0, U0,
         U0, U0, U0, U0, U0;

    JMat <<  0, J2, J2, J4, J4,
         J2,  0, J3, J1, J1,
         J2, J3,  0, J1, J1,
         J4, J1, J1,  0, J1,
         J4, J1, J1, J1,  0;

    UMat = UMat-2*JMat;

    ifroot std::cout << "U: "<<U <<", J: "<<J<<"\n";
    ifroot std::cout << "UMat:\n" << UMat <<"\nJMat:\n" << JMat <<"\n";


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
                    if(i!=j and k!=l) { //c^+_i c^+_j  = c_k c_l = 0.
                        if ( ispin== kspin and jspin==lspin) {
                            //density-density terms
                            if(i==k and j==l) {
                                // c^+_is1 c^+_js2 c_js2 c_is1 =  n_is1 n_js2 ,
                                coeff = UMat(iorb,jorb);
                                if(ispin == jspin) coeff -= JMat(iorb,jorb);
                            }
                        }
                    }
                    if(coeff >= 1e-5) {
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
    ifroot std::cout << "rot U matrix\n";
    Eigen::MatrixXcd SolverBasis_adj = SolverBasis.adjoint();
    int length = (int) Utensor.size();
    int nH=0;


//    ifroot std::cout << "new(rot) U matrix\n";
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
//                        if(
//                            Uindex[nH](0) >= n_spinorb or
//                            Uindex[nH](1) >= n_spinorb or
//                            Uindex[nH](2) >= n_spinorb or
//                            Uindex[nH](3) >= n_spinorb)
//                        {
//                            std::cout <<  Uindex[nH](0)   << "\n";
//                            std::cout <<  Uindex[nH](1)   << "\n";
//                            std::cout <<  Uindex[nH](2)   << "\n";
//                            std::cout <<  Uindex[nH](3)   << "\n";
//                        }
                        interaction +=  SolverBasis_adj(i, Uindex[nH](0)) * SolverBasis_adj(j, Uindex[nH](1)) *
                                        Utensor[nH] *
                                        SolverBasis(Uindex[nH](2), k) *  SolverBasis(Uindex[nH](3),l);
                        nH++;
                    }//nH
                    if(std::abs(interaction) > 1e-4) {
                        rotUindex.push_back(temp);
                        rotUtensor.push_back(interaction);
                    }
                    //                    }
                }//l
            }//k
        }//j
    }///i
    ifroot  std::cout << "rotU size: " <<  (int) rotUtensor.size() <<"n";
    //    write_Uijkl(index, Utensor, Uindex, indexC);
}


void proj_Uijkl(
    std::vector<cmplx > & Utensor, std::vector<Eigen::VectorXi>  & Uindex,
    std::vector<cmplx > & projUtensor, std::vector<Eigen::VectorXi>  & projUindex,
    std::vector<int> & sub2full, int n_spinorb
) {
    ifroot std::cout << "proj U matrix size: ";
    int length = (int) Utensor.size();
    ifroot std::cout <<  length << "\n";
    int nH=0;


    ifroot std::cout << "new(proj) U matrix\n";
    for(int i=0; i<n_spinorb; i++) {
        for(int j=0; j<n_spinorb; j++) {
            Eigen::VectorXi temp(4);
            temp(0)=i;
            temp(1)=j;
            temp(2)=i;
            temp(3)=j;
            int nH=0;
            cmplx  interaction=0;
            while(nH < length) {
                if (  ( Uindex[nH](0) == sub2full[temp(0)] and
                        Uindex[nH](1) == sub2full[temp(1)] and
                        Uindex[nH](2) == sub2full[temp(2)] and
                        Uindex[nH](3) == sub2full[temp(3)]     )

                   ) {
                    interaction  +=  Utensor[nH];
                }
                if (  ( Uindex[nH](0) == sub2full[temp(0)] and
                        Uindex[nH](1) == sub2full[temp(1)] and
                        Uindex[nH](2) == sub2full[temp(3)] and
                        Uindex[nH](3) == sub2full[temp(2)]     )

                   ) {
                    interaction  -=  Utensor[nH];
                }
                nH++;
            }
            ifroot std::cout << i <<" " <<j<<" " <<i<<" " <<j <<" : " << interaction<<"\n";
            if(i!=j) {
                projUindex.push_back(temp);
                projUtensor.push_back(interaction);
            }
        }//j
    }///i
    //    write_Uijkl(index, Utensor, Uindex, indexC);
}
