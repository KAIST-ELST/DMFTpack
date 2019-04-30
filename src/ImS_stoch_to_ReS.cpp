#include "tight_common.h"
#include "mpi.h"
#include "TB.h"

void rewrite_retardedSw() {
    int NumCorrOrbit = NSpinOrbit_per_atom * NumCorrAtom;
    char chr1[50];
    char chr2[50];
    std::vector<std::vector<std::vector<double> > >  ImS(NumCorrOrbit);
    std::vector<std::vector<std::vector<double> > >  ReS(NumCorrOrbit);
    double S_left_realPart,  S_left_imagPart;
    double S_right_realPart, S_right_imagPart;
    for(int orbital1=0; orbital1<NumCorrOrbit; orbital1++) {
        ImS[orbital1].resize(NumCorrOrbit);
        ReS[orbital1].resize(NumCorrOrbit);
        for(int orbital2=0; orbital2<NumCorrOrbit; orbital2++) {
            ImS[orbital1][orbital2].resize(Spectral_EnergyGrid+1);
            ReS[orbital1][orbital2].resize(Spectral_EnergyGrid+1);
            ImS[orbital1][orbital2][Spectral_EnergyGrid] = magicNumber;
            ReS[orbital1][orbital2][Spectral_EnergyGrid] = magicNumber;
        }
    }


    for(int orbital1=0; orbital1<NumCorrOrbit; orbital1++) {
        for(int orbital2=0; orbital2<NumCorrOrbit; orbital2++) {

            sprintf(chr1, "%d",orbital1+1);
            sprintf(chr2, "%d",orbital2+1);
            ifroot std::cout << "Read " << (std::string("realFreq_Sw.dat_")+chr1+"_"+chr2).c_str()<<  "\n";
            FILE * dataIN  = fopen((std::string("realFreq_Sw.dat_")+chr1+"_"+chr2).c_str(), "r");

            if(dataIN!=NULL) {
                double freq_left =  0;
                double freq_right ;
                fscanf(dataIN,"%lf  %lf %lf\n",
                       &freq_right,&S_right_realPart, &S_right_imagPart);
                int k=0;
                assert(ImS[orbital1][orbital2][Spectral_EnergyGrid] == magicNumber);
                assert(ReS[orbital1][orbital2][Spectral_EnergyGrid] == magicNumber);
                while ( k<Spectral_EnergyGrid and real(E0+k*dE) < freq_right ) {
                    ReS[orbital1][orbital2][k] = 0;
                    ImS[orbital1][orbital2][k] = 0;
                    k++;
                }
                assert(ImS[orbital1][orbital2][Spectral_EnergyGrid] == magicNumber);
                assert(ReS[orbital1][orbital2][Spectral_EnergyGrid] == magicNumber);

                freq_left = freq_right;
                S_left_realPart = S_right_realPart;
                S_left_imagPart = S_right_imagPart;
                fscanf(dataIN,"%lf  %lf %lf\n",&freq_right,&S_right_realPart, &S_right_imagPart);

                while (!feof(dataIN) ) {
                    while ( k<Spectral_EnergyGrid and real(E0+k*dE) < freq_right ) {
                        double Ek = real(E0+k*dE);
                        ReS[orbital1][orbital2][k] =    ( (Ek - freq_left)*S_right_realPart + (freq_right-Ek)*S_left_realPart)    /(freq_right-freq_left);
                        ImS[orbital1][orbital2][k] =    ( (Ek - freq_left)*S_right_imagPart + (freq_right-Ek)*S_left_imagPart)    /(freq_right-freq_left);
                        k++;
                    }
                    freq_left = freq_right;
                    S_left_realPart = S_right_realPart;
                    S_left_imagPart = S_right_imagPart;
                    fscanf(dataIN,"%lf  %lf %lf\n",&freq_right,&S_right_realPart, &S_right_imagPart);
                }

                assert(ImS[orbital1][orbital2][Spectral_EnergyGrid] == magicNumber);
                assert(ReS[orbital1][orbital2][Spectral_EnergyGrid] == magicNumber);

                while ( k<Spectral_EnergyGrid ) {
                    ReS[orbital1][orbital2][k] = 0;
                    ImS[orbital1][orbital2][k] = 0;
                    k++;
                }
                fclose(dataIN);
            }
            else {
                ifroot std::cout << "##Warning: We can't read realFreq_Sw.dat_a_b##\n" ;
                int k=0;
                while ( k<Spectral_EnergyGrid ) {
                    ReS[orbital1][orbital2][k] = 0;
                    ImS[orbital1][orbital2][k] = 0;
                    k++;
                }
            }
            assert(ImS[orbital1][orbital2][Spectral_EnergyGrid] == magicNumber);
            assert(ReS[orbital1][orbital2][Spectral_EnergyGrid] == magicNumber);
        }
    }//orbital
    ifroot{
        FILE * datap1 = fopen("realFreq_Sw.dat", "w");
        for( int n=0; n<Spectral_EnergyGrid; n++) {
            fprintf(datap1, "%+0.3f   ", real(E0+n*dE));
            for(int orbital1=0; orbital1<NumCorrOrbit; orbital1++) {
                fprintf(datap1, "%+0.8f  %0.8f    ", ReS[orbital1][orbital1][n],  ImS[orbital1][orbital1][n]    );
            }
            fprintf(datap1, "\n");
        } //n
        fclose(datap1);


        datap1 = fopen("realFreq_Sw.full.dat", "w");
        for( int n=0; n<Spectral_EnergyGrid; n++) {
            for(int orbital1=0; orbital1<NumCorrOrbit; orbital1++) {
                for(int orbital2=0; orbital2<NumCorrOrbit; orbital2++) {
                    fprintf(datap1, "%d %d %d  %+0.8f  %0.8f\n",n, orbital1, orbital2, ReS[orbital1][orbital2][n],  ImS[orbital1][orbital2][n]    );
                }
            }
        } //n
        fclose(datap1);
        sleep(5);
    }
    MPI_Barrier(MPI_COMM_WORLD);
}
