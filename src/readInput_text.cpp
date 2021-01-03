#include "tight_common.h"
#include <fstream>

bool iequals(const std::string& a, const std::string& b)
{

    if (a==b) {
        return true;
    }
    else return false;
}


double read_double(const std::string &Inputfile, const std::string &keyword, bool defaultBool,  double dft) {
    double result;
    if(mpi_rank==0) {
        std::string IsitKeyword;
        std::ifstream input(Inputfile.c_str());
        int read=0;

        while(!input.eof()) {
            input >> IsitKeyword;

            if( iequals(IsitKeyword,keyword) ) {
                input >> IsitKeyword; // read "="
                input >> result;
                std::cout<<"Reading "<<keyword<<" : "<<result<<std::endl;
                read=1;
//            return result;
                break;
            }
        }
        if(read==0 and  defaultBool == false  ) {
            std::cout <<"InputFileError\n********************\n"<< keyword << "\n********************\n";
            assert(0);
        }
        if(read==0 and defaultBool == true ) {
            std::cout<<"Reading "<<keyword<<" : "<<dft<< "(default)" << std::endl;
            result = dft;
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&result, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    return result;
}


//double read_double(const std::string &Inputfile1,const std::string &Inputfile2, const std::string &keyword, double dft) {
//    double result;
//    result = read_double(Inputfile1, keyword, -100*dft);
//    if(result == -100*dft)  result =  read_double(Inputfile2, keyword, dft);
//    return  result;
//}

int read_int(const std::string &Inputfile, const std::string &keyword, bool defaultBool, int dft) {
    int result;
    if(mpi_rank==0) {
        std::string IsitKeyword;
        std::ifstream input(Inputfile.c_str());
        int read=0;

        while(!input.eof()) {
            input >> IsitKeyword;

            if( iequals(IsitKeyword, keyword) ) {
                input >> IsitKeyword; // read "="
                input >> result;
                if(mpi_rank==0) std::cout<<"Reading "<<keyword<<" : "<<result<<std::endl;
                read=1;
//            return result;
                break;
            }
        }
        if(read==0 and defaultBool==false) {
            std::cout <<"InputFileError\n********************\n"<< keyword << "\n********************\n";
            assert(0);
        }
        if(read==0 and defaultBool==true) {
            std::cout<<"Reading "<<keyword<<" : "<<dft<< "(default)" << std::endl;
            result = dft;
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&result, 1, MPI_INT, 0, MPI_COMM_WORLD);
    return result;
}
int read_int(const std::string &Inputfile1,const std::string &Inputfile2, const std::string &keyword, int dft) {
    int result;
    result = read_int(Inputfile1, keyword, -9999999999*dft);
    if(result == -9999999999*dft)  result =  read_int(Inputfile2, keyword, dft);
    return  result;
}




void read_int_array(const std::string &Inputfile, const std::string &keyword, std::vector<int>  & var,  int length, bool defaultBool, int dft) {
    if(mpi_rank==0) {
        std::string IsitKeyword;


        std::ifstream input(Inputfile.c_str());
        int read=0;
        FILE * test = fopen(Inputfile.c_str(),"r");
        if(test ==NULL) {
            std::cout << "Fail to open " << Inputfile.c_str() <<"\n";
            return;
        }
        fclose(test);
        while(!input.eof()) {
            input >> IsitKeyword;
            if( iequals(IsitKeyword,keyword) ) {
                input >> IsitKeyword; // read "="
                if(mpi_rank==0) std::cout<<"Reading "<<keyword<<": ";
                for (int i =0; i<length; i++) {
                    input >>  var[i];
                    if(mpi_rank==0) std::cout<< var[i] <<" ";
                    if(input.eof()) {
                        std::cout <<"Please check input "<<keyword<<std::endl;
                        exit(1);
                    }
                }
                if(mpi_rank==0) std::cout <<"\n";
                read=1;
            }
        }
        if(read==0 and defaultBool == false) {
            std::cout <<"InputFileError\n********************\n"<< keyword << "\n********************\n";
            exit(1);
        } else if(read==0 and defaultBool == true) {
            ifroot std::cout <<"We use default value for " << keyword <<"\n";
        }
        for (int i =0; i<length; i++) {
            if(mpi_rank==0) std::cout<< var[i] <<" ";
        }
        if(mpi_rank==0) std::cout <<"\n";
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&var[0], length, MPI_INT, 0, MPI_COMM_WORLD);
}



void read_double_array(const std::string &Inputfile, const std::string &keyword, std::vector<double> & var,  int length, bool required, double dft) {

    if(mpi_rank==0) {
        std::string IsitKeyword;
        std::ifstream input(Inputfile.c_str());
        int read=0;
        FILE * test = fopen(Inputfile.c_str(),"r");
        if(test ==NULL) {
            std::cout << "Fail to open " << Inputfile.c_str() <<"\n";
            return;
        }
        fclose(test);

        while(!input.eof()) {
            input >> IsitKeyword;

            if( iequals(IsitKeyword,keyword) ) {
                input >> IsitKeyword; // read "="
                if(mpi_rank==0) std::cout<<"Reading "<<keyword<<": ";
                for (int i =0; i<length; i++) {
                    input >>  var[i];
                    if(mpi_rank==0) std::cout<< var[i] <<" ";
                    if(input.eof()) {
                        std::cout <<"Please check input "<<keyword<<std::endl;
                        exit(1);
                    }
                }
                if(mpi_rank==0) std::cout <<"\n";
                read=1;
            }
        }
        if(read==0 and required==false) {
            ifroot std::cout <<"We don't read " << keyword <<", use default value\n";
            for (int i =0; i<length; i++) {
                var[i] = dft;
            }
        }
        else if(read==0 and required==true) {
            std::cout <<"InputFileError\n********************\n"<< keyword << "\n********************\n";
            exit(1);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&var[0], length, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}





std::string read_string(const std::string &Inputfile, const std::string &keyword, bool defaultBool, std::string defaultVal  ) {
    char  ctemp[100];
    std::string  result;

    if(mpi_rank==0) {
        std::string IsitKeyword;
        std::ifstream input(Inputfile.c_str());
        int read=0;
        while(!input.eof()) {
            input >> IsitKeyword;

            if( iequals(IsitKeyword, keyword) ) {
                input >> IsitKeyword; // read "="
                input >> result;
                if(mpi_rank==0) std::cout<<"Reading "<<keyword<<" : "<<result<<std::endl;
                read=1;
//                return result;
                break;
            }
        }
        if(read==0 and defaultBool==false) {
            std::cout <<"InputFileError\n********************\n"<< keyword << "\n********************\n";
            assert(read==1);
        }
        if(read==0 and defaultBool==true) {
            std::cout<<"Reading "<<keyword<<" : "<<defaultVal <<"(default)" <<std::endl;
            result  = defaultVal;
        }
        std::strcpy(ctemp,  result.c_str());
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(ctemp, 100, MPI_CHAR, 0, MPI_COMM_WORLD);

    result = ctemp;
    return result;
}
