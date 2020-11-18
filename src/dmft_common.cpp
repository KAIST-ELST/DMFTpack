#include "ImgTimeFtn.h"
#include <cstring>
#include <stdlib.h>
#include "tight_common.h"
#ifdef READ_HDF5_dmft
#include <H5Cpp.h>
#endif


void Wait_Run( char fileName[100], int checkTime, int mpi_rank,int maxTime  ) {
//    ifroot system( (std::string("rm ") + fileName).c_str() );
    int solvercheck=0;
    int STOP_message =0;
    while (1) {
        if(mpi_rank==0) {
            if( 0!= access(fileName, F_OK) ) {
                solvercheck++;
            }
            else {
                STOP_message = 1;
            }
        }/*ifroot*/
        sleep(checkTime);
//        if( (mpi_rank==0 and STOP_message==1) or mpi_rank!=0) MPI_Bcast(&STOP_message, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&STOP_message, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if(STOP_message==1) break;
    }
    sleep(10);
    ifroot  std::cout <<"imp.sol. Finished..\n";
}



std::string intToString(int n) {
    std::stringstream s;
    s << n;
    return s.str();
}
