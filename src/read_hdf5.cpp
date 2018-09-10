#ifdef READ_HDF5_dmft
#include "tight_common.h"
#include <iostream>
#include <string>
#include <complex>
#ifndef H5_NO_NAMESPACE
#ifndef H5_NO_STD
using std::cout;
using std::endl;
#endif // H5_NO_STD
#endif
//#include "/home/users1/jhsim/miniconda2/include/H5Cpp.h"
#include <H5Cpp.h>
#ifndef H5_NO_NAMESPACE
using namespace H5;
#endif

const H5std_string FILE_NAME ("input.out.h5");
const H5std_string GROUP_NAME_what ("gf");
const H5std_string DATASET_NAME( "data" );

//const int    NX = 500;        // output buffer dimensions
//const int    NY = 6;
//const int    NZ = 6;
//const int    NW = 2;
//const int    RANK_OUT = 4;



int read_hdf5 (int ATOM){

    sleep(10);
    /*
     * Output buffer initialization.
     */
    std::complex<double>  data_out[N_freq][NSpinOrbit_per_atom][NSpinOrbit_per_atom]; /* output buffer */
    for (int i = 0; i < N_freq; i++) {
        for (int j = 0; j < NSpinOrbit_per_atom; j++) {
            for (int k = 0; k < NSpinOrbit_per_atom ; k++) {
                data_out[i][j][k] = 0;
            }
        }
    }
    /*
     * Try block to detect exceptions raised by any of the calls inside it
     */
    try
    {
        /*
         * Turn off the auto-printing when failure occurs so that we can
         * handle the errors appropriately
         */
        Exception::dontPrint();
        /*
         * Open the specified file and the specified dataset in the file.
         */
        H5File file( FILE_NAME, H5F_ACC_RDONLY );
        Group what= file.openGroup( GROUP_NAME_what );
        DataSet dataset = what.openDataSet( DATASET_NAME );

        /*
         * Get the class of the datatype that is used by the dataset.
         */
        H5T_class_t type_class = dataset.getTypeClass();
        /*
         * Get class of datatype and print message if it's an integer(JHS:complex),.
         */
        if( type_class == H5T_IEEE_F64LE )
        {
            cout << "Data set has INTEGER type" << endl;
            /*
            * Get the integer datatype
             */
            IntType intype = dataset.getIntType();
            /*
             * Get order of datatype and print message if it's a little endian.
             */
            H5std_string order_string;
            H5T_order_t order = intype.getOrder( order_string );
            cout << order_string << endl;
            /*
             * Get size of the data element stored in file and print it.
             */
            size_t size = intype.getSize();
            cout << "Data size is " << size << endl;
        }
        /*
         * Get dataspace of the dataset.
         */
        DataSpace dataspace = dataset.getSpace();
        /*
         * Get the number of dimensions in the dataspace.
         */
        int rank = dataspace.getSimpleExtentNdims();
        /*
         * Get the dimension size of each dimension in the dataspace and
         * display them.
         */
        hsize_t dims_out[4];
        int ndims = dataspace.getSimpleExtentDims( dims_out, NULL);
      ifroot  cout << "rank " << rank << ", dimensions " <<
             (unsigned long)(dims_out[0]) << " x " <<
             (unsigned long)(dims_out[1]) << " x " <<
             (unsigned long)(dims_out[2]) << " x " <<
             (unsigned long)(dims_out[3]) << endl;

        dataset.read( data_out, H5T_IEEE_F64LE);
        if(mpi_rank==0) {
           std::stringstream ss;
            ss << ATOM;

            FILE *datap4 = fopen( (std::string("Gw_imp.full.dat")+ ss.str()).c_str()  , "w");
            for(int w=0; w<N_freq; w++) {
                for(int m=0; m< NSpinOrbit_per_atom; m++) {
                    for(int n=0; n< NSpinOrbit_per_atom; n++) {
                        fprintf(datap4, "%d %d %d %+0.8f %+0.8f\n",w,m,n,real(data_out[w][m][n]), imag(data_out[w][m][n]));
                    }//n
                }//m
            }//w
            fclose(datap4);
        }
        sleep(1);
    }  // end of try block
    // catch failure caused by the H5File operations
    catch( FileIException error )
    {
        cout <<"ee1\n";
        error.printError();
        return -1;
    }
    // catch failure caused by the DataSet operations
    catch( DataSetIException error )
    {
        cout <<"ee2\n";
        error.printError();
        return -1;
    }
    // catch failure caused by the DataSpace operations
    catch( DataSpaceIException error )
    {
        cout <<"ee3\n";
        error.printError();
        return -1;
    }
    // catch failure caused by the DataSpace operations
    catch( DataTypeIException error )
    {
        cout <<"ee4\n";
        error.printError();
        return -1;
    }
    return 0;  // successfully terminated
}
#endif
