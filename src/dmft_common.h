#ifndef dmft_common_INCLUDE
#define dmft_common_INCLUDE
#include <stdio.h>
#include <stdlib.h>
#include <complex>
#include <iostream>
#include <sstream>
#include <vector>



#define ifroot   if(mpi_rank==0)
//Compile option
//#define DEBERG 1 

//mathematical constant
#define     cmplx    std::complex<double>
#define   dimension_3d             3
#define   Spin_DegreeOfFreedom     2

//physical constant
//Energy=eV, length=\AA, time =  hbar / eV == ut, k = 1/(unit_vector)   
#define hbar =   1 
#define mass_e  =  13.123437306590125   //mass = eV*(ut/\AA)^2 = eV^-1 * (hbar/\AA)^2 
//  e [mass] 
//= e [kg]/[mass] * [mass]
//= e [kg] /1/eV*(hbar/A)^2
//= 9.1093821*10e-31/(1/(1.6021764*10e-19)* 1.054571*10e-34/(1.0*10e-10)* 1.054571*10e-34/(1.0*10e-10)  )
#define charge_e =  1.    // e
#define conductivity_fromUnit_to_invOhmcm  0.41082359052809814

// computational options
#define magicNumber  -999


//Lattice information
#define     ax               pi          //lattice const
#define     ay               pi          //lattice const
#define     az               pi          //lattice const


static cmplx I(0,1);
const static  double pi = 3.1415926535897932384626;
extern int mpi_numprocs, mpi_rank , mpi_namelen, node_local_numprocs, node_local_rank ;



std::string intToString(int n);






#endif
