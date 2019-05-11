#include "tight_common.h"
#include "model.h"
#include <fstream>
#include "ImgTimeFtn.h"

//LongRangeOrder =
//0 1 2 3
//2 3 0 1
//HartreeAtoms =    index from "1"
//HartreeOrbitlas = index from "1"



//HartreeOrbital_idx = starting from '0'
Eigen::Vector3d cross_product(Eigen::Vector3d a, Eigen::Vector3d b) {
    Eigen::Vector3d c;
    c(0) =   a(1)*b(2)-a(2)*b(1);
    c(1) =   a(2)*b(0)-a(0)*b(2);
    c(2) =   a(0)*b(1)-a(1)*b(0);
    return c;
}


void analysis_example(std::string  scfout_file);


//mpi
int mpi_numprocs, mpi_rank, mpi_namelen, node_local_rank, node_local_numprocs;
char processor_name[100];

//TB.cpp
int NumSubSpace;
int *FromOrbitalToAtom;
int *FromOrbitalToLocalOrbital_DFT;
int *FromOrbitalToAtom_model;
int **HartrRange;
//int **HartrRange_DFT;
//std::vector<Eigen::MatrixXcd> densityMatDFT;


//downfolding
double   lower_model_window, upper_model_window;
double   lower_spectrum_window, upper_spectrum_window;
int downfolding;


//TB
int Nkpath;
Eigen::Vector3d RecipUnitVector_b1;
Eigen::Vector3d RecipUnitVector_b2;
Eigen::Vector3d RecipUnitVector_b3;
double KpathPoint[30*3];
int num_subshell;
Eigen::VectorXi subshell;
Eigen::VectorXi Rydberg_set;
Eigen::VectorXi rot_sym;



int  k_pointx, k_pointy, k_pointz,  k_grid=0, mu_adjust, NumHartrOrbit, NumCorrAtom, BraLatt_x,BraLatt_y,BraLatt_z, Measure_nn, NumAtom, NumOrbit;
int NumAtom_per_cluster, NumCluster, NumHartrOrbit_per_cluster;
int H0_from_OpenMX;
//int  SOCCal,  interOrbitHop, currentIt=0, DFTIt=0;
int  DFTIt=0;

double beta,  EnergyUnit, doublecounting, NumberOfElectron, Band_renorm, TotalCorrEle_in_Imp;
Eigen::Matrix2cd * Zeeman_field_spin;
int  N_freq, N_tau,  magnetism;
unsigned long long Num_MC_steps;
std::string mode;
std::string solver_bin;
double infinitesimal; // EWIN;
double nominal_charge;
//std::string delta_file;
std::string system_name;
std::string SOLVERexe, SOLVERdir;
std::string dctype, localOrbitalType, SOLVERtype, Lowlevel_SOLVERtype, Coulombtype;
//int * LongRangeOrder
//int * LongRangeOrder_Hart, *LongRangeOrder_DFT;
//int *accumulated_Num_SpinOrbital ;
int * KS2Hartr;

int  maxDmftIt, maxDFTIt,maxTime=60,restart, NSpinOrbit_per_atom, N_peratom_HartrOrbit;
double  UHubb, Uprime, JHund,mixing;

int knum, knum_mpiGlobal, myksta,mykend, NsBath=-1 ;
//int SOLVERtype;
int mixingType, mixingStart, mixingFtn;
cmplx dE, E0;
int Spectral_EnergyGrid;
cmplx dw, w0;
//cmplx ** delta_w ;
//cmplx ** delta_t ;
Eigen::MatrixXcd impurity_site_Hamiltonian;
//cmplx **  Sw_inf;
//cmplx  ** Sw_Hartree, ** Sw_doublecounting;
//cmplx  **  NumMatrix;
//Eigen::MatrixXcd Sw_Hartree, Sw_doublecounting;
Eigen::MatrixXcd  NumMatrix;
//double   **  UMatrix;
int * isOrbitalCorr, *isOrbitalCorrinHart,  * isOrbitalHartr, * isOrbitalHartrDFT, *CorrIndex,  *HartrIndex, *HartrIndex_inDFT, *Hart2Corr;
//*CorrToHartr,
Eigen::MatrixXi CorrToHartr;
int impurityBasisSwitch;

//Eigen::MatrixXi Uindex;
//Eigen::VectorXcd Utensor;



std::vector<Eigen::VectorXi> Uindex;
std::vector<cmplx > Utensor;
std::vector<Eigen::VectorXi> Uindex_stronglyCorr;
std::vector<cmplx > Utensor_stronglyCorr;


//ImgFreqFtn doublecounting_w(0);


void read_Uijkl(   std::vector<cmplx >  & Utensor,  std::vector<Eigen::VectorXi>  & Uindex) ;

void read_inputFile(const std::string &hamiltonian) {





    downfolding = 1;
    k_grid=0;
//ED outdated
//input var. with  dependence

    infinitesimal  = read_double(std::string("input.parm"), std::string("Infinitesimal"), true, 0.05)  ;  //NOTE : EWIN / Spectral_E should be smaller than infinitesimal
    Spectral_EnergyGrid             = read_double(std::string("input.parm"), std::string("SPECTRAL_ENERGY_GRID"), true, 1000) ;



//Computational options
    system_name=  read_string(std::string("input.parm"), std::string("SYSTEM_NAME"), false)   ;
    mu_adjust =   read_int(std::string("input.parm"), std::string("mu_adjust"), true, 1)           ;   //0=fixed mu, 1=Adjusting mu to Filling Factor
    maxDFTIt    =read_int(std::string("input.parm"), std::string("MAX_DFT_ITER"),true, 1);
    maxDmftIt    =read_int(std::string("input.parm"), std::string("MAX_DMFT_ITER"),true, 20);
    restart =     read_int(std::string("input.parm"), std::string("RESTART"), false)  ;
    mixing =      read_double(std::string("input.parm"), std::string("MIXING"), true, 0.8 )  ;
    mixingFtn =  read_int(std::string("input.parm"), std::string("MIXING_FUNCTION"), true ,1)  ;    //0=hybridization,  [1=self-energy]


    /*solver option*/
    SOLVERtype  = read_string(std::string("input.parm"), std::string("SOLVER_TYPE"), false);//TB, ALPS_CTSEG, ALPS_CTHYB, RUTGERS_CTSEG, RUTGERS_CTHYB, SC2PT, IPT, 2PT
    std::string solverexe_default, solverdir_default;
    if(SOLVERtype.find(std::string("ALPS"))  != std::string::npos  ) {
        solverexe_default = std::string("hybridization");
    }
    else if(SOLVERtype.find(std::string("RUTGERS"))  != std::string::npos  ) {
        solverexe_default = std::string("ctqmc_rutgers");
    }
    solverdir_default = std::string("~/bin/");
    SOLVERexe   =read_string(std::string("input.parm"), std::string("SOLVER_EXE"), true, solverexe_default);
    SOLVERdir   =read_string(std::string("input.parm"), std::string("SOLVER_DIR"), true, solverdir_default);
    Lowlevel_SOLVERtype  = read_string(std::string("input.parm"), std::string("Lowlevel_SOLVERTYPE"),true, std::string("HF"));  //HF, 2PT


    N_freq    =    read_int(std::string("input.solver"),std::string("input.solver1"), std::string("N_MATSUBARA"),false);
    N_tau =        read_int(std::string("input.solver"),std::string("input.solver1"), std::string("N_TAU"),false)  ;
    maxTime =     60;

//
//Lattice information
//
    H0_from_OpenMX =   read_int(std::string("input.parm"), std::string("H0_FROM_OPENMX"),true, 0);
    if(H0_from_OpenMX!=0  ) {
        ifroot analysis_example(system_name+std::string(".scfout"));
        MPI_Barrier(MPI_COMM_WORLD);
        sleep(5);
    }
    NumAtom_per_cluster    =   read_int(std::string("input.parm"), std::string("N_ATOMS_CLUSTER"), true , 1)   ;   //num atoms for each cluster. Thus, num cluster = NumAtom/NumAtom_per_cluster
    NumAtom          =   read_int(std::string("input.parm"), std::string("N_ATOMS"), false )   ;   //total Num of atoms
    NumCorrAtom      =   read_int(std::string("input.parm"), std::string("N_CORRELATED_ATOMS"), false) ; //Num of atoms, which inlude correlated and/or HF orbitals
    NumberOfElectron = read_double(std::string("input.parm"), std::string("N_ELECTRONS"), false, -1) ;    // Number of electron per unit cell = NumOrbit * FillingFactor
    EnergyUnit =    read_double(std::string("input.parm"), std::string("EnergyUnit"), true, 1)    ;    //Input (Hopping Hamiltonian) -> eV

    std::vector<int> k_point(3);
    read_int_array(std::string("input.parm"), std::string("K_POINTS"),     k_point,     3, false, -1) ;
    k_pointx = k_point[0];
    k_pointy = k_point[1];
    k_pointz = k_point[2];

    magnetism   =         read_int(std::string("input.parm"), std::string("MAGNETISM"),true,   2)    ; //0=para, 1=col, 2=noncol


//Impurity information
    UHubb =                   read_double(std::string("input.solver"), std::string("U"), true, -2)    ;
    Uprime =                  read_double(std::string("input.solver"), std::string("U'"),true, UHubb)    ;
    JHund =                   read_double(std::string("input.solver"), std::string("J"), true, -2)    ;
    Coulombtype =             read_string(std::string("input.solver"), std::string("COULOMB_TYPE"), true, "dd")    ; //, dd (density-density) k(Kanamori), r (read Uijkl_Hart.dat)
    beta                    = read_double(std::string("input.solver"), std::string("BETA"), false, -1) ;
    NSpinOrbit_per_atom     =    read_int(std::string("input.solver"), std::string("N_ORBITALS"),false)  ; //correlated orbital per Correlated atom
    N_peratom_HartrOrbit    =    read_int(std::string("input.solver"), std::string("N_HARTREE_ORBITALS"), true, NSpinOrbit_per_atom)   ; //Hartree orbital per correlated atom
    N_peratom_HartrOrbit += NSpinOrbit_per_atom;
    impurityBasisSwitch  =       read_int(std::string("input.solver"), std::string("impurityBasisSwitch"), true , 0);

//etc..
    Band_renorm =         read_double(std::string("input.parm"), std::string("BAND_RENORM"), true,  1)    ;


//projector
    dctype           =   read_string(std::string("input.parm"), std::string("DC_TYPE"), true,std::string("fll")); //fll,  nominal, fll_Uprime
    if(dctype == std::string("none")) {
        doublecounting =     0;
    }
    else {
        doublecounting =     1;
        if ( dctype.find(std::string("nominal")) != std::string::npos      ) {
            nominal_charge   =   read_double(std::string("input.parm"), std::string("N_d"), false,  -1);
        }
    }
    localOrbitalType = read_string(std::string("input.parm"), std::string("LOCAL_ORBITAL_TYPE"), true, std::string("nao_direct")); //nao_direct, nao_recip, nao_hyb, nao_sc,  and... _F

    upper_model_window = read_double(std::string("input.parm"), std::string("MODEL_WINDOW_U"), true,  10 ) ;
    lower_model_window = read_double(std::string("input.parm"), std::string("MODEL_WINDOW_D"), true, -10) ;

    upper_spectrum_window = read_double(std::string("input.parm"), std::string("Spectrum_window_u"), true,  10 ) ;
    lower_spectrum_window = read_double(std::string("input.parm"), std::string("Spectrum_window_d"), true, -10) ;
    if(upper_model_window > upper_spectrum_window)  upper_spectrum_window =  std::abs(upper_model_window*1.5);
    if(lower_model_window < lower_spectrum_window)  lower_spectrum_window = -std::abs(lower_model_window*1.5);


    E0 = lower_spectrum_window;
    dE =  (upper_spectrum_window-lower_spectrum_window) / (Spectral_EnergyGrid-1);
    assert(Spectral_EnergyGrid>1);


    NumHartrOrbit  =          N_peratom_HartrOrbit * NumCorrAtom;
    assert(N_peratom_HartrOrbit>=NSpinOrbit_per_atom);
    assert( NumAtom>=NumCorrAtom );



//other variables...
    NumCluster = NumCorrAtom / NumAtom_per_cluster;
    NumHartrOrbit_per_cluster = N_peratom_HartrOrbit * NumAtom_per_cluster;
    ifroot std::cout << "We have " << NumCluster <<" clusters with "<< NumHartrOrbit_per_cluster <<" orbitals for each cluster\n";



//DOS & Bandoption
    //real or imaginary time
    dw = I*(2*pi)/beta;
    w0 = I*pi/beta;




//Lattice information
    FILE *dataIN = fopen(hamiltonian.c_str(), "r");
    std::vector<int> accumulated_Num_SpinOrbital_local(NumAtom+1);
    accumulated_Num_SpinOrbital_local[0]=0;
    if (dataIN == NULL) {
        printf("Cannot open Hopping data file...\n");
        exit(1);
    }
    NumOrbit=0;
    for(int i=0; i<NumAtom; i++) {
        int temp1=0;
        fscanf(dataIN,"Norbital of Atom = %d\n", &temp1);
        accumulated_Num_SpinOrbital_local[i+1]=temp1+accumulated_Num_SpinOrbital_local[i];
        NumOrbit += temp1;
    }
    rewind(dataIN);
    fclose(dataIN);
    ifroot std::cout<<"Reading Num of orbit : "<<NumOrbit <<"\n";


    std::vector<int> num_subshell_forAtom(NumAtom);
    for(int i=0; i<NumAtom; i++)  num_subshell_forAtom[i] = 1;
    read_int_array(std::string("input.parm"), std::string("num_subshell"),      num_subshell_forAtom,     NumAtom,true,  1) ;


    num_subshell=0;
    for(int i=0; i<NumAtom; i++) {
        num_subshell += num_subshell_forAtom[i];
    }
    subshell.setZero(num_subshell+1);
    Rydberg_set.setZero(num_subshell);
    std::vector<int>    subshell_temp(num_subshell);
    std::vector<int> Rydberg_set_temp(num_subshell);
    std::vector<int>     rot_sym_temp(num_subshell);
//    if(num_subshell == NumAtom) {
    for(int i=0; i<num_subshell; i++) {
        subshell_temp[i] = accumulated_Num_SpinOrbital_local[i+1] - accumulated_Num_SpinOrbital_local[i];
        Rydberg_set_temp[i]  = 1;
    }
//   }
    if(num_subshell != NumAtom) {
        read_int_array(std::string("input.parm"), std::string("subshell"),        subshell_temp,        num_subshell,true,-1) ;
        read_int_array(std::string("input.parm"), std::string("Rydberg_set"),     Rydberg_set_temp,     num_subshell,true,-1) ;
    }
    for(int i=0; i<num_subshell; i++) rot_sym_temp[i] = i;
    read_int_array(std::string("input.parm"), std::string("rot_sym"),     rot_sym_temp,     num_subshell,true,-1) ;
    rot_sym.setZero(num_subshell);
    int test=0;
    subshell(0)=0;
    for(int i=0; i<num_subshell; i++) {
        Rydberg_set[i] = Rydberg_set_temp[i];
        rot_sym[i] = rot_sym_temp[i];
        subshell[i+1] = subshell_temp[i] + subshell[i];
        test += subshell_temp[i] ;
    }
    int temp=0;
    for(int i=0; i<NumAtom; i++) {
        for(int j=temp; j<num_subshell_forAtom[i]+temp; j++) {
            rot_sym[j] += i*1000 ;
        }
        temp += num_subshell_forAtom[i];
    }
    for(int i=0; i<num_subshell; i++) {
        ifroot std::cout <<"SubShell "<<i<<",\tdim:" << subshell_temp[i]<<"\t rot_sym:" << rot_sym[i] <<"\n";
    }
    ifroot std::cout <<"total num of states nlm = " << test<<"\n";
    assert(test==NumOrbit);



//Memory allocation
    int i=0;
    FromOrbitalToAtom =new int [NumOrbit];
    FromOrbitalToLocalOrbital_DFT=new int [NumOrbit];
    for (int ob=0; ob< accumulated_Num_SpinOrbital_local[NumAtom]; ob++) {
        if (ob >= accumulated_Num_SpinOrbital_local[i+1]) i++;
        FromOrbitalToAtom[ob]=i;
        FromOrbitalToLocalOrbital_DFT[ob]=ob-accumulated_Num_SpinOrbital_local[i];
    }
    KS2Hartr       = new int    [NumOrbit];

    NumMatrix.setZero(N_peratom_HartrOrbit*NumCorrAtom, N_peratom_HartrOrbit*NumCorrAtom);
/////////////////////////////////////////////
//And read additional information
/////////////////////////////////////////////
//DFT+DMFT
//NumMatDFT = new NumMatDFT[k_pointx*k_pointy*k_pointz];
//Initial occupation & UMatrix
    /*set Hartree space index*/
    std::vector<int>  HartreeAtom_idx(NumCorrAtom);
    read_int_array(std::string("input.parm"), std::string("HARTREE_ATOMS"),     HartreeAtom_idx,     NumCorrAtom, false,-1) ;
    for(int i=0; i<NumCorrAtom; i++) {
        HartreeAtom_idx[i]-=1;   //[0,1,2,...NumCorrAtom-1]   // number of element = NumCorrAtom
    }

    std::vector<int>HTO_idx_temp(2);
    read_int_array(std::string("input.parm"), std::string("HARTREE_ORBITALS_RANGE"),  HTO_idx_temp,  2, false, -1) ;
    HTO_idx_temp[0]-=1;   //[0,N_peratom_HartrOrbit]   //used in for..   i=HTO_idx[0];i<HTO_idx[1];i++

    std::vector<int>HartreeOrbital_idx(NumCorrAtom * 2);
    for(int i=0; i<NumCorrAtom; i++) {
        HartreeOrbital_idx[i*2+0] = HTO_idx_temp[0] + accumulated_Num_SpinOrbital_local[HartreeAtom_idx[i]];
        HartreeOrbital_idx[i*2+1] = HTO_idx_temp[1] + accumulated_Num_SpinOrbital_local[HartreeAtom_idx[i]];
    }
    //1s= 1:2  2s=3:4
    //1p= 5:10 2p=11:16
    //1d=17:26

    std::vector<double>  Zeeman_spin_read(NumCorrAtom * 4);
    Zeeman_field_spin = new Eigen::Matrix2cd [NumAtom];
    for(int at=0; at<NumAtom; at++) {
        Zeeman_field_spin[at].setZero(2,2);
    }
    read_double_array(std::string("input.parm"), std::string("Zeeman_spin"), Zeeman_spin_read, NumCorrAtom * 4, false, 0.0);
    Eigen::Matrix2cd  sigmax, sigmay, sigmaz;
    sigmax << 0.,1.,1.,0.;
    sigmay << 0.,-1*I,I,0.;
    sigmaz << 1.,0.,0.,-1.;
    for(int at=0; at<NumCorrAtom; at++) {
        double unitVectorNorm = std::sqrt( std::pow(Zeeman_spin_read[at*4+1],2) +   std::pow(Zeeman_spin_read[at*4+2],2) + std::pow(Zeeman_spin_read[at*4+3],2) );
        Zeeman_spin_read[at*4+1] /= (unitVectorNorm +1e20);
        Zeeman_spin_read[at*4+2] /= (unitVectorNorm +1e20);
        Zeeman_spin_read[at*4+3] /= (unitVectorNorm +1e20);
        Zeeman_field_spin[HartreeAtom_idx[at]] =
            Zeeman_spin_read[at*4+0]/2* ( Zeeman_spin_read[at*4+1] * sigmax + Zeeman_spin_read[at*4+2] * sigmay + Zeeman_spin_read[at*4+3] * sigmaz);
    }


    /*
    * Define Correlated space
    * */
    isOrbitalCorrinHart = new int [N_peratom_HartrOrbit *NumCorrAtom];
    Hart2Corr =           new int [N_peratom_HartrOrbit * NumCorrAtom];
    HartrIndex =  new int [N_peratom_HartrOrbit * NumCorrAtom];
    CorrToHartr.setZero(NumCorrAtom, NSpinOrbit_per_atom);
    CorrIndex = new int   [NSpinOrbit_per_atom * NumCorrAtom];

    isOrbitalCorr =       new int [NumOrbit];
    isOrbitalHartr= new int  [NumOrbit];
    HartrRange = new int * [NumCorrAtom];
//    HartrRange_DFT = new int * [NumCorrAtom];
    for(int i=0; i<NumCorrAtom; i++) {
        HartrRange[i] = new int [2];
//        HartrRange_DFT[i] = new int [2];
    }
    FromOrbitalToAtom_model=new int [NumOrbit];

    setCorrelatedSpaceIndex(HartreeOrbital_idx, NumCorrAtom);

    isOrbitalHartrDFT= new int  [NumOrbit];
    HartrIndex_inDFT = new int [N_peratom_HartrOrbit    * NumCorrAtom];
    for(int i=0; i<NumOrbit; i++) {
        isOrbitalHartrDFT[i] = isOrbitalHartr[i];
    }
    for(int at=0; at<NumCorrAtom; at++) {
        for (int i=0; i<N_peratom_HartrOrbit; i++) {
            HartrIndex_inDFT[at*N_peratom_HartrOrbit + i] = HartrIndex[at * N_peratom_HartrOrbit +i];
        }
    }

    /*Construct Utensor */
    std::vector<Eigen::VectorXi> Uindex_atom;
    std::vector<cmplx > Utensor_atom;
    if (Coulombtype == std::string("r")     ) {
        read_Uijkl(   Utensor_atom,  Uindex_atom) ;
    }
//    else {}
//        if ( SOLVERtype.find(std::string("SEG")) != std::string::npos ) {}
    else if ( Coulombtype == "dd" ) {
        ifroot std::cout << "gen_Uijkl_dd\n";
        gen_Uijkl_density_density(N_peratom_HartrOrbit, UHubb, Uprime, JHund, Utensor_atom, Uindex_atom);
    }
    else if (Coulombtype =="k")  gen_Uijkl(N_peratom_HartrOrbit, UHubb, Uprime, JHund, Utensor_atom, Uindex_atom);
//    {}


    for(int  cl =0; cl < NumAtom_per_cluster; cl++) {
        for(int idx1=0;  idx1 <  Utensor_atom.size(); idx1++) {
            Eigen::VectorXi temp(4);
            temp(0) = Uindex_atom[idx1](0)  + cl *N_peratom_HartrOrbit;
            temp(1) = Uindex_atom[idx1](1)  + cl *N_peratom_HartrOrbit;
            temp(2) = Uindex_atom[idx1](2)  + cl *N_peratom_HartrOrbit;
            temp(3) = Uindex_atom[idx1](3)  + cl *N_peratom_HartrOrbit;
            Uindex.push_back(temp);
            Utensor.push_back(Utensor_atom[idx1]);
        }
    }

    for(int idx1=0;  idx1 <  Utensor_atom.size(); idx1++) {
        if (    isOrbitalCorrinHart[Uindex_atom[idx1](0)] and
                isOrbitalCorrinHart[Uindex_atom[idx1](1)] and
                isOrbitalCorrinHart[Uindex_atom[idx1](2)] and
                isOrbitalCorrinHart[Uindex_atom[idx1](3)] ) {

            Eigen::VectorXi temp(4);
            temp(0) = Hart2Corr[Uindex_atom[idx1](0)];
            temp(1) = Hart2Corr[Uindex_atom[idx1](1)];
            temp(2) = Hart2Corr[Uindex_atom[idx1](2)];
            temp(3) = Hart2Corr[Uindex_atom[idx1](3)];
            Uindex_stronglyCorr.push_back(temp);
            Utensor_stronglyCorr.push_back(Utensor_atom[idx1]);
        }
    }

    /*input dependency parameters*/
    if (SOLVERtype==std::string("TB")) {
        mode =      read_string(std::string("input.parm"), std::string("MODE"),false)  ; //dos, qsdos, band, qsband
        if(mode.find(std::string("band")) != std::string::npos) {
            Nkpath =   read_int   (std::string("input.parm"), std::string("N_K_PATH"),false)   ;
            assert(Nkpath<9);
            std::vector<double> KpathPoint_vec((Nkpath+1)*3);
            read_double_array (std::string("input.parm"), std::string("K_PATH"),KpathPoint_vec, (Nkpath+1)*3,  true   ) ;
            for(int i=0; i<(Nkpath+1)*3; i++) KpathPoint[i] = KpathPoint_vec[i];
            k_grid =   read_int   (std::string("input.parm"), std::string("K_GRID_BAND"),40)   ; //band k grid per band line
            std::vector<double> Rn(9);
            read_double_array(std::string("input.parm"), std::string("UNIT_VECTORS"),     Rn,     9,true) ;
            Eigen::Vector3d UnitVector_a1(Rn[0], Rn[1],Rn[2]);
            Eigen::Vector3d UnitVector_a2(Rn[3], Rn[4],Rn[5]);
            Eigen::Vector3d UnitVector_a3(Rn[6], Rn[7],Rn[8]);
            ifroot std::cout <<"lattice vectors: "
                             << UnitVector_a1<<"\n"
                             << UnitVector_a2<<"\n"
                             << UnitVector_a3<<"\n";
            Eigen::Vector3d temp =  cross_product(UnitVector_a2, UnitVector_a3);
            double volume = std::abs(UnitVector_a1.dot( temp));

            Eigen::Vector3d temp1 =  cross_product(UnitVector_a2, UnitVector_a3);
            Eigen::Vector3d temp2 =  cross_product(UnitVector_a3, UnitVector_a1);
            Eigen::Vector3d temp3 =  cross_product(UnitVector_a1, UnitVector_a2);
            RecipUnitVector_b1 = 2*pi * temp1 /volume;
            RecipUnitVector_b2 = 2*pi * temp2 /volume;
            RecipUnitVector_b3 = 2*pi * temp3 /volume;
            ifroot std::cout <<"Reciprocal vectors: "
                             << RecipUnitVector_b1<<"\n"
                             << RecipUnitVector_b2<<"\n"
                             << RecipUnitVector_b3<<"\n";
        }
        mu_adjust = 0;                                                                      //chemical pot. is given by previous DMFT cal
        mixingType=0;
    }
}/*read input*/




void setCorrelatedSpaceIndex( std::vector<int> HartreeOrbital_idx, int NumCorrAtom   ) {
    int test=0;
    /*Strong correlation space (Hartree space)*/
    for(int i=0; i<NumOrbit; i++) {
        isOrbitalHartr[i] = 0;
        FromOrbitalToAtom_model[i] = 0;
    }
    for(int i=0; i<NumCorrAtom; i++) {
        for(int j=HartreeOrbital_idx[i*2+0]; j<HartreeOrbital_idx[i*2+1]; j++) {
            isOrbitalHartr[ j ] = 1;
            FromOrbitalToAtom_model[j] = i;
        }
        HartrRange[i][0] = HartreeOrbital_idx[i*2+0];
        HartrRange[i][1] = HartreeOrbital_idx[i*2+1];
    }
    /*HartrIdex*/
    test=0;
    for(int at=0; at<NumCorrAtom; at++) {
        for(int i=HartrRange[at][0]; i<HartrRange[at][1]; i++) {
            HartrIndex[test] = i ;
            test++;
        }
    }

    int k =0;
    for(int at=0; at<NumCorrAtom; at++) {
        for(int j=HartreeOrbital_idx[at*2+0]; j<HartreeOrbital_idx[at*2+1]; j++) {
            KS2Hartr[j] = k;
            k++;
        }
    }


    //Dynamical correlation space
    std::vector<int> isOrbitalCorrinHart_vec(N_peratom_HartrOrbit*NumCorrAtom);
    for (int i=0; i<N_peratom_HartrOrbit*NumCorrAtom; i++) isOrbitalCorrinHart_vec[i] = 1;
    read_int_array(std::string("input.parm"), std::string("DYNAMIC_ORBITALS"),  isOrbitalCorrinHart_vec, N_peratom_HartrOrbit * NumCorrAtom, true, 1);
    for (int i=0; i<N_peratom_HartrOrbit*NumCorrAtom; i++) isOrbitalCorrinHart[i] = isOrbitalCorrinHart_vec[i];

    ifroot std::cout << "is Dynamical orbital\n";
    for(int i=0; i<NumOrbit; i++) isOrbitalCorr[i] = 0;
    int temp=0;
    for(int i=0; i<NumOrbit; i++) {
        if(isOrbitalHartr[i] == 1) {
            ifroot std::cout << "correlated orbitals: " <<i<<"\n";
            if (isOrbitalCorrinHart[temp]==1) isOrbitalCorr[i] = 1;
            temp++;
        }
    }

    int        test2=0;
    for(int at=0; at<NumCorrAtom; at++) {
        int test=0; //test=0,1,...  NSpinOrbit_per_atom;
        ifroot std::cout << "Hartree orbital range for atom "<<at<< ": from " << HartreeOrbital_idx[at*2+0]  <<"to" << HartreeOrbital_idx[at*2+1] <<"\n";
        for(int i=HartreeOrbital_idx[at*2+0]; i<HartreeOrbital_idx[at*2+1]; i++) {
            int j= i-(HartreeOrbital_idx[at*2+0]);
            if(isOrbitalCorr[i]) {
                CorrToHartr(at,test) = at*N_peratom_HartrOrbit + j ;
                Hart2Corr[at*N_peratom_HartrOrbit+j] = test2;
                test++;
                test2++;
            }
            else {
                Hart2Corr[j] = -999999;
            }
        }
    }
    if(not( test2 == NumCorrAtom*NSpinOrbit_per_atom)) {
        ifroot std::cout
                << "Please check input,  total number of Correlated orbital !=  NumCorrAtom*NSpinOrbit_per_atom\n"
                << test2<<", " <<NumCorrAtom*NSpinOrbit_per_atom ;
        MPI_Barrier(MPI_COMM_WORLD);
        exit(1);
    }
    test2=0;
    for(int at=0; at<NumCorrAtom; at++) {
        test=0;
        for(int i=HartreeOrbital_idx[at*2+0]; i<HartreeOrbital_idx[at*2+1]; i++) {
            if(isOrbitalCorr[i]) {
                CorrIndex[test2] = HartrIndex[CorrToHartr(at,test)] ;
                test++;
                test2++;
            }
        }
    }
    assert(test2== NSpinOrbit_per_atom * NumCorrAtom);
    MPI_Barrier(MPI_COMM_WORLD);
    ifroot std::cout << "CorrelatedSpace index.\n";
}//setCorrelatedSpaceIndex;





void on_the_fly_control() {
    mixing = read_double(std::string("input.parm"), std::string("MIXING"), true, 0.8)  ;
    maxTime =     60;
    maxDmftIt    =   read_int(std::string("input.parm"), std::string("MAX_DMFT_ITER"),true, 20);
    maxDFTIt    =read_int(std::string("input.parm"), std::string("MAX_DFT_ITER"),true,1);
}

int isSameAtom(int o1, int o2, int model) {
    assert(isOrbitalHartr[o1] and isOrbitalHartr[o2]);
    if(model==0) {
        if (FromOrbitalToAtom[o1] == FromOrbitalToAtom[o2]) return 1;
        else return 0 ;
    }
    else {
        if (FromOrbitalToAtom_model[o1] == FromOrbitalToAtom_model[o2]) return 1;
        else return 0 ;
    }
}

//int nonNeg(int & a) {
//    assert (a>=0);
//    return a;
//}


cmplx operator*( int a, cmplx b)
{
    return (double) a* b;
}
cmplx operator*( cmplx a, int b)
{
    return  a*(double) b;
}
cmplx operator/( cmplx  a, int b)
{
    return  a/(double) b;
}
bool operator!=( cmplx  a, int b)
{
    return  a != (double) b;
}
bool operator==( cmplx  a, int b)
{
    return  a == (double) b;
}



// basis:
// default = truncated (low-energy) space = valence
// DFT = (R,\alpha) index from DFT
// KS  = band index from DFT = valence
// valence = dp space ( = truncated space)


//NSpinOrbit_per_atom   = # of Corr orbitals per atom
//N_peratom_HartrOrbit  = # of Hartree orbitals per atom

//isOrbitalHartr[NumOrbit]  = 1 for hartree, 0 else...
//isOrbitalHartrDFT[NumOrbit]  = 1 for hartree, 0 else...
//isOrbitalCorr[NumOrbit]
//isOrbitalCorrinHart = new int [N_peratom_HartrOrbit *NumCorrAtom];

//KS2Hartr     = from KS to  Hartr
//CorrIndex    = from corr     index to whole KS orbital index
//HartrIndex   = from hartreee index to whole KS orbital index

//CorrToHartr  = from corr to hartree
//Hart2Corr
//HartrRange




//DF_CorrBase  = <lowenergy_projected_orbal |  low_energy_KS_band >
//KS_eigenVectors_orthoBasis  = <orthogonal basis |  KS_band >
//transformMatrix_k =  basis transformation from orthogonal to non-orthogonal basis, e.g.,  //\ket{ortho_j} = \ket{DFT_AO_i} T_{ij}
