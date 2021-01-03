TARGET_dmft=dmft
DEPFILE = dependency

##################
##osx
#################################
#
##CC=/Users/jae-hoonsim/local/bin/mpic++
##CCFLAGS=   -stdlib=libc++
#
#
#CCFLAGS=   -std=c++11  
#
#
#OFLAGS = -O3      # -g  -ggdb  
#
#DEFINES= -DMPI_ENABLED   \
#         -DEIGEN_INITIALIZE_MATRICES_BY_ZERO  
#

#
#
#
#########################
##gcc
#########################

#
#CC=mpicxx
##OFLAGS = -O0 -std=c++11       -g  -ggdb   -fstack-protector-all
#CCFLAGS=  -O3 -std=c++11     #-fstack-protector-all
#OFLAGS =      -lblas -llapack 
#
#DEFINES= -DMPI_ENABLED   \
#         -DEIGEN_INITIALIZE_MATRICES_BY_ZERO    # -DEIGEN_USE_LAPACKE    -DEIGEN_USE_BLAS

##########################
### gnu+mkl
###########################

CC = mpicxx
MKL = ${MKLROOT}

CCFLAGS=  -O3   -std=c++11   -I${MKL}/include       -finline-functions
OFLAGS = -O3

#CCFLAGS=  -O0   -std=c++11   -I${MKL}/include       -finline-functions   -g  -ggdb  
#OFLAGS = -O0       -g  -ggdb  

DEFINES= -DMPI_ENABLED   \
         -L${MKL}/lib/intel64 -lmkl_intel_lp64  -lmkl_core -lifcore   -lmkl_sequential\
         -DEIGEN_USE_MKL_ALL -DEIGEN_INITIALIZE_MATRICES_BY_ZERO 




###########################
###Intel MKL
###########################
#
#CC = mpiicpc   #-I/usr/local/include 
#MKL = ${MKLROOT}
#CCFLAGS= -O3    -std=c++11   -I${MKL}/include -ip  -inline-level=2   -inline-forceinline   -w2  #-check-pointers=rw
##CCFLAGS= -O0    -std=c++11   -I${MKL}/include   #-check-pointers=rw
#OFLAGS =  ""
##OFLAGS =  " -fPIC"
#DEFINES= -L${MKL}/lib/intel64  -mkl  \
#         -DEIGEN_USE_MKL_ALL -DEIGEN_INITIALIZE_MATRICES_BY_ZERO  # -DMPI_ENABLED   \
#
#

#######################
#EIGEN
#######################
#EIGEN_INC = -I$(CURDIR)/../../EIGEN3.3.7/
EIGEN_INC = -I/home/jaehoon.sim/src_nobkup/EIGEN3.3.7

SOURCE =    band.cpp   FourierTransform.cpp  tight_common.cpp Construct_hyb_delta.cpp GreenFtn.cpp  \
            dos.cpp  readInput_text.cpp mpi_interface.cpp dmft_main.cpp  ImgTimeFtn.cpp   SOLVER.cpp  \
            TB.cpp               ConstructModelHamiltonian.cpp  ImS_stoch_to_ReS.cpp \
            upfolding.cpp  qs_spectra.cpp gen_Uijkl.cpp    dmft_common.cpp pulay.cpp   \
            SCGF2_fullLoop.cpp  ctqmc_rutgers.cpp  read_TB_parm.cpp downfolding.cpp  get_weight.cpp SpreadFtn.cpp \
            read_scfout.cpp


OBJS_dmft = $(SOURCE:%.cpp=%.o) 

.SUFFIXES : .cpp .o
.cpp.o: 
	$(CC)  -c     $(CCFLAGS)   $(DEFINES)    $(EIGEN_INC)    $<


$(TARGET_dmft) : $(OBJS_dmft)  
	$(CC) $(CCFLAGS) $(OFLAGS)   -o $@ $(OBJS_dmft)    $(DEFINES)     
	cp $(TARGET_dmft) ../example
all : style dep $(TARGET_dmft)
	cp $(TARGET_dmft)  ../example
dep : 
	$(CC) $(CCFLAGS) -MM $(EIGEN_INC)    $(SOURCE) >  $(DEPFILE)
style : 
	astyle $(SOURCE)
clean : 
	rm $(TARGET_dmft) $(OBJS_dmft)  *.cpp.orig  $(DEPFILE)
-include $(DEPFILE)
