#//Here $@ means  target i.e. $(TARGET) at this rule
#See http://developinghappiness.com/?p=28
##http://developinghappiness.com/?p=174
#see icc warning : http://stackoverflow.com/questions/27343360/how-to-turn-on-icc-icpc-warnings


TARGET_dmft=dmft
DEPFILE = dependency




CC = mpic++

##########################
###gcc
##########################
##CCFLAGS=   -std=c++11     -fstack-protector-all
##OFLAGS = -O0      -g  -ggdb    

#CCFLAGS=   -std=c++11  
#OFLAGS = -O3      

#DEFINES= -DMPI_ENABLED   \
#         -DEIGEN_INITIALIZE_MATRICES_BY_ZERO  
##


##########################
###Intel MKL
##########################

INTEL = /opt/intel/intel_2017_beta
MKL = ${INTEL}/mkl

CCFLAGS=   -std=c++11   -I${MKL}/include -ip  -inline-level=2   -inline-forceinline   -w2 # -check-pointers=rw
OFLAGS = -O3       -g  -ggdb  

#OFLAGS = -O3
#CCFLAGS=   -std=c++11   -I${MKL}/include -ip  -inline-level=2   -inline-forceinline   -w2  

DEFINES= -DMPI_ENABLED   \
         -L${MKL}/lib/intel64 -lmkl_intel_lp64  -lmkl_core -lifcore   -lmkl_sequential\
         -DEIGEN_USE_MKL_ALL -DEIGEN_INITIALIZE_MATRICES_BY_ZERO 


#gcc





#EIGEN
EIGEN_INC = -I$(CURDIR)/../../eigen-git-mirror/







SOURCE =    band.cpp   FourierTransform.cpp  tight_common.cpp Construct_hyb_delta.cpp GreenFtn.cpp  \
            dos.cpp  readInput_text.cpp mpi_interface.cpp dmft_main.cpp  ImgTimeFtn.cpp   SOLVER.cpp  \
            TB.cpp     diag_EIGEN.cpp          ConstructModelHamiltonian.cpp  ImS_stoch_to_ReS.cpp \
            upfolding.cpp  qs_spectra.cpp gen_Uijkl.cpp    dmft_common.cpp pulay.cpp   \
            SCGF2_fullLoop.cpp  ctqmc_rutgers.cpp  read_TB_parm.cpp downfolding.cpp  get_weight.cpp SpreadFtn.cpp \
            read_scfout.cpp


OBJS_dmft = $(SOURCE:%.cpp=%.o) 

.SUFFIXES : .cpp .o
.cpp.o: 
	$(CC)  -c  $(OFLAGS)    $(CCFLAGS)   $(DEFINES)    $(EIGEN_INC)    $<


$(TARGET_dmft) : $(OBJS_dmft)  
	$(CC) $(CCFLAGS) $(OFLAGS)  -fPIC    -o $@ $(OBJS_dmft)    $(DEFINES)     
	cp $(TARGET_dmft) ../example
all : style dep $(TARGET_dmft)
	cp $(TARGET_dmft)  ../example
dep : 
	$(CC) $(CCFLAGS) -MM $(EIGEN_INC)    $(SOURCE) >  $(DEPFILE)
style : 
	astyle $(SOURCE)
clean : 
	rm $(TARGET_dmft) $(OBJS_dmft)  *.cpp.orig
-include $(DEPFILE)
