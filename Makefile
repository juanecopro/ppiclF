FC               = mpiifort # Valid MPI Fortran compiler
CC               = mpiicc  # Valid MPI C compiler
#FFLAGS          += -O2 -g -cpp
#CPFLAGS         += -O2 -g -E
#FFLAGS          += -qoffload-option,mic,compiler, "-O3" -O3 -fpe0 -axCORE-AVX2 -xSSE4.2 -cpp
#CPFLAGS         += -qoffload-option,mic,compiler, "-O3" -O3 -axCORE-AVX2 -xSSE4.2 -E
#OPTFLAGS        = -xmic-avx512
#OPTFLAGS        = -xCORE-AVX2
#OPTFLAGS        = -axMIC-AVX512,CORE-AVX2
OPTFLAGS        = -g -traceback
FFLAGS          += $(OPTFLAGS) -cpp
CPFLAGS         += $(OPTFLAGS) -E
FFLAGS          += -I${HOME}/Zoltan/build/include 
 	           
INSTALL_LOCATION = .

####################
# DO NOT TOUCH BELOW
####################

# PPICLF LIBRARY
SOURCE_ROOT_PPICLF=$(INSTALL_LOCATION)/source
PPICLF_IFLAGS+=-I$(SOURCE_ROOT_PPICLF)

# GSLIB LIBRARY
SOURCE_ROOT_GSLIB=$(INSTALL_LOCATION)/3rd_party/gslib
SOURCE_ROOT_GSLIB_OBJ=$(SOURCE_ROOT_GSLIB)/gslib/src
GSLIB_IFLAGS+=-I$(SOURCE_ROOT_GSLIB)/include

SRC = $(SOURCE_ROOT_PPICLF)/ppiclf_user.f     \
      $(SOURCE_ROOT_PPICLF)/ppiclf_comm.f     \
      $(SOURCE_ROOT_PPICLF)/ppiclf_op.f       \
      $(SOURCE_ROOT_PPICLF)/ppiclf_io.f       \
      $(SOURCE_ROOT_PPICLF)/ppiclf_solve.f    \
      $(SOURCE_ROOT_PPICLF)/ppiclf_mxm.f      
OBJ = $(SOURCE_ROOT_PPICLF)/ppiclf_user.o     \
      $(SOURCE_ROOT_PPICLF)/ppiclf_comm.o     \
      $(SOURCE_ROOT_PPICLF)/ppiclf_op.o       \
      $(SOURCE_ROOT_PPICLF)/ppiclf_io.o       \
      $(SOURCE_ROOT_PPICLF)/ppiclf_solve.o    \
      $(SOURCE_ROOT_PPICLF)/ppiclf_mxm.o      
PRE = $(SOURCE_ROOT_PPICLF)/PPICLF_SOLN.h     \
      $(SOURCE_ROOT_PPICLF)/PPICLF_GRID.h     \
      $(SOURCE_ROOT_PPICLF)/PPICLF_OPT.h      \
      $(SOURCE_ROOT_PPICLF)/PPICLF_PARALLEL.h \
      $(SOURCE_ROOT_PPICLF)/PPICLF_GEOM.h       
POS = $(SOURCE_ROOT_PPICLF)/PPICLF
      

# Make commands
default: makeThird preProcess zoltanmod getObjs libObjs

zoltanmod: $(SOURCE_ROOT_PPICLF)/zoltanRCB.f90
	$(FC) -O2 $(FFLAGS) -c $(SOURCE_ROOT_PPICLF)/zoltanRCB.f90

libObjs: $(SOURCE_ROOT_PPICLF)/ppiclf.o $(SOURCE_ROOT_PPICLF)/zoltanRCB.o
	@ar crv $(SOURCE_ROOT_PPICLF)/libppiclF.a $(SOURCE_ROOT_PPICLF)/*.o $(SOURCE_ROOT_GSLIB_OBJ)/*.o
	@echo "                       "
	@echo "***********************"
	@echo "*** LIBRARY SUCCESS ***"
	@echo "***********************"
	@echo "                       "

getObjs: $(SOURCE_ROOT_PPICLF)/ppiclf.f
	$(FC) -O2 $(FFLAGS) -c $(SOURCE_ROOT_PPICLF)/ppiclf.f $(GSLIB_IFLAGS) 
	mv *.o $(SOURCE_ROOT_PPICLF)
	@echo "                              "
	@echo "******************************"
	@echo "*** LIBRARY OBJECT SUCCESS ***"
	@echo "******************************"
	@echo "                              "

preProcess: $(PRE)
	$(CC) $(CPFLAGS) $(PRE) > $(POS)_tmp
	grep "^[^#]" $(POS)_tmp > $(POS)
	rm $(POS)_tmp
	echo '#include "PPICLF_USER.h"' > $(SOURCE_ROOT_PPICLF)/ppiclf.f
	echo '#include "PPICLF_STD.h"' >> $(SOURCE_ROOT_PPICLF)/ppiclf.f
	cat $(SRC)  >> $(SOURCE_ROOT_PPICLF)/ppiclf.f
	@echo "                       "
	@echo "********************"
	@echo "*** PREPROCESSED ***"
	@echo "********************"
	@echo "                       "
	

makeThird: $(SOURCE_ROOT_GSLIB)/install
	cd $(SOURCE_ROOT_GSLIB); \
	CFLAGS="$(OPTFLAGS)" ./install $(CC) $(FC)
	@echo "                       "
	@echo "***********************"
	@echo "*** INSTALLED GSLIB ***"
	@echo "***********************"
	@echo "                       "

clean:
	 rm -r $(SOURCE_ROOT_PPICLF)/*.o                        
	 rm -r $(SOURCE_ROOT_PPICLF)/ppiclf.f
	 rm -r $(SOURCE_ROOT_PPICLF)/PPICLF
	 rm -r $(SOURCE_ROOT_PPICLF)/libppiclF.a                
	 rm -r $(SOURCE_ROOT_GSLIB)/gslib                       
	 rm -r $(SOURCE_ROOT_GSLIB)/lib                         
	 rm -r $(SOURCE_ROOT_GSLIB)/include                     
	 rm -r $(SOURCE_ROOT_GSLIB)/*.tar.gz                    
	 rm -r $(INSTALL_LOCATION)/short_tests/*.pyc            
	 rm -r $(INSTALL_LOCATION)/short_tests/lib/*.pyc        
	 rm -r $(INSTALL_LOCATION)/short_tests/test1/test1.log*
	 rm -r $(INSTALL_LOCATION)/short_tests/test1/*.vtu
	 cd $(INSTALL_LOCATION)/short_tests/test1/; \
	make clean

