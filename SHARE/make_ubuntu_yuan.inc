#######################################################################
#            Define Basic Utilities
#######################################################################
  SHELL = /bin/sh
  PWD1 = `pwd`
  MYHOME = $(HOME)/code/bin
  PRECOMP:= /lib/cpp -traditional -DLINUX
  COMPILE = $(FC)
  COMPILE_FREE = $(FC) -ffree-form -ffree-line-length-none -ffixed-line-length-none
  LINK    = ld $(FLAGS) -o
  LINK_AR = ar -ruv
  LINK_C  = $(CC) -shared -Wl,-z-defs

#######################################################################
#            Define Compiler Flags
#######################################################################
  FLAGS_R = -O2 -fexternal-blas -fbacktrace -fcheck=all,no-array-temps \
			-Wno-argument-mismatch
  FLAGS_D = -g -O0 -fexternal-blas -fbacktrace -fcheck=all \
			-Wno-argument-mismatch
			
  LIBS    = -L/usr/lib -lblas -llapack -lscalapack-openmpi
  
# Need to allow argument mismatch for gfortran version 10 and above
  GFORTRAN_VERSION := $(shell $(FC) -dumpversion)
  ifeq "$(shell [ $(GFORTRAN_VERSION) -gt 9 ] && echo true)" "true"
	 FLAGS_R += -fallow-argument-mismatch
     FLAGS_D += -fallow-argument-mismatch
  endif

#######################################################################
#            MPI Options
#######################################################################
  LMPI    = T
  MPI_COMPILE = mpif90.openmpi
  MPI_COMPILE_FREE = mpif90 -ffree-form \
                     -ffree-line-length-none -ffixed-line-length-none
  MPI_COMPILE_C = mpicc.openmpi
  MPI_LINK = mpif90.openmpi -shared -fPIC
  MPI_RUN = mpiexec
  MPI_RUN_OPTS = -np 4
#  MPI_RUN_OPTS_SM = -np 16
#  MPI_RUN_OPTS_MD = -np 64
#  MPI_RUN_OPTS_LG = -np 256

#######################################################################
#            NAG Options
#######################################################################
  LNAG = F
  NAG_LIB = -L$(NAG_ROOT)/lib -lnag_nag

#######################################################################
#            NETCDF Options
#######################################################################
  LNETCDF = T
  NETCDF_INC = $(shell nc-config --fflags)
  NETCDF_LIB = $(shell nc-config --flibs)

#######################################################################
#            NTCC Options
#######################################################################
  LNTCC = F
  NTCC_INC = -I/usr/include
  NTCC_LIB = -L/usr/lib -laladdinsub -lr8slatec -ladpak\
             -lcppsub -lcomput -lpspline -lportlib -lezcdf -lmds_sub \
             -lmdstransp -lvaxonly

#######################################################################
#            HDF5 Options
#######################################################################
  LHDF5 = T
  HDF5_INC = -I/usr/include/hdf5/openmpi
  HDF5_LIB = -L/usr/lib/x86_64-linux-gnu/hdf5/openmpi \
             -lhdf5 -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran

#######################################################################
#             PGPLOT Options
#######################################################################
  LPGPLOT = F
  PGPLOT_INC = -I$(PGPLOT_DIR)
  PGPLOT_LIB = -L$(PGPLOT_DIR) -lpgplot -L/usr/lib/x86_64-linux-gnu -lX11

#######################################################################
#             SILO Options
#######################################################################
  LSILO = F
  SILO_INC = -I$(SILOHOME)/include
  SILO_LIB = -L$(SILOHOME)/lib/x86_64-linux-gnu -lsiloh5

#######################################################################
#            FFTW3 Options
#######################################################################
  LFFTW3 = F
  FFTW3_INC = 
  FFTW3_LIB = 

#######################################################################
#            DKES/NEO Options
#######################################################################
  LDKES = T
  LNEO  = T

#######################################################################
#            GENE Options
#######################################################################
  LGENE = F
  GENE_INC = -I$(GENE_PATH)
  GENE_DIR = $(GENE_PATH)
  LIB_GENE = libgene.a
  GENE_LIB = $(GENE_DIR)/$(LIB_GENE) \
             -L/u/slazerso/src/GENE17_2016/external/pppl_cluster/futils/src -lfutils \
             -L$(FFTWHOME)/lib -lfftw3 \
             -L$(SLEPC_DIR)/$(PETSC_ARCH)/lib -lslepc \
             -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -lpetsc -lX11

#######################################################################
#            COILOPT++ Options
#######################################################################
  LCOILOPT = F
  COILOPT_INC = -I$(COILOPT_PATH)
  COILOPTPP_DIR = $(COILOPT_PATH)
  LIB_COILOPTPP = libcoilopt++.a
  COILOPT_LIB = $(COILOPT_PATH)/$(LIB_COILOPTPP) \
                -L$(GSLHOME)/lib/x86_64-linux-gnu -lgsl -lgslcblas -lstdc++ -lmpi_cxx

#######################################################################
#            REGCOIL Options
#######################################################################
  LREGCOIL= F
  REGCOIL_DIR = $(REGCOIL_PATH)
  REGCOIL_INC = -I$(REGCOIL_DIR) 
  LIB_REGCOIL = libregcoil.a
  REGCOIL_LIB = $(REGCOIL_DIR)/$(LIB_REGCOIL) -fopenmp

#######################################################################
#            TERPSICHORE Options
#######################################################################
  LTERPSICHORE= F
  TERPSICHORE_INC = -I$(TERPSICHORE_PATH)
  TERPSICHORE_DIR = $(TERPSICHORE_PATH)
  LIB_TERPSICHORE = libterpsichore.a
  TERPSICHORE_LIB = $(TERPSICHORE_DIR)/$(LIB_TERPSICHORE)

#######################################################################
#            TRAVIS Options
#######################################################################
  LTRAVIS= F
  TRAVIS_DIR = $(TRAVIS_PATH)
  LIB_TRAVIS = libtravis64_sopt.a
  LIB_MCONF  = libmconf64.a
  TRAVIS_LIB = $(TRAVIS_DIR)/lib/$(LIB_TRAVIS) \
               $(TRAVIS_DIR)/mag_conf/lib/$(LIB_MCONF) -lstdc++

#######################################################################
#            Available Energy Options
#######################################################################
  LAEOPT= F
  AEOPT_DIR = $(AEOPT_PATH)
  AEOPT_INC = -I$(AEOPT_DIR) 
  LIB_AEOPT = libtrapAE.a
  AEOPT_LIB = $(AEOPT_PATH)/$(LIB_AEOPT)

#######################################################################
#            MANGO Options
#######################################################################
  LMANGO = T
  MANGO_DIR = /home/yuan/code/mango
  LIB_MANGO = libmango.a 
	
# If LMANGO=T, any flags specified in the next 2 lines will be added to 
# MANGO_F_COMPILE_FLAGS and MANGO_F_LINK_FLAGS respectively
# as defined in $(MANGO_DIR)/lib/mangoMakeVariables

  MANGO_INC = -I$(MANGO_DIR)/include
  MANGO_LIB = $(MANGO_DIR)/lib/$(LIB_MANGO) -fopenmp

#######################################################################
#            System-specific path and binary definitions
#######################################################################
  FC = gfortran
  BLASHOME = /usr/lib/x86_64-linux-gnu
  BLACS_HOME = /usr/lib/x86_64-linux-gnu
  SCALAPACK_HOME = /usr/lib/x86_64-linux-gnu
  MPIHOME = /usr/lib/x86_64-linux-gnu
  NETCDF_HOME = /usr/lib/x86_64-linux-gnu
  HDF5_HOME = /usr/lib/x86_64-linux-gnu
  PGPLOT_DIR = /usr/lib/x86_64-linux-gnu
  SILOHOME = /usr/lib/x86_64-linux-gnu

#######################################################################
#            LIBSTELL Shared Options
#######################################################################
LIB_SHARE = -lc -lgfortran -lstdc++ -lmpi -lmpi_mpifh -lz -lc -lm -lpthread $(LIBS) -lc

# LIB_SHARE = $(BLASHOME)/lib/libblas.so \
          $(SCALAPACK_HOME)/lib/libscalapack-openmpi.so \
          $(BLACS_HOME)/lib/libblacs-openmpi.so  $(BLACS_HOME)/lib/libblacsCinit-openmpi.so $(BLACS_HOME)/lib/libblacsF77init-openmpi.so \
          $(HDF5_HOME)/lib/x86_64-linux-gnu/libhdf5_hl.so $(HDF5_HOME)/lib/x86_64-linux-gnu/libhdf5_fortran.so $(HDF5_HOME)/lib/x86_64-linux-gnu/libhdf5hl_fortran.so $(HDF5_HOME)/lib/x86_64-linux-gnu/libhdf5.so \
          $(HDF5_HOME)/lib/x86_64-linux-gnu/libmpi_usempif08.so $(HDF5_HOME)/lib/x86_64-linux-gnu/libhdf5_openmpi_fortran.so \
          $(NETCDF_HOME)/lib/x86_64-linux-gnu/libnetcdf.so $(NETCDF_HOME)/lib/x86_64-linux-gnu/libnetcdff.so $(NETCDF_HOME)/lib/x86_64-linux-gnu/libnetcdf_c++.so \
          $(SILOHOME)/lib/x86_64-linux-gnu/libsiloh5.so \
          $(GSLHOME)/lib/x86_64-linux-gnu/libgsl.so \
          $(GCC6_HOME)/libgfortran.so $(GCC6_HOME)/libstdc++.so \
          $(MPIHOME)/lib/x86_64-linux-gnu/libmpi.so $(MPIHOME)/lib/x86_64-linux-gnu/libmpi_mpifh.so \
          /usr/lib/x86_64-linux-gnu/libm.so /usr/lib/liblapack.so /usr/lib/x86_64-linux-gnu/libdl.so
# LIB_SHARE = $(BLASHOME)/lib/libblas.so \
          $(SCALAPACK_HOME)/lib/libscalapack-openmpi.so.1 \
          $(BLACS_HOME)/lib/libblacs-openmpi.so.1 \
          $(BLACS_HOME)/lib/libblacsCinit-openmpi.so.1 \
          $(BLACS_HOME)/lib/libblacsF77init-openmpi.so.1 \
          $(HDF5_HOME)/lib/x86_64-linux-gnu/libhdf5_openmpi_fortran.so \
          $(NETCDF_HOME)/lib/libnetcdf.so \
          $(NETCDF_HOME)/lib/libnetcdff.so \
          $(NETCDF_HOME)/lib/libnetcdf_c++.so \
          $(SILOHOME)/lib/x86_64-linux-gnu/libsiloh5.so \
          $(GSLHOME)/lib/x86_64-linux-gnu/libgsl.so \
          $(GCC49_HOME)/libgfortran.so \
          $(GCC49_HOME)/libstdc++.so \
          $(MPIHOME)/lib/libmpi.so \
          $(MPIHOME)/lib/libmpif77.so \
          /usr/lib/x86_64-linux-gnu/libm.so \
          /usr/lib/liblapack.so \
          /usr/lib/x86_64-linux-gnu/libdl.so

#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#   This is universal to all build systems, do not edit unless your
#   adding a NEW code or library to the build.
#   New means NEW not new version or new system.
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
# Some Initial stuff
# Home bin specifies where the build scripts are located.
HOME_BIN = ../../
MOD1_PATH = -I$(MYHOME)/libstell_dir
MOD_PATH = 
CLEAN_ARCH = 
LIB_DIR = ../../LIBSTELL
LIB = libstell.a

# Handle Which Code is being CALLED
ifeq ($(MAKECMDGOALS),xstelloptv2)
  PRECOMP += -DBNORM -DBEAMS3D_OPT -DJINV_OPT -DTXPORT_OPT -DDIAGNO_OPT -DSKS -DSKS2
  VMEC_DIR ?= ../../VMEC2000
  LIB_VMEC= libvmec.a
  MOD_PATH+= -I$(VMEC_DIR)/$(LOCTYPE)
  BEAMS3D_DIR = ../../BEAMS3D
  LIB_BEAMS3D = libbeams3d.a
  MOD_PATH+= -I$(BEAMS3D_DIR)/$(LOCTYPE)
  BOOTSJ_DIR = ../../BOOTSJ
  LIB_BOOTSJ = libbootsj.a
  MOD_PATH+= -I$(BOOTSJ_DIR)/$(LOCTYPE)
  BNORM_DIR = ../../BNORM
  LIB_BNORM = libbnorm.a
  MOD_PATH+= -I$(BNORM_DIR)/$(LOCTYPE)
  BOOZ_DIR = ../../BOOZ_XFORM
  LIB_BOOZ= libbooz.a
  MOD_PATH+= -I$(BOOZ_DIR)/$(LOCTYPE)
  COBRA_DIR = ../../COBRAVMEC
  LIB_COBRA = libcobra.a
  MOD_PATH+= -I$(COBRA_DIR)/$(LOCTYPE)
  DIAGNO_DIR = ../../DIAGNO
  LIB_DIAGNO= libdiagno.a
  MOD_PATH+= -I$(DIAGNO_DIR)/$(LOCTYPE)
  FIELDLINES_DIR = ../../FIELDLINES
  LIB_FIELDLINES = libfieldlines.a
  MOD_PATH+= -I$(FIELDLINES_DIR)/$(LOCTYPE)
  JINV_DIR = ../../J_INVARIANT
  LIB_JINV = libjinv.a
  MOD_PATH+= -I$(JINV_DIR)/$(LOCTYPE)
  MGRID_DIR= ../../MAKEGRID
  LIB_MGRID= libmakegrid.a
  MOD_PATH+= -I$(MGRID_DIR)/$(LOCTYPE)
  LIB_LINK=  $(VMEC_DIR)/$(LOCTYPE)/$(LIB_VMEC)  \
	$(BEAMS3D_DIR)/$(LOCTYPE)/$(LIB_BEAMS3D) \
	$(BOOTSJ_DIR)/$(LOCTYPE)/$(LIB_BOOTSJ) \
	$(BNORM_DIR)/$(LOCTYPE)/$(LIB_BNORM) \
	$(BOOZ_DIR)/$(LOCTYPE)/$(LIB_BOOZ) \
	$(COBRA_DIR)/$(LOCTYPE)/$(LIB_COBRA) \
	$(DIAGNO_DIR)/$(LOCTYPE)/$(LIB_DIAGNO) \
	$(FIELDLINES_DIR)/$(LOCTYPE)/$(LIB_FIELDLINES) \
	$(JINV_DIR)/$(LOCTYPE)/$(LIB_JINV) \
	$(MGRID_DIR)/$(LOCTYPE)/$(LIB_MGRID) 
  ifeq ($(LDKES),T)
    DKES_DIR = ../../DKES
    LIB_DKES = libdkes.a
    MOD_PATH += -I$(DKES_DIR)/$(LOCTYPE)
    PRECOMP += -DDKES_OPT
    LIB_LINK += $(DKES_DIR)/$(LOCTYPE)/$(LIB_DKES)
  else
    DKES_DIR = 
    LIB_DKES = 
  endif
  ifeq ($(LNEO),T)
    NEO_DIR = ../../NEO
    LIB_NEO = libneo.a
    MOD_PATH += -I$(NEO_DIR)/$(LOCTYPE)
    PRECOMP += -DNEO_OPT
    LIB_LINK += $(NEO_DIR)/$(LOCTYPE)/$(LIB_NEO)
  else
    NEO_DIR = 
    LIB_NEO = 
  endif
else ifeq ($(MAKECMDGOALS),xthrift)
  PRECOMP += -DSKS  -DSKS2
  VMEC_DIR = ../../VMEC2000
  LIB_VMEC= libvmec.a
  MOD_PATH+= -I$(VMEC_DIR)/$(LOCTYPE)
  BOOTSJ_DIR = ../../BOOTSJ
  LIB_BOOTSJ = libbootsj.a
  MOD_PATH+= -I$(BOOTSJ_DIR)/$(LOCTYPE)
  BOOZ_DIR = ../../BOOZ_XFORM
  LIB_BOOZ= libbooz.a
  MOD_PATH+= -I$(BOOZ_DIR)/$(LOCTYPE)
  DIAGNO_DIR = ../../DIAGNO
  LIB_DIAGNO= libdiagno.a
  MOD_PATH+= -I$(DIAGNO_DIR)/$(LOCTYPE)
  LIB_LINK=  $(VMEC_DIR)/$(LOCTYPE)/$(LIB_VMEC)  \
  $(BOOTSJ_DIR)/$(LOCTYPE)/$(LIB_BOOTSJ) \
  $(BOOZ_DIR)/$(LOCTYPE)/$(LIB_BOOZ) \
  $(DIAGNO_DIR)/$(LOCTYPE)/$(LIB_DIAGNO)
  LGENE = F
  LCOILOPT = F
  LSFINCS = F
  LREGCOIL= F
  LMANGO = F
  LTERPSICHORE = F
else ifeq ($(MAKECMDGOALS),help)
else ifeq ($(MAKECMDGOALS),test_stellopt)
else ifeq ($(MAKECMDGOALS),test_stellopt_coilopt)
else ifeq ($(MAKECMDGOALS),test_stellopt_gene)
else ifeq ($(MAKECMDGOALS),test_stellopt_regcoil)
else ifeq ($(MAKECMDGOALS),test_stellopt_sfincs)
else ifeq ($(MAKECMDGOALS),test_stellopt_terpsichore)
else ifeq ($(MAKECMDGOALS),test_stellopt_travis)
else
  LGENE = F
  LCOILOPT = F
  LSFINCS = F
  LREGCOIL= F
  LMANGO = F
  LTERPSICHORE = F
  LTRAVIS = F
			 
endif
ifeq ($(MAKECMDGOALS),xanimec)
  PRECOMP += -D_ANIMEC
endif
ifeq ($(MAKECMDGOALS),xparvmec)
  PRECOMP += -DSKS  -DSKS2
endif
ifeq ($(MAKECMDGOALS),xsiesta)
  PRECOMP += -DSKS  -DSKS2
endif
ifeq ($(MAKECMDGOALS),xgtovmi)
  ifeq ($(LPGPLOT),T)
    MOD1_PATH += $(PGPLOT_INC)
    LIB_LINK += $(PGPLOT_LIB)
  endif
endif
ifeq ($(MAKECMDGOALS),xdes_plot)
  ifeq ($(LPGPLOT),T)
    MOD1_PATH += $(PGPLOT_INC)
    LIB_LINK += $(PGPLOT_LIB)
  endif
endif

# Other Codes
ifeq ($(LGENE),T)
  MOD1_PATH += $(GENE_INC)
  PRECOMP += -DGENE
  MOD_PATH += -I$(GENE_DIR)
else
  GENE_LIB = 
  GENE_DIR = 
  LIB_GENE = 
endif
ifeq ($(LCOILOPT),T)
  MOD1_PATH += $(COILOPT_INC)
  PRECOMP += -DCOILOPTPP
  MOD_PATH += -I$(COILOPTPP_DIR)
else
  COILOPT_LIB = 
  COILOPTPP_DIR = 
  LIB_COILOPTPP = 
endif
ifeq ($(LREGCOIL),T)
  MOD1_PATH += $(REGCOIL_INC)
  PRECOMP += -DREGCOIL
  MOD_PATH += -I$(REGCOIL_DIR)
else
  REGCOIL_LIB = 
  REGCOIL_DIR = 
  LIB_REGCOIL = 
endif
ifeq ($(LTERPSICHORE),T)
  MOD1_PATH += $(TERPSICHORE_INC)
  PRECOMP += -DTERPSICHORE
  MOD_PATH += -I$(TERPSICHORE_DIR)
else
  TERPSICHORE_LIB = 
  TERPSICHORE_DIR = 
  LIB_TERPSICHORE = 
endif
ifeq ($(LTRAVIS),T)
  PRECOMP += -DTRAVIS
else
  TRAVIS_LIB = 
  TRAVIS_DIR = 
  LIB_TRAVIS = 
endif
ifeq ($(LSFINCS),T)
  MOD1_PATH += $(SFINCS_INC)
  PRECOMP += -DSFINCS
  MOD_PATH += $(SFINCS_DIR)
else
  SFINCS_LIB =
  SFINCS_DIR =
  LIB_SFINCS =
endif
ifeq ($(LAEOPT),T)
  MOD1_PATH += $(AEOPT_INC)
  PRECOMP += -DAEOPT
  MOD_PATH += -I$(AEOPT_DIR)
else
  AEOPT_LIB = 
  AEOPT_DIR = 
  LIB_AEOPT = 
endif
ifeq ($(LMANGO),T)
  # This next line will define the variables MANGO_C_COMPILE_FLAGS, MANGO_F_COMPILE_FLAGS, MANGO_C_LINK_FLAGS, and MANGO_F_LINK_FLAGS
  include $(MANGO_DIR)/lib/mangoMakeVariables
  # Add any flags defined by the individual stellopt make_*.inc files:
  MANGO_INC += $(MANGO_F_COMPILE_FLAGS)
  MANGO_LIB += $(MANGO_F_LINK_FLAGS)
  MOD1_PATH += $(MANGO_INC)
  PRECOMP += -DMANGO
  MOD_PATH += $(MANGO_INC)
else
  MANGO_LIB =
  MANGO_DIR =
  LIB_MANGO =
endif

# Handle other codes
LIB_LINK += $(GENE_LIB) $(COILOPT_LIB) $(REGCOIL_LIB) $(SFINCS_LIB) $(TERPSICHORE_LIB) $(TRAVIS_LIB) $(AEOPT_LIB) $(MANGO_LIB) 

# Add libraries to LIB_LINK
# LIB_LINK+= $(MYHOME)/libstell.a $(LIBS)
# This version avoids referencing the bin directory for the symbolic link to lib.
# Which we can't reference since we now have the awk file in the main directory.
LIB_LINK+= $(LIB_DIR)/$(LOCTYPE)/$(LIB) $(LIBS)

# 3rd Party Libs
ifeq ($(LMPI),T)
  PRECOMP += -DMPI_OPT
  COMPILE = $(MPI_COMPILE)
  COMPILE_FREE = $(MPI_COMPILE_FREE)
  LINK    = $(MPI_COMPILE) $(FLAGS) -o 
  LINK_C  = $(MPI_LINK)
endif
ifeq ($(LNETCDF),T)
  MOD1_PATH += $(NETCDF_INC)
  PRECOMP += -DNETCDF
  LIB_LINK += $(NETCDF_LIB)
  MOD_PATH += $(NETCDF_INC)
  LIB_SHARE += $(NETCDF_LIB)
endif
ifeq ($(LFFTW3),T)
  MOD1_PATH += $(FFTW3_INC)
  PRECOMP += -DFFTW3
  LIB_LINK += $(FFTW3_LIB)
  MOD_PATH += $(FFTW3_INC)
  LIB_SHARE += $(FFTW3_LIB)
endif
ifeq ($(LHDF5),T)
  MOD1_PATH += $(HDF5_INC)
  PRECOMP += -DLHDF5 -DH5_USE_16_API
  LIB_LINK += $(HDF5_LIB)
  MOD_PATH += $(HDF5_INC)
  LIB_SHARE += $(HDF5_LIB)
endif
ifeq ($(LNAG),T)
  PRECOMP += -DNAG
  LIB_LINK += $(NAG_LIB)
  LIB_SHARE += $(NAG_LIB)
endif
ifeq ($(LNTCC),T)
  MOD1_PATH += $(NTCC_INC)
  PRECOMP += -DLADAS -DNTCC
  LIB_LINK += $(NTCC_LIB)
  MOD_PATH += $(NTCC_INC)
  LIB_SHARE += $(NTCC_LIB)
endif
ifeq ($(LSILO),T)
  MOD1_PATH += $(SILO_INC)
  PRECOMP += -DSILO_AVAI
  LIB_LINK += $(SILO_LIB)
  LIB_SHARE += $(SILO_LIB)
endif

# This is hear because of MPI above resetting LINK
ifeq ($(MAKECMDGOALS),libstell.a)
  LINK = $(LINK_AR)
  COMPILE += -fPIC
  COMPILE_FREE += -fPIC
  INC_PATH = 
endif
ifeq ($(MAKECMDGOALS),libstell.so)
  LINK = $(LINK_C)
  LIB_LINK = $(LIB_SHARE)
  COMPILE += -fPIC
  COMPILE_FREE += -fPIC
  INC_PATH = 
endif

# Pass Git information to code
GIT_REPO    = $(shell git ls-remote --get-url)
GIT_BRANCH  = $(shell git rev-parse --abbrev-ref HEAD)
GIT_VERSION = $(shell git describe --abbrev=4 --dirty --always --tags)
GIT_HASH    = $(shell git show -s --format=%H)
TIME        = $(shell date +"%d.%m.%Y %H:%M:%S")
PRECOMP += -DGIT_VERSION_EXT="'$(GIT_VERSION)'"
PRECOMP += -DGIT_HASH_EXT="'$(GIT_HASH)'"
PRECOMP += -DGIT_BRANCH_EXT="'$(GIT_BRANCH)'"
PRECOMP += -DGIT_REPO_EXT="'$(GIT_REPO)'"
PRECOMP += -DBUILT_ON_EXT="'$(TIME)'"

# This must come here
MOD_PATH+= -I
