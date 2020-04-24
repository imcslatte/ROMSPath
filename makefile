#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#                                                                       :::
#                          ROMSPath Makefile                              :::
#                                                                       :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


#==========================================================================
#                         USER-DEFINED OPTIONS                            =
#==========================================================================

#CPPFLAGS  = -DGROWTH -DWETDRY  -DSTOKES
#CPPFLAGS  =  -DWETDRY

#------------------------------------------------
#    Set compiler and flags
#------------------------------------------------
#
#    Turn one of the following on:
IFORT:= 
GFORTRAN :=    
PGI  := 

ifdef IFORT 
  FC = ifort
#  NETCDF_INCDIR = /home/rjdave/local/include
#  NETCDF_LIBDIR = /home/rjdave/local/lib 
  CPPFLAGS  += -DIFORT
  FFLAGS = -fp-model strict -mcmodel=medium  -O3 -fpp  $(CPPFLAGS) -I$(NETCDF_INCDIR)
#   FFLAGS = -g -O0 -traceback -check all -check bounds  -cpp  $(CPPFLAGS) -I$(NETCDF_INCDIR)
endif

ifdef GFORTRAN
  FC = gfortran
#  NETCDF_INCDIR = /home/rjdave/local/include
#  NETCDF_LIBDIR = /home/rjdave/local/lib
  CPPFLAGS  += -DGFORTRAN
  FFLAGS = -march=k8 -ffast-math -fno-cx-limited-range -O3 -funroll-loops --param max-unroll-times=4 -ffree-line-length-none -cpp  $(CPPFLAGS) -I$(NETCDF_INCDIR)
#   FFLAGS = -g -O0   -ffree-line-length-none -cpp  $(CPPFLAGS) -I$(NETCDF_INCDIR)
endif

ifdef PGI
  FC = pgf90
  CPPFLAGS  += -DPGI
#  NETCDF_INCDIR = /home/rjdave/local/include
#  NETCDF_LIBDIR = /home/rjdave/local/lib
  FFLAGS := -g -cpp  $(CPPFLAGS) -I$(NETCDF_INCDIR) 
endif

#------------------------------------------------
#    Set NetCDF Library Locations.
#    If NetCDF was compiled with HDF5, set:
#        HDF5 := on
#    Otherwise, leave blank:
#        HDF5 :=
#------------------------------------------------

HDF5 := on
NFCONFIG := on
#==========================================================================
# End of user-defined options. Nothing should be changed below this point =
#==========================================================================

OBJS          = parameter_module.o grid_module.o random_module.o   \
				hydrodynamic_module.o interpolation_module.o boundary_module.o\
				pdf_module.o hor_turb_module.o	pppack.o ver_turb_module.o advection_module.o\
				growth_module.o behavior_module.o

ifdef NFCONFIG        
	NF_CONFIG ?= nf-config
#	ifdef GFORTRAN
#		NF_CONFIG ?= /opt/sw/apps/gcc-7.3.0/netcdf/4.6.1/bin/nf-config
#	endif
#	ifdef IFORT
#		NF_CONFIG ?= /opt/sw/apps/intel-18.0.1/netcdf/4.6.1/bin/nf-config
#	endif
    NETCDF_INCDIR ?= $(shell $(NF_CONFIG) --prefix)/include
    LIBS := $(shell $(NF_CONFIG) --flibs)
else
	ifdef HDF5
		ifdef PGI_USGS
			LIBS      = -L$(NETCDF_LIBDIR) -lnetcdf -lnetcdff -L/share/apps/hdf5/lib -lhdf5_hl -lhdf5 -lz -lm -L/share/apps/szip/lib -lsz -lcurl
		else
			LIBS      = -L$(NETCDF_LIBDIR) -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lcurl -lz
		endif
	else
		LIBS      = -L$(NETCDF_LIBDIR) -lnetcdff -lnetcdf
	endif

endif

ROMSPath : $(OBJS)

	@echo "  Compiling ROMSPath.f90 "
	@echo " Using $(FC) "
	@echo " $(FC) $(FFLAGS) -o ROMSPath.exe ROMSPath.f90 $(OBJS) $(LIBS)  -save-temps "
	@$(FC) $(FFLAGS) -o ROMSPath.exe ROMSPath.f90 $(OBJS) $(LIBS)  -save-temps
	@echo "  "
	@echo "  Compilation Successfully Completed"
	@echo "  "

%.o: %.f90
	@echo "  Compiling $<"
	@$(FC) $(FFLAGS) -cpp  $(CPPFLAGS) -save-temps -c $<

clean:
	\rm *.o *.mod *.i90 *.s *.exe

