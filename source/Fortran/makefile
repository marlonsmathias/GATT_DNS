FCOMP = mpifort
ARGS = -ffree-line-length-0

SOURCEFOLDER = ../../source/Fortran

include makefile_extra

include $(DECOMPDIR)/src/Makefile.inc

INCLUDE = -I$(MATROOT)/extern/include -I$(SOURCEFOLDER) -I$(DECOMPDIR)/include
LINKS = -L$(MATROOT)/bin/glnxa64 -lmat -lmx -Wl,-rpath=$(MATROOT)/bin/glnxa64/ -L$(DECOMPDIR)/lib -l2decomp_fft
OBJECTS = readWriteMat.o derivs.o nseq.o timeIntegration.o


main: disturbances.o $(OBJECTS)
	$(FCOMP) $(ARGS) -J. $(OBJECTS) disturbances.o $(SOURCEFOLDER)/main.F90 $(LINKS) -o main $(INCLUDE)
	
disturbances.o:
	$(FCOMP) $(ARGS) -J. -c disturbances.F90 -o disturbances.o
	
%.o: $(SOURCEFOLDER)/%.F90
	$(FCOMP) $(ARGS) -c -J. $< $(INCLUDE) -o $@