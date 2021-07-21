#IGNORE:
CMPLR = gfortran
FFLAGS= -Wno-tabs -ffree-line-length-none   -Wuninitialized -fbounds-check -fdefault-real-8 
Compil = gfortran


LIBDIR = lib
SRCDIR = src
BINDIR = bin
BLDDIR = build

VPATH =  %.f90 $(SRCDIR) $(LIBDIR)

SRC = sub_standard.f90 GrainThermalFluctuations.f90  main.f90 
OBJ = sub_standard.o GrainThermalFluctuations.o main.o
LIBSRC = nrtype.f90 nrutil.f90

LIBOBJ = nrtype.o nrutil.o


run.x:  $(LIBOBJ:%.o=$(BLDDIR)/%.o) $(OBJ:%.o=$(BLDDIR)/%.o) 
	$(CMPLR) $(FFLAGS)  $^ -o $(BINDIR)/$@ ;\


$(BLDDIR)/%.o: 	%.f90 | $(BLDDIR)
	mkdir -p $(BLDDIR)
	mkdir -p $(BINDIR)
	$(CMPLR) $(FFLAGS) -c $< -o $@ -J$(BLDDIR)

clean:  
	rm -f $(OBJ)
	rm -f *.mod
run:    
	cd ./bin; time ./run.x 
