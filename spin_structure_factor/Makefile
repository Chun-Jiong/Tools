
FC=gfortran

F90=\
	SSF.f90

EXEC=SSF

all: $(F90)
	$(FC) -O1 -cpp -o $(EXEC) $(F90)

clean:
	-rm -rf $(EXEC) *.mod
