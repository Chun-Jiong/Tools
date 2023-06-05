
FCC=gfortran

F90=\
	SSF.f90

EXEC=SSF
FC=-O1 -cpp
#FC=-O0 -g -cpp -DGFORT -fbounds-check -Wuninitialized

all: $(F90)
	$(FCC) $(FC) -o $(EXEC) $(F90)

clean:
	-rm -rf $(EXEC) *.mod
