FCC=ifort
FCC=gfortran
FC=-O1 -cpp

#F90=entropy_E.f90
#EXEC=entropy_E
F90=entropy_Cv.f90
EXEC=entropy_Cv

all:
	$(FCC) $(FC) $(F90) -o $(EXEC)
	-rm -f *.mod
