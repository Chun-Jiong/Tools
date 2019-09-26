# Note: we can compare the times with option -O2 and -O3 and if they are huge
# different that means the code can be optimized very much
#
FCC=gfortran
#FCC=ifort

FC=-O0 -fbounds-check -Wuninitialized
#FC=

SRC=DynamicalExponent.f90
EXEC=a.out  # the dynamics of QSI with pairing

# FFTW
LIB=$(HOME)/opt/FFTW-3.3.8/lib/libfftw3.a

$(EXEC): $(SRC)
	$(FCC) $(FC) $(SRC) $(LIB) -o $(EXEC)
	rm -f *.mod

run:
	./$(EXEC)

ed:
	vim $(SRC)

clean:
	rm -f *.mod $(EXEC)
