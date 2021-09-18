# Note: we can compare the times with option -O2 and -O3 and if they are huge
# different that means the code can be optimized very much
#
FCC=gfortran
#FCC=ifort

FC=-O0 -cpp -fbounds-check -Wuninitialized
#FC=

SRC=DynamicalExponent.f90
EXEC=dyz  # the dynamics of QSI with pairing

# FFTW
MY_INCLUDE_PATH=-I$(HOME)/opt/fftw-gcc/include
MY_LIBRARY_PATH=-L$(HOME)/opt/fftw-gcc/lib
LIBS=-lfftw3

$(EXEC): $(SRC)
	$(FCC) $(FC) $(SRC) -o $(EXEC) $(MY_INCLUDE_PATH) $(MY_LIBRARY_PATH) $(LIBS)
	rm -f *.mod

run:
	./$(EXEC)

ed:
	vim $(SRC)

clean:
	rm -f *.mod $(EXEC)
