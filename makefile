FC = gfortran
FFLAGS = -O2 -Wall

SRC = Precision.f90 CIUtils.f90 Constants.f90 Input.f90 operators.f90 \
      Integrals.f90 BuildHam.f90 Lanczos.f90 Main.f90

EXE = Xcode

all: $(EXE)

$(EXE): $(SRC)
	$(FC) $(FFLAGS) $(SRC) -framework Accelerate -o $(EXE)

clean:
	rm -f *.o *.mod $(EXE)
