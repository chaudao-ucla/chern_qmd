##----------------------------------------
## taken from Makefile.inc...

PLAT = linux
CPP = g++


LAPACK_LIB   = -llapack
BLAS_LIB     = -lblas
ARMA_LIB     = -larmadillo
ARPACK_LIB   = -larpack


CPP_WARNINGS = -Wall -pedantic-errors
CPP_DEBUG    = -g
CPP_OPTIM    = -O2

CXXFLAGS    = -std=c++11 $(CPP_DEBUG) $(CPP_WARNINGS) $(CPP_OPTIM)

ALL_LIBS     = $(LAPACK_LIB) $(BLAS_LIB)  $(ARPACK_LIB) $(ARMA_LIB)

## --- to use the latest version of armadillo in MPIPKS
INC = -I/usr/local/armadillo/usr/local/include

DEPS = util.hpp model.hpp vec3.hpp analysis.hpp test.hpp

##----------------------------------------
## taken from examples/product/simple/Makefile...


gvn: GvNMD.o util.o analysis.o tdga.o qmd.o lattice.o test.o
	$(CPP) $(INC) $(CXXFLAGS) -o gvn GvNMD.o util.o analysis.o tdga.o qmd.o lattice.o test.o $(ALL_LIBS)

GvNMD.o: GvNMD.cpp $(DEPS)
	$(CPP) $(INC) $(CXXFLAGS) -c $<

qmd.o: qmd.cpp $(DEPS)
	$(CPP) $(INC) $(CXXFLAGS) -c $<

tdga.o: tdga.cpp $(DEPS)
	$(CPP) $(INC) $(CXXFLAGS) -c $<

lattice.o: lattice.cpp $(DEPS)
	$(CPP) $(INC) $(CXXFLAGS) -c $<

%.o: %.cpp %.hpp $(DEPS)
	$(CPP) $(INC) $(CXXFLAGS) -c $<
	


.PHONY: clean
clean:
	rm -f *~ *.o core gvn



