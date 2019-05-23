FC=gfortran
LHAPDFCONFIG=lhapdf-config
FFLAGS=-O2
CXXFLAGS=`$(LHAPDFCONFIG) --cflags` -O2
LDFLAGS=`$(LHAPDFCONFIG) --ldflags` -lgfortran -lfmt
all:gridgen gridtest
gridgen:gridgen.o dgauss.o ElectroweakFlux.o gridcalc.o
	g++ $^ -o $@ $(CXXFLAGS) $(LDFLAGS)
ElectroweakFlux.o:ElectroweakFlux.f ElectroweakFlux.inc
gridtest:gridtest.o dfint.o kerset.o ewa_pp.o gridcalc.o dgauss.o ElectroweakFlux.o
	g++ $^ -o $@ $(CXXFLAGS) $(LDFLAGS)
ewa_pp.f:gridgen
	./gridgen
