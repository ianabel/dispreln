# 
# Copyright Ian G Abel 2018
#

# If boost is installed site-wide, then don't set anything and walk on
# If you have installed a local boost install, set BOOSTDIR to the location
# of the install, or you can set BOOSTINCDIR and BOOSTLIBDIR manually

ifdef BOOSTDIR
BOOSTINCDIR = $(BOOSTDIR)/include
BOOSTLIBDIR = $(BOOSTDIR)/lib
endif

ifdef BOOSTINCDIR
CXXFLAGS += -I$(BOOSTINCDIR) -L$(BOOSTLIBDIR) -Wl,-rpath-link=$(BOOSTLIBDIR)  
endif

CXXFLAGS += -Wall -pedantic -std=c++17
CXXDEBUGFLAGS +=  -Og -g
CXXTESTFLAGS +=   -DBOOST_TEST_DYN_LINK -lboost_unit_test_framework
CXXRELEASEFLAGS += -Ofast


all: DispSolver

DispSolver: DispSolver.cpp Solver.cpp RootFinder.o RootFinder.h Faddeeva.o Faddeeva.hh DispReln.cpp Config.cpp
	$(CXX) $(CXXFLAGS) -Ofast -lboost_program_options -o $@ DispSolver.cpp RootFinder.cpp Solver.cpp Faddeeva.o DispReln.cpp Config.cpp

test: RootFindingTests DispRelnTests
	LD_LIBRARY_PATH=$(BOOSTLIBDIR) ./RootFindingTests
	LD_LIBRARY_PATH=$(BOOSTLIBDIR) ./DispRelnTests


RootFindingTests: RootFindingTests.cpp RootFinder.h RootFinder.cpp Solver.cpp
	$(CXX) $(CXXFLAGS) $(CXXDEBUGFLAGS) $(CXXTESTFLAGS) -o $@ RootFindingTests.cpp RootFinder.cpp Solver.cpp

DispRelnTests: DispRelnTests.cpp RootFinder.h DispReln.cpp DispReln.h RootFinder.cpp Solver.cpp  Faddeeva.o Faddeeva.hh Config.h Config.cpp
	$(CXX) $(CXXFLAGS) $(CXXDEBUGFLAGS) $(CXXTESTFLAGS) -o $@ DispRelnTests.cpp RootFinder.cpp Solver.cpp Faddeeva.o DispReln.cpp Config.cpp

# We never need a debug-enabled version of this.
Faddeeva.o: Faddeeva.cc Faddeeva.hh
	$(CXX) $(CXXFLAGS) $(CXXRELEASEFLAGS) -c -o $@ Faddeeva.cc


clean:
	rm -f DispRelnTests RootFindingTests *.o

.PHONY: clean test
