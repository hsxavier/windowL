# Variable definitions (change things here):
COMP    = g++
HEALDIR = /home/skems/prog/Healpix_3.11
BIN     = .
OBJ     = .
SRC     = .
CXXOMP  = -fopenmp
LIBS    = -lgsl -lgslcblas -lhealpix_cxx -lcxxsupport -lcfitsio -lsharp -lfftpack -lc_utils -lgomp


# Options for Splinter:
ifeq ($(HOSTNAME), splinter-login.local)
HEALDIR = /home/hsxavier/Healpix_3.11
CXXFITS = -I/usr/include/cfitsio
endif
# Healpix files:
HEALDATA = $(HEALDIR)/data
CXXHEAL  = -I$(HEALDIR)/src/cxx/generic_gcc/include
LDHEAL   = -L$(HEALDIR)/src/cxx/generic_gcc/lib


# General instructions:
all: $(BIN)/windowL

clean:
	rm -f $(BIN)/windowL
	rm -f $(OBJ)/windowL-made.cpp
	rm -f $(OBJ)/*.o


# Compiling and etc:
$(BIN)/windowL: $(OBJ)/windowL.o $(OBJ)/Utilities.o
	$(COMP) $(LDHEAL) $(OBJ)/windowL.o $(OBJ)/Utilities.o -o $@ $(LIBS)

$(OBJ)/windowL.o: $(SRC)/windowL.cpp $(SRC)/Utilities.hpp
	sed 's|#define HEALPIX_DATA .*|#define HEALPIX_DATA \"$(HEALDATA)\"|g' $(SRC)/windowL.cpp > $(OBJ)/windowL-made.cpp
	$(COMP) $(CXXHEAL) -c $(OBJ)/windowL-made.cpp -o $@ $(CXXOMP)

$(OBJ)/Utilities.o: $(SRC)/Utilities.cpp $(SRC)/Utilities.hpp
	$(COMP) $(CXXHEAL) -c $(SRC)/Utilities.cpp -o $@
