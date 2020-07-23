SHELL = /bin/sh

# fortran compiler 
FC = ifort #

FFLAGS_D = -fPIC -shared

FFLAGS_S = -c

INC_DIR = .include/

LIB_DIR = lib/

all: $(LIB_DIR)libforsys.a  

system_mod.o: src/system_mod.f90 

	$(FC) $(FFLAGS_S) src/system_mod.f90  -module include/

constants_mod.o: src/constants_mod.f90

	$(FC) $(FFLAGS_S) src/constants_mod.f90  -module include/

$(LIB_DIR)libforsys.a: system_mod.o constants_mod.o  

	ar -cvq $(LIB_DIR)libforsys.a system_mod.o constants_mod.o 

.PHONY: all clean count

clean:

	rm *.o 

