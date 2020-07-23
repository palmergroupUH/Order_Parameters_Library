SHELL = /bin/sh

# intel fortran compiler 
FC = ifort #

FFLAGS_D = -fPIC -shared

FFLAGS_S = -c 

INC_DIR = .include/

LIB_DIR = lib/

linker_syslib = -L../general/lib/ -lforsys 

sys_include = -I ../general/include

local_include = -module include/

system_file = ../general/src/system_mod.f90 

constatns_file = ../general/src/constants_mod.f90  

include  = -I include

all: factorial_mod.o spherical_harmonics_mod.o $(LIB_DIR)libsph.so $(LIB_DIR)libformath.a

factorial_mod.o: src/factorial_mod.f90

	$(FC) $(FFLAGS_S) $(linker_syslib) src/factorial_mod.f90 $(local_include) $(sys_include)

spherical_harmonics_mod.o: src/spherical_harmonics_mod.f90 

	$(FC) $(FFLAGS_S) $(linker_syslib) src/spherical_harmonics_mod.f90 $(local_include) $(sys_include) 

$(LIB_DIR)libsph.so: src/factorial_mod.f90  src/spherical_harmonics_mod.f90 

	$(FC) $(FFLAGS_D) $(sys_include) src/factorial_mod.f90  src/spherical_harmonics_mod.f90 -o $(LIB_DIR)libsph.so $(local_include)

$(LIB_DIR)libformath.a:  

	ar -cvq $(LIB_DIR)libformath.a factorial_mod.o spherical_harmonics_mod.o  

.PHONY: all clean count

clean:

	rm *.o 

