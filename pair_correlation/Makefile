SHELL = /bin/sh

# intel fortran compiler
FC = ifort #

FFLAGS_D = -fPIC -shared

LIB_DIR = lib/

INC_DIR = include/

sys_include = -I ../general/include

local_include = -module include/

mathlib = ../mathlib/src/

all: $(LIB_DIR)lib_rdf.so $(LIB_DIR)lib_ssf.so  

$(LIB_DIR)lib_rdf.so: src/compute_rdf_mod.f90 

	$(FC) $(FFLAGS_D) $(sys_include) src/compute_rdf_mod.f90 -o $(LIB_DIR)lib_rdf.so $(local_include) 

$(LIB_DIR)lib_ssf.so: src/compute_sf_mod.f90

	$(FC) $(FFLAGS_D) $(sys_include) $(mathlib)distribution.f90 src/compute_sf_mod.f90 -o $(LIB_DIR)lib_ssf.so $(local_include) 

.PHONY: all clean count

clean:

	rm *.o 

