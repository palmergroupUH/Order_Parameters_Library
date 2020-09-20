SHELL = /bin/sh

# intel fortran compiler 
FC = ifort #

FFLAGS_D = -fPIC -shared

LIB_DIR = lib/

INC_DIR = include/

sys_include = -I ../general/include

local_include = -module include/

math_lib_sph = ../mathlib/src/spherical_harmonics_mod.f90

math_lib_wign3j = ../mathlib/src/wigner3j_symobol_mod.f90

math_lib_factorial = ../mathlib/src/factorial_mod.f90

sorting_lib  = ../algorithm/src/sort_mod.f90

cluster_lib = ../algorithm/src/single_linkage_cluster.f90

op_lib = ../order_parameter/src/order_parameter_mod.f90 

tetra_lib = ../order_parameter/src/tetrahedral_op_mod.f90 

pair_lib = ../pair_correlation/src/compute_rdf_mod.f90 

all: $(LIB_DIR)lib_bhatia_thornton.so 

$(LIB_DIR)lib_bhatia_thornton.so: src/bhatia_thornton_analysis.f90 

	$(FC) $(FFLAGS_D) $(sys_include) $(math_lib_factorial) $(math_lib_sph) $(math_lib_wign3j) $(sorting_lib) $(op_lib) $(tetra_lib) $(pair_lib) src/bhatia_thornton_analysis.f90 -o $(LIB_DIR)lib_bhatia_thornton.so $(local_include) 

.PHONY: all clean count

clean:

	rm *.o 

