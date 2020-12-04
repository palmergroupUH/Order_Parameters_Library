SHELL = /bin/sh

# intel fortran compiler 
FC = ifort #

# shared library and O3 optimization level
FFLAGS_D = -fPIC -shared -O3

LIB_DIR = lib/

INC_DIR = include/

sys_lib = ../general/src/system_mod.f90

const_include = -I ../general/include

local_include = -module include/

math_lib_sph = ../mathlib/src/spherical_harmonics_mod.f90

math_lib_wign3j = ../mathlib/src/wigner3j_symobol_mod.f90

math_lib_factorial = ../mathlib/src/factorial_mod.f90

sorting_lib  = ../algorithm/src/sort_mod.f90

cluster_lib = ../algorithm/src/single_linkage_cluster.f90

sph_lib = src/sph_op_mod.f90

nb_list_lib = src/nb_list_mod.f90 

all:  $(LIB_DIR)lib_HMC.so $(LIB_DIR)lib_HMC_q3.so 

# ----------------------------------------------------------------------------
#                             HMC
# ----------------------------------------------------------------------------

$(LIB_DIR)lib_HMC.so: $(cluster_lib) $(sys_lib) $(math_lib_factorial) $(math_lib_sph) $(math_lib_wign3j) $(sorting_lib) $(sph_lib) $(nb_list_lib) src/HMC_mod.f90

	@ mkdir -p $(LIB_DIR)
	@ mkdir -p $(INC_DIR)

	$(FC) $(FFLAGS_D) $(const_include) $(sys_lib) $(math_lib_factorial) $(math_lib_sph) $(math_lib_wign3j) $(sorting_lib) $(cluster_lib) $(sph_lib) $(nb_list_lib) src/HMC_mod.f90 -o $(LIB_DIR)lib_HMC.so $(local_include)

$(LIB_DIR)lib_HMC_q3.so: src/oplib_q3_only.f90  

	$(FC) $(FFLAGS_D) src/oplib_q3_only.f90 -module include/ -o $(LIB_DIR)lib_HMC_q3.so 

.PHONY: all clean count

clean:

	rm *.o 

