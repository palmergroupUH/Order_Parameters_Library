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

all: $(LIB_DIR)lib_RussoRomanoTanaka.so $(LIB_DIR)lib_CHILL.so $(LIB_DIR)lib_Ql_Wl.so $(LIB_DIR)lib_Li_et_al.so $(LIB_DIR)lib_q_tetra.so 

$(LIB_DIR)lib_q_tetra.so: src/tetrahedral_op_mod.f90 

	@ mkdir -p $(LIB_DIR)
	@ mkdir -p $(INC_DIR)

	$(FC) $(FFLAGS_D) $(sys_include) $(sorting_lib) src/tetrahedral_op_mod.f90 -o $(LIB_DIR)lib_q_tetra.so $(local_include)

$(LIB_DIR)lib_RussoRomanoTanaka.so: $(math_lib_factorial) $(math_lib_sph) $(math_lib_wign3j) src/order_parameter_mod.f90 src/Russo_Romano_Tanaka_mod.f90

	@ mkdir -p $(LIB_DIR)
	@ mkdir -p $(INC_DIR)

	$(FC) $(FFLAGS_D) $(sys_include) $(math_lib_factorial) $(math_lib_sph) $(math_lib_wign3j) $(sorting_lib) $(cluster_lib) src/order_parameter_mod.f90 src/Russo_Romano_Tanaka_mod.f90 -o $(LIB_DIR)lib_RussoRomanoTanaka.so $(local_include)

$(LIB_DIR)lib_CHILL.so: $(math_lib_factorial) $(math_lib_sph) $(cluster_lib) $(math_lib_wign3j) src/order_parameter_mod.f90 src/CHILL_mod.f90

	@ mkdir -p $(LIB_DIR)
	@ mkdir -p $(INC_DIR)

	$(FC) $(FFLAGS_D) $(sys_include) $(math_lib_factorial) $(math_lib_sph) $(math_lib_wign3j) $(sorting_lib) $(cluster_lib) src/order_parameter_mod.f90 src/CHILL_mod.f90 -o $(LIB_DIR)lib_CHILL.so $(local_include)

$(LIB_DIR)lib_Ql_Wl.so: $(math_lib_factorial) $(math_lib_sph) $(math_lib_wign3j) $(cluster_lib) src/order_parameter_mod.f90 src/test_Ql_Wl.f90

	@ mkdir -p $(LIB_DIR)
	@ mkdir -p $(INC_DIR)

	$(FC) $(FFLAGS_D) $(sys_include) $(math_lib_factorial) $(math_lib_sph) $(math_lib_wign3j) $(sorting_lib) $(cluster_lib) src/order_parameter_mod.f90 src/test_Ql_Wl.f90 -o $(LIB_DIR)lib_Ql_Wl.so $(local_include)

$(LIB_DIR)lib_Li_et_al.so: $(math_lib_factorial) $(math_lib_sph) $(math_lib_wign3j) $(cluster_lib) src/order_parameter_mod.f90 src/Li_et_al_mod.f90

	@ mkdir -p $(LIB_DIR)
	@ mkdir -p $(INC_DIR)

	$(FC) $(FFLAGS_D) $(sys_include) $(math_lib_factorial) $(math_lib_sph) $(math_lib_wign3j) $(sorting_lib) $(cluster_lib) src/order_parameter_mod.f90 src/Li_et_al_mod.f90 -o $(LIB_DIR)lib_Li_et_al.so $(local_include)
.PHONY: all clean count

clean:

	rm *.o 

