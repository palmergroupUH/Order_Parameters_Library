SHELL = /bin/sh

# intel fortran compiler 
FC = ifort #

FFLAGS_D = -fPIC -shared

LIB_DIR = lib/

sys_include = -I ../general/include

local_include = -module include/

math_lib_sph = ../mathlib/src/spherical_harmonics_mod.f90

math_lib_factorial = ../mathlib/src/factorial_mod.f90

all: $(LIB_DIR)lib_RussoRomanoTanaka.so 

$(LIB_DIR)lib_RussoRomanoTanaka.so: $(math_lib_factorial) $(math_lib_sph) src/order_parameter_mod.f90 src/Russo_Romano_Tanaka_mod.f90

	$(FC) $(FFLAGS_D) $(sys_include) $(math_lib_factorial) $(math_lib_sph) src/order_parameter_mod.f90 src/Russo_Romano_Tanaka_mod.f90 -o $(LIB_DIR)lib_RussoRomanoTanaka.so $(local_include)

.PHONY: all clean count

clean:

	rm *.o 

