SHELL = /bin/sh

# intel fortran compiler 
FC = ifort #

FFLAGS_D = -fPIC -shared

FFLAGS_S = -c 

INC_DIR = .include/

LIB_DIR = lib/

linker_syslib = -L../general/lib/ -lforsys 

sys_include = -module ../general/include

all: sort_mod.o $(LIB_DIR)libforalgorithm.a 

sort_mod.o: src/sort_mod.f90

	$(FC) $(FFLAGS_S) $(linker_syslib) src/sort_mod.f90 -module include/ $(sys_include) 

$(LIB_DIR)libformath.a:  

	ar -cvq $(LIB_DIR)libforalgorithm.a sort_mod.o 

.PHONY: all clean count

clean:

	rm *.o 

