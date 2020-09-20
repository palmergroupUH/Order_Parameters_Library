SHELL = /bin/sh

# intel fortran compiler 
FC = ifort #

# compile dynamic library:
FFLAGS_D = -fPIC -shared

# compile static library
FFLAGS_S = -c 

INC_DIR = .include/

LIB_DIR = lib/

linker_syslib = -L../general/lib/ -lforsys 

sys_include = -I ../general/include

all: single_linkage_cluster.o sort_mod.o linked_list_mod.o linspace_mod.o $(LIB_DIR)libforalgorithm.a  

sort_mod.o: src/sort_mod.f90

	$(FC) $(FFLAGS_S) $(linker_syslib) src/sort_mod.f90 -module include/ $(sys_include) 

single_linkage_cluster.o: src/single_linkage_cluster.f90

	$(FC) $(FFLAGS_S) $(linker_syslib) src/single_linkage_cluster.f90 -module include/ $(sys_include) 

linked_list_mod.o: src/linked_list_mod.f90

	$(FC) $(FFLAGS_S) $(linker_syslib) src/linked_list_mod.f90 -module include/ $(sys_include) 

linspace_mod.o: src/linspace_mod.f90 

	$(FC) $(FFLAGS_S) $(linker_syslib) src/linspace_mod.f90 -module include/ $(sys_include)  

$(LIB_DIR)libforalgorithm.a: sort_mod.o linked_list_mod.o linspace_mod.o  

	ar -cvq $(LIB_DIR)libforalgorithm.a sort_mod.o linked_list_mod.o linspace_mod.o  

.PHONY: all clean count

clean:

	rm *.o 

