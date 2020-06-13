# Python standard library: 
import numpy as np  
import os 

# Local library: 

# fortran API: 
import 
# Third-party libraries: 


fortranlib_address =  os.path.join(os.path.dirname(IO.__file__) , "fortran") 

potential_lib = CDLL(os.path.join(fortranlib_address,"fortran_lj_pot.so")) 


def call_lj_potential_due_to_i( ):  

    return None 

