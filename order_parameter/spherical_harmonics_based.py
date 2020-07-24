# Python standard library
import numpy as np
import os
import sys
import itertools
import logging
from ctypes import CDLL, c_int, c_double, c_float, byref

# Local library:
import IO.reader
from IO.type_conversion import string_to_ctypes_string,\
                               int_to_ctypes_int,\
                               np_to_ctypes_array
import order_parameter 

# determine the library address  
fortranlib_address = os.path.join(os.path.dirname(order_parameter.__file__), "lib")

# Load RussoRomanoTanaka library 
RRTanaka_lib = CDLL(os.path.join(fortranlib_address, "lib_RussoRomanoTanaka.so"))

# Load CHILL library 
CHILL_lib = CDLL(os.path.join(fortranlib_address, "lib_CHILL.so"))


# ----------------------------------------------------------------------------
#                           RussoRomanoTanaka 
# ----------------------------------------------------------------------------

def intialize_RussoRomanoTanaka(l):

    l_dim = 2*l + 1 

    sph_const = np.ctypeslib.as_ctypes(np.zeros(l_dim, dtype=np.float64))
    
    Plm_const = np.ctypeslib.as_ctypes(np.zeros((l+1)*4, dtype=np.float64))

    l = c_int(l)

    RRTanaka_lib.initialize_RussoRomanoTanaka(byref(l), sph_const, Plm_const)

    return sph_const, Plm_const 

def call_RussoRomanoTanaka(sph_const, Plm_const, total_atoms, l, xyz, box, maxnb, nnb, cutoff_sqr):

    cij = np.ctypeslib.as_ctypes(np.zeros(nnb*total_atoms, dtype=np.float64))    

    cutoff_sqr = c_double(cutoff_sqr)

    l = c_int(l)

    total_atoms = c_int(total_atoms)

    maxnb = c_int(maxnb)

    nnb = c_int(nnb) 

    RRTanaka_lib.call_RussoRomanoTanaka(sph_const,
                                        Plm_const,
                                        byref(total_atoms),
                                        byref(maxnb),
                                        byref(nnb),  
                                        byref(l),
                                        byref(cutoff_sqr),
                                        box,
                                        xyz,
                                        cij)

    cij = np.ctypeslib.as_array(cij)

    return cij

# ----------------------------------------------------------------------------
#                                 CHILL 
# ----------------------------------------------------------------------------

def intialize_CHILL(l):

    l_dim = 2*l + 1 

    sph_const = np.ctypeslib.as_ctypes(np.zeros(l_dim, dtype=np.float64))
    
    Plm_const = np.ctypeslib.as_ctypes(np.zeros((l+1)*4, dtype=np.float64))

    l = c_int(l)

    CHILL_lib.initialize_CHILL(byref(l), sph_const, Plm_const)

    return sph_const, Plm_const

def call_CHILL(keyword, sph_const, Plm_const, total_atoms, l, xyz, box, maxnb, nnb, cutoff_sqr):
    
    CHILL_keyword, strlength = string_to_ctypes_string(keyword) 
    
    cij = np.ctypeslib.as_ctypes(np.zeros(nnb*total_atoms, dtype=np.float64))    
    
    chill_ID_list = np.ctypeslib.as_ctypes(np.zeros(total_atoms, dtype=np.int32))    

    cutoff_sqr = c_double(cutoff_sqr)

    l = c_int(l)

    total_atoms = c_int(total_atoms)

    maxnb = c_int(maxnb)

    nnb = c_int(nnb)

    CHILL_lib.call_CHILL(sph_const,
                         Plm_const,
                         CHILL_keyword,
                         byref(strlength),
                         byref(total_atoms),
                         byref(maxnb),
                         byref(nnb),  
                         byref(l),
                         byref(cutoff_sqr),
                         box,
                         xyz,
                         cij,
                         chill_ID_list)

    return np.ctypeslib.as_array(cij), np.ctypeslib.as_array(chill_ID_list)
