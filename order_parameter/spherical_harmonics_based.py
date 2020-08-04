# Python standard library
import numpy as np
import os
import sys
import itertools
import logging
from ctypes import CDLL, c_int, c_double, c_float, byref, POINTER

# Local library:
import IO.reader
from IO.type_conversion import string_to_ctypes_string,\
                               int_to_ctypes_int,\
                               np_to_ctypes_array

from mathlib.math_API import compute_wigner_3j
import order_parameter 

# determine the library address  
fortranlib_address = os.path.join(os.path.dirname(order_parameter.__file__), "lib")

# Load RussoRomanoTanaka library 
RRTanaka_lib = CDLL(os.path.join(fortranlib_address, "lib_RussoRomanoTanaka.so"))

# Load CHILL library 
CHILL_lib = CDLL(os.path.join(fortranlib_address, "lib_CHILL.so"))

# Load general Ql and Wl library for testing
Ql_Wl_lib = CDLL(os.path.join(fortranlib_address, "lib_Ql_Wl.so"))

# Load Li et al 
Li_lib = CDLL(os.path.join(fortranlib_address, "lib_Li_et_al.so"))

# ----------------------------------------------------------------------------
#                             General Ql and Wl 
# ----------------------------------------------------------------------------

def generate_wigner3j_symobol_mat(l):

    wigner3j_symobl_vals = []

    m_index = []

    num_pairs = 0
    
    for m1 in range(-l, l+1):

        for m2 in range(-l, l+1):

            for m3 in range(-l, l+1):

                if ((m1 + m2 + m3) == 0):

                    num_pairs += 1
                    
                    matrix = np.array([[l, l, l], [m1, m2, m3]], dtype=np.int32)
                               
                    wigner3j = compute_wigner_3j(matrix) 
 
                    wigner3j_symobl_vals.append(wigner3j)  

                    m_index.append([m1, m2, m3]) 
 
    return num_pairs, np.array(wigner3j_symobl_vals), np.array(m_index) 

def initialize_Ql_Wl(l):

    num_pairs_w3j, wigner3j_symobl_ary, m_index_ary = generate_wigner3j_symobol_mat(l)

    l_dim = 2*l + 1

    sph_const = np.ctypeslib.as_ctypes(np.zeros(l_dim, dtype=np.float64))

    Plm_const = np.ctypeslib.as_ctypes(np.zeros((l+1)*4, dtype=np.float64))

    l = c_int(l)

    Ql_Wl_lib.initialize_test_Wl_Ql(byref(l), sph_const, Plm_const)

    return num_pairs_w3j, wigner3j_symobl_ary, m_index_ary, sph_const, Plm_const

def call_Ql_Wl(num_pairs_w3j,
               wigner3j_symobl_ary,
               m_index_ary,
               sph_const,
               Plm_const,
               total_atoms,
               l,
               xyz,
               box,
               maxnb,
               nnb,
               cutoff_sqr):

    num_pairs_w3j = c_int(num_pairs_w3j) 
    
    wigner3j_symobl_ary = np.ctypeslib.as_ctypes(wigner3j_symobl_ary.astype(np.float64)) 

    m_index_ary = np.ctypeslib.as_ctypes(m_index_ary.astype(np.int32).reshape(1,m_index_ary.size))

    cutoff_sqr = c_double(cutoff_sqr)

    l = c_int(l)

    total_atoms = c_int(total_atoms)

    maxnb = c_int(maxnb)

    nnb = c_int(nnb)

    cluster_id_list = POINTER(c_int)() 

    num_ice = c_int()  

    Ql_Wl_lib.call_test_Wl_Ql(sph_const,
                              Plm_const,
                              m_index_ary, 
                              wigner3j_symobl_ary,
                              byref(num_pairs_w3j),
                              byref(total_atoms),
                              byref(maxnb),
                              byref(nnb),
                              byref(l),
                              byref(cutoff_sqr),
                              box,
                              xyz)

    return None 


# ----------------------------------------------------------------------------
#                           RussoRomanoTanaka 
# ----------------------------------------------------------------------------


def intialize_RussoRomanoTanaka(l):

    num_pairs_w3j, wigner3j_symobl_ary, m_index_ary = generate_wigner3j_symobol_mat(l)

    l_dim = 2*l + 1 

    sph_const = np.ctypeslib.as_ctypes(np.zeros(l_dim, dtype=np.float64))
    
    Plm_const = np.ctypeslib.as_ctypes(np.zeros((l+1)*4, dtype=np.float64))

    l = c_int(l)

    RRTanaka_lib.initialize_RussoRomanoTanaka(byref(l), sph_const, Plm_const)

    return num_pairs_w3j, wigner3j_symobl_ary, m_index_ary, sph_const, Plm_const 

def call_RussoRomanoTanaka(sph_const, Plm_const, total_atoms, l, xyz, box, maxnb, nnb, cutoff_sqr, crys_cut):

    cij = np.ctypeslib.as_ctypes(np.zeros(nnb*total_atoms, dtype=np.float64))    
    
    count_bonds = np.ctypeslib.as_ctypes(np.zeros(total_atoms, dtype=np.float64))

    cutoff_sqr = c_double(cutoff_sqr)

    crys_cut = c_double(crys_cut)

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
                                        byref(crys_cut),
                                        box,
                                        xyz,
                                        cij,    
                                        count_bonds)

    cij = np.ctypeslib.as_array(cij)

    count_bonds = np.ctypeslib.as_array(count_bonds) 
    
    return cij, count_bonds 

# ----------------------------------------------------------------------------
#                                 CHILL/CHILL+ 
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

    cij_hist = np.ctypeslib.as_ctypes(np.zeros(20, dtype=np.float64))

    chill_ID_list = np.ctypeslib.as_ctypes(np.zeros(total_atoms, dtype=np.int32))    

    cutoff_sqr = c_double(cutoff_sqr)

    l = c_int(l)

    total_atoms = c_int(total_atoms)

    maxnb = c_int(maxnb)

    nnb = c_int(nnb)

    num_ice = c_int()  

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
                         cij_hist,
                         chill_ID_list,
                         byref(num_ice))

    return np.ctypeslib.as_array(cij_hist), np.ctypeslib.as_array(chill_ID_list)

# ----------------------------------------------------------------------------
#                              Q6 FFS by Li et al. 
# ----------------------------------------------------------------------------

def intialize_Li_et_al(l):

    l_dim = 2*l + 1

    sph_const = np.ctypeslib.as_ctypes(np.zeros(l_dim, dtype=np.float64))
    
    Plm_const = np.ctypeslib.as_ctypes(np.zeros((l+1)*4, dtype=np.float64))

    l = c_int(l)

    Li_lib.initialize_Li_et_al(byref(l), sph_const, Plm_const)

    return sph_const, Plm_const

def call_Li_et_al(sph_const, Plm_const, total_atoms, maxnb, nnb, l, cutoff_sqr,  box, xyz): 

    Iij = np.ctypeslib.as_ctypes(np.zeros(total_atoms, dtype=np.float64))    

    cutoff_sqr = c_double(cutoff_sqr)

    l = c_int(l)

    total_atoms = c_int(total_atoms)

    maxnb = c_int(maxnb)

    nnb = c_int(nnb)

    num_ice = c_int()

    Li_lib.call_Li_et_al(sph_const,
                         Plm_const,
                         byref(total_atoms),
                         byref(maxnb),
                         byref(nnb),
                         byref(l),
                         byref(cutoff_sqr),
                         box,
                         xyz,
                         Iij)

    return np.ctypeslib.as_array(Iij) 



