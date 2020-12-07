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
RRTanaka_lib = CDLL(os.path.join(fortranlib_address, "lib_RussoTanaka.so"))

# Load CHILL library 
CHILL_lib = CDLL(os.path.join(fortranlib_address, "lib_CHILL.so"))

# Load general Ql and Wl library for testing
Ql_Wl_lib = CDLL(os.path.join(fortranlib_address, "lib_Ql_Wl.so"))

# Load Li et al 
Li_lib = CDLL(os.path.join(fortranlib_address, "lib_Li_et_al.so"))

# Load HMC library
HMC_lib = CDLL(os.path.join(fortranlib_address, "lib_HMC.so"))

# Load HMC q3 library
HMC_q3_only_lib = CDLL(os.path.join(fortranlib_address, "lib_HMC_q3.so"))

# ----------------------------------------------------------------------------
#                             HMC
# ----------------------------------------------------------------------------

def calc_q3_cluster_HMC(n_atoms,n_q3_neigh,q3_cutoff,n_bonds,box,x): 

    n_atoms = c_int(n_atoms)

    n_q3_neigh = c_int(n_q3_neigh)

    q3_cutoff = c_double(q3_cutoff)
    
    n_bonds = c_int(n_bonds)

    largest_cluster = c_int() 

    box = np_to_ctypes_array(box)  

    x = np_to_ctypes_array(x) 

    HMC_q3_only_lib.call_calc_q3_cluster(byref(n_atoms),
                                         byref(n_q3_neigh),
                                         byref(q3_cutoff),
                                         byref(n_bonds),
                                         box,
                                         x,
                                         byref(largest_cluster))

    return largest_cluster.value 


def calc_ql_cluster(style,
                    total_atoms,
                    box,
                    xyz,
                    maxnb,
                    nnb,
                    l,
                    ql_cutoff,
                    cluster_cutoff,
                    n_bonds):

    style_key, strlength = string_to_ctypes_string(style)

    total_atoms = c_int(total_atoms)

    maxnb = c_int(maxnb)

    nnb = c_int(nnb)

    l = c_int(l)
        
    ql_cutoff_sqr = c_double(ql_cutoff * ql_cutoff)

    cluster_cutoff = c_double(cluster_cutoff)

    n_bonds = c_int(n_bonds)

    nxtl = c_int()

    lcl = c_int()
    
    HMC_lib.call_calc_ql_cluster(style_key,  
                                 byref(strlength), 
                                 byref(total_atoms),
                                 box,
                                 xyz, 
                                 byref(maxnb),
                                 byref(nnb),
                                 byref(l),
                                 byref(ql_cutoff_sqr),
                                 byref(cluster_cutoff),
                                 byref(n_bonds),
                                 byref(nxtl),
                                 byref(lcl))


    return nxtl.value, lcl.value  



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


def call_RussoRomanoTanaka(total_atoms,
                           xyz,
                           box,
                           maxnb,
                           nnb,
                           l,
                           ql_cutoff,
                           n_bonds): 
                           

    label = np.ctypeslib.as_ctypes(np.zeros(total_atoms, dtype=np.int32))  

    Ql_invar = np.ctypeslib.as_ctypes(np.zeros(total_atoms, dtype=np.float64))  

    Wl_invar = np.ctypeslib.as_ctypes(np.zeros(total_atoms, dtype=np.float64)) 

    ql_cutoff_sqr = c_double(ql_cutoff * ql_cutoff)
    
    total_atoms = c_int(total_atoms)

    l = c_int(l)

    maxnb = c_int(maxnb)

    nnb = c_int(nnb) 

    n_bonds = c_int(n_bonds) 

    nxtl = c_int()

    RRTanaka_lib.call_Russo_Romano_Tanaka(byref(total_atoms),  
                                          box,
                                          xyz, 
                                          byref(maxnb),
                                          byref(nnb),  
                                          byref(l), 
                                          byref(ql_cutoff_sqr),
                                          byref(n_bonds), 
                                          byref(nxtl), 
                                          label,
                                          Ql_invar, 
                                          Wl_invar)
        
    label = np.ctypeslib.as_array(label)
    
    Ql_invar = np.ctypeslib.as_array(Ql_invar)
    
    Wl_invar = np.ctypeslib.as_array(Wl_invar)

    num_crys = nxtl.value 

    return num_crys, label, Ql_invar, Wl_invar 

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


