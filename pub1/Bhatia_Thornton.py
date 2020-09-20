# Python standard library
import numpy as np
import os
import sys
import itertools
import logging
from ctypes import CDLL, c_int, c_double, c_float, byref, POINTER

# Local library:
import mathlib
import pair_correlation 
import IO.reader
from IO.type_conversion import string_to_ctypes_string,\
                               int_to_ctypes_int,\
                               np_to_ctypes_array

import pub1.Bhatia_Thornton

# determine the library address

fortranlib_address = os.path.join(os.path.dirname(pub1.__file__), "lib") 

hist_lib_addr = os.path.join(os.path.dirname(mathlib.__file__), "lib") 

pair_corr_addr = os.path.join(os.path.dirname(pair_correlation.__file__), "lib") 

# Load Bhatia_Thornton library:

bt_lib = CDLL(os.path.join(fortranlib_address, "lib_bhatia_thornton.so"))  

# Load the pair_correlatikon library

pair_correl_lib = CDLL(os.path.join(pair_corr_addr, "lib_rdf.so"))

# Load the histogram library: 
hist_lib = CDLL(os.path.join(hist_lib_addr, "libforhist.so"))

def initialize_hist(num_bins, lower, upper): 

    hist = np.ctypeslib.as_ctypes(np.zeros(num_bins, dtype=np.float64)) 

    r_mid = np.ctypeslib.as_ctypes(np.zeros(num_bins, dtype=np.float64))

    num_bins = c_int(num_bins) 

    lower = c_double(lower) 

    upper = c_double(upper) 

    interval = c_double() 

    hist_lib.call_init_hist(byref(lower), byref(upper), byref(num_bins), hist, byref(interval), r_mid)

    return hist, interval, r_mid  

def update_hist(N, array, interval, lower, upper, num_bins, hist):  

    num_bins = c_int(num_bins) 

    lower = c_double(lower) 

    upper = c_double(upper) 

    N = c_int(N) 

    hist_lib.call_array_into_hist(byref(N), array, byref(interval), byref(upper), byref(lower), byref(num_bins), hist)

    return hist

def normalize_hist(r_mid, num_bins, interval, q_dist):

    num_bins = c_int(num_bins) 

    hist_lib.call_normalize_hist(byref(num_bins), byref(interval), q_dist)

    return np.ctypeslib.as_array(r_mid), np.ctypeslib.as_array(q_dist) 

def compute_q_tetra(total_atoms, xyz, box):

    sorted_q_indx = np.ctypeslib.as_ctypes(np.zeros(total_atoms, dtype=np.int32)) 

    q_tetra = np.ctypeslib.as_ctypes(np.zeros(total_atoms, dtype=np.float64)) 

    total_atoms = int_to_ctypes_int(total_atoms) 

    maxnb = c_int(35)

    nnb = c_int(4)

    cutoff_sqr = c_double(5.8*5.8)
    
    bt_lib.call_label_high_q_low_q(byref(total_atoms), byref(maxnb), byref(nnb), xyz, box, byref(cutoff_sqr), sorted_q_indx, q_tetra)

    return np.ctypeslib.as_array(sorted_q_indx), q_tetra  

def get_num_species(keyword, total_atoms):

    if (keyword =="HDL_HDL" or keyword == "LDL_LDL" or keyword == "HDL_LDL"):

        return int(total_atoms/2) 

    elif (keyword =="All"): 

        return total_atoms 

    else: 

        sys.exit("keyword not recognized ! Choose: 'HDL-HDL', 'LDL-LDL', 'HDL-LDL', and 'All'")

def compute_pair_correl(keyword, num_bins, total_atoms, xyz, box, sorted_q_indx):

    RDF_hist = np.ctypeslib.as_ctypes(np.zeros(num_bins, dtype=np.float64))  

    num_bins = c_int(num_bins) 

    if (keyword =="HDL_HDL"): 
        
        n_atoms = get_num_species("HDL_HDL", total_atoms)   
    
        norm_atom = n_atoms

        HDL_index = np.ctypeslib.as_ctypes(sorted_q_indx[0:n_atoms])  

        n_atoms = c_int(n_atoms) 

        bt_lib.call_homo_q_pair_corrl(byref(n_atoms), HDL_index, byref(num_bins), xyz, box, RDF_hist)

    elif (keyword=="LDL_LDL"): 

        n_atoms = get_num_species("LDL_LDL", total_atoms)  

        norm_atom = n_atoms

        LDL_index = np.ctypeslib.as_ctypes(sorted_q_indx[n_atoms:]) 

        n_atoms = c_int(n_atoms) 

        bt_lib.call_homo_q_pair_corrl(byref(n_atoms), LDL_index, byref(num_bins), xyz, box, RDF_hist)

    elif (keyword=="HDL_LDL"): 
    
        total_atoms_A = get_num_species("HDL_LDL", total_atoms) 
    
        norm_atom = total_atoms_A

        # equal molar binary mixture
        total_atoms_B = total_atoms_A

        HDL_index = np.ctypeslib.as_ctypes(sorted_q_indx[0:total_atoms_A])  

        LDL_index = np.ctypeslib.as_ctypes(sorted_q_indx[total_atoms_A:])  

        total_atoms = c_int(total_atoms) 

        bt_lib.call_hetero_q_pair_corrl(byref(total_atoms), byref(total_atoms_A), byref(total_atoms_B), HDL_index, LDL_index, byref(num_bins), xyz, box, RDF_hist) 

    elif (keyword == "All"):

        all_index = np.ctypeslib.as_ctypes(np.arange(1,total_atoms+1, dtype=np.int32))  
        
        norm_atom = total_atoms

        total_atoms = c_int(total_atoms) 

        bt_lib.call_homo_q_pair_corrl(byref(total_atoms), all_index, byref(num_bins), xyz, box, RDF_hist)

    else: 

        sys.exit("keyword not recognized ! Choose: 'HDL-HDL', 'LDL-LDL', 'HDL-LDL', and 'All'")

    return RDF_hist, norm_atom

def update_pair_correl_hist(num_bins, new_hist, RDF_hist): 

    num_bins = c_int(num_bins)

    hist_lib.call_combine_hist(byref(num_bins), new_hist, RDF_hist)  

    return None   

def normalize_pair_correl_hist(rdf_hist, num_bins, cutoff, natoms, num_configs, box):  
    
    density = c_double( natoms / np.prod(np.ctypeslib.as_array(box)/10) )  

    gr = np.ctypeslib.as_ctypes(np.zeros(num_bins, dtype=np.float64))
    
    r2hr = np.ctypeslib.as_ctypes(np.zeros(num_bins, dtype=np.float64)) 

    num_bins = c_int(num_bins) 
    
    cutoff = c_double(cutoff) 

    natoms = c_int(natoms) 

    num_configs = c_int(num_configs)
    
    pair_correl_lib.call_normalize_hist(rdf_hist, byref(num_bins), byref(cutoff), byref(natoms), byref(num_configs), byref(density), gr, r2hr)

    return np.ctypeslib.as_array(r2hr), np.ctypeslib.as_array(gr)  


def write_hist_file(keyword, total_atoms, box, num_configs, interval, num_bins, hist):
   
    num_atoms = get_num_species(keyword, total_atoms) 

    # all atoms: total_atoms * num_configs
    first_row = num_atoms * num_configs 
    
    # all volumes: volume per frame * num_configs
    second_row = np.prod(np.ctypeslib.as_array(box)/10.0) * num_configs
    
    # all atoms except reference particle at the center: (N-1) * num_configs
    third_row = (num_atoms - 1) * num_configs

    # number of configurations
    fourth_row = num_configs 

    # interval: cutoff / num_bins 
    fifth_row = interval.value 
    
    filename = "RDF_HIST_" + "%s" % keyword + "_%d" % (num_configs) + ".txt" 

    hist = np.prod(np.ctypeslib.as_array(hist)) 

    np.savetxt(filename, np.c_[ [first_row, second_row, third_row, fourth_row, fifth_row]])

    #for i in range(num_bins): 

    return None 





