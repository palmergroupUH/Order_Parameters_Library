# Python standard library
import numpy as np
import os
import sys
import itertools
import logging
from ctypes import CDLL, c_int, c_double, c_float, byref, POINTER, c_bool

# Local library:
import mathlib
import pair_correlation 
import IO.reader
from IO.type_conversion import string_to_ctypes_string,\
                               int_to_ctypes_int,\
                               double_to_ctypes_double,\
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

# Load the static structure factor library: 

ssf_lib = CDLL(os.path.join(pair_corr_addr, "lib_ssf.so"))

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

        total_atoms_A = c_int(total_atoms_A) 

        total_atoms_B = c_int(total_atoms_B) 

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

def normalize_from_hist_file(rdf_hist_data): 

    num_configs = rdf_hist_data[3]

    dr = rdf_hist_data[4]

    rdf_hist = np.copy(rdf_hist_data[5:])  

    num_bins =  rdf_hist.size  

    natoms_all_configs = rdf_hist_data[0]   
     
    volume_all_configs = rdf_hist_data[1] 

    cutoff = dr * num_bins 

    r_mid = np.zeros(num_bins) 

    for i in range(num_bins):  

        r_mid[i] = 0.5 * dr + dr * i  
    
    # Convert Python data type into C types:

    gr = np.ctypeslib.as_ctypes(np.zeros(num_bins, dtype=np.float64))
    
    r2hr = np.ctypeslib.as_ctypes(np.zeros(num_bins, dtype=np.float64)) 

    num_bins = c_int(num_bins)

    natoms = c_int(int(natoms_all_configs / num_configs))

    rho = c_double(natoms_all_configs / volume_all_configs)

    cutoff = c_double(cutoff)

    num_configs = c_int(int(num_configs))

    rdf_hist = np.ctypeslib.as_ctypes(rdf_hist)
    
    pair_correl_lib.call_normalize_hist(rdf_hist, byref(num_bins), byref(cutoff), byref(natoms), byref(num_configs), byref(rho), gr, r2hr)

    return natoms.value, rho.value, num_bins.value, r_mid, np.ctypeslib.as_array(gr), np.ctypeslib.as_array(r2hr)


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

    hist = np.ctypeslib.as_array(hist) 

    with open(filename, "w") as output: 

        np.savetxt(output, np.c_[ [first_row, second_row, third_row, fourth_row, fifth_row]])
        
        np.savetxt(output, np.c_[hist]) 

    return None 

def load_multiple_hist_file(filename_lst):

    for i, histfile in enumerate(filename_lst):

        data = IO.reader.np_loadtxt(histfile) 

        all_rdf_data = np.zeros((data[:,0].size, 3))

        if ("LDL_LDL" in histfile):

            all_rdf_data[:,0] = data[:,0]

        elif ("HDL_HDL" in histfile): 

            all_rdf_data[:,1] = data[:,0]

        elif ("HDL_LDL" in histfile):

            all_rdf_data[:,2] =  data[:,0]

        else: 

            sys.exit("No keyword 'HDL_HDL', 'LDL_LDL', 'HDL_LDL' found in the loaded RDF histogram file")

    return all_rdf_data  


def compute_pair_correl_ssf(keyword,
                            corrected,
                            partial,
                            same_species,
                            N1N2,
                            x_spec_1_2, 
                            N,
                            R_ranges,
                            Q_ranges,
                            num_bins_q,
                            num_bins_gr,
                            rho,
                            r_mid,
                            gr):

    corrected = c_bool(corrected)

    partial = c_bool(partial)

    same_species = c_bool(same_species)

    N1N2 = double_to_ctypes_double(N1N2)

    x_spec_1_2 = double_to_ctypes_double(x_spec_1_2)

    N = double_to_ctypes_double(N)  

    R_size = int_to_ctypes_int(R_ranges.size) 

    R_ranges = np_to_ctypes_array(R_ranges.astype(np.float64)) 
    
    Q_ranges = np_to_ctypes_array(Q_ranges.astype(np.float64))   

    num_bins_q = int_to_ctypes_int(num_bins_q) 

    num_bins_gr = int_to_ctypes_int(num_bins_gr) 
    
    rho = double_to_ctypes_double(rho) 
    
    r_mid = np_to_ctypes_array(r_mid) 

    gr = np_to_ctypes_array(gr) 

    if (keyword == "pair_correl"): 

        r2hr = np.ctypeslib.as_ctypes(np.zeros(num_bins_gr.value, dtype = np.float64))

        ssf_lib.call_hr_fs_correct(byref(corrected),
                                   byref(partial),
                                   byref(same_species),  
                                   byref(N1N2),
                                   byref(x_spec_1_2), 
                                   byref(N),
                                   byref(R_size),
                                   R_ranges,
                                   Q_ranges,
                                   byref(num_bins_q),
                                   byref(num_bins_gr),
                                   byref(rho),
                                   r_mid, 
                                   gr, 
                                   r2hr)

        return np.ctypeslib.as_array(r_mid), np.ctypeslib.as_array(r2hr)

    elif (keyword == "structure_factor"):

        SQ = np.ctypeslib.as_ctypes(np.zeros(num_bins_q.value, dtype=np.float64)) 
    
        Q_mid = np.ctypeslib.as_ctypes(np.zeros(num_bins_q.value, dtype=np.float64)) 
       
        S0_fs = c_double() 
 
        ssf_lib.call_ssf_fs_correct(byref(corrected), 
                                   byref(partial), 
                                   byref(same_species), 
                                   byref(x_spec_1_2),  
                                   byref(N), 
                                   byref(R_size),  
                                   R_ranges, 
                                   Q_ranges,  
                                   byref(num_bins_q), 
                                   byref(num_bins_gr), 
                                   byref(rho),  
                                   r_mid, 
                                   gr, 
                                   byref(S0_fs), 
                                   Q_mid,
                                   SQ)
        
        return S0_fs.value, np.ctypeslib.as_array(Q_mid), np.ctypeslib.as_array(SQ)  


def compute_pair_correl_from_hist_file(file_lst,
                                       ssf_or_pair,
                                       corrected,
                                       R_ranges,
                                       Q_ranges,
                                       num_bins_q):
    r_or_Q = [] 

    pair_correl_or_ssf = []

    S0_fs_lst = [] 

    for rdf_file in file_lst:
        
        RDF_hist = IO.reader.np_loadtxt(rdf_file)
        
        N, rho, num_bins_gr, r_mid, gr, r2hr = normalize_from_hist_file(RDF_hist)
        
        N1N2 = N * N 

        if ("LDL_LDL" in rdf_file):
        
            partial = True  
    
            same_species = True
    
            # the density from histogram file 
            # is already equal to rho * sqrt(x1*x2)  
            x_spec_1_2 = 1   

        elif ("HDL_HDL" in rdf_file):

            partial = True 

            same_species = True 

            # the density from histogram file 
            # is already equal to rho * sqrt(x1*x2)
            x_spec_1_2 = 1  

        elif ("HDL_LDL" in rdf_file):

            partial = True

            same_species = False

            # the density from histogram file 
            # is already equal to rho * sqrt(x1*x2)  
            x_spec_1_2 = 1  

        elif ("All" in rdf_file):

            partial = False

            same_species = True

            # the density from histogram file 
            # is already equal to rho * sqrt(x1*x2)  
            x_spec_1_2 = 1 

        else:

            sys.exit("keyword not found in %s file" % rdf_file)

        if (ssf_or_pair == "pair_correl"): 

            r_mid, r2hr = compute_pair_correl_ssf(
                                                  "pair_correl",
                                                  corrected,
                                                  partial,
                                                  same_species,
                                                  N1N2,
                                                  x_spec_1_2, 
                                                  N,
                                                  R_ranges, 
                                                  Q_ranges,
                                                  num_bins_q,
                                                  num_bins_gr,
                                                  rho,
                                                  r_mid,
                                                  gr)

            r_or_Q.append(r_mid) 

            pair_correl_or_ssf.append(r2hr)  

        elif (ssf_or_pair == "structure_factor"):
            
            S0_fs, Q_mid, S_Q = compute_pair_correl_ssf(
                                                 "structure_factor",
                                                  corrected,
                                                  partial,
                                                  same_species, 
                                                  N1N2,
                                                  x_spec_1_2,  
                                                  N,
                                                  R_ranges, 
                                                  Q_ranges,
                                                  num_bins_q,
                                                  num_bins_gr,
                                                  rho,
                                                  r_mid,
                                                  gr)

            
            S0_fs_lst.append(S0_fs) 

            r_or_Q.append(Q_mid) 

            pair_correl_or_ssf.append(S_Q)  

    if (ssf_or_pair == "pair_correl"):

        return r_or_Q, pair_correl_or_ssf 

    elif (ssf_or_pair == "structure_factor"):

        return S0_fs_lst, r_or_Q, pair_correl_or_ssf


def compute_static_sf(file_lst, R_ranges, Q_ranges, num_bins_q):

    corrected = True

    num_bins_q = int_to_ctypes_int(num_bins_q)     

    ssf_or_pair = "structure_factor" 
 
    r_or_Q, pair_correl_or_ssf = compute_pair_correl_from_hist_file(file_lst,
                                                                    ssf_or_pair,
                                                                    corrected,
                                                                    R_ranges,
                                                                    Q_ranges) 
    
    return r_or_Q[0], pair_correl_or_ssf[0]
    

def compute_Bhatia_Thornton_ssf(file_lst, R_ranges, Q_ranges, num_bins_q, x_H_frac, x_L_frac):

    corrected = True 

    ssf_or_pair = "structure_factor"  

    x_H_frac = c_double(x_H_frac) 

    x_L_frac = c_double(x_L_frac) 

    SNN = np.ctypeslib.as_ctypes(np.zeros(num_bins_q, dtype = np.float64)) 

    SCC = np.ctypeslib.as_ctypes(np.zeros(num_bins_q, dtype = np.float64)) 

    SNC = np.ctypeslib.as_ctypes(np.zeros(num_bins_q, dtype = np.float64)) 

    theta = np.ctypeslib.as_ctypes(np.zeros(num_bins_q, dtype = np.float64)) 

    SA = np.ctypeslib.as_ctypes(np.zeros(num_bins_q, dtype = np.float64)) 

    num_bins_q = int_to_ctypes_int(num_bins_q)     
 
    S0_fs_lst, r_or_Q, pair_correl_or_ssf = compute_pair_correl_from_hist_file(file_lst,
                                                                               ssf_or_pair,
                                                                               corrected,
                                                                               R_ranges,
                                                                               Q_ranges,
                                                                               num_bins_q)
    SHH, SLL, SHL = pair_correl_or_ssf

    SHH = np.ctypeslib.as_ctypes(SHH) 

    SLL = np.ctypeslib.as_ctypes(SLL) 

    SHL = np.ctypeslib.as_ctypes(SHL) 

    ssf_lib.call_Bhatia_Thornton(byref(num_bins_q),
                                 byref(x_H_frac),
                                 byref(x_L_frac),
                                 SHH,
                                 SLL,
                                 SHL,
                                 SNN,
                                 SCC,
                                 SNC,
                                 theta,
                                 SA)
 
    return (np.ctypeslib.as_array(r_or_Q[0]),  
            np.ctypeslib.as_array(SNN),
            np.ctypeslib.as_array(SCC),
            np.ctypeslib.as_array(SNC),
            np.ctypeslib.as_array(theta),
            np.ctypeslib.as_array(SA))


def compute_the_lorenzian(SA_ranges, SCC_ranges, Q_mid, SA, SCC):

    num_SA = c_int(SA_ranges.size) 

    num_SCC = c_int(SCC_ranges.size) 

    Q_SA_selected = Q_mid[SA_ranges]

    Q_SCC_selected = Q_mid[SCC_ranges]

    correl_SA = c_double()

    correl_SCC = c_double()

    SA = np_to_ctypes_array(SA[SA_ranges].astype(np.float64))

    SCC = np_to_ctypes_array(SCC[SCC_ranges].astype(np.float64))

    Q_SA_selected = np_to_ctypes_array(Q_SA_selected)

    Q_SCC_selected = np_to_ctypes_array(Q_SCC_selected)

    bt_lib.call_SA_SCC_lorentzian(byref(num_SCC),
                                  Q_SCC_selected,
                                  SCC,
                                  byref(num_SA),
                                  Q_SA_selected,
                                  SA,
                                  byref(correl_SA),
                                  byref(correl_SCC))  

    return correl_SA.value, correl_SCC.value 


