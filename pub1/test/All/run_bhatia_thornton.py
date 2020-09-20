# Python standard library
import os
import sys
from ctypes import CDLL, c_int, c_double, c_long, c_float, byref

# Local library:
import IO.reader
from IO.type_conversion import string_to_ctypes_string,\
                               int_to_ctypes_int,\
                               np_to_ctypes_array

from pub1.Bhatia_Thornton import compute_q_tetra, initialize_hist,\
                                 update_pair_correl_hist, normalize_hist, \
                                 compute_pair_correl, update_hist, \
                                 normalize_pair_correl_hist  

# Third-party:
import numpy as np
import matplotlib.pyplot as plt



def compute_high_low_q_correlation(): 

    dcdfile = "/project/palmer/Jingxiang/Trajectories/Publication_mWAC/3375_run2.dcd"    

    total_atoms, total_frames = IO.reader.call_read_header(dcdfile) 

    keyword = "All"

    # convert Angstrom to nm
    An_to_nm = 10

    num_bins_q_dist = 300

    num_bins_rdf = 250

    # initialize q tetrahedral histogram 
    q_tetra_hist, q_interval, q_range = initialize_hist(num_bins_q_dist, -3.0, 1.0)    
    
    # get box size
    xyz, box = IO.reader.call_read_dcd_xyz_box(dcdfile, 1, total_atoms, return_numpy=True) 

    half_box_size  = box.min()/2.0/An_to_nm 
    
    # initialize pair_correlation histogram
    RDF_hist, rdf_interval, rdf_range = initialize_hist(num_bins_rdf, 0.0, half_box_size)  

    frame_counter = 0

    for i in range(total_frames): 
    
        # read xyz coordinates
        xyz, box = IO.reader.call_read_dcd_xyz_box(dcdfile, i+1, total_atoms, return_numpy=False) 
       
        # compute q tetrahedral order parameter 
        q_sorted_indx, q_tetra = compute_q_tetra(total_atoms, xyz, box)

        # update q distribution
        q_tetra_hist = update_hist(total_atoms, q_tetra, q_interval, -3.0, 1.0, num_bins_q_dist, q_tetra_hist)

        # compute the pair correlation between species
        new_hist, norm_atom = compute_pair_correl(keyword, num_bins_rdf, total_atoms, xyz, box, q_sorted_indx)
    
        # update pair correlation distribution
        update_pair_correl_hist(num_bins_rdf, new_hist, RDF_hist)
    
        frame_counter += 1  

    r_mid, q_norm = normalize_hist(q_range, num_bins_q_dist, q_interval, q_tetra_hist)

    r2hr, gr = normalize_pair_correl_hist(RDF_hist, num_bins_rdf, half_box_size, norm_atom, frame_counter, box)

    np.savetxt("r2hr.txt", np.c_[r2hr])

    np.savetxt("q_normalized_3375.txt", np.c_[r_mid, q_norm])

    return None 

compute_high_low_q_correlation() 
