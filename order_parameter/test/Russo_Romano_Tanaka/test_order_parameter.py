# Python standard library
import os
import sys
from ctypes import CDLL, c_int, c_double, c_long, c_float, byref

# Local library:
import IO.reader
from IO.type_conversion import string_to_ctypes_string,\
                               int_to_ctypes_int,\
                               np_to_ctypes_array

# import RussoRomanoTanaka related functions
from order_parameter.spherical_harmonics_based import intialize_RussoRomanoTanaka  
from order_parameter.spherical_harmonics_based import call_RussoRomanoTanaka 

# Third-party:  
import numpy as np
import matplotlib.pyplot as plt  

def histogram(bins_range, num_bins):

    interval = (bins_range[-1] - bins_range[0])/num_bins

    bins_range_hist = np.array([bins_range[0] - 0.5*interval, bins_range[-1] + 0.5*interval])

    return np.zeros(num_bins+1), bins_range_hist

def test_RussoRomanoTanaka():

    #dcdfile = "../data/RussoRomenoTanaka/Ih/187K/traj.dcd"
    #dcdfile = "/project/palmer/Jingxiang/custom/order_parameter/data/RussoRomenoTanaka/liquid/traj.dcd"
    dcdfile = "/project/palmer/Jingxiang/custom/order_parameter/data/RussoRomenoTanaka/Ih/187K/traj.dcd"
    #dcdfile = "/project/palmer/Jingxiang/ours_optimization/tutorial/Tutorial_04_preparation/ReferenceData/rdf/mW_300K_1bar/traj.dcd"

    total_atoms, total_frames = IO.reader.call_read_header(dcdfile)

    nnb = 16

    maxnb = 40

    cutoff = 5.8
    
    connect_cut = 0.75
    
    crys_cut = 12

    cutoff_sqr = cutoff * cutoff

    num_pairs_w3j, wigner3j_symobl_vals, m_index_ary = intialize_RussoRomanoTanaka()

    num_bins = 80

    hist, bins_range = histogram([0, 1], num_bins)

    for i in range(2500,total_frames):

        xyz, box = IO.reader.call_read_dcd_xyz_box(dcdfile, i+1, total_atoms, return_numpy=False)

        cij = call_RussoRomanoTanaka(total_atoms, xyz, box, maxnb, nnb, cutoff_sqr, connect_cut, crys_cut)
    
        hist_each, bin_edges = np.histogram(cij, bins=num_bins+1, range=(bins_range[0], bins_range[-1]))
        
        hist += hist_each  

    hist = hist/sum(hist)/(bin_edges[1]-bin_edges[0])

    mid = (bin_edges[1:] + bin_edges[:-1])/2.0
    
    np.savetxt("Ih_dot_product.txt", np.c_[mid, hist])

    return None 
    
test_RussoRomanoTanaka()
