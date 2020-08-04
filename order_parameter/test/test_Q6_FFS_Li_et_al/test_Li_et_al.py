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
from order_parameter.spherical_harmonics_based import intialize_Li_et_al
from order_parameter.spherical_harmonics_based import call_Li_et_al

# Third-party:
import numpy as np
import matplotlib.pyplot as plt  

def histogram(bins_range, num_bins):

    interval = (bins_range[-1] - bins_range[0])/num_bins

    bins_range_hist = np.array([bins_range[0] - 0.5*interval, bins_range[-1] + 0.5*interval])

    return np.zeros(num_bins+1), bins_range_hist

def test_Li_et_al():

    dcdfile = "/project/palmer/Jingxiang/custom/order_parameter/data/Q6_FFS_Li_et_al/liquid_240K/traj.dcd"

    #dcdfile = "/project/palmer/Jingxiang/custom/order_parameter/data/CHILL/CHILL_plus/Ih/traj.dcd"

    total_atoms, total_frames = IO.reader.call_read_header(dcdfile)

    # order of spherical harmonics:

    l = 6

    nnb = 0

    maxnb = 40

    cutoff = 3.4

    cutoff_sqr = cutoff * cutoff

    sph_const, Plm_const = intialize_Li_et_al(l)

    num_bins = 80

    hist, bins_range = histogram([-1, 1], num_bins)
    
    for i in range(500, total_frames):
        
        xyz, box = IO.reader.call_read_dcd_xyz_box(dcdfile, i+1, total_atoms, return_numpy=False)
    
        Iij = call_Li_et_al(sph_const, Plm_const, total_atoms, maxnb, nnb, l, cutoff_sqr,  box, xyz)
        
        hist_each, bin_edges = np.histogram(Iij, bins=num_bins+1, range=(bins_range[0], bins_range[-1]))      
        
        hist += hist_each  
    
    hist = hist/sum(hist)/(bin_edges[1]-bin_edges[0])

    mid = (bin_edges[1:] + bin_edges[:-1])/2.0
    
    np.savetxt("Iij_test_liquid.txt", np.c_[mid, hist])

    #plt.plot(mid, hist)  

    #plt.show()

    return None 

test_Li_et_al()
