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
from order_parameter.spherical_harmonics_based import intialize_CHILL
from order_parameter.spherical_harmonics_based import call_CHILL

# Third-party:  
import numpy as np
import matplotlib.pyplot as plt  

def test_CHILL():

    dcdfile = "/project/palmer/Jingxiang/custom/order_parameter/data/CHILL/traj_T096_32000_11sq2.dcd"

    #dcdfile = "/project/palmer/Jingxiang/ours_optimization/tutorial/Tutorial_04_preparation/ReferenceData/rdf/mW_300K_1bar/traj.dcd"

    total_atoms, total_frames = IO.reader.call_read_header(dcdfile)

    # order of spherical harmonics:

    l = 3

    nnb = 4

    maxnb = 40

    cutoff = 4.8

    cutoff_sqr = cutoff * cutoff

    sph_const, Plm_const = intialize_CHILL(l)

    num_bins = 200

    hist = np.zeros(num_bins)

    CHILL_keyword = "CHILL"

    for i in range(3000):

        xyz, box = IO.reader.call_read_dcd_xyz_box(dcdfile, i+1, total_atoms, return_numpy=False)

        cij, chill_id_list = call_CHILL(CHILL_keyword, sph_const, Plm_const, total_atoms, l, xyz, box, maxnb, nnb, cutoff_sqr)
   
        print( np.count_nonzero(chill_id_list !=4))

        hist_each, bin_edges = np.histogram(cij, bins=num_bins, range=(0,1))

        hist += hist_each 

    hist = hist/np.sum(hist)/(bin_edges[1] - bin_edges[0])

    mid_point = (bin_edges[:-1] + bin_edges[1:])/2.0

    #plt.semilogy(mid_point, hist)

    #plt.ylim([0.03,100])

    #plt.xlim([0,1])

    #plt.show()

    return None 

test_CHILL()
