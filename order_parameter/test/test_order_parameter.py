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

def test_RussoRomanoTanaka():

    #dcdfile = "../data/RussoRomenoTanaka/Ih/187K/traj.dcd"
    dcdfile = "/project/palmer/Jingxiang/custom/order_parameter/data/RussoRomenoTanaka/traj.dcd"
    #dcdfile = "/project/palmer/Jingxiang/ours_optimization/tutorial/Tutorial_04_preparation/ReferenceData/rdf/mW_300K_1bar/traj.dcd"

    total_atoms, total_frames = IO.reader.call_read_header(dcdfile)

    # order of spherical harmonics:

    l = 12

    nnb = 16

    maxnb = 40

    cutoff = 5.8
    
    crys_cut = 0.75

    crys_bond = 12

    cutoff_sqr = cutoff * cutoff

    sph_const, Plm_const = intialize_RussoRomanoTanaka(l)

    num_bins = 200

    cij_hist = np.zeros(num_bins)

    connect_hist = np.zeros(nnb+1)

    for i in range(5000,10000):

        xyz, box = IO.reader.call_read_dcd_xyz_box(dcdfile, i+1, total_atoms, return_numpy=False)

        cij, count_bonds = call_RussoRomanoTanaka(sph_const, Plm_const, total_atoms, l, xyz, box, maxnb, nnb, cutoff_sqr, crys_cut)

        cij_each, bin_edges = np.histogram(cij, bins=num_bins, range=(0,1))
    
        connect_each, bin_edges = np.histogram(count_bonds, bins=nnb+1, range=(-0.5,nnb+0.5))

        print (bin_edges) 

        cij_hist += cij_each 

        connect_hist += connect_each 

    cij_hist_norm = cij_hist/np.sum(cij_hist)/(bin_edges[1] - bin_edges[0])

    connect_hist_norm = connect_hist/np.sum(connect_hist)/(bin_edges[1] - bin_edges[0])

    #for i in cij_hist_norm: 

        #print (i)

    for i in connect_hist_norm: 
        
        print (i) 

    #mid_point = (bin_edges[:-1] + bin_edges[1:])/2.0

    #plt.semilogy(mid_point, hist)

    #plt.ylim([0.03,100]) 

    #plt.xlim([0,1]) 

    #plt.show()

    return None 

test_RussoRomanoTanaka()
