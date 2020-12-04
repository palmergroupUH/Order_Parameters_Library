# Python standard library
import os
import sys
#from oplib import * 
from oplib_q12_q6 import * 
from ctypes import CDLL, c_int, c_double, c_long, c_float, byref

# Local library:
import IO.reader
from IO.type_conversion import string_to_ctypes_string,\
                               int_to_ctypes_int,\
                               np_to_ctypes_array

# import HMC related functions:
from order_parameter.op_lib import calc_ql_cluster 
from order_parameter.op_lib import calc_q3_cluster_HMC 

# Third-party:
import numpy as np

def test_HMC_ReinhardtDoye():

    dcdfile = "test_traj.dcd"

    total_atoms, total_frames = IO.reader.call_read_header(dcdfile)

    maxnb = 40

    nnb = 0

    l = 3

    ql_cutoff = 3.6

    cluster_cutoff = ql_cutoff 

    n_bonds = 3

    # The bond criterion must be modified in the source code
    style = "ReinhardtDoye"

    for i in range(total_frames):

        xyz, box = IO.reader.call_read_dcd_xyz_box(dcdfile, i+1, total_atoms, return_numpy=False)
        
        nxtl, largest_cluster = calc_ql_cluster(style, total_atoms, box, xyz, maxnb, nnb, l, ql_cutoff, cluster_cutoff, n_bonds)
        
        box = np.ctypeslib.as_array(box) 
        
        xyz = np.ctypeslib.as_array(xyz)
        
        nstar = calc_q3_cluster_HMC(n_atoms=total_atoms, n_q3_neigh=0,q3_cutoff=3.60,n_bonds=3,box=box,x=xyz) 

        assert largest_cluster == nstar

    return None 

def test_HMC_RussoTanaka():
    
    dcdfile = "test_traj.dcd"

    total_atoms, total_frames = IO.reader.call_read_header(dcdfile)

    maxnb = 60

    nnb = 16

    l = 12

    ql_cutoff = 6 

    cluster_cutoff = 3.6 

    n_bonds = 12

    # The bond criterion must be modified in the source code
    bond_crit = "RussoTanaka"

    for i in range(total_frames): 
        
        xyz, box = IO.reader.call_read_dcd_xyz_box(dcdfile, i+1, total_atoms, return_numpy=False)
        
        nxtl, largest_cluster_lib = calc_ql_cluster(bond_crit, total_atoms, box, xyz, maxnb, nnb, l, ql_cutoff, cluster_cutoff, n_bonds)
    
        box = np.ctypeslib.as_array(box) 
        
        xyz = np.ctypeslib.as_array(xyz)

        nxtl_new, largest_cluster_q12 = calc_q12_cluster(n_q12_neigh=16, q12_cutoff=6, crys_cutoff=3.6, n_bonds=12, box=box, x=xyz, n_atoms=total_atoms)  

        print (largest_cluster_lib) 

        assert largest_cluster_lib == largest_cluster_q12  

    return None 

test_HMC_ReinhardtDoye()

test_HMC_RussoTanaka()

