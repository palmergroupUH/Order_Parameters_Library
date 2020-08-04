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
from order_parameter.spherical_harmonics_based import initialize_Ql_Wl
from order_parameter.spherical_harmonics_based import call_Ql_Wl

# Third-party:
import numpy as np
import matplotlib.pyplot as plt  

def test_Ql_Wl(xyzfile):

    total_atoms, total_frames = IO.reader.call_read_header(xyzfile)
    
    # order of spherical harmonics:

    l = 4

    nnb = 0

    maxnb = 60

    cutoff = 1.0

    cutoff_sqr = cutoff * cutoff

    num_pairs_w3j, wigner3j_symobl_ary, m_index_ary, sph_const, Plm_const = initialize_Ql_Wl(l)
    
    for i in range(1):

        xyz, box = IO.reader.call_read_traj(xyzfile, i+1, total_atoms, return_numpy=False, xyz_keyword={"read_box": True})
       
        call_Ql_Wl(num_pairs_w3j,
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
                   cutoff_sqr)
   
    return None 

# --------------------------- My Q4, Q6, W4, W6 ----------------------
#           Q4                    Q6                      W4                    W6 
# fcc 0.190940653956492    0.574524259714069      -0.159317373133081     -1.316060073064693E-002 
# bcc 0.08832 (nnb=12)     0.538048050351382 (nnb=12) 0.159317373133081 (nnb=14)  +1.316060073064693E-002 (nnb=14 or nnb=8) 
# sc  0.763762615825974    0.353553390593275          0.159317373133081   1.316060073064693E-002

# References: 

# [1]:  P. J. Steinhardt, D. R. Nelson, and M. Ronchetti.
# Bond-orientational order inliquids and glasses.
# Phys. Rev. B, 28(2):784, 1983 


# [2]:  Y.  Wang,  S.  Teitel,  and  C.  Dellago.
# Melting  of  icosahedral  gold  nanoclus-ters from molecular dynamics simulations.
# The Journal of Chemical Physics,122(21):214722, 2005


# [3]: I. Saika-Voivod, P. H. Poole, and R. K. Bowles. 
# Test of classical nucleation theoryon deeply supercooled high-pressure simulated silica.
# The Journal of ChemicalPhysics, 124(22):224709, 2006. 

test_Ql_Wl("bcc.xyz" )
