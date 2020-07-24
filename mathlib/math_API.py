# Python standard library
import numpy as np
import multiprocessing as mp
import os
import sys
import itertools
import logging
from ctypes import CDLL, c_int, c_double, c_long, c_float, byref

# Local library:
import mathlib
from IO.type_conversion import string_to_ctypes_string,\
                               int_to_ctypes_int,\
                               np_to_ctypes_array

# Third-party library:

fortranlib_address = os.path.join(os.path.dirname(mathlib.__file__), "lib")

sph_lib = CDLL(os.path.join(fortranlib_address, "libsph.so"))

def compute_associated_legendre_poly(l, m, cos_theta):

    preconst = np.ctypeslib.as_ctypes(np.zeros((l+1)*4, dtype=np.float64))

    Plm = np.ctypeslib.as_ctypes(np.zeros((2*l+1), dtype=np.float64))

    l = c_int(l)

    sph_lib.call_associated_legendre_const(byref(l), preconst)

    cos_theta = c_double(cos_theta) 
    
    sph_lib.call_associated_legendre_poly(preconst,byref(cos_theta),byref(l),Plm)  

    return None 

def compute_spherical_harmonics(l,cos_theta, cos_phi, sin_phi, m=None):

    # 2 is for real and imaginary part a complex vector
   
    l_dim = 2 * l + 1 

    sph_const = np.ctypeslib.as_ctypes(np.zeros(l_dim, dtype=np.float64))

    plm_const  = np.ctypeslib.as_ctypes(np.zeros((l+1)*4, dtype=np.float64))

    l = c_int(l)

    sph_lib.call_spherical_harmonics_const(byref(l), sph_const)

    sph_lib.call_associated_legendre_const(byref(l), plm_const)

    cos_theta = c_double(cos_theta)

    cos_phi = c_double(cos_phi)  

    sin_phi = c_double(sin_phi)

    Ylm_complex = np.ctypeslib.as_ctypes(np.zeros(l_dim*2, dtype=np.float64))
    
    sph_lib.call_spherical_harmonics(sph_const,
                                     plm_const,
                                     byref(l),
                                     byref(cos_theta),
                                     byref(cos_phi),
                                     byref(sin_phi),
                                     Ylm_complex)  

    Ylm_complex = np.ctypeslib.as_array(Ylm_complex).reshape(2,l_dim).T

    # if m is provided, return m value
    if (m is not None):

        index = m + l.value
    
        # unpack 2d arrays into real and img part
        return complex(*Ylm_complex[index])

    # if m is not provided, return Ylm 
    else:

        return Ylm_complex

def compute_optimized_Y12(sin_theta, cos_theta, cos_phi, sin_phi, m=None):

    # 2 is for real and imaginary part a complex vector
   
    cos_theta = c_double(cos_theta)

    sin_theta = c_double(sin_theta)

    cos_phi = c_double(cos_phi) 

    sin_phi = c_double(sin_phi) 

    Ylm_complex = np.ctypeslib.as_ctypes(np.zeros(25*2, dtype=np.float64))

    sph_lib.call_optimized_Y12(byref(sin_theta),
                               byref(cos_theta),
                               byref(cos_phi), 
                               byref(sin_phi), 
                               Ylm_complex)

    Ylm_complex = np.ctypeslib.as_array(Ylm_complex).reshape(2,25).T
    
    # if m is provided, return m value
    if (m is not None):

        index = m + l.value
    
        # unpack 2d arrays into real and img part
        return complex(*Ylm_complex[index])

    # if m is not provided, return Ylm 
    else:

        return Ylm_complex  

