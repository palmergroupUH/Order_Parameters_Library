# Standard library:
import time
import numpy as np

# Local library:
from mathlib.math_API import (compute_spherical_harmonics, 
                              compute_optimized_Y12)

# Third-party library
import pytest


def test_any_spherical_harmonics():

    # test range for.
    #l = [1,3,4,6,8,10,12]
    l = [12] 

    mathmatica_results  = [complex(0.089495, -0.13938),
                           0.428789,
                           complex(-0.089495, -0.13938)]
    cos_theta = np.cos(0.5) 

    sin_theta = np.sqrt((1-cos_theta**2)) 

    phi = 1 

    cos_phi = np.cos(phi) 

    sin_phi = np.sin(phi)

    counter = 0 

    for order in l: 
        
        m_list = [-order, 0, order]

        for m in m_list:
            
            results = compute_optimized_Y12(sin_theta, cos_theta, cos_phi, sin_phi) 
            print (results) 
            results = compute_spherical_harmonics(order, cos_theta, cos_phi, sin_phi, m) 
            print (results)  
            #val = results - mathmatica_results[counter]
            
            #print (val.real < 10**-9 and val.imag < 10**-9)  

            #assert (val.real < 10**-6 and val.imag < 10**-6)

            counter += 1 

    return None

def test_optimized_Y12():

    theta = 0

    phi = 1 

    results = compute_optimized_Y12(np.sin(theta), np.cos(theta), np.cos(phi), np.sin(phi))

    #print (results) 
    theta = 1 

    results = compute_optimized_Y12(np.sin(theta), np.cos(theta), np.cos(phi), np.sin(phi))
    #print (results) 
    phi = 0 

    return None

test_any_spherical_harmonics()

#test_optimized_Y12()
