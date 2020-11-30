# Standard library:
import time
import numpy as np

# Local library:
from mathlib.math_API import (compute_spherical_harmonics, 
                              compute_optimized_Y12,
                              compute_optimized_Y8,
                              compute_optimized_Y6, 
                              compute_optimized_Y4)

# Third-party library
from pytest import approx
from scipy.special import sph_harm

def scipy_spherical_harmonics(l,phi, theta):
    # Note: 
    # theta  is  [0, Pi] 
    # phi is [0, 2pi] 

    sph_array = np.zeros(2*l+1, dtype=np.complex128)

    real_img = np.zeros((2*l+1,2), dtype=np.float64)

    for indx, m in enumerate(range(-l, l+1)): 
    
        sph_array[indx] = sph_harm(m,l,phi,theta)  

    real_img[:,0] = sph_array.real 

    real_img[:,1] = sph_array.imag 

    return real_img 

def test_any_spherical_harmonics():

    # test order range for.
    l = [1,3,4,6,8,10,12]

    for order in l: 
        
        # theta  is  [0, Pi] 
        theta_ary = np.array([0,0.25,  0.5, 1, 2.5, 3.14], dtype=np.float64) 

        # phi  is  [0, 2Pi]
        phi_ary = np.array([0, 0.25, 0.5, 1, 2.5, 3.14, 5, 6.28], dtype=np.float64)

        for theta in theta_ary: 

            for phi in phi_ary:  

                results = compute_spherical_harmonics(order, np.cos(theta), np.cos(phi), np.sin(phi)) 
        
                if (order == 12): 
            
                    optimized_12_results = compute_optimized_Y12(np.sin(theta), np.cos(theta), np.cos(phi), np.sin(phi)) 
                    
                    assert approx(results, 10**-24)  == optimized_12_results

                elif (order == 4):

                    optimized_4_results = compute_optimized_Y4(np.sin(theta), np.cos(theta), np.cos(phi), np.sin(phi))  
                    
                    assert approx(results, 10**-24)  == optimized_4_results 

                else: 

                    scipy_results = scipy_spherical_harmonics(order, phi, theta)     

                    assert approx(results, 10**-24)  == scipy_results

    return None

def test_optimized_Y12():

    # theta  is  [0, Pi] 
    theta_ary = np.array([0,0.25,  0.5, 1, 2.5, 3.14], dtype=np.float64) 

    # phi  is  [0, 2Pi]
    phi_ary = np.array([0, 0.25, 0.5, 1, 2.5, 3.14, 5, 6.28], dtype=np.float64)

    for theta in theta_ary: 

        for phi in phi_ary:  

            results = compute_optimized_Y12(np.sin(theta), np.cos(theta), np.cos(phi), np.sin(phi))

            scipy_results = scipy_spherical_harmonics(12, phi, theta)
            
            assert approx(results, 10**-24)  == scipy_results
    
    return None

def test_optimized_Y8():

    # theta  is  [0, Pi] 
    theta_ary = np.array([0, 0.25, 0.5, 1, 2.5, 3.14], dtype=np.float64) 

    # phi  is  [0, 2Pi] 
    phi_ary = np.array([0, 0.25, 0.5, 1, 2.5, 3.14, 5, 6.28], dtype=np.float64)

    for theta in theta_ary: 

        for phi in phi_ary:  

            results = compute_optimized_Y8(np.sin(theta), np.cos(theta), np.cos(phi), np.sin(phi))

            scipy_results = scipy_spherical_harmonics(8, phi, theta)
            
            assert approx(results, 10**-24)  == scipy_results
    
    return None 

def test_optimized_Y4():

    # theta  is  [0, Pi] 
    theta_ary = np.array([0, 0.25, 0.5, 1, 2.5, 3.14], dtype=np.float64) 

    # phi  is  [0, 2Pi] 
    phi_ary = np.array([0, 0.25, 0.5, 1, 2.5, 3.14, 5, 6.28], dtype=np.float64)

    for theta in theta_ary: 

        for phi in phi_ary:  

            results = compute_optimized_Y4(np.sin(theta), np.cos(theta), np.cos(phi), np.sin(phi))

            scipy_results = scipy_spherical_harmonics(4, phi, theta)
            
            assert approx(results, 10**-24)  == scipy_results
    
    return None 

def test_optimized_Y6():

    # theta  is  [0, Pi] 
    theta_ary = np.array([0, 0.25, 0.5, 1, 2.5, 3.14], dtype=np.float64) 

    # phi  is  [0, 2Pi] 
    phi_ary = np.array([0, 0.25, 0.5, 1, 2.5, 3.14, 5, 6.28], dtype=np.float64)

    for theta in theta_ary: 

        for phi in phi_ary:  

            results = compute_optimized_Y6(np.sin(theta), np.cos(theta), np.cos(phi), np.sin(phi))

            scipy_results = scipy_spherical_harmonics(6, phi, theta)

            assert approx(results, 10**-24)  == scipy_results
    
    return None 




test_optimized_Y12()
test_optimized_Y4()
test_optimized_Y6()
test_optimized_Y8()
test_any_spherical_harmonics()
