# Standard library:
import time
import numpy as np

# Local library:
from mathlib.math_API import compute_wigner_3j 

# Third-party library
from pytest import approx 

def compute_num_paris(l): 

    num_pairs = 0 

    m_lst = [] 

    for i in range(-l, l+1):

        for j in range(-l, l+1): 

            for k in range(-l, l+1): 
    
                if (i+j+k ==0):
        
                    num_pairs = num_pairs + 1  
                  
                    m_lst.append([i,j,k]) 
                     
    return num_pairs, m_lst  

def test_wigner_3j_order4(): 

    num_pairs, m_lst = compute_num_paris(l=4)

    selected = [[-4,4,0], [1,2,-3], [0,-1,1], [0, 0, 0], [-1,0,1], [3,1,-4]]

    # Mathmatica results: 
    # For example, ThreeJSymbol[{4,0}, {4,0}, {4,0}]
     
    Ref_mathmatica = [0.104297703129124,
                      -0.06232979933389715,
                      -0.06704852344015112,
                      0.1340970468803022,
                      -0.06704852344015112,
                      -0.1649091483060512    
                      ]

    for i, pairs in enumerate(selected):

        matrix = np.array([[4, 4, 4], pairs], dtype=np.int32)

        wigner3j = compute_wigner_3j(matrix) 

        assert approx(Ref_mathmatica[i], 10**-24)  == wigner3j
        
    return None 

def test_wigner_3j(): 

    l_range = [3, 6] 

    Ref_mathmatica = [[-0.1543033499620919,0,-0.1543033499620919], [0.05118272511620992, -0.0930595002112908, 0.04652975010564537 ]]

    selected = [[-2,2,0],[0, 0, 0], [-1,0,1]] 

    # Mathmatica results: 
    # For example, ThreeJSymbol[{4,0}, {4,0}, {4,0}]

    for i_l, l in enumerate(l_range):
     
        for i, pairs in enumerate(selected):

            matrix = np.array([[l, l, l], pairs], dtype=np.int32)

            wigner3j = compute_wigner_3j(matrix) 
            
            assert approx(Ref_mathmatica[i_l][i], 10**-24)  == wigner3j

    return None 


def test_wigner_3j_test(): 

    l_range = [4] 

    Ref_mathmatica = [[-0.1543033499620919,0,-0.1543033499620919], [0.05118272511620992, -0.0930595002112908, 0.04652975010564537 ]]

    selected = [[2,-4,2],[3, 1, -4], [0,-1,1]] 

    # Mathmatica results: 
    # For example, ThreeJSymbol[{4,0}, {4,0}, {4,0}]

    for i_l, l in enumerate(l_range):
     
        for i, pairs in enumerate(selected):

            matrix = np.array([[l, l, l], pairs], dtype=np.int32)

            wigner3j = compute_wigner_3j(matrix) 
            
            print (wigner3j)  

    return None 

test_wigner_3j_test()
