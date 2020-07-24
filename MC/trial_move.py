# Python standard library: 
import numpy as np  
import math

# Local library: 

# fortran API: 


# Third-party libraries: 

def translate(particle_coord, max_disp): 

    # give the coordinate a random displacement  
    # factor "2" is make the translate to be within (-1,1) 
    # trial translate vector is 3 dimensions 
    perturb_rand_vec_3d = 2*(np.random.rand(3)-0.5)*max_disp

    # apply periodical bounardy condition and minimum image convention  

    new_coord = perturb_rand_vec_3d - box*(perturb_rand_vec_3d/box) 

    return  old_coord, new_coord 

def delta_potential(xyz,box,cutoff,old_coord,new_coord): 

    # compute energy at old configuration  

    pe_old = compute_pe_due_to_one(xyz,box,old_coord,cutoff) 
     
    # compute energy at new configuration  

    pe_new = compute_pe_due_to_one(xyz,box,new_coord,cutoff) 

    return pe_new - pe_old 

def translate_chance(chance):  

    if (np.random.uniform() < chance):  

        return True  

    else: 

        return False 
 

