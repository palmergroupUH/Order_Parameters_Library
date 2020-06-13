# Python standard library: 
import numpy as np  
import math

# Local library: 

# fortran API: 

from potential import potential_lj 

# Third-party libraries: 

def randomly_translate_a_particle(xyz,total_atoms,max_disp): 

    # pick a random particle 

    selected_particle = int(np.random.uniform()*total_atoms) 

    # find its current coordinate 

    old_coord  = xyz[selected_particle,:] 

    # give the coordinate a random displacement  
    # factor "2" is make the translate to be within (-1,1) 
    # trial translate vector is 3 dimensions 
    perturb_rand_vec_3d = 2*(np.random.rand(3)-0.5)*max_disp

    return selected_particle, old_coord, perturb_rand_vec_3d + old_coord  

def change_of_potential(xyz,box,cutoff,old_coord,new_coord): 

    # compute energy at old configuration  

    pe_old = compute_pe_due_to_one(xyz,box,old_coord,cutoff) 
     
    # compute energy at new configuration  

    pe_new = compute_pe_due_to_one(xyz,box,new_coord,cutoff) 

    return pe_new - pe_old 

def perform_trial_translate(chance):  

    if (np.random.uniform() < chance):  

        return True  

    else: 

        return False 
 

