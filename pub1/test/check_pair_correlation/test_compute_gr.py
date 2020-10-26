# Python standard library
import os
import sys
from ctypes import CDLL, c_int, c_double, c_long, c_float, byref

# Local library:
import IO.reader
from IO.type_conversion import string_to_ctypes_string,\
                               int_to_ctypes_int,\
                               np_to_ctypes_array

# Perform Bhatia Thornton anlysis
from pub1.Bhatia_Thornton import load_multiple_hist_file,\
                                 compute_static_sf,\
                                 compute_pair_correl_from_hist_file,\
                                 compute_Bhatia_Thornton_ssf,\
                                 compute_the_lorenzian 

# Third-party:
import numpy as np
import matplotlib.pyplot as plt



Temperature = [3375, 3400, 3425, 3450, 3500, 3550, 3600, 3700, 3800] 

subdir = ["LDL_LDL", "HDL_HDL", "HDL_LDL"] 
#subdir = ["All"]

R_ranges = np.linspace(4.28, 4.32, 5) 

Q_ranges = np.array([0, 50])

num_bins_q = 800 

x_H_frac = 0.5 

x_L_frac = 0.5 

S0_fs_lst_all_T = [] 

for temp in Temperature: 

    data_lst = [] 

    for folder in subdir: 

        directory = os.path.join("%s_K" % temp, folder)  

        all_files = list(os.walk(directory))

        feature = "RDF_HIST_" + folder

        for ifile in all_files[0][2]:

            if (feature in ifile):

               file_path = os.path.join(directory, ifile)
               
               data_lst.append(file_path)  

    """ 
    r_mid_lst, r2hr_lst = compute_pair_correl_from_hist_file(data_lst,
                                                             "pair_correl",
                                                             corrected=True,
                                                             R_ranges=R_ranges,
                                                             Q_ranges=Q_ranges,
                                                             num_bins_q = num_bins_q)

    
    S0_fs_lst, Q_lst, SQ_lst = compute_pair_correl_from_hist_file(data_lst,
                                                       "structure_factor",
                                                       corrected=True,
                                                       R_ranges=R_ranges,
                                                       Q_ranges=Q_ranges,
                                                       num_bins_q = num_bins_q)
  
    S0_fs_lst_all_T.extend(S0_fs_lst) 

    #np.savetxt("total_pair_correl/S_tot_%d_K.txt"%temp, np.c_[Q_lst[0], SQ_lst[0]])
    #np.savetxt("total_pair_correl/r2hr_all_%d_K.txt"%temp, np.c_[r_mid_lst[0], r2hr_lst[0]])
    #np.savetxt("total_pair_correl/inset_r2hr_all_%d_K.txt"%temp, np.c_[r_mid_lst[0][201:], r2hr_lst[0][201:]])
    """
    
    #np.savetxt("total_pair_correl/S0_tot.txt", np.c_[Temperature, S0_fs_lst_all_T])

    Q_mid, SNN, SCC, SNC, theta, SA = compute_Bhatia_Thornton_ssf(data_lst, R_ranges, Q_ranges, num_bins_q, x_H_frac, x_L_frac) 
    np.savetxt("bhatia/SNN_%d_K.txt"%temp, np.c_[Q_mid, SNN]) 
    np.savetxt("bhatia/SCC_%d_K.txt"%temp, np.c_[Q_mid, SCC]) 
    np.savetxt("bhatia/SNC_%d_K.txt"%temp, np.c_[Q_mid, SNC]) 
    np.savetxt("bhatia/theta_%d_K.txt"%temp, np.c_[Q_mid, theta]) 
    np.savetxt("bhatia/SA_%d_K.txt"%temp, np.c_[Q_mid, SA]) 
     
    Q_mid_sqr = Q_mid**2

    QCC = np.arange(23,28) 
    
    QA = np.arange(23,30)  

    correl_SA, correl_SCC =compute_the_lorenzian(QA, QCC, Q_mid, SA, SCC)     

    print (correl_SA, correl_SCC) 

