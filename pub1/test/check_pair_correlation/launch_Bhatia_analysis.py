# Python standard library
import os
import sys
import time 

# Local


# Third-Party libraries


# list of temperature (K) to be analyzed:
Temperature = [3375, 3400, 3425, 3450, 3500, 3550, 3600, 3700, 3800] 
#Temperature = [3375] 

# list of trajectories:
traj = ["3375_run2.dcd", "3400_run2.dcd", "3425_run2.dcd", "3450_run2.dcd", "3500_run2.dcd", "3550.dcd", "3600.dcd", "3700.dcd", "3800.dcd"]
#traj = ["3375_run2.dcd"]

# trajectories path:
traj_path = '/project/palmer/Jingxiang/Trajectories/Publication_mWAC/traj/'

keyword = '"All"'

HOME = os.getcwd() 

modified_file = "run_bhatia_thornton.py"

for temp, dcdfile in zip(Temperature, traj) : 

    trajpath = os.path.join(traj_path, dcdfile)

    trajpath = '"' + trajpath + '"'

    shell_cmd_1 = "sed -i '89c T = %d' " % temp 

    shell_cmd_2 = "sed -i '92c dcdfile = %s'" % trajpath 

    shell_cmd_3 = "sed -i '97c keyword = %s'" % keyword 

    # first level of directory
    wk_folder = "%d_K" % temp

    # second level
    wk_folder = os.path.join(wk_folder, keyword.replace('"', ''))  
    
    copy_cmd = "cp template/* %s" % wk_folder 
   
    if (not os.path.isdir(wk_folder)):

        os.makedirs(wk_folder)

    os.system(copy_cmd)  

    os.chdir(wk_folder)

    os.system(shell_cmd_1 + " " + modified_file) 

    os.system(shell_cmd_2 + " " + modified_file) 

    os.system(shell_cmd_3 + " " + modified_file) 

    os.system("sbatch run_script")

    os.chdir(HOME) 

    time.sleep(1)  

