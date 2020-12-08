from ctypes import *
import sys,random,math
from lammps import lammps
#from oplib_q12_q6 import * 
from order_parameter.op_lib import calc_ql_cluster, calc_Ql 

#from oplib import * 

#-------------------- Run Parameters --------------------#

# Thermodynamic state -  configured for LAMMPS' "real" units
T = 185.0             # Temperature [K]
P = 1579.077           # Pressure [atm]

# Umbrella sampling parameters
nxtl_max = 30
krho = 5000.0 
rhoc = 0.93

# MC parameters
Nsweeps  = 1000000      # Number of MC Sweeps (1 Sweep = 1 HMC or Volume Move)
Nsteps = 25              # Number of MD Steps per HMC trial
dt = 10.0                 # time step [fs]
prob_hmc = 0.8           # Probability of selecting HMC move (vs Volume Move); >= 1.0 for NVT ensemble       
max_lnvol = 0.040        # Maximum log volume displacement
rseed = 1234             # Random number seed (fixed for testing)

# File output frequency
freq_thermo = 50         #  Thermodynamic output
freq_traj = 5000         # XYZ trajectory
freq_restart = 5000      # Restart file (LAMMPS write_data)
freq_flush = 500         # Flush files

#-------------------- Physical Constants and Conversion Factors --------------------#
# Note: These should changes if LAMMPS' units change

# Constants.
kB = 1.380648520000e-23  # Boltzmann's constant [m^2 kg s^-2 K^-1]
Nav = 6.02214090000e23   # Avagadro's number [molecules mol^-1]
R = kB * Nav / 1000.0    # Gas constant [kJ mol^-1 K^-1]

# Thermal energy 
kTs = R*T                # SI units: [kJ mol^-1]
kTL = kTs/4.184          # LAMMPS units [real units, kcal mol^-1] 

# Velocity prefactor 
vf = 1e-4*R*T            # LAMMPS units [real units, g A^2 * fs^-2 mol^-1]            

# Pressure prefactor for volume change move
Pb = P*1.01325           # Pressure [bar]
Pc = kB*1.0e30*1.0e-5    # Conversion factor [bar A^3 K^-1]
Pf = Pb/(Pc*T)           # Prefactor for Metropolis criterion [A^-3]

# Density prefactor
df = 1.0/(Nav*1e-24)           # [ mol molecule^-1 A^3 cm^-3]

#-------------------- Initialization --------------------#

#Seed random number generator
random.seed(rseed)

# Initialize LAMMPS

lmp = lammps(name="",cmdargs=["-log","none","-screen","none"])

# Load script/data
lmp.file("in.nve")

# Set the integration time step
lmp.command("timestep %f" % dt)

# Define compute for kinetic energy and virial pressure
lmp.command("compute thermo_ke all ke")

# Get initial system properties
natoms = lmp.extract_global("natoms",0)
mass = lmp.extract_atom("mass",2)
atomid = lmp.gather_atoms("id",0,1)
atype = lmp.gather_atoms("type",0,1)

# Compute the total mass of the system, update density prefactor
mass_total = 0.0
for i in range(natoms):
  mass_total += mass[atype[i]]
df *= mass_total  # [g A^3/cm^3] 

# Allocate coordinte and velocity arrays 
x=(3*natoms*c_double)()
x_new = (3*natoms*c_double)()
v=(3*natoms*c_double)()

# Initialize properties
pe = 0.0
ke = 0.0
etot = 0.0
box = 0.0
vol = 0.0

# Initialize counters
n_acc_hmc = 0.0
n_try_hmc = 0.0
n_acc_vol = 0.0
n_try_vol = 0.0

# Get initial position and velociteies
x = lmp.gather_atoms("x",1,3)
v = lmp.gather_atoms("v",1,3)


# Compute initial PE [dimensionless]

pe = lmp.extract_compute("thermo_pe",0,0)/kTL

# Compute box length and volume (assumes cubic!)
boxlo = lmp.extract_global("boxxlo",1)
boxhi = lmp.extract_global("boxxhi",1)
box = boxhi- boxlo
vol = math.pow(box,3)


# Compute initial OPs 

#nxtl, largest_cluster = calc_q12_cluster(n_q12_neigh=16,q12_cutoff=6,crys_cutoff=3.5, n_bonds=12,box=box,x=x,n_atoms=natoms)
#print (nxtl, largest_cluster) 
    
nxtl, largest_cluster = calc_ql_cluster(bond_crit="RussoTanaka", n_atoms=natoms, box=box, x=x, maxnb=60, n_ql_neigh=16, l=12, ql_cutoff=6, cluster_cutoff=3.5, n_bonds=12) 


rho = df/vol

# Open files for writing
thermo = open('thermo_test_3.6.dat', 'w')
traj = open('traj.dat', 'w')

#-------------------- Support Functions --------------------#

# Velocity initialization
# -Draw initial velocities from Maxwell-Boltzmann distribution
# -Send velocities to LAMMPS
# -"Run 0" to set velocities internally
# WARNING: THIS ROUTINE ONLY WORKS FOR POINT PARTICLES.  DO NOT USE FOR RIGID BODIES.
def init_vel():
  for i in range(natoms):
    indx = 3*i
    sigma = math.sqrt(vf/mass[atype[i]])
    v[indx] =  random.gauss(0.0,sigma) 
    v[indx+1] = random.gauss(0.0,sigma)
    v[indx+2] = random.gauss(0.0,sigma)
  lmp.scatter_atoms("v",1,3,v)
  lmp.command("run 0")
  return


def calc_bias(move):
  bias = 0.0
  if(move == 0):
    if(nxtl_new > nxtl_max):
      bias = 600.0
  elif (move == 1): 
    bias = 0.5*krho*((rho_new - rhoc)*(rho_new - rhoc) - (rho - rhoc)*(rho - rhoc))
    if(nxtl_new > nxtl_max):
      bias = 600.0
  return bias

# Computes MC move acceptance ratio
def acc_ratio(acc,trys):
  if(trys == 0.0):
    ratio = 0.0
  else:
    ratio = acc/trys
  return ratio

# Debugging function - check for discrepancy between current MC energy and energy calculated from scratch
def check_eng(pe):
  lmp.scatter_atoms("x",1,3,x)
  lmp.command("run 0")
  pe_lmps = lmp.extract_compute("thermo_pe",0,0)/kTL
  err = math.fabs((pe-pe_lmps)/pe_lmps)
  if(err >= 1e-14): print (err, pe,pe_lmps)


#-------------------- MC --------------------#

for isweep in range(Nsweeps):
  
  # HMC Move
  if(random.random() <= prob_hmc): 

    # Update number of trial 
    n_try_hmc += 1.0

    # Scatter coordinates
    lmp.scatter_atoms("x",1,3,x)
    
    # Generate initial velocities; compute KE [dimensionless]
    init_vel()
    ke = lmp.extract_compute("thermo_ke",0,0)/kTL
    etot = pe + ke
   
    # Run Nsteps MD steps
    lmp.command("run %d" % Nsteps)
   
    # Compute new PE, KE, and total energy [dimensionless]
    pe_new = lmp.extract_compute("thermo_pe",0,0)/kTL
    ke_new = lmp.extract_compute("thermo_ke",0,0)/kTL
    etot_new = pe_new + ke_new

    # Compute the argument [dimensionless]
    arg = (etot_new - etot)

    # Compute OP and umbrella bias
    x_new = lmp.gather_atoms("x",1,3)
    #nxtl_new,largest_cluster_new =calc_q12_cluster(n_q12_neigh=16,q12_cutoff=6,crys_cutoff=3.5, n_bonds=12,box=box,x=x_new,n_atoms=natoms) 

    nxtl_new, largest_cluster_new = calc_ql_cluster(bond_crit="RussoTanaka", n_atoms=natoms, box=box, x=x_new, maxnb=60, n_ql_neigh=16, l=12, ql_cutoff=6, cluster_cutoff=3.5, n_bonds=12) 
    arg += calc_bias(0) # unbiased  
    
    # Apply Metropolis acceptance criterion
    if random.random() <= math.exp(-arg):
      n_acc_hmc += 1
      pe = pe_new
      nxtl = nxtl_new
      largest_cluster = largest_cluster_new
      for i in range(3*natoms):
        x[i] = x_new[i]
#    else: # Reject, restore old state
#      check_eng(pe) 
 
  # Volume MC Move
  else:

    # Update number of trials
    n_try_vol += 1.0

    # Scatter coordinates (not necessary here)
    lmp.scatter_atoms("x",1,3,x)

    # Compute random displacement in ln(V)
    lnvol = math.log(vol) + (random.random() - 0.5)*max_lnvol

    # Calculate new box volue, size and scale factor
    vol_new = math.exp(lnvol)
    box_new  = math.pow(vol_new, 1.0/3.0)
    lmp.command("change_box all x final 0.0 %.10f y final 0.0 %.10f z final 0.0 %.10f units box" % (box_new, box_new, box_new))

    # Scale the coordinates and send to LAMMPS
    scalef = box_new/box
    for i in range(3*natoms):   
      x_new[i] = scalef*x[i]
    lmp.scatter_atoms("x",1,3,x_new)
    lmp.command("run 0")

    # Compute the new PE [dimensionless]
    pe_new = lmp.extract_compute("thermo_pe",0,0)/kTL

    # Calculate argument factor [dimensionless]
    arg = (pe_new-pe) + Pf*(vol_new-vol) - (float(natoms) + 1.0)*math.log(vol_new/vol)

    # Compute OP and umbrella bias
    rho_new = df/vol_new
    #nxtl_new,largest_cluster_new =calc_q12_cluster(n_q12_neigh=16,q12_cutoff=6,crys_cutoff=3.5, n_bonds=12,box=box_new,x=x_new,n_atoms=natoms) 

    nxtl_new, largest_cluster_new = calc_ql_cluster(bond_crit="RussoTanaka", n_atoms=natoms, box=box_new, x=x_new, maxnb=60, n_ql_neigh=16, l=12, ql_cutoff=6, cluster_cutoff=3.5, n_bonds=12) 
    arg += calc_bias(1) # unbiased  

    # Apply Metropolis acceptance criterion
    if random.random() <= math.exp(-arg):
      n_acc_vol += 1.0
      pe = pe_new
      vol = vol_new
      box = box_new
      rho = rho_new
      nxtl = nxtl_new
      largest_cluster = largest_cluster_new
      for i in range(3*natoms):
        x[i] = x_new[i]
    else: # Reject, restore the old state
      lmp.command("change_box all x final 0.0 %.10f y final 0.0 %.10f z final 0.0 %.10f units box" % (box, box, box))      
#     check_eng(pe)

 
  # File Input and Output
  if((isweep + 1) % freq_thermo == 0):  # Thermodynamic data
    hmc_acc = acc_ratio(n_acc_hmc, n_try_hmc)
    vol_acc = acc_ratio(n_acc_vol, n_try_vol)
    virial = lmp.extract_compute("thermo_press",0,0)  # Get virial pressure
    #Q6 = calc_q6_optimize(n_q6_neigh=4,q6_cutoff=5.5,box=box,x=x,n_atoms=natoms)
    #Q6_old = calc_q6(n_q6_neigh=4,q6_cutoff=5.5,box=box,x=x,n_atoms=natoms) 
    Q6 = calc_Ql(l=6, n_ql_neigh=4, ql_cutoff=5.5, box=box, x=x, n_atoms=natoms, maxnb=60)
    
    thermo.write("%d %f %f %f %f %f %f %f %f %f\n" % (isweep + 1, kTL*pe, virial, vol, hmc_acc, vol_acc, nxtl,largest_cluster,Q6,rho))
  if((isweep + 1) % freq_traj == 0):    # Trajectory
    traj.write("%d \n" % natoms)
    traj.write("%.10f %.10f %.10f \n" % (box,box,box))
    for i in range(natoms):
      indx = 3*i
      traj.write("%d %.10f %.10f %.10f \n" % (atomid[i], x[indx], x[indx+1], x[indx+2]))
  if((isweep + 1) % freq_restart == freq_restart/2): lmp.command("write_data restart_a.dat") # Alternate restart files for redundancy
  if((isweep + 1) % freq_restart == 0): lmp.command("write_data restart_b.dat")
  if((isweep + 1) % freq_flush == 0):  # Force flush
    thermo.flush()
    traj.flush()

