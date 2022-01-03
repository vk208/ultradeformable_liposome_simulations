##########
# CREDIT #
##########

# Written by: Jiaming Xu 
# Edited by: Vyshnavi Karra 
# University: Northeastern University 
# Advisor: Francisco Hung 

# Last Updated: 01-03-2022 

# For any questions, please create an issue on GitHub or email Jiaming Xu (xu.jiam@northeastern.edu)

##########
# SET UP #
##########

import os
import MDAnalysis as mda
import numpy as np
from matplotlib import pyplot as plt

############# 
# FUNCTIONS #
#############

def LonlatCreator (array, bins = 20):
    
    x_min = np.min(array[:,0])
    x_max = np.max(array[:,0])
    y_min = np.min(array[:,1])
    y_max = np.max(array[:,1])
    lon = np.linspace(x_min, x_max, bins+3)
    lat = np.linspace(y_min, y_max, bins+3)
    
    return lon, lat

def bounceZBack(array, is_solvent= True):
    '''
    
    This function is designed to move the whole system by centering the lipid average position in z direction. This
    is because the lipid membrane moves in z direction in NAMD simulation.
    
    This function first abstract the average z position of all atoms in lipids, then move the outside part after
    substraction (both above and below original pbc) accross the z boundary.
    
    Return new positions as an array.
        
    '''
    lipid_ave_z = np.average(u.select_atoms('resname DYPC or resname DDPC or resname DOPC').positions[:,2])
    centered_pos = array - [0, 0, lipid_ave_z]
    
    if is_solvent:
        z_range = np.ptp(array[:,2])
        

        if lipid_ave_z < 0:
            z_max = np.max( array[:,2])
            centered_pos[centered_pos[:,2] > z_max ] -= [0, 0, z_range]
        else:
            z_min = np.min( array[:,2])
            centered_pos[centered_pos[:,2] < z_min ] += [0, 0, z_range]

    return centered_pos

    def aveNearTwo(array):
    a = array[:-1]
    b = array[1:]
    return ( a + b ) /2
    
    
def WaterDensMap(universe, resname=False, natoms=3):
    '''
    
    Three inputs.
    1. Positions array
    2. Residue name for selection
    3. Atom name for selection
    
    upper and lower density map as output.
    '''
    

    upper_exp_dens = []
    lower_exp_dens = []
    xdim = []


    for i in u.trajectory[:]:

        pos_water_ini = u.select_atoms('resname %s'%resname).wrap()

        # Determine if water is separated by lipids. in other words, if lipids accross z boundary.
        if np.ptp(pos_water_ini, axis=0)[2] > 63:
            dens, zedges = np.histogram(pos_water_ini[:,2], bins= 20)
            center_of_blank = np.average(aveNearTwo(zedges)[dens == 0])    
            pos_water_ini = u.select_atoms('resname %s'%resname).translate([0,0,-center_of_blank]).wrap()

        
        center_line = np.average(pos_water_ini, axis=0)[2]
    
        
        xdim.append(i.dimensions[0])

        upper_pos = pos_water_ini[pos_water_ini[:,2] <= center_line]
        lower_pos = pos_water_ini[pos_water_ini[:,2] > center_line]   # Reverse due to pbc
        # Calculating heat map for upper and lower water for each frame
        upper_dens, _, _ = np.histogram2d(upper_pos[:,0], upper_pos[:,1], bins = 21)
        lower_dens, _, _ = np.histogram2d(lower_pos[:,0], lower_pos[:,1], bins = 21)    

        # Adding heat maps to summary
        upper_exp_dens.append(upper_dens)
        lower_exp_dens.append(lower_dens)
        
    # Convert to nm^2
    xdim = np.average(np.array(xdim))/10 ** 2
    
    # Averaging over frames
    upper_exp_dens = np.average(np.array(upper_exp_dens), axis=0) / xdim / natoms
    lower_exp_dens = np.average(np.array(lower_exp_dens), axis=0) / xdim / natoms
    
    return upper_exp_dens, lower_exp_dens

    def LipidDensMap(universe, resname, natoms = 1):
    upper_exp_dens = []
    lower_exp_dens = []
    xdim = []
    
    for i in u.trajectory[:]:
        xdim.append(i.dimensions[0])
        u.atoms.wrap()
        upper_lipid_pos = u.select_atoms('resid 1:100 and resname %s'%resname).positions
        lower_lipid_pos = u.select_atoms('resid 101:200 and resname %s'%resname).positions
        #lon_bins, lat_bins = LonlatCreator(lipid)
        upper_dens, _, _ = np.histogram2d(upper_lipid_pos[:,0], upper_lipid_pos[:,1], bins=21)
        lower_dens, _, _ = np.histogram2d(lower_lipid_pos[:,0], lower_lipid_pos[:,1], bins=21)    

        # Adding heat maps to summary
        upper_exp_dens.append(upper_dens)
        lower_exp_dens.append(lower_dens)

        # Averaging over frames
    xdim = np.average(np.array(xdim))/10 ** 2
    upper_exp_dens = np.average(np.array(upper_exp_dens), axis=0) / xdim /natoms
    lower_exp_dens = np.average(np.array(lower_exp_dens), axis=0) / xdim /natoms
    
    return upper_exp_dens, lower_exp_dens

#################
# SYSTEM SET UP # 
#################
path = ''
coord_file = '' # GRO/PSF/TPR 
traj_file = ''  # XTC/TRR/DCD

# Change Directory
os.chdir(path) 

# Load Trajectory 
u = mda.Universe(coord_file,traj_file) 

################
# RUN ANALYSIS # 
################

# IMPORTANT: Change *RESNAME* to select different groups
waterup, waterlo = WaterDensMap(u, resname='PW') 

# IMPORTANT: Change *RESNAME* and *NUMBER* (number of atoms/beads in residue)
ddpcup, ddpclo = LipidDensMap(u, 'DTPC', 8) 
dopcup, dopclo = LipidDensMap(u, 'DOPC', 12)

############
# PLOTTING # 
############

fig, axs = plt.subplots(nrows = 3, ncols = 2, figsize=(8,12), constrained_layout=True)

vmax = np.max(np.vstack((waterup, waterlo)))
im1 = axs[0,0].imshow(waterup, vmax=vmax)
axs[0,0].set_title('Water Density above membrane')
im2 = axs[0,1].imshow(waterlo, vmax=vmax)
axs[0,1].set_title('Water Density below membrane')

vmax = np.max(np.vstack((dopclo, dopcup)))#, ddpclo, ddpcup)))
im3 = axs[1,0].imshow(ddpcup, vmax=vmax)
axs[1,0].set_title('DDPC Density in Upper Leaflet')
im4 = axs[1,1].imshow(ddpclo, vmax=vmax)
axs[1,1].set_title('DDPC Density in Lower Leaflet')

im5 = axs[2,0].imshow(dopcup, vmax=vmax)
axs[2,0].set_title('DOPC Density in Upper Leaflet')
im6 = axs[2,1].imshow(dopclo, vmax=vmax)
axs[2,1].set_title('DOPC Density in Lower Leaflet')


fig.colorbar(im1, ax=axs[0,1],shrink=0.9)
fig.colorbar(im5, ax=axs[1,1],shrink=0.9)
fig.colorbar(im5, ax=axs[2,1],shrink=0.9)

#IMPORTANT: Change the *TITLE* 
fig.suptitle('75% DOPC, and Water Density Map at 0 mN/m Surface Tension')# 0 mN/m surface tension')
plt.show()
