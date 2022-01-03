##########
# CREDIT #
##########

# Written by: Vyshnavi Karra 
# Edited by: Vyshnavi Karra 
# University: Northeastern University 
# Advisor: Francisco Hung 

# Last Updated: 01-03-2022 

# For any questions, please create an issue on GitHub

##########
# SET UP #
##########

import numpy as np
import scipy as sp
import glob
import matplotlib.pyplot as plt 
import pandas as pd
import re 
import seaborn as sns
import pandas as pd
import os

from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import MaxNLocator


sns.set_style('ticks')

#############
# FUNCTIONS #
#############

def read_file(path): 
    """
    This function takes in the local path of the energy .xvg file that was generated using the command 
    "echo 13 14 15 0 | gmx_mpi energy -f <filename>.edr -s <filename>.tpr -o <system>_areaperlipid.xvg" in 
    Discovery. Remember that openmpi and gromacs should be loaded before running that command in Slurm, 
    can check typing "module list" in Terminal/Command Prompt. 
    
    This .xvg file should have 3 main variables: 
    1. time (in ps)
    2. box_x (in nm)
    3. box_y (in nm)
    4. box_z (in nm)
    """
    
    data = []
    time = []        #in picoseconds (ps)
    box_x = []       #in nm, box vector elements
    box_y = []       #in nm, box vector elements
    box_z = []       #in nm, box vector elements 
    
    #open file and put all of the variable information into a list
    with open(path,'r') as f:
         for line in f:
            if (not (line.lstrip().startswith('#')) and not ((line.lstrip().startswith('@')))):
                cols = line.split() 
                if len(cols)==(4): 
                    for i in range(len(cols)): 
                        data.append(cols[i])
                else: 
                    print('This .xvg file is missing a variable (box x, box y, box z, density)')
                    
    #convert the list into an array and reshape it
    data = np.reshape(np.array(data),[int(len(data)/4),4])  
    
    for i in range(len(data)): 
        time.append(float(data[i][0]))
        box_x.append(float(data[i][1]))
        box_y.append(float(data[i][2]))
        box_z.append(float(data[i][3]))
    
    
    return data, time, box_x, box_y, box_z

def read_density_file(path): 
    """
    This function takes in the local path of the energy .xvg file that was generated using the command 
    "echo 2 2 | gmx_mpi density -f <filename>.xtc -s <filename>.tpr -center 
    -o <system>_<lipid if necessary>density.xvg" in Discovery. For the mixed lipid systems, also do 
    "echo 3 3 | gmx_mpi density -f <filename>.xtc -s <filename>.tpr -center 
    -o <system>_<lipid if necessary>density.xvg" in Discovery. Continue to do that as many times as there are 
    lipids in the systems. Remember to rename the files between each one if the -o flag isn't added. 
    
    This .xvg file should have 2 main variables: 
    1. relative position from the center of the simulation box (in nm)
    2. density (in kg/m^3)
    """
    
    data = []
    position = []
    density = [] 
    
    #open file and put all of the variable information into a list
    with open(path,'r') as f:
         for line in f:
            if (not (line.lstrip().startswith('#')) and not ((line.lstrip().startswith('@')))):
                cols = line.split() 
                if len(cols)==(2): 
                    for i in range(len(cols)): 
                        data.append(cols[i])
                else: 
                    print('This .xvg file is missing a variable (position, density)')
                    
    #convert the list into an array and reshape it
    data = np.reshape(np.array(data),[int(len(data)/2),2])  
    
    for i in range(len(data)): 
        position.append(float(data[i][0]))
        density.append(float(data[i][1]))
    
    
    return data, position, density 


def areal_strain(box_x1,box_y1,box_x2,box_y2): 
    """
    This function will calculate the surface area of the bilayers to then calculate the compressibility modulus
    as defined in this paper: https://doi.org/10.1038/s41598-019-44318-9 
    
    Ka is defined as: 
    
    """
    area1 = []
    for i in range(len(box_x1)): 
        if len(box_x1)==len(box_y1):
            area1.append(box_x1[i]*box_y1[i])
        else: 
            print('Error: The xy dimensions of the 1st system are not equal in length. x is {} and y is {} long'.format(len(box_x1),len(box_y1)))
    
    area2 = []
    for i in range(len(box_x2)): 
        if len(box_x2)==len(box_y2):
            area2.append(box_x2[i]*box_y2[i])
        else: 
            print('Error: The xy dimensions of the 2nd system are not equal in length. x is {} and y is {} long'.format(len(box_x2),len(box_y2)))
    
    epsilonA = []
    for i in range(len(area1)): 
        if len(area1)==len(area2): 
            epsilonA.append((area1[i]/area2[i]) - 1)
        else: 
            print('Error: The lengths of the area arrays of the 2 systems are different lengths. System 1 is {} and System 2 is {} long'.format(len(area1),len(area2)))
            

    return epsilonA  

def bilayer_thickness(position,density): 
    """
    This function will calculate the bending modulus as defined in this 
    paper: https://doi.org/10.1038/s41598-019-44318-9 
    
    Kc is defined as: 
    
    """
    
    new_density = []
    new_position = []
    for i in range(len(density)): 
        if density[i] != 0: 
            new_density.append(density[i])
            new_position.append(position[i])
    
    return (np.max(new_position) + abs(np.min(new_position)))

############
# GET DATA #
############

# System 1 # 
path = ''
data, time, box_x, box_y, box_z = read_file(path)

# IMPORTANT: Do this ^ for all systems

#####################
# CALCULATE EPSILON #
#####################

# System 1, 2 
epsilonA1 = areal_strain(box_x2,box_y2,box_x0,box_y0)
avg_epsilonA1 = np.sum(epsilonA1)/len(epsilonA1)

# IMPORTANT: Do this ^ for all systems

###############
# PLOT GRAPHS #
###############

# System Array 1 
gamma = np.array([0, 2, 5])
epsilonA = np.array([0, avg_epsilonA1, avg_epsilonA2])

plt.figure(figsize=[10,5])
plt.plot(epsilonA, gamma, 'o', color='cyan', markersize=15)
plt.plot(epsilonA, gamma, '--', color='cyan', linewidth=3)
plt.title('Coarse Grained DOPC System',fontsize=20)
plt.xlabel('Average Areal Strain',fontsize=18)
plt.ylabel('Surface Tension (mN/m)',fontsize=18)
plt.tick_params(axis='both',which='major',labelsize=16)
#plt.tight_layout()
plt.savefig('filename.png')
plt.show()

# IMPORTANT: Do this ^ for all systems, remember to change the *TITLE*, and the save fig *FILENAME*

##########
# GET KA #
##########

# System Array 1 
Ka = []
for i in range(len(gamma)-1): 
    if len(gamma)==len(epsilonA):
        Ka.append((gamma[i+1]-gamma[i])/(epsilonA[i+1]-epsilonA[i]))
    else: 
        print('gamma and epsilon are not the same length. gamma is {} and epsilon is {} long'.format(len(gamma),len(epsilonA)))

avg_Ka = np.sum(Ka)/len(Ka)
print('The compressibility modulus for DOPC is {:3f}'.format(avg_Ka))

# IMPORTANT: Do this ^ for all systems

########################################
# GET DATA FROM DOPC IN MIXED BILAYERS #
########################################

#~ h0 ~# 

# System 1 
path = ''
data, position, density = read_density_file(path)
thickness = bilayer_thickness(position,density)

# IMPORTANT: Do this ^ for all systems

################
# CALCULATE Kc #
################

# System 1, 2
Kc1 = (avg_Ka * (thickness2 - thickness0))/24

# IMPORTANT: Do this ^ for all systems

###############
# PLOT GRAPHS #
###############

# System Array 1 
gamma = np.array([0, 2, 5])
Kc = np.array([0, Kc1, Kc2])

plt.figure(figsize=[10,5])
plt.plot(gamma, Kc, 'o', color='cyan', markersize=15)
plt.plot(gamma, Kc, '--', color='cyan', linewidth=3)
plt.title('Coarse Grained DOPC System',fontsize=20)
plt.xlabel('Surface Tension (mN/m)',fontsize=18)
plt.ylabel('Bending Modulus',fontsize=18)
plt.tick_params(axis='both',which='major',labelsize=16)
#plt.tight_layout()
plt.savefig('filename.png')
plt.show()

# IMPORTANT: Do this ^ for all systems, remember to change the *TITLE*, and the save fig *FILENAME*





