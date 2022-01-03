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

colorlist1 = ['#F8B195', '#F67280', '#C06C84', '#6C5B7B', '355C7D']
colorlist = ['#A8E6CE','#DCEDC2','#FFD385','#FFAAA6','#FF8C94']

from freud.locality import Voronoi
import freud
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os 
import MDAnalysis as mda

# Local paths 
path = ''
coord_file = ''
traj_file = ''

# Forcefield for MDAnalysis 
forcefield = 'martini'

# Change directory
os.chdir(path)

# Load trajectory
u = mda.Universe(coord_file,traj_file, type=forcefield)

####################
# VORONOI ANALYSIS #
####################

filenames = []

for ii in u.trajectory[::1]:
  # IMPORTANT: Change the *RESNAME*, the *BEAD NAME*, and the *RESID* for each system that is analyzed!!! 

   # Process positions of all lipids for Voronoi plot. Upper leaflet
   upper_dxpcpos = u.select_atoms('resname D*PC and name PO4 and resid 1:1520').wrap()
   upper_dxpcpos[:,2] = 0

   # Load positions of each lipid for coloring them separately.
   upper_dopcpos = u.select_atoms('resname DOPC and name PO4 and resid 1:1520').wrap()
   upper_ddpcpos = u.select_atoms('resname DTPC and name PO4 and resid 1:1520').wrap()

   # Find the index of lipid
   upper_lipidA_index = []
   for i in upper_dopcpos:
       same_pos = np.where(upper_dxpcpos == [i])[0] # Comparing positions in one lipid to overall lipids.
       values, counts = np.unique(same_pos, return_counts=True)
       upper_lipidA_index.append(values[counts != 1])

   upper_lipidB_index = []
   for i in upper_ddpcpos:
       same_pos = np.where(upper_dxpcpos == [i])[0] # Comparing positions in one lipid to overall lipids.
       values, counts = np.unique(same_pos, return_counts=True)
       upper_lipidB_index.append(values[counts != 1])

   # Process positions of all lipids for Voronoi plot. Upper leaflet
   lower_dxpcpos = u.select_atoms('resname D*PC and name PO4 and resid 1521:3040').wrap()
   lower_dxpcpos[:,2] = 0

   # Load positions of each lipid for coloring them separately.
   lower_dopcpos = u.select_atoms('resname DOPC and name PO4 and resid 1521:3040').wrap()
   lower_ddpcpos = u.select_atoms('resname DTPC and name PO4 and resid 1521:3040').wrap()

   # Find the index of lipid
   lower_lipidA_index = []
   for i in lower_dopcpos:
       same_pos = np.where(lower_dxpcpos == [i])[0] # Comparing positions in one lipid to overall lipids.
       values, counts = np.unique(same_pos, return_counts=True)
       lower_lipidA_index.append(values[counts != 1])

   lower_lipidB_index = []
   for i in lower_ddpcpos:
       same_pos = np.where(lower_dxpcpos == [i])[0] # Comparing positions in one lipid to overall lipids.
       values, counts = np.unique(same_pos, return_counts=True)
       lower_lipidB_index.append(values[counts != 1])

   box = freud.box.Box(Lx=ii.dimensions[0], Ly=ii.dimensions[1], is2D=True)
   voro = freud.locality.Voronoi()


   ## Plotting ##

   fig, axs = plt.subplots(nrows = 1, ncols = 2, figsize=(16,8), constrained_layout=True)

   XLIMS = [-10,ii.dimensions[0]+10]
   #tensions3 = ['-7','0','7','15']

   # Plot for upper leaflet
   ax = axs[0]
   a = np.array(voro.compute((box, upper_dxpcpos)).polytopes, dtype=object)
   for i in range(len(a)):
       c = np.append(a[i],[a[i][0]], axis=0)
       x = c[:,0]
       y = c[:,1]
       ax.plot(x,y, color='gray', linewidth=.5)
       
       # IMPORTANT: Change the *LABEL* as needed!!! 
       if i in upper_lipidA_index:
           ax.fill(x,y, color=colorlist[0], label="DOPC" if i == upper_lipidA_index[-1] else "_nolegend_")

       elif i in upper_lipidB_index:
           ax.fill(x,y, color=colorlist[-2], label="DTPC" if i == upper_lipidB_index[-1] else "_nolegend_")

       ax.scatter(upper_dxpcpos[:,0], upper_dxpcpos[:,1], s=1, c='gray', zorder=2)
       ax.set_title('Upper Leafleat')
       ax.set_xlim(XLIMS)
       ax.set_ylim(XLIMS)

   # Plot for lower leaflet
   ax = axs[1]
   b = np.array(voro.compute((box, lower_dxpcpos)).polytopes, dtype=object)
   for i in range(len(b)):
       c = np.append(b[i],[b[i][0]], axis=0)
       x = c[:,0]
       y = c[:,1]
       ax.plot(x,y, color='gray', linewidth=.5)

       if i in lower_lipidA_index:
           ax.fill(x,y, color=colorlist[0])
       elif i in lower_lipidB_index:
           ax.fill(x,y, color=colorlist[-2])
       ax.scatter(lower_dxpcpos[:,0], lower_dxpcpos[:,1], s=1, c='gray', zorder=2)    
       ax.set_title('Lower Leafleat')
       ax.set_xlim(XLIMS)
       ax.set_ylim(XLIMS)

   # IMPORTANT: Change the *TITLE* as needed!!! 
   fig.suptitle('Voronoi Plot of 75% DOPC 25% DDPC 30 nm X 30 nm Membrane')
   fig.legend(loc='upper center', bbox_to_anchor=(0.5, 0.95), ncol=3, fancybox=True, shadow=True)

   #IMPORTANT: You can change the *FILENAME* here if desired!!
   filename = 'f%i.png'%ii.frame
   filenames.append(filename)
   plt.savefig(filename, dpi=400)
   plt.close(fig)

#############
# BUILD GIF # 
#############

import imageio

# IMPORTANT: You can export to mp4 file for better interaction. but have to install ffmpeg.
mp4_name = ''

with imageio.get_writer(mp4_name, mode='I') as writer: 
    for filename in filenames:
        image = imageio.imread(filename)
        writer.append_data(imageio.imread(filename))
        

