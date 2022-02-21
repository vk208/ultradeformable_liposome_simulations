##########
# CREDIT #
##########

# Written by: Vyshnavi Karra (with help on classes from David Farina Jr)
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

sns.set()

###########
# CLASSES # 
###########

class Bead():
    
    def __init__(self, index=-1, residue='', bead_type='', coords=None):
        
        self.residue = residue
        self.index = int(index)
        self.bead_type = bead_type
        self.coords = coords
        
    @property
    def res_id(self):
        return int(re.search('\d+',self.residue).group(0))
    
    @property
    def res_name(self):
        return re.search('\D+',self.residue).group(0)   
    
    def __repr__(self):
        
        return f'<Bead {self.index} {self.bead_type}>'

class Ensemble():
    
    def __init__(self, name='',beads=None, box_dimensions=None):
        
        self.name = name
        self._beads = beads
        self.box_dimensions = box_dimensions
        
    @property
    def beads(self):
        return self._beads
    
    @beads.setter
    def beads(self,beads):
        for bead in beads:
            assert isinstance(bead,Bead)
        self._beads = beads
        
    @property
    def residues(self):
        rez = []
        for bead in self.beads:
            if bead.res_name not in rez:
                rez.append(bead.res_name)
        return rez
    
        
    def __repr__(self):
        
        return f'<Ensemble {self.name}>'
    
    def read_gro_file(self,path):
        
        with open(path,'r') as f:
            lines = f.readlines()
        
        self.name = lines[0].strip()
        n_beads = int(lines[1].strip())
        beads = []
        
        for line in lines[2:-1]:
            residue,bead_type,index,x,y,z = line.split()
            beads.append(Bead(index=int(index), residue=residue, bead_type=bead_type, 
                             coords=(float(x),float(y),float(z))))

        box_dim = tuple(map(float,lines[-1].split()))
        self._beads = beads
        self.box_dimensions = box_dim
        #assert len(self.beads) == n_beads
        
    def to_gro_file(self,path):
        
        if os.path.dirname(path):
            os.makedirs(os.path.dirname(path),exist_ok=True)
        
        with open(path,'w') as f:
            f.write(f"{self.name}\n")
            f.write(f" {len(self.beads)}\n")
            for b in self.beads:
                line = ("{0:>5d}{1:<5s}{2:>5s}{3:>5d}{4:>8.3f}{5:>8.3f}{6:>8.3f}\n".format(
                    b.res_id, b.res_name, b.bead_type, b.index, b.coords[0], b.coords[1], b.coords[2]))
                f.write(line)
            f.write("  {0:.4f}  {1:.4f}  {2:.4f}\n".format(*self.box_dimensions))

    
    def get_beads(self, index=np.array([]), residue=np.array([]), res_name=np.array([]), res_id=np.array([]), bead_type=np.array([])):
        
        query_beads = []
        for bead in self.beads:
            if index:
                if bead.index not in index:
                    continue
            if residue:
                if bead.residue not in residue:
                    continue
            if bead_type:
                if bead.bead_type not in bead_type:
                    continue
            if res_name:
                if bead.res_name not in res_name:
                    continue
            if res_id:
                if bead.res_id not in res_id:
                    continue
            query_beads.append(bead)
        
        return query_beads
    
    def sort(self):
        self.beads.sort(key=lambda x: x.index)
        
    def get_positions(self):
    
        positions = np.zeros((len(self.beads),3))
        
        for i,bead in enumerate(self.beads):
            positions[i,:] = bead.coords
            
        return positions
        
    def plot(self):
        
        res = self.residues
        color_i = np.random.randint(0,len(cm.colors.BASE_COLORS),len(res))
        res_colors = {}
        for i,r in zip(color_i,res):
            res_colors[r] = list(cm.colors.BASE_COLORS.values())[i]
            
        colors = [res_colors[bead.res_name] for bead in self.beads]
        
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        
        xyzs = self.get_positions()
        
        ax.scatter(xyzs[:,0],xyzs[:,1],xyzs[:,2],s=1,c=colors)
        
        return ax
    
    def to_dataframe(self):
        
        return pd.DataFrame([b.__dict__ for b in ensemble.beads])

################        
# GET WF BEADS # 
################

# System 1 
path = ''
system = Ensemble()
system.read_gro_file(path=path)
W_beads = system.get_beads(bead_type=['W'])
for i in range(0,int(len(W_beads)/10)): 
    W_beads[i].bead_type='WF'
    W_beads[i].residue+='F'
    
new_path = ''
system.to_gro_file(path= new_path)

WF_beads = int(len(W_beads)/10)
print(WF_beads)
print(int(len(W_beads))-WF_beads)      


