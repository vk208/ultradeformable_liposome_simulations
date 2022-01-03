
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
import math 
import scipy as sp
import glob 
import re
from collections import defaultdict

##########
# GROUPS # 
##########

# Surface Tensions 
tensions = [0,2,5,10,20,30,40,50,60,70]

################
# ENERGY PLOTS #
################

# System 1

path = '' # working directory 

for i in range(len(tensions)): 
    with open('S1_energy{}.submit'.format(int(tensions[i])), 'w') as energy: 
        energy.write('#!/bin/bash\n')
        energy.write('#SBATCH -N 1\n')
        energy.write('#SBATCH --tasks-per-node=1\n')
        energy.write('#SBATCH --cpus-per-task=4\n')
        energy.write('#SBATCH --output=energy{}.out\n'.format(int(tensions[i])))
        energy.write('#SBATCH --error=energy{}.err\n'.format(int(tensions[i])))
        energy.write('#SBATCH -J energy{}\n'.format(int(tensions[i])))
        energy.write('#SBATCH -p {}\n'.format(partition))
        energy.write('#SBATCH -t 4:00:00\n')
        energy.write('#SBATCH --chdir={}\n'.format(path))
        energy.write('\n')
        energy.write('module purge \n')
        energy.write('module load openmpi/3.1.1\n')
        energy.write('module load gromacs/2018.4-cpu\n')
        energy.write('\n')
        energy.write('source /shared/centos7/gromacs/2018.4-cpu/bin/GMXRC\n')
        energy.write('\n')
        energy.write('echo 13 15 17 22 23 0 | gmx_mpi energy -f dmpc_md_gamma{}.edr -s dmpc{}.tpr -o dmpc_gamma{}_energy.xvg\n'.format(int(tensions[i]),int(tensions[i]),int(tensions[i])))
energy.close()

with open('S1energy.bash','w') as bash: 
    bash.write('#!/bin/bash\n')
    for i in range(len(tensions)): 
        bash.write('sbatch S1_energy{}.submit\n'.format(int(tensions[i])))
    bash.close()

#################################
# TRJCONV FOR SMALLER XTC FILES # 
#################################

# System 1

path = '' # working directory 
partition = ''

for i in range(len(tensions)): 
    with open('S1_trjconv{}.submit'.format(int(tensions[i])), 'w') as trjconv: 
        trjconv.write('#!/bin/bash\n')
        trjconv.write('#SBATCH -N 1\n')
        trjconv.write('#SBATCH --tasks-per-node=1\n')
        trjconv.write('#SBATCH --cpus-per-task=4\n')
        trjconv.write('#SBATCH --output=trj{}.out\n'.format(int(tensions[i])))
        trjconv.write('#SBATCH --error=trj{}.err\n'.format(int(tensions[i])))
        trjconv.write('#SBATCH -J trj{}\n'.format(int(tensions[i])))
        trjconv.write('#SBATCH -p {}\n'.format(partition))
        trjconv.write('#SBATCH -t 1-00:00:00\n')
        trjconv.write('#SBATCH --chdir={}\n'.format(path))
        trjconv.write('\n')
        trjconv.write('module purge \n')
        trjconv.write('module load openmpi/3.1.1\n')
        trjconv.write('module load gromacs/2018.4-cpu\n')
        trjconv.write('\n')
        trjconv.write('source /shared/centos7/gromacs/2018.4-cpu/bin/GMXRC\n')
        trjconv.write('\n')
        trjconv.write('gmx_mpi trjconv -f dmpc_6eq_gamma{}.xtc -o dmpc_gamma{}_md_video.xtc -skip 1000\n'.format(int(tensions[i]),int(tensions[i])))
trjconv.close()

with open('S1trj.bash','w') as bash: 
    bash.write('#!/bin/bash\n')
    for i in range(len(tensions)): 
        bash.write('sbatch S1_trjconv{}.submit\n'.format(int(tensions[i])))
    bash.close()

#########################
# SURFACE TENSION PLOTS #
#########################

# System 1 

path = '' # working directory 
partition = '' 

for i in range(len(tensions)): 
    with open('S1_surf{}.submit'.format(int(tensions[i])), 'w') as surf: 
        surf.write('#!/bin/bash\n')
        surf.write('#SBATCH -N 1\n')
        surf.write('#SBATCH --tasks-per-node=1\n')
        surf.write('#SBATCH --cpus-per-task=4\n')
        surf.write('#SBATCH --output=energy{}.out\n'.format(int(tensions[i])))
        surf.write('#SBATCH --error=energy{}.err\n'.format(int(tensions[i])))
        surf.write('#SBATCH -J energy{}\n'.format(int(tensions[i])))
        surf.write('#SBATCH -p {}\n'.format(partition))
        surf.write('#SBATCH -t 4:00:00\n')
        surf.write('#SBATCH --chdir={}\n'.format(path))
        surf.write('\n')
        surf.write('module purge \n')
        surf.write('module load openmpi/3.1.1\n')
        surf.write('module load gromacs/2018.4-cpu\n')
        surf.write('\n')
        surf.write('source /shared/centos7/gromacs/2018.4-cpu/bin/GMXRC\n')
        surf.write('\n')
        surf.write('echo 21 35 39 43 0 | gmx_mpi energy -f dmpcddpc_md_gamma{}.edr -s dmpcddpc{}.tpr -o dmpcddpc_gamma{}_surfacetension.xvg\n'.format(int(tensions[i]),int(tensions[i]),int(tensions[i])))
surf.close()

with open('S1surf.bash','w') as bash: 
    bash.write('#!/bin/bash\n')
    for i in range(len(tensions)): 
        bash.write('sbatch S1_surf{}.submit\n'.format(int(tensions[i])))
    bash.close()

#############
# APL PLOTS #
#############

# System 1 

path = '' # working directory 
partition = '' 

for i in range(len(tensions)): 
    with open('S1_apl{}.submit'.format(int(tensions[i])), 'w') as apl: 
        apl.write('#!/bin/bash\n')
        apl.write('#SBATCH -N 1\n')
        apl.write('#SBATCH --tasks-per-node=1\n')
        apl.write('#SBATCH --cpus-per-task=4\n')
        apl.write('#SBATCH --output=energy{}.out\n'.format(int(tensions[i])))
        apl.write('#SBATCH --error=energy{}.err\n'.format(int(tensions[i])))
        apl.write('#SBATCH -J energy{}\n'.format(int(tensions[i])))
        apl.write('#SBATCH -p {}}\n'.format(partition))
        apl.write('#SBATCH -t 4:00:00\n')
        apl.write('#SBATCH --chdir={}}\n'.format(path))
        apl.write('\n')
        apl.write('module purge \n')
        apl.write('module load openmpi/3.1.1\n')
        apl.write('module load gromacs/2018.4-cpu\n')
        apl.write('\n')
        apl.write('source /shared/centos7/gromacs/2018.4-cpu/bin/GMXRC\n')
        apl.write('\n')
        apl.write('echo 19 20 21 0 | gmx_mpi energy -f dmpc_md_gamma{}.edr -s dmpc{}.tpr -o dmpc_gamma{}_areaperlipid.xvg\n'.format(int(tensions[i]),int(tensions[i]),int(tensions[i])))
apl.close()

with open('S1apl.bash','w') as bash: 
    bash.write('#!/bin/bash\n')
    for i in range(len(tensions)): 
        bash.write('sbatch S1_apl{}.submit\n'.format(int(tensions[i])))
    bash.close()

#################
# DENSITY PLOTS #
#################

# System 1 

path = '' # working directory 
partition = '' 

for i in range(len(tensions)): 
    with open('S1_density{}.submit'.format(int(tensions[i])), 'w') as density: 
        density.write('#!/bin/bash\n')
        density.write('#SBATCH -N 1\n')
        density.write('#SBATCH --tasks-per-node=1\n')
        density.write('#SBATCH --cpus-per-task=4\n')
        density.write('#SBATCH --output=energy{}.out\n'.format(int(tensions[i])))
        density.write('#SBATCH --error=energy{}.err\n'.format(int(tensions[i])))
        density.write('#SBATCH -J energy{}\n'.format(int(tensions[i])))
        density.write('#SBATCH -p {}}\n'.format(partition))
        density.write('#SBATCH -t 4:00:00\n')
        density.write('#SBATCH --chdir={}\n'.format(path))
        density.write('\n')
        density.write('module purge \n')
        density.write('module load openmpi/3.1.1\n')
        density.write('module load gromacs/2018.4-cpu\n')
        density.write('\n')
        density.write('source /shared/centos7/gromacs/2018.4-cpu/bin/GMXRC\n')
        density.write('\n')
        density.write('echo 3 3 0 | gmx_mpi density -f dmpc_6eq_gamma{}.xtc -s dmpc{}.tpr -o dmpc_gamma{}_density.xvg\n'.format(int(tensions[i]),int(tensions[i]),int(tensions[i])))
density.close()

with open('S1density.bash','w') as bash: 
    bash.write('#!/bin/bash\n')
    for i in range(len(tensions)): 
        bash.write('sbatch S1_density{}.submit\n'.format(int(tensions[i])))
    bash.close()

