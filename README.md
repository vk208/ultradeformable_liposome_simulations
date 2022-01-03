# Ultradeformable Liposome Simulations
This is for the simulation work done on ultradeformable liposomes by Dr. Francisco Hung's computational lab at Northeastern University (Boston, MA). 

This was created by Vyshnavi Karra (GitHub: vk208) and Jiaming Xu (GitHub: JP-Xu).

For more information, please create an issue on GitHub or email Jiaming Xu (xu.jiam@northeastern.edu) 

# Codes 
VoroAnalysis.py takes coordinate and trajectory files from your lipid bilayer simulation, creates a Voronoi diagram for each timestep, and compiles those diagrams into an mp4 gif. 

DensityHeatMap.py takes coordinate and trajectory files from your lipid bilayer simulation and creates water and lipid density heat maps for both the lower and upper leaflets. 

CompressibilityModulus.py reads the xvg files, calculates the area compressibility modulus and the bending modulus based on this paper (DOI: https://doi.org/10.1038/s41598-019-44318-9), and plots them. 
