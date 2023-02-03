# cg-hair-surface-models
*Coarse-grained molecular models of the surface of hair* [1]  
E Weiand, JP Ewen, PH Koenig, Y Roiter, SH Page, S Angioletti-Uberti and D Dini

**Maintainer / Corresponding author:**  
Erik Weiand  
Tribology Group, Imperial College London  
erik.weiand19@imperial.ac.uk  

This repository contains the relevant files for studying wetting and surface structure of the 18-MEA monolayer on hair by means of coarse-grained molecular simulations (CGMD) in LAMMPS. The corresponding Soft Matter paper can be found at [https://doi.org/10.1039/D1SM01720A](https://doi.org/10.1039/D1SM01720A "https://doi.org/10.1039/D1SM01720A") (last updated 2022/02/04). Please cite this repository using the following DOI:


[![DOI](https://zenodo.org/badge/421338855.svg)](https://zenodo.org/badge/latestdoi/421338855)


## Content
The following subdirectories are available at the moment:
- 01_CGMD_Simulations/ - Coarse-grained molecular dynamics initial configurations and LAMMPS input files for monolayer simulations in vacuum and with water/n-hexadecane droplets.
- 02_AAMD_Simulations/ - All atom (AA) reference simulation initial topology files and LAMMPS input scripts of fully-functionalized monolayers in vacuum and with water droplets. The OPLS force field and SPC/E water are used in agreement with [2].
- 03_Initial_Configurations/ - Initial configurations and corresponding workflow descriptions for hair surfaces and liquid bulk systems.
- 04_Post_Processsing/ - MATLAB post-processing scripts, primarily for contact angle analysis from LAMMPS trajectory files (.traj).

Further workflow descriptions are added in *README.txt* files in the subdirectories and will be merged into this document in the future. For making full use of this repository, the following tools should be accessible to the user:
- LAMMPS (https://www.lammps.org/)
- Moltemplate (https://www.moltemplate.org/)
- packmol (https://m3g.github.io/packmol/)
- MATLAB (https://de.mathworks.com/products/matlab.html)
- A molecular visualization tool, such as VMD (https://www.ks.uiuc.edu/Research/vmd/)

Open-source accessibility for the post-processing routines is currently not planned but should not be a great challenge to achieve by transferring the MATLAB scripts to Python/NumPy.

## References
[1] E Weiand, JP Ewen, PH Koenig, Y Roiter, SH Page, S Angioletti-Uberti and D Dini, "Coarse-grained molecular models of the surface of hair," *Soft Matter* (2022).  
[2] DW Cheong, FCH Lim, L Zhang, "Insights into the Structure of Covalently Bound Fatty Acid Monolayers on a Simpliﬁed Model of the Hair Epicuticle from Molecular Dynamics Simulations," *Langmuir* , **28**, 13008-13017 (2012)
