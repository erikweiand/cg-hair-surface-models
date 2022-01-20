Content of this subdirectory:

	- BulkWater/: Subdirectory for bulk water equilibration
	- BulkWater/Moltemplate/: Subdirectory to set up initial configuration for
	bulk water.
        - run.in.droplet: LAMMPS input file to extract a hemispherical droplet
        - system.in.*: LAMMPS setting files

Workflow for creating hemispherical droplet:

	1) Create a cubic bulk water system with sufficient box length. The box length can be smaller than the droplet diameter, as the box can
	be replicated in x,y,z before extracting the droplet. However, it needs to be large for sufficient decay of correlation functions.
	Use the "BulkWater/Moltemplate/" folder to set up a domain with randomly distributed polarizable MARTINI water beads:
		a) With "packmol < packWater.inp" create the initial coordinates.
		b) With "moltemplate.sh -xyz system.xyz -nocheck water_box.lt" create the LAMMPS data file.
		c) Copy the "water_box.data" file from "BulkWater/Moltemplate/" to "BulkWater/".
	2) Equilibrate the water box in LAMMPS in "BulkWater/" at ambient pressure:
		
		lmp -i run.in.min
		lmp -i run.in.npt

	3) Copy the equilibrated box "system_after_npt.data" from "BulkWater/" to the current directory. Extract the hemispherical droplet from 
	the equilibrated box using "lmp -i run.in.droplet". Enable tiling if the box is smaller than the droplet radius by uncommenting the
	corresponding line in "run.in.droplet".
