Content of this subdirectory:

	- BulkHexadecane/: Subdirectory for bulk hexadecane equilibration
	- BulkHexadecane/Moltemplate/: Subdirectory to set up initial configuration for
	bulk hexadecane.
	- run.in.droplet: LAMMPS input file to extract a hemispherical droplet
	- system.in.*: LAMMPS setting files

Workflow for creating hemispherical droplet:

	1) Create a cubic bulk hexadecane system with sufficient box length. The box length can be smaller than the droplet diameter, as the box can
	be replicated in x,y,z before extracting the droplet. However, it needs to be large for sufficient decay of correlation functions.
	Use the "BulkHexadecane/Moltemplate/" folder to set up a domain with randomly distributed MARTINI hexadecane molecules:
		a) With "packmol < packHex.inp" create the initial coordinates.
		b) With "moltemplate.sh -xyz system.xyz -nocheck hexadecane_box.lt" create the LAMMPS data file.
		c) Copy the "hexadecane_box.data" file from "BulkHexadecane/Moltemplate/" to "BulkHexadecane/".
	2) Equilibrate the hexadecane box in LAMMPS in "BulkHexadecane/" at ambient pressure:
		
		lmp -i run.in.min
		lmp -i run.in.npt

	Note that the system.in.* files could be optimized to exclude any Coulombic and thus long-range interactions, to significantly speed up the
	bulk simulations. The system.in.* files provided here are derived from the contact angle simulations for consistency and therefore include
	long-range electrostatics by default.
	3) Copy the equilibrated box "system_after_npt.data" from "BulkHexadecane/" to the current directory. Extract the hemispherical droplet from 
	the equilibrated box using "lmp -i run.in.droplet". Enable tiling if the box is smaller than the droplet radius by uncommenting the
	corresponding line in "run.in.droplet".
