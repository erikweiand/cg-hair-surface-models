Content of this subdirectory:

	- graphene.data: Graphene sheet for periodic box in x/y. Created using the Nanobuilder in VMD.
	- populateHexBeads.m: MATLAB script to initialize lipids at a grafting distance of 0.65nm for 
	all relevant damage ratios.
	- run.in.combine: LAMMPS script for merging the graphene sheet with the lipids.
	- system.in.*: Setting files required for running LAMMPS
	- Counterions/: Subdirectory to create a thin sheet of counterions to balance the number of 
	charged beads on the hair surface, using packmol and Moltemplate.

Workflow for creating a (damaged) hair surface:

	1) Create a graphene base (already available here) with a tool of your choice (Moltemplate, VMD, Python etc.)
	2) Create a lipid monolayer with "populateHexBeads.m" in MATLAB. The routine could be easily adapted for Python applications.
	3) Create a counterion sheet when using damaged surfaces. The number of charged beads in the monolayer can be obtained in MATLAB 
	with "sum(types == 11)" in MATLAB, after running "populateHexBeads.m". This should be the number of ion beads to be used in packmol/Moltemplate.
	Inside "Counterions", first launch packmol using "packmol < packSodium.inp", then Moltemplate using "moltemplate -xyz system.xyz -nocheck sodium.lt",
	assuming both packmol and Moltemplate are correctly installed on your machine. Please refer to the corresponding manuals for that. The "-nocheck" flag
	for Moltemplate is required because the .lt file directly specifies the bead type as 7, rather than letting Moltemplate assign a number (starting from 1).
	4) Copy the ion data file from "Counterions/" to this directory and merge ion, graphene and lipids using "run.in.combine" with LAMMPS.
