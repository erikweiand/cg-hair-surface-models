In order to run a simulation, simply copy the following files from "runfiles" into the desired subdirectory:

	1) run.in.min - Energy minimization
	2) run.in.nvt - Production run in the NVT ensemble
	3) system.in.init - General settings
	4) system.in.settings - Force field settings/atom groups

Each subdirectory already includes the initial configuration topology consisting of fully-functionalized 
surfaces. Note that these are all atom simulations using OPLS! Simulations should be run using LAMMPS:

	lmp -i run.in.min
	lmp -i run.in.nvt

Paper results were obtained using parallel runs with 48 cores on a single node each.
