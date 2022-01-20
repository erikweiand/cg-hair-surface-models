In order to run a simulation, simply copy the following files from "runfiles" into the desired subdirectory:

	1) run.in.min - Energy minimization
	2) run.in.nvt - Production run in the NVT ensemble
	3) system.in.init - General settings
	4) system.in.settings - Force field settings/atom groups

Each subdirectory already includes the initial configuration topology consisting of (damaged) surfaces
with deposited hemispherical droplets. Counterions have been equilibrated beforehand and are
located closely to the charged surfaces where applicable. Simulations should be run using LAMMPS:

	lmp -i run.in.min
	lmp -i run.in.nvt

Paper results were obtained using parallel runs with 32 to 48 cores on a single node each.
