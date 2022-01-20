This directory contains post-processing tools for LAMMPS hair surface simulations, primarily for contact angle
simulations:

	- readTraj.m: MATLAB script to read trajectory files. No bond or angle information
	is included.
	- combineMat.m: MATLAB script to combine trajectory files in .mat format.
	- plotDensityProfiles: MATLAB script to obtain time-averaged density profiles from
	trajectory files. Alternatively, use LAMMPS runtime averaging and binning.
	- calcContactAngle.m: MATLAB script to estimate a single, time-averaged contact 
	angle from trajectories.
	- calcContactAngle_withErrors.m. MATLAB script to estimate time-averaged contact
	angles with error estimates from trajectories. This is the preferred method for final
	contact angle values and has been used in the paper as well.
	
