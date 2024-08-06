# DensityofModes
Scripts for Calculating the Vibrational Density of Modes from voltage data or velocity data of granular matter. 

File Name: getDOMvoltage:

This file takes a string of voltages as a vector as input (can be multiple strings in rows as a matrix for multiple entries) and it should output the velocities of particles, Velocity Autocorr and Vibrational Density of Modes individually for each device and also do an average. 

Input Expected: Columns of voltages proportional to particle accelerations (no time required). 

File Name: getDOMvelocity:

This file takes a string of velocities and gives you the Autocorr and DOM invidually for each device and also does an average. 

Both files above pad data to double length and do a velocity autocorr in fourier space. They also include a 60hz filter which is on by default. You need the SIGNAL PROCESSING TOOLBOX to use these files.  

Input Expected: Columns of particle velocities (no time required). 
