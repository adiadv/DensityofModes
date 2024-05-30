# DensityofModes
Scripts for Calculating the Vibrational Density of Modes from voltage data or velocity data of granular matter. 

File Name: getDOMvoltage:

This file takes a string of voltages as a vector as input (can be multiple strings in rows as a matrix for multiple entries) and it should output the velocities of particles, Velocity Autocorr and Vibrational Density of Modes individually for each device and also do an average. 

File Name: getDOMvelocity:

This file takes a string of velocities and gives you the Autocorr and DOM invidually for each device and also does an average. 

Both files above can pad the data to 2^n, but do not do this in general. You need to uncomment and change those lines manually. 
