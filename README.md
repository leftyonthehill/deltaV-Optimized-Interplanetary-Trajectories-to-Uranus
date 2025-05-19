# DeltaV-Optimized-Interplanetary-Trajectories-to-Uranus

To successfuly run this project, first open mgaGenticCall.m. This is the parent script that calls all other MATLAB functions in the folder.
If interested in running my code, I would consider looking at the user variables at the top of the script. When the mission string is set
equal to "Uranus", then the script will try and find a route to Uranus and plot the trajectory. If the mission type is set equal to 
"Galileo", then the script will try to recreate Galileo's flight data and produce a series of 4 plots comparing the two data sets 
together. Feel free to run the script as much as you wish to create new trajectories to Uranus, just be aware that by not running the code in 
parallel, see line 264 to see its status, it will take a few minutes to run just one. On the other hand, if you are interested in seeing 
the code perform its verification, the code is set to only produce one trajectory and one set of graphs because we are looking to verify 
rather than discover.

If you are interested in reviewing my data, first load the data file you would like to review (either bestVerifiedGalileo.mat or 
uranusData60.mat). Then update the value of mission and only run the first section at the top of the script. Then, finally, go to the last
section, at the bottom of the script, and run. This will get you dates and plots of the mission you are interested in.

Unfortunately, I was not able to include all of the necesary ephemeride files. You will need to go to JPL's Navigation Ancillary Information 
Facility (NAIF) to download the remaining files. 
  - The latest leap seconds: https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/latest_leapseconds.tls
  - JPL Planetary Ephemerides: https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440.bsp

Additionally, to be able to utilize these files, you will also need to download the JPL's SPICE (Spacecraft, Planet, Instrament, C-matrix, Events) 
toolkit. You can get the tool box for MATLAB here: https://naif.jpl.nasa.gov/naif/toolkit_MATLAB.html. Within the mice.zip download, you will need 
to retreive the following files:
  - mice/src/mice/cspice_furnsh.m
  - mice/src/mice/cspice_kclear.m
  - mice/src/mice/cspice_spkezr.m
  - mice/src/mice/cspice_str2et.m
  - mice/lib/mice.mexw64

If you wish to generate a comparison between my work and NASA's Galileo mission, you will need to download the follwing files from JPL's NAIF
  - Galileo's Interplanetary Cruise (s970311a.bsp): https://naif.jpl.nasa.gov/pub/naif/GLL/kernels/spk/s970311a.bsp
  - Galileo's Primary Tour (s980326a.bsp): https://naif.jpl.nasa.gov/pub/naif/GLL/kernels/spk/s980326a.bsp


If there are any questions about running the provided MATLAB scripts, please reach out to me at geneluevano@gmail.com.

If you would like to see me give my 15 minute final presentation, you can find the recording here:
  - https://youtu.be/2_SKgY0i2k8
