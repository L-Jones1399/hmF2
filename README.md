config.dat - config file that gives the location of where it saves/accesses files, along with a handful of other user changed variables

scripting.py - python script that acccesses the GIRO website and saves the information that we specify into a singular .dat file

Cosmic_Analysis.py - python script that compares every COSMIC observation, trying to find a match in the GIRO observations that is within a specified radius of it, as well as the closest in time. Can generate a range of graphs based on what function is run. Currently it generates a scatter plot of each station broken into 6-hr bins, as well as two polymesh plots of latitude vs Time and Geomagnetic activity respectfully

StationList.dat - The file that CleanUp.py should give you after running it

COSMIC2023.dat - A file containing the COSMIC observations that are used
