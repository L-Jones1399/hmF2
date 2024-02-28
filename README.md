THINGS TO BE DONE
1. Ensure that the file StationCodes.txt is in the save_directory, where the other files will be generated
2. Ensure that all the COSMIC files (the files that you specify as COSMIC_Files in the config file) is also in the save_directory, where the output files will be generated

config.dat - config file that gives the location of where it saves/accesses files, along with a handful of other user changed variables

scripting.py - python script that acccesses the GIRO website and saves the information that we specify into a singular .dat file

Cosmic_Analysis.py - python script that compares every COSMIC observation, trying to find a match in the GIRO observations that is within a specified radius of it, as well as the closest in time. Can generate a range of graphs based on what function is run. Currently it generates a scatter plot of each station broken into 6-hr bins, as well as two polymesh plots of latitude vs Time and Geomagnetic activity respectfully

StationList_2023.03.01_to_2023.04.30.dat - The file that scripting.py should give you if you're running it over these two months for hmF2. This file is not called, it is just here for reference.

COSMIC2023.dat - A file containing the COSMIC observations that are used for my graphs. This file, and all the other COSMIC ones you specified in COSMIC_Files (in the config file) are called.

StationCodes.txt - A text file containing all GIRO station's codes, which are fed into scripting.py. This file is called
