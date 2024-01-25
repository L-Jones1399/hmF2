config.dat - config file that gives the location of where it saves/accesses files.

scripting.py - python script that acccesses the GIRO website and saves into text files the information that we specify for all stations

CleanUp.py - python script that takes the information previously generated and cleans up the text files to include only relevant information. Then saves all the text files into one .dat file for ease of access

Cosmic_Analysis.py - python script that compares every COSMIC observation, trying to find a match in the GIRO observations that is within a specified radius of it, as well as the closest in time. Also creates a scatter plot of the two

StationList.dat - The file that CleanUp.py should give you after running it

COSMIC2023.dat - A file containing the COSMIC observations that are used
