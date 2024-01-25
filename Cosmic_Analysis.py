#Name: Cosmic_Analysis.py
#Author: Liam Jones
#Purpose: Analyse COSMIC observations and compare with the relevant ground stations
#Last Modified: 25/01/2024
#04/01/2024: Code cleans up .dat files that we need
#24/01/2024: Added station data .dat file in, as well as performing spatial and temporal analysis
           # to determine which point is closer in both space and time. Also added basic scatter plot
#25/01/2024: Changed scatter plot to contain discretised bins, which the points are coloured differently
           # depending on how far away they are.

#Imports
import configparser
from datetime import datetime
import numpy as np
from math import radians, sin, cos, sqrt, atan2, pi
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from matplotlib import colors


#Setting up the save_directory
config = configparser.ConfigParser()
config.read('config.dat')
#Sets the save directory based on the config file
save_directory = config['Settings']['save_directory']
#Set up the max distance that still results in a match
max_distance = float(config['Settings']['max_distance'])
print(max_distance)

#Defining functions

#Name: Haversine
#Purpose: Calculate distance between 2 lat and lon values on the surface of Earth
#Inputs: lat and lon of Primary point (degrees), lat and lon of Secondary point (degrees)
#Outputs: Distance between the two points (km)
def haversine(lat1, lon1, lat2, lon2):
    # Convert latitude and longitude from degrees to radians
    lat1, lon1, lat2, lon2 = map(radians, [lat1, lon1, lat2, lon2])

    # Adjust longitude values to the -180 to 180 range (but in radians)
    lon1 = (lon1 + pi) % (2 * pi) - pi
    lon2 = (lon2 + pi) % (2 * pi) - pi

    # Haversine formula
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = sin(dlat / 2)**2 + cos(lat1) * cos(lat2) * sin(dlon / 2)**2
    c = 2 * atan2(sqrt(a), sqrt(1 - a))

    # Radius of the Earth in kilometers (mean value)
    radius = 6378.137

    # Calculate the distance
    distance = radius * c

    return distance

#Name: is_within_radius
#Purpose: Run Haversine and return the distance between the two points
#Inputs: lat and lon of Primary point (degrees), lat and lon of Secondary point (degrees), desired radius of acceptance
#Outputs: Distance between the two points (km), boolean whether the point is acceptable
def is_within_radius(primary_lat, primary_lon, check_lat, check_lon, radius):
    # Check if the distance is within the specified radius
    distance = haversine(primary_lat, primary_lon, check_lat, check_lon)
    #print(distance)
    return distance <= radius, distance

#Name: ClosestTime
#Purpose: Find the GIRO data point that is the closest in time to the COSMIC date, for a given ground station
#Inputs: StationName (str) - the station of interest. COSMIC_Date (datetime) - the primary date of interest
#Outputs: Closest_Date (datetime) - the closest matching GIRO date. Closest_Date_Index (int) - the index of the date
         #Closest_Date_hmF2 (float) - the hmF2 value that corresponds to the date. StationName (str) - the name of the station
def ClosestTime(StationName, COSMIC_Date):
    #Grab all the Observations for that station
    #Arrange it into a 2d array with columns of Date and hmF2 reading
    StationData = np.array(StationConfig[StationName]["data"].split()).reshape(-1,2)
    #print(StationData[0][0]) #Date
    #print(StationData[0][1]) #hmF2
    #Extract the date from it
    #Need to ignore the time zone info, was UTC but now none
    Date_vals = [datetime.fromisoformat(entry).replace(tzinfo=None) for entry in StationData[:, 0]]
    #print(Date_vals[0])

    #Comparing the Station date to the COSMIC date to get the closest time
    Closest_Date = min(Date_vals, key = lambda date: abs(date-COSMIC_Date))
    #Gets the index value of the closest date
    Closest_Date_Index = Date_vals.index(Closest_Date)
    Closest_Date_hmF2 = float(StationData[Closest_Date_Index][1])
    # print(Closest_Date_Index)
    return Closest_Date, Closest_Date_Index, Closest_Date_hmF2, StationName




#Reading in COSMIC data
'''
COSMIC data - read one line at a time. First priority is location, so compare the lat and lon from both.
Can be more than one valid point - use closest or average?
Then match time - closest.

START TIGHT AND EXPAND
DO SCATTER PLOT AND SEE RELATIONSHIP, THEN RELAX RULES
get correlation coefficient, RMSE, R
'''
#Read COSMIC data .dat file
filepath = save_directory + "\\" + "COSMIC2023.dat"
file = open(filepath, 'r')
contents = file.readlines()
file.close()

#Deconstructing each line into each data point
contents = [element.split() for element in contents]

# Delete the irrelevant columns
Cleaned_array = np.delete(contents, [5,6,7], axis=1)
Cleaned_array = (Cleaned_array.astype(float)) #ALL ELEMENTS MUST BE SAME DATA TYPE
#print(Cleaned_array[5]) #Dont round :(
#Delete the rows that fall outside the wanted timeframe
for i, values in enumerate(Cleaned_array):
    if Cleaned_array[i][1] == 3 and Cleaned_array[i][2] == 1:
        Cleaned_array = Cleaned_array[i:, :]
        break

#Do the same with the desired end date, but making sure that you have all the datapoints
for i, values in enumerate(Cleaned_array):
    if (Cleaned_array[i][1] == 5 and Cleaned_array[i][2] == 1) and (Cleaned_array[i+1][1] == 5 and Cleaned_array[i+1][2] == 1):
        Cleaned_array = Cleaned_array[:i, :]
        break




#Cleaned_array2 = np.hstack((Cleaned_array2, Cleaned_array[:, -3:]))
#print(Cleaned_array2[:, :5])
#print(type(Cleaned_array2[0][1]))

#Change the date into a format that can be processed
#Cannot convert to date, all elements MUST BE SAME DATA TYPE 
# for i, rows in enumerate(Cleaned_array):
#     print(Cleaned_array[i][0])
#     print(type(int(Cleaned_array[i][0])))
#     #Replace the first column with this new date time, later we will remove all the old ones
#     Cleaned_array[i][0] = datetime((Cleaned_array[i][0].astype(int)), (Cleaned_array[i][1].astype(int)), (Cleaned_array[i][2].astype(int)), (Cleaned_array[i][3].astype(int)), (Cleaned_array[i][4].astype(int)))

'''
Importing the Station data from the .dat file
'''
StationFilepath = save_directory + '\\StationList.dat'
StationConfig = configparser.ConfigParser()
StationConfig.read(StationFilepath)

#Grab all section headers (individual stations)
Station_Names = StationConfig.sections()
#COSMIC_Date = (datetime((2023), (3), (1), (0), (30))).strftime("%Y-%m-%dT%H:%M")
#COSMIC_Date = datetime.fromisoformat(COSMIC_Date)
#print(COSMIC_Date)
#Closest_Date, Closest_Date_Index, Closest_Date_hmF2 = ClosestTime("(AL945) ALPENA", COSMIC_Date)

'''
Comparing both point's lat and lon to find match using Haversine's formula
'''

COSMIC_Match_hmF2 = []
GIRO_Match_hmF2 = []
Distances_List = []


#Go through every line in the COSMIC data
for i, lines in enumerate(Cleaned_array):

    MatchFound = False
    Smallest_Distance = 20040 #1/2 Circumference of Earth, nothing should be larger than this
    

    #Get the lat and lon for that point
    primary_point = (Cleaned_array[i][6], Cleaned_array[i][7])  # COSMIC Data point

    #Go through every station to find the one that is closest
    for index, stations in enumerate(Station_Names):
        check_point = (float(StationConfig[stations]['lat']), float(StationConfig[stations]['lon']))  # Station Data point

        #DistBool - Boolean, whether the point is within the acceptable range. Distance (km) - distance between the two points
        DistBool, Distance = is_within_radius(primary_point[0], primary_point[1], check_point[0], check_point[1], max_distance)

        #If there are two or more points within the acceptable range
        #If the distance is in the acceptable region, AND the distance is smaller than the one already stored
        if DistBool and Distance < Smallest_Distance:
            #print("The Ground Station is within the specified radius of the COSMIC Observation.")
            Smallest_Distance = Distance
            Closest_Station = stations
            MatchFound = True
    if MatchFound == True: 
        print("Iteration No. " + str(i+1))
        #print("Closest Station is: " + Closest_Station + " at " + str(Smallest_Distance))
        #Convert the COSMIC date into a readable date format
        COSMIC_Date = datetime((Cleaned_array[i][0].astype(int)), (Cleaned_array[i][1].astype(int)), (Cleaned_array[i][2].astype(int)), (Cleaned_array[i][3].astype(int)), (Cleaned_array[i][4].astype(int))).strftime("%Y-%m-%dT%H:%M")
        COSMIC_Date = datetime.fromisoformat(COSMIC_Date)
        #Call function ClosesetTime, which finds and returns the closest time from a given station, and the hmF2 value at that time
        Closest_Date, Closest_Date_Index, Closest_Date_hmF2, StationName = ClosestTime(Closest_Station, COSMIC_Date)
        COSMIC_Match_hmF2.append(Cleaned_array[i][5])
        GIRO_Match_hmF2.append(Closest_Date_hmF2)
        Distances_List.append(Smallest_Distance)
        #Now to break it down in time - find the closest ground station observation.
        #print("GIRO: " + str(Closest_Date_hmF2) + ' COSMIC: ' + str(Cleaned_array[i][5]))
        #print("GIRO: " + str(Closest_Date) + ' COSMIC: ' + str(COSMIC_Date))
        #print(StationName)

#Create a scatter plot of COSMIC and GIRO observations
correlation_coefficient, _ = pearsonr(COSMIC_Match_hmF2, GIRO_Match_hmF2)

#Creating the colour-coding
bin_edges = [0,20, 40, 60, 80, 100]
NoOfBins = np.linspace(1, 5, 5)
#cmap2 = colormaps['hsv']()
colormap = plt.get_cmap('RdYlBu')
bin_indices = np.digitize(Distances_List, bin_edges, right=True)
norm = colors.BoundaryNorm(np.arange(0.5, 6, 1), colormap.N)

#print(bin_indices)
print(Distances_List[0:10])
print(bin_indices[0:10])

scatter = plt.scatter(COSMIC_Match_hmF2, GIRO_Match_hmF2, c=bin_indices, cmap=colormap, norm=norm, edgecolors='black', linewidth=0.5, label = "Data Points")
#Add colour bar
cbar = plt.colorbar(ticks=NoOfBins)
#cbar = plt.colorbar(scatter, ticks=np.arange(0, len(bin_edges)), label='Distance Bins (km)')
#midpoints = (NoOfBins[:-1] + NoOfBins[1:]) / 2

# Set ticks to the midpoints
#cbar.set_ticks(midpoints)

#TickLabels = [f'{bin_edges[i-1]}-{bin_edges[i]}' for i in range(1, len(bin_edges))]
TickLabels = ['0-20', '20-40', '40-60', '60-80', '80-100']
cbar.set_ticklabels(TickLabels)
cbar.set_label("Distance between Observations (km)")

#cbar.set_ticklabels([f'{bin_edges[i]}' for i in range(0, len(bin_edges))])
#cbar.ax.set_yticklabels([f'{bin_edges[i]}' for i in range(0, len(bin_edges))])
#Adding the 1:1 line
plt.plot([min(COSMIC_Match_hmF2), max(COSMIC_Match_hmF2)], [min(COSMIC_Match_hmF2), max(COSMIC_Match_hmF2)], color='red', linestyle='--', label='1:1 Line')
# Annotate the correlation coefficient on the plot
plt.annotate(
    f'Correlation: {correlation_coefficient:.2f}',
    xy=(0.5, 0.95),
    xycoords='axes fraction',
    ha='center',
    fontsize=10,
    color='blue'
)
plt.xlabel("COSMIC Data")
plt.ylabel("GIRO Data")
plt.legend()
plt.title("Scatter Plot of GIRO and COSMIC data within " + str(max_distance) + " km")
plt.show()