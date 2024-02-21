#Name: Cosmic_Analysis.py
#Author: Liam Jones
#Purpose: Analyse COSMIC observations and compare with the relevant ground stations
#Last Modified: 09/02/2024
#04/01/2024: Code cleans up .dat files that we need
#24/01/2024: Added station data .dat file in, as well as performing spatial and temporal analysis
           # to determine which point is closer in both space and time. Also added basic scatter plot
#25/01/2024: Changed scatter plot to contain discretised bins, which the points are coloured differently
           # depending on how far away they are.
#08/02/2024: Changed the way COSMIC data is read, so it now can read in multiple files. Also added a filter
           # to only select the COSMIC data within a specific date range
#15/02/2024: Added another function called SingleStation that plots GIRO and COSMIC hmF2 for one station only
#21/02/2024: Added a function that plots the data from SingleStation as multiple graphs over the local time

#Imports
import configparser
from datetime import datetime, timedelta
import numpy as np
from math import radians, sin, cos, sqrt, atan2, pi
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from matplotlib import colors
import ast
from scipy import stats
from sklearn import linear_model, datasets

def Initialise():

    #Setting up the save_directory
    config = configparser.ConfigParser()
    config.read('config.dat')
    #Sets the save directory based on the config file
    save_directory = config['Settings']['save_directory']
    #Set up the max distance that still results in a match
    max_distance = float(config['Settings']['max_distance'])
    #Get the file location of the GIRO measurements
    StationList = config['Settings']['StationList']
    #Get the names of the COSMIC data files
    COSMIC_Files = ast.literal_eval(config['Settings']['COSMIC_Files'])
    #Get the start date
    start_date = ast.literal_eval(config['Settings']['start_date'])
    end_date = ast.literal_eval(config['Settings']['end_date'])
    #print(max_distance)
    return save_directory, max_distance, StationList, COSMIC_Files, start_date, end_date

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
def ClosestTime(StationName, COSMIC_Date, StationConfig):
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
#Name: ReadCOSMIC
#Purpose: Read in the multiple COSMIC.dat files to be analysed
#Inputs: save_directory (str) - location of save directory, COSMIC_Files (list) - names of COSMIC files
#Output: Cleaned_array (np.array) - array of all COSMIC measurements, including date, location and hmF2 value
def ReadCOSMIC(save_directory, COSMIC_Files):
    #Read COSMIC data .dat files
    for i, files in enumerate(COSMIC_Files):

        filepath = save_directory + "\\" + files
        file = open(filepath, 'r')
        contents = file.readlines()
        file.close()

        #Deconstructing each line into each data point
        contents = [element.split() for element in contents]
        if i == 0:
            # Delete the irrelevant columns
            Cleaned_array = np.delete(contents, [5,6,7], axis=1)
            Cleaned_array = (Cleaned_array.astype(float)) #ALL ELEMENTS MUST BE SAME DATA TYPE
        else:
            Cleaned_array2 = np.delete(contents, [5,6,7], axis=1)
            Cleaned_array2 = (Cleaned_array2.astype(float)) #ALL ELEMENTS MUST BE SAME DATA TYPE
            Cleaned_array = np.vstack((Cleaned_array, Cleaned_array2))

    return Cleaned_array



'''
Importing the Station data from the .dat file
'''
#Name: ImportStationData
#Purpose: Import the station hmF2 data from its .dat fiel
#Inputs: save_directory (str), Station_filename (str) - filename of the .dat file
#Outputs: StationConfig (ConfigParser) - Config file with all GIRO measurements, StationNames - list of all stations
def ImportStationData(save_directory, StationList):
    StationFilepath = save_directory + '\\' + StationList
    StationConfig = configparser.ConfigParser()
    StationConfig.read(StationFilepath)

    #Grab all section headers (individual stations)
    Station_Names = StationConfig.sections()
    return StationConfig, Station_Names

'''
Comparing both point's lat and lon to find match using Haversine's formula
'''
#Name: FindMatches
#Purpose: Find the closest GIRO station to each COSMIC measurement, then the closest hmF2 measurement in time
#Inputs: save_directory (str), Cleaned_array (ndarray), max_distance (int)
#Outputs: COSMIC_match_hmF2 (list) - list of all hmF2 from COSMIC that had a match. GIRO_match_hmF2 (list) - list of corresponding matches from GIRO
        # Distances_List (list) - List containing the distance each point was from each other
def FindMatches(save_directory, Cleaned_array, max_distance, StationList):
    StationConfig, Station_Names = ImportStationData(save_directory, StationList)
    
    #Creating blank lists that get appended to further down
    COSMIC_Match_hmF2 = []
    GIRO_Match_hmF2 = []
    Distances_List = []
    print_index = 0

    #Go through every line in the COSMIC data
    for i, lines in enumerate(Cleaned_array):
        """
        Do it for 1 day, then loop
        """

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
                #Erase the previous match and replace it with this one
                Smallest_Distance = Distance
                Closest_Station = stations
                MatchFound = True
        #If there is a point in the specified radius
        if MatchFound == True: 
            
            

            #Convert the COSMIC date into a readable date format
            COSMIC_Date = datetime((Cleaned_array[i][0].astype(int)), (Cleaned_array[i][1].astype(int)), (Cleaned_array[i][2].astype(int)), (Cleaned_array[i][3].astype(int)), (Cleaned_array[i][4].astype(int))).strftime("%Y-%m-%dT%H:%M")
            COSMIC_Date = datetime.fromisoformat(COSMIC_Date)
            
            #Should only allow every 100th match to be printed, to save on terminal space
            print_index += 1
            if print_index % 100 ==0:  
                print("Current Date: " + str(COSMIC_Date) + f" ({str(i)})")

            #Call function ClosesetTime, which finds and returns the closest time from a given station, and the hmF2 value at that time
            Closest_Date, Closest_Date_Index, Closest_Date_hmF2, StationName = ClosestTime(Closest_Station, COSMIC_Date, StationConfig)

            #Append the hmF2 and distance data to the lists to be graphed
            COSMIC_Match_hmF2.append(Cleaned_array[i][5])
            GIRO_Match_hmF2.append(Closest_Date_hmF2)
            Distances_List.append(Smallest_Distance)

    return COSMIC_Match_hmF2, GIRO_Match_hmF2, Distances_List


'''
Create a scatter plot of COSMIC and GIRO observations
'''
#Name: ColourScatter
#Purpose: Plot a scatter plot of COSMIC and GIRO hmF2 against each other. It also colour codes the points based on their proximity to each other
#Inputs: COSMIC_Match_hmF2 (list), GIRO_Match (list), max_distance (int), Distances_list (list)
#Outputs: Scatter plot
def ColourScatter(COSMIC_Match_hmF2, GIRO_Match_hmF2, max_distance, Distances_List):

    #calculate the correlation coefficient
    correlation_coefficient, _ = pearsonr(COSMIC_Match_hmF2, GIRO_Match_hmF2)

    #Creating the colour-coding
    #bin_edges = [0,20, 40, 60, 80, 100]
    #Sets the edges of the bin based on the max distance
    bin_edges = np.arange(0,max_distance + 1,20)

    #Sets the number of bins
    NoOfBins = np.linspace(1, 5, 5)

    #Find an appropriate colour map
    colormap = plt.get_cmap('RdYlBu')
    #Discretise the distance into the bins specified above
    bin_indices = np.digitize(Distances_List, bin_edges, right=True)
    #Discretise the colours from the map
    norm = colors.BoundaryNorm(np.arange(0.5, max_distance, 1), colormap.N)

    #Create the scatter plot
    scatter = plt.scatter(COSMIC_Match_hmF2, GIRO_Match_hmF2, c=bin_indices, cmap=colormap, norm=norm, edgecolors='black', linewidth=0.5, label = "Data Points")
    print(len(COSMIC_Match_hmF2))
    #Add colour bar
    cbar = plt.colorbar(ticks=NoOfBins)

    #Change the colour bar labels to their respective distances
    TickLabels = [f'{int(bin_edges[i-1])}-{int(bin_edges[i])}' for i in range(1, len(bin_edges))]
    #TickLabels = ['0-20', '20-40', '40-60', '60-80', '80-100']
    cbar.set_ticklabels(TickLabels)
    cbar.set_label("Distance between Observations (km)")

    #Adding the 1:1 line
    plt.plot([0, 600], [0, 600], color='red', linestyle='--', label='1:1 Line')

    # Annotate the correlation coefficient on the plot
    plt.annotate(
        f'Correlation: {correlation_coefficient:.2f}',
        xy=(0.5, 0.95),
        xycoords='axes fraction',
        ha='center',
        fontsize=10,
        color='blue'
    )

    #Formatting the scatter plot
    plt.xlabel("COSMIC Data")
    plt.ylabel("GIRO Data")
    plt.xlim(0, 600)
    plt.ylim(0, 600)
    plt.legend()
    plt.title("Scatter Plot of GIRO and COSMIC data within " + str(max_distance) + " km")
    plt.show()

'''
This section creates scatter plots for different bins, e.g 0-20, 0-40 etc.
'''
#Name: SpecificScatter
#Purpose: Create scatter plots for smaller sections of the data
#Inputs: COSMIC_Match_hmF2 (list), GIRO_Match (list), max_distance (int), Distances_list (list)
#Outputs: Scatter plot of smaller distances
def SpecificScatter(COSMIC_Match_hmF2, GIRO_Match_hmF2, Distances_List, max_distance):
    specified_bin_index = [1,2,3,4,5]

    bin_edges = np.linspace(0, max_distance, 6).astype(int)
    bin_indices = np.digitize(Distances_List, bin_edges, right=True)

    COSMIC_Match_hmF2 = np.asarray(COSMIC_Match_hmF2)
    GIRO_Match_hmF2 = np.asarray(GIRO_Match_hmF2)
    Distances_List = np.asarray(Distances_List)

    selected_COSMIC = COSMIC_Match_hmF2[np.isin(bin_indices, specified_bin_index)]
    selected_GIRO = GIRO_Match_hmF2[np.isin(bin_indices, specified_bin_index)]
    selected_bin_indices = bin_indices[np.isin(bin_indices, specified_bin_index)]

    correlation_coefficient, _ = pearsonr(selected_COSMIC, selected_GIRO)


    # Create the scatter plot
    for bin_index in specified_bin_index:
        bin_COSMIC = selected_COSMIC[selected_bin_indices == bin_index]
        bin_GIRO = selected_GIRO[selected_bin_indices == bin_index]
        plt.scatter(bin_COSMIC, bin_GIRO, label=f'Bin {bin_index}')


    # Customize the plot
    # Annotate the correlation coefficient on the plot
    plt.annotate(
        f'Correlation: {correlation_coefficient:.2f}',
        xy=(0.5, 0.95),
        xycoords='axes fraction',
        ha='center',
        fontsize=10,
        color='blue'
    )
    plt.plot([0, 600], [0, 600], color='red', linestyle='--', label='1:1 Line')
    plt.title('Scatter Plot for 0-100 km')
    plt.xlabel('COSMIC hmF2')
    plt.ylabel('GIRO hmF2')
    plt.legend()
    plt.xlim(0, 600)
    plt.ylim(0, 600)
    plt.grid(True)

    # Show the plot
    plt.show()

#Name: RMSE_Scatter
#Purpose: Generate a plot of Correlation Coefficient and RMSE vs distance
#Inputs: COSMIC_Match_hmF2 (list), GIRO_Match (list), max_distance (int), Distances_list (list)
#Outputs: Plot of RMSE vs distance, Plot of Correlation Coeffient vs distance
def RMSE_Scatter(COSMIC_Match_hmF2, GIRO_Match_hmF2, Distances_List, max_distance):
    print("RMSE")
    #Do RMSE from 0-20, 20-40 etc
    BinLimits = np.arange(0,max_distance + 1,20)
    bin_indices = np.digitize(Distances_List, BinLimits, right=True)
    COSMIC_Match_hmF2 = np.asarray(COSMIC_Match_hmF2)
    GIRO_Match_hmF2 = np.asarray(GIRO_Match_hmF2)

    rmse = []
    distance = []
    Correlation_Co = []
    #Loop it for all different bins
    for i in range(1, len(BinLimits)):
        specified_bin_index = i
        #Select all those that have that particular bin
        selected_COSMIC = COSMIC_Match_hmF2[bin_indices == specified_bin_index]
        selected_GIRO = GIRO_Match_hmF2[bin_indices == specified_bin_index]
        #Calculate RMSE
        # Calculate squared differences
        squared_diff = (np.asarray(selected_COSMIC) - np.asarray(selected_GIRO)) ** 2

        # Calculate mean squared difference
        mse = np.mean(squared_diff)

        # Calculate RMSE
        rmse.append(np.sqrt(mse))
        distance.append(f'{int(BinLimits[i-1])}-{int(BinLimits[i])}')

        #Calculate the correlation coefficient
        #calculate the correlation coefficient
        if len(selected_COSMIC) >= 2:
            
            correlation_coefficient, _ = pearsonr(selected_COSMIC, selected_GIRO)
            Correlation_Co.append(correlation_coefficient)
        else:
            #Checking if there are enough values in the selected bin to calculate
            Correlation_Co.append(0)

    plt.plot(distance, Correlation_Co)
    plt.savefig(save_directory + "\\CorrCo.png")
    plt.xlabel("Bins")
    plt.ylabel("Correlation Coefficient")
    plt.show()
    plt.plot(distance, rmse)
    plt.savefig(save_directory + "\\RMSE.png")
    plt.xlabel("Bins")
    plt.ylabel("RMSE")
    plt.show()

def SingleStation(save_directory, max_distance, StationList, Cleaned_array):
    print("Single station")
    '''
    In the two month for one station, go through every point and find the closest COSMIC spatially and temporally
    '''
    #Importing the station data
    StationConfig, Station_Names = ImportStationData(save_directory, StationList)
    #Taking the first station's data
    StationData = np.array(StationConfig[Station_Names[11]]["data"].split()).reshape(-1,2)

    #Add a Column to the COSIMC data that contains their row index and zeros
    num_rows = Cleaned_array.shape[0]
    # Create an array of row indices
    row_indices = np.arange(num_rows).reshape(-1, 1)
    Cleaned_array = np.hstack((Cleaned_array, row_indices))
    num_rows = Cleaned_array.shape[0]
    # Create an array of zeros with the same number of rows as arr_with_indices
    zeros_column = np.zeros((num_rows, 1), dtype=Cleaned_array.dtype)
    # Add the zeros column to the right of arr_with_indices
    Cleaned_array = np.hstack((Cleaned_array, zeros_column))

    #Filter these points so only the ones within 100 km are kept
    GIRO_Lat = float(StationConfig[Station_Names[11]]["lat"])
    GIRO_Lon = float(StationConfig[Station_Names[11]]["lon"])
    ValidPoints = []
    ValidPointsDistance = []
    GIRO_hmF2 = []

    #Filter all COSMIC points to keep only ones within 100 km
    for x in range(len(Cleaned_array)):
        COSMIC_Lat = Cleaned_array[x][6]
        COSMIC_Lon = Cleaned_array[x][7]
        DistBool, distance = is_within_radius(GIRO_Lat, GIRO_Lon, COSMIC_Lat, COSMIC_Lon, max_distance)
        if DistBool:
            ValidPoints.append(Cleaned_array[x])
            ValidPointsDistance.append(distance)
    #print(StationData[:, 1])
    ValidPoints = np.asarray(ValidPoints) #All (clos enough) COSMIC point's full 10 rows
    Weighted_hmF2 = []
    Weighted_Distances = []
    ValidGIRODates = []
    #Looping through all GIRO observations
    for i in range(len(StationData)):
        #print(StationData[i])
        #Get the COSMIC data
        #Get all points that are within the 30 min window either side
        Time_Filtered, date_filtered, distance_filtered = Time_Window(ValidPoints, StationData[i][0], ValidPointsDistance)
        if len(Time_Filtered) == 0:
            continue
        #print(date_filtered)
        #print(StationData[i][0])
        MP_Date = date_filtered #All the matching COSMIC point's dates
        Time_Filtered = np.asarray(Time_Filtered) #All Matching COSMIC point's full 10 rows
        #print(Time_Filtered[:, 8])
        MP_hmF2 = Time_Filtered[:, 5] #All the Matching COSMIC point's hmF2 values
        GIRO_hmF2.append(float(StationData[i][1])) #The GIRO Point's hmF2 value

        '''
        Oversampling control
        '''
        #Adding 1 to the duplicate flag
        for b in Time_Filtered:
            Cleaned_array[int(b[8])][9] = Cleaned_array[int(b[8])][9] + 1
        
        #print(len(Time_Filtered))
        #print(Time_Filtered[0])
        #print(Time_Filtered[-1])
        
        #Filter these points so only the ones within 100 km are kept
        # GIRO_Lat = float(StationConfig[Station_Names[0]]["lat"])
        # GIRO_Lon = float(StationConfig[Station_Names[0]]["lon"])

        # #Iterate through all filtered times to see which are within the specified boundary
        # Matching_PointsDistance = []
        # MP_hmF2 = []
        # MP_Date = []
        # for y in range(len(Time_Filtered)):
        #     COSMIC_Lat = Time_Filtered[y][6]
        #     COSMIC_Lon = Time_Filtered[y][7]
        #     DistBool, distance = is_within_radius(GIRO_Lat, GIRO_Lon, COSMIC_Lat, COSMIC_Lon, max_distance)
        #     if DistBool == True:
        #         Matching_PointsDistance.append(distance)
        #         MP_hmF2.append(Time_Filtered[y][5])
        #         MP_Date.append(date_filtered[y])
        weighted_average, weighted_dist_avg = calculate_weight(MP_Date, distance_filtered, MP_hmF2, StationData[i][0])
        if len(MP_Date) >= 2:
            print(len(MP_Date))
        Weighted_hmF2.append(weighted_average)
        Weighted_Distances.append(weighted_dist_avg)
        ValidGIRODates.append(StationData[i][0])
        #print(weighted_average)

        #Then either weight the remaining, or choose the best one
    #After going through all data points, plot a graph of Weighted COSMIC vs GIRO
    GIRO_hmF2 = np.asarray(GIRO_hmF2)

    '''
    Do Some Oversampling Control - go through weighted hmF2. If more than 1 same, choose one that is smallest distance and delete rest
    '''
    # Create a dictionary to store the minimum distance for each unique value
    min_distance_dict = {}

    for value, distance, GIROHMF2, GIRO_Dates in zip(Weighted_hmF2, Weighted_Distances, GIRO_hmF2, ValidGIRODates):
        if value in min_distance_dict:
            # Update minimum distance if current distance is smaller
            if min(min_distance_dict[value][0], distance):
                min_distance_dict[value] = [abs(min_distance_dict[value][0]- distance), GIROHMF2, GIRO_Dates]
        else:
            # Add value to dictionary if it's encountered for the first time
            min_distance_dict[value] = [distance, GIROHMF2, GIRO_Dates]

    filtered_values = np.asarray(list(min_distance_dict.keys()))

    print(len(set(Weighted_hmF2)))
    print(len(set(filtered_values)))
    filtered_GIRO_hmF2 = np.asarray(list(min_distance_dict.values()))
    filtered_GIRO_Dates = filtered_GIRO_hmF2[:, 2]
    filtered_GIRO_hmF2 = filtered_GIRO_hmF2[:, 1].astype(float)
    



    '''
    Removing outliers, stdev method
    '''
    slope, intercept, r_value, p_value, std_err = stats.linregress(filtered_GIRO_hmF2, filtered_values)
    line = slope * filtered_GIRO_hmF2 + intercept
    # Calculate residuals
    residuals = filtered_values - line

    # Define threshold for outliers (e.g., beyond 2 standard deviations)
    threshold = 2 * np.std(residuals)

    x_filtered2 = []
    y_filtered2 = []
    Final_Date_filtered = []
    # Remove outliers
    #Keep dates attached to the GIRO points as well
    for i in range(len(filtered_GIRO_hmF2)):
        if np.abs(residuals[i]) < threshold:
            x_filtered2.append(filtered_GIRO_hmF2[i])
            y_filtered2.append(filtered_values[i])
            Final_Date_filtered.append(filtered_GIRO_Dates[i])
    x_filtered = filtered_GIRO_hmF2[np.abs(residuals) < threshold]
    y_filtered = filtered_values[np.abs(residuals) < threshold]
    #filtered_data = PointstoGraph[(PointstoGraph[:, 0]< threshold) & (PointstoGraph[:, 1] <= threshold)]

    TimeBins(x_filtered2, y_filtered2, Final_Date_filtered, GIRO_Lat, GIRO_Lon)


    #Recalculate the line of best fit
    slope_cleaned, intercept_cleaned, _, _, _ = stats.linregress(x_filtered, y_filtered)
    line_of_best_fit_cleaned = slope_cleaned * x_filtered + intercept_cleaned

    '''
    Using Sklearn's RANSAC algorithm
    '''
    # ransac = linear_model.RANSACRegressor()
    # print(y_filtered)
    # #print(x_filtered.reshape(-1, 1))
    # ransac.fit(x_filtered.reshape(-1, 1), y_filtered)
    # line_X = np.arange(x_filtered.min(), x_filtered.max())[:, np.newaxis]
    # line_y_ransac = ransac.predict(line_X)
    # RANSAC_slope = ransac.estimator_.coef_[0]
    # RANSAC_intercept = ransac.estimator_.intercept_

    correlation_coefficient, _ = pearsonr(x_filtered, y_filtered)

    rmse = np.sqrt(np.mean((y_filtered - line_of_best_fit_cleaned)**2))

    plt.scatter(x_filtered, y_filtered)
    plt.plot(x_filtered, line_of_best_fit_cleaned, color='red', label='Line of Best Fit')
    #plt.plot(line_X, line_y_ransac, color='green', label='RANSAC')
    print()


    # Add text annotations for correlation coefficient and RMSE
    plt.figtext(.13, .74, f'Correlation Coefficient: {correlation_coefficient:.2f}')
    plt.figtext(0.13, .70, f'RMSE: {rmse:.2f}')
    plt.figtext(0.13, .66, f'Line Equation: {str(round(slope_cleaned, 2)) + "*x + " + str(round(intercept_cleaned, 2))}')
    #plt.figtext(0.13, .62, f'RANSAC Line Equation: {str(round(RANSAC_slope, 2)) + "*x + " + str(round(RANSAC_intercept, 2))}')

    # Add labels and legend
    plt.xlabel('GIRO hmF2')
    plt.ylabel('Weighted COSMIC hmF2')
    plt.title(f"Scatter plot for {Station_Names[11]}")
    plt.legend()

    # Show plot
    plt.show()

    

def Time_Window(COSMIC_Data, GIRO_date, distance, window = 20):
    #Get the datetime of the GIRO main point
    known_datetime = datetime.fromisoformat(GIRO_date).replace(tzinfo=None)
    #known_datetime = datetime.strptime(GIRO_date, "%Y-%m-%d %H:%M:%S")  # Convert known date to datetime object
    start_time = known_datetime - timedelta(minutes=window)  # Calculate start time of time window
    end_time = known_datetime + timedelta(minutes=window)    # Calculate end time of time window
    
    filtered_data = []
    filtered_date = []
    filtered_distance = []
    #Iterate through every COSMIC data point, only keeping those that fall within the 1 hr period
    for i, row in enumerate(COSMIC_Data):
        #row_datetime = datetime.strptime(row[0], "%Y-%m-%d %H:%M:%S")  # Assuming datetime is in the first column
        row_datetime = datetime((row[0].astype(int)), (row[1].astype(int)), (row[2].astype(int)), (row[3].astype(int)), (row[4].astype(int))).strftime("%Y-%m-%dT%H:%M")
        row_datetime = datetime.fromisoformat(row_datetime)
        if start_time <= row_datetime <= end_time:
            filtered_data.append(row)
            filtered_date.append(row_datetime)
            filtered_distance.append(distance[i])
    
    return filtered_data, filtered_date, filtered_distance

def calculate_weight(date_list, distance_list, value_list, target_date):
    weights = []
    distance_weights = []
    for date, distance in zip(date_list, distance_list):
        #date2 = datetime.fromisoformat(target_date)
        date_diff = abs(date - datetime.fromisoformat(target_date).replace(tzinfo=None))
        distance_weight = 1 / (distance + 1)
        time_weight = (1 / (date_diff.total_seconds() + 1))
        weight = distance_weight*time_weight
        weights.append(weight)
        distance_weights.append(distance_weight)
    
    weighted_sum = sum(weight * value for weight, value in zip(weights, value_list))
    total_weight = sum(weights)
    
    weight_dist_avg = sum(weight * value for weight, value in zip(distance_weights, distance_list))/sum(distance_weights)

    if total_weight == 0:
        return 0  # Avoid division by zero error
    
    weighted_average = weighted_sum / total_weight
    return weighted_average, weight_dist_avg

#Name: TimeBins
#Purpose: Break the single station's GIRO points down by local time of day
#Inputs: x_filtered (list), y_filtered (list), date_Filtered (list), GIRO_Lat (float), GIRO_Lon (float)
#Outputs: 4 scatter plots of GIRO vs COSMIC hmF2 readings
def TimeBins(x_filtered, y_filtered, date_filtered, GIRO_Lat, GIRO_Lon):
    print("TimeBins")
    '''
    Divide into local time sectors, 6hrs. I want date of GIRO point
    '''
    Date_vals = [datetime.fromisoformat(entry) for entry in date_filtered]

    #Converting to local time, the 15 degree way
    #Convert the lon to -180 to 180
    GIRO_Lon = (GIRO_Lon + 180) % (2 * 180) - 180
    Time_diff = round(GIRO_Lon/15)

    modified_dates_list = [date - timedelta(hours=Time_diff) for date in Date_vals]
    #Take the hours from this, and divide the x and y into bins based on this
    hour_list = [date.hour for date in modified_dates_list]

    bin_edges = [-1, 6, 12, 18, 24]
    bin_indices = np.digitize(hour_list, bin_edges, right=True)
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    for i, ax in enumerate(axes.flat):
        #Iterate through x_filtered, if bin_indicies = i+1 add it to new lists
        bin_x = []
        bin_y = []
        for b in range(len(x_filtered)):
            if bin_indices[b] == (i+1):
                bin_x.append(x_filtered[b])
                bin_y.append(y_filtered[b])

        ax.scatter(bin_x, bin_y)
        if i == 0:
            ax.set_title(f"Hours of {bin_edges[i]+1} to {bin_edges[i+1]}")
        else:
            ax.set_title(f"Hours of {bin_edges[i]} to {bin_edges[i+1]}")
        ax.set_xlabel("GIRO hmF2")
        ax.set_ylabel("COSMIC hmF2")

    plt.tight_layout()
    plt.title("Breakdown of time for Darwin")
    plt.show()

    print("DDD")




if __name__ == "__main__": 
    save_directory, max_distance, StationList, COSMIC_Files, start_date, end_date = Initialise()
    Cleaned_array = ReadCOSMIC(save_directory, COSMIC_Files)
    '''
    Select 1 day's data, then call
    '''
    #Only selects the dates in the COSMIC data that fall within the set range
    # Convert dates to integers for comparison
    start_date_int = start_date[0] * 10000 + start_date[1] * 100 + start_date[2]
    end_date_int = end_date[0] * 10000 + end_date[1] * 100 + end_date[2]

    # Convert date columns in data to integer format for comparison
    dates_int = Cleaned_array[:, 0] * 10000 + Cleaned_array[:, 1] * 100 + Cleaned_array[:, 2]

    # Filter the rows based on the date range criteria
    filtered_rows = Cleaned_array[(dates_int >= start_date_int) & (dates_int <= end_date_int)]

    SingleStation(save_directory, max_distance, StationList, filtered_rows)

    #COSMIC_Match_hmF2, GIRO_Match_hmF2, Distances_List = FindMatches(save_directory, filtered_rows, max_distance, StationList)
    #ColourScatter(COSMIC_Match_hmF2, GIRO_Match_hmF2, max_distance, Distances_List)
    #RMSE_Scatter(COSMIC_Match_hmF2, GIRO_Match_hmF2, Distances_List, max_distance)