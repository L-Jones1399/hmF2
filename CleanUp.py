#Task: Develop a routine to extract hmF2 from the GIRO website - Code to work on Linux system and can be used on the Bureau system for future applications
#Name: CleanUp.py
#Author: Liam Jones
#Purpose: Take text files previously generated and reformat them, finally saving all observations into a single .dat file
#Last Modified: 24/01/2024
#04/01/2024: Now the program saves all station data to a single .dat file, rather than individual .txt files. 
#24/01/2024: Edited slightly so now it removes the N or E from the lat and lon values

'''
If needed, I can combine this with scripting.py to not even save the .txt files, just write to the .dat
'''

#imports
import matplotlib.pyplot as plt
import os
import configparser
import re

'''
Grab each file ending in txt from same folder and format it correctly, discarding those with no data??
If it has no data the program just skips it, meaning it won't be altered and saved again
'''
#Gets the name of every text file in the specified folder
#Reads the config file
config = configparser.ConfigParser()
config.read('config.dat')
#Sets the save directory based on the config file
save_directory = config['Settings']['save_directory']
#Gets a list of every text file in the save directory
station_files = [f for f in os.listdir(save_directory) if f.endswith('.txt')] #station_files is now a list

#Iterate through each text file and format it appropriately
for i, txtfiles in enumerate(station_files):
    #Each line is saved as a new entry into the contents data list
    txtfilepath = save_directory + "\\" + txtfiles
    data = open(txtfilepath, 'r')
    #Reading each text file
    contents = data.readlines()
    data.close()
    #print(contents)

    #Checks if these is actually data in the file
    if len(contents) < 25:
        continue #Maybe here write a blank text file with Lat Long and nothing else


    #Removes unnecessary lines of data
    #Locating the position of the data we want to keep (Lat, Long, and hmF2)
    start_indices = [i for i, elem in enumerate(contents) if "# Location" in str(elem)]
    end_indices = [i for i, elem in enumerate(contents) if "#Time" in str(elem)]
    #Deleting the rest of the unnecessary information
    #print(txtfiles)
    del contents[start_indices[0] + 1 :end_indices[-1]]
    del contents[0:start_indices[0]]
    #Discard the final line
    del contents[-1]

    #removes excess spaces from header
    excess_space = contents[1].split()
    contents[1] = " ".join(excess_space)

    #Remove b' and \\n from each line
    contents = [element.lstrip("b'").lstrip("#").replace("\\n'", "") for element in contents]
    #print(contents[0:2])
    #print(contents[2])

    #Split all strings into their specific data types, 2d array
    data_table = [element.split() for element in contents]
    #print(data_table[1:-1])
    hmF2_vals = [float(entry[2]) for entry in data_table[2:]]
    #print(hmF2_vals[0:10])
    #print(type(hmF2_vals[0]))

    #remove the unnecessary data in columns 2 and 4, while keeping the lat and lon only
    data_table = [data_table[0][2:4]] + [[row[i] for i in range(len(row)) if i not in (1, 3)] for row in data_table[1:]]

    #Save each file
    #txtSavePath = save_directory + "\\Clean" + txtfiles
    #cleanedText = open(txtSavePath , 'w')
    #for row in data_table:
    #    row_str = " ".join(map(str, row))
    #    cleanedText.write(row_str + '\n')
    #cleanedText.close()

    #Append each file into the same .dat file
    filepath = save_directory + '\\StationList.dat'
    config = configparser.ConfigParser()
    config.read(filepath) #Multiline strings need to be tabbed.

    #Each section name is the station name
    section_name = txtfiles[:-4]
    #Retrieving the Lat and Lon values, disregarding the N or E at the end
    LatVal = re.sub(r'[^0-9.-]', '',data_table[0][0])
    LonVal = re.sub(r'[^0-9.-]', '',data_table[0][1])

    #Writing data array to a string to add to the .dat file
    array_string = '\n'.join([' '.join(map(str, row)) for row in data_table[2:]])

    #Defining what we want to add to the .dat file
    section_content = {
        'Lat': LatVal,
        'Lon': LonVal,
        'Data': array_string
    }

    #Writing relevant information to the .dat file
    config[section_name] = section_content
    with open(filepath, 'w') as config_file:
        config.write(config_file)


    '''
    The output .dat file will now contain all the measurements for each station.
    Each different section is a different ground station
    Each section contains the lat and lon, written in 0 to 360 degree lon
    The hmF2 data is tied to the timestamp, which is written in ISO 8601 standard.
    '''

    #plotting data
    # plt.plot(hmF2_vals)
    # plt.xlabel("Data points from March 1 to April 30 2023")
    # plt.ylabel("hmF2 (km)")
    # plt.show()
    print(section_name)
