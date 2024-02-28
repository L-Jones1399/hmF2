#Name: scripting.py
#Author: Liam Jones
#Purpose: Get data from the GIRO website and download it, ready to be analysed
#Last Modified: 28/02/2024


#TO INSTALL PACKAGES : py -m pip install _____

#Imports
import ast
import configparser
import urllib.request, urllib.error
import re

'''
TODO - Stop using Selenium, use API calls, e.g. https://lgdc.uml.edu/fastchar/getbest?ursiCode=AL945&charName=hmF2&fromDate=2023.03.01&toDate=2023.03.02, will get more
instructions later
'''

def CleanData(save_directory, contents, StationCode):
    #print("Clean")
    #Removes unnecessary lines of data
    #Locating the position of the data we want to keep (Lat, Long, and hmF2)
    start_indices = [i for i, elem in enumerate(contents) if "# Location" in str(elem)]
    end_indices = [i for i, elem in enumerate(contents) if "# Time" in str(elem)]
    #Deleting the rest of the unnecessary information
    #print(txtfiles)
    del contents[start_indices[0] + 1 :end_indices[-1]]
    del contents[0:start_indices[0]]
    #Discard the final line
    del contents[-1]

    #removes excess spaces from header
    excess_space = contents[1].split()
    contents[1] = " ".join(excess_space)
    #Grab the Station Name from the header before deleting it
    pattern = r'URSI-Code\s+(.*)'
    match = re.search(pattern, contents[0])
    # Get the matched substring after "URSI-Code"
    StationName = match.group(1)
    StationName = StationName[7:]
    #print(StationName)
    #Split all strings into their specific data types, 2d array
    data_table = [element.split() for element in contents]
    #print('yep')
    #remove the unnecessary data in columns 2 and 4, while keeping the lat and lon only
    data_table = [data_table[0][4:5] + data_table[0][6:7]] + [[row[i] for i in range(len(row)) if i not in (1, 3)] for row in data_table[2:]]
    #print('done')
    SaveData(save_directory, data_table, StationName, StationCode)

def SaveData(save_directory, data_table, StationName, StationCode):
    #print("Save")
    filepath = save_directory + f'\\StationList_{formatted_start_date}_to_{formatted_end_date}.dat'
    config = configparser.ConfigParser()
    config.read(filepath) #Multiline strings need to be tabbed.

    section_name = f'({StationCode}) {StationName}'

    #Writing data array to a string to add to the .dat file
    array_string = '\n'.join([' '.join(map(str, row)) for row in data_table[1:]])

    #Defining what we want to add to the .dat file
    section_content = {
        'Lat': data_table[0][0],
        'Lon': data_table[0][1],
        'Data': array_string
    }

    #Writing relevant information to the .dat file
    config[section_name] = section_content
    with open(filepath, 'w') as config_file:
        config.write(config_file)



if __name__ == "__main__": 
    print("Running the API caller")

    #Reading the config file using the configparser library
    config = configparser.ConfigParser()
    config.read('config.dat')
    #Sets the save directory based on the config file
    save_directory = config['Settings']['save_directory']
    #Gets the start and end date from the config file
    start_date = ast.literal_eval(config['Settings']['start_date'])
    formatted_start_date = f"{start_date[0]:04d}.{start_date[1]:02d}.{start_date[2]:02d}"
    end_date = ast.literal_eval(config['Settings']['end_date'])
    formatted_end_date = f"{end_date[0]:04d}.{end_date[1]:02d}.{end_date[2]:02d}"

    #Get the wanted variable name from the config file (e.g hmF2)
    Wanted_Variable = config['Settings']['Relevant_Variable']

    '''
    This loads a premade text file that contains every station code. 
    '''
    #Load txt file with all station codes
    filepath = save_directory + "\\StationCodes.txt"
    file = open(filepath, 'r')
    StationCodes = file.readlines()
    #Contains a list of all Station codes
    StationCodes = [string.rstrip('\n') for string in StationCodes]
    file.close()

    #Formatting the url
    url = "https://lgdc.uml.edu/fastchar/getbest?ursiCode=AL945&charName=hmF2&fromDate=2023.03.01&toDate=2023.03.02"
    #Modify the url
    pattern_url = r'ursiCode=([^&]+)'
    pattern_sdate = r'fromDate=([^&]+)'
    pattern_edate = r'toDate=\d{4}\.\d{2}\.\d{2}'
    #Change what variable we want to pull
    start_index = url.find('charName=') + len('charName=')
    end_index = url.find('&fromDate')
    #sDate = '2020.01.01'
    #eDate = '2023.06.06'

    #Modify end and start dates outside loop because only need to do once
    new_url = re.sub(pattern_sdate, f'fromDate={formatted_start_date}', url )
    new_url = re.sub(pattern_edate, f'toDate={formatted_end_date}', new_url )
    new_url = new_url[:start_index] + Wanted_Variable + new_url[end_index:]
    for i, code in enumerate(StationCodes):
        #print(code)
        #modify the url
        new_url = re.sub(pattern_url, f'ursiCode={code}', new_url)

        #Using urllib to open the url
        #Error handling - will keep trying to download the data until it reaches three attempts
        attempt = 1
        max_attempts = 3
        while attempt <= max_attempts:
            try:
                page = urllib.request.urlopen(new_url)
                break
            except urllib.error.URLError as e:
                print(e)
                print('\n On station: ' + code + ' retrying')
                attempt += 1
        #page = urllib.request.urlopen(new_url)
        html_content = page.read()
        decoded_content = html_content.decode('utf-8')
        lines = decoded_content.split('\n')
        '''
        Now we need to clean up the data so it can be saved in the .dat file
        '''
        if len(lines) < 25:
            continue
        CleanData(save_directory, lines, code)
    print("Done")
