#Name: scripting.py
#Author: Liam Jones
#Purpose: Get data from the GIRO website and download it, ready to be analysed
#Last Modified: 03/01/2024


#TO INSTALL PACKAGES : py -m pip install _____

#Imports
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import Select
import configparser

#User defined variables (start date, end date)
StartDate = "2023-03-01 00:00"
EndDate = '2023-04-30 23:59'
RelevantData = "hmF2 -- Peak height F2-layer"

#Reading the config file using the configparser library
config = configparser.ConfigParser()
config.read('config.dat')
#Sets the save directory based on the config file
save_directory = config['Settings']['save_directory']


# Using Chrome to access web
driver = webdriver.Chrome()
# Open the website
driver.get('https://giro.uml.edu/didbase/scaled.php')

#Find the drop down list of all station names
location_menu = driver.find_element(By.NAME, "location")
#Select the drop down list of station names
selectLocation = Select(location_menu)

#Store each location in a list for later use
Location_List = []
for option in selectLocation.options:
    Location_List.append(option.text)

#Iterates through the list of station names
for StationLocation in Location_List:

    #Website retrieved again because we leave this page when the data is generated
    driver.get('https://giro.uml.edu/didbase/scaled.php')

    #select the station that you want (iterates)
    #Find the location drop down list again
    location_menu = driver.find_element(By.NAME, "location")
    #Select the drop down list
    selectLocation = Select(location_menu)
    #Select the station that corresponds to the iteration that we are on
    selectLocation.select_by_visible_text(StationLocation)

    #Find and replace the startdate with your own
    #Find the startdate text field
    startdate = driver.find_element(By.NAME, 'date_start')
    #Clear the text in it
    startdate.clear()
    #Paste in our own date
    startdate.send_keys(StartDate)

    #Find and replace the enddate with your own
    #Find the enddate text field
    enddate = driver.find_element(By.NAME, 'date_end')
    #Clear the text in it
    enddate.clear()
    #Paste in our own date
    enddate.send_keys(EndDate)

    #Manipulating the data checkbox (for hmF2)
    #Find the Data drop down list
    hmF2_menu = driver.find_element(By.NAME, "chosenchars[]")
    #Select it
    select = Select(hmF2_menu)
    #Select the data that we want (hmF2)
    select.select_by_visible_text(RelevantData)

    #Click the button to retrieve the data
    submit_button = driver.find_element(By.NAME, 'query_submit')
    # Click button
    submit_button.click()

    #Save webpage to a text file
    #Get the url of the new page we are on (it has the data we need)
    data_page_url = driver.current_url
    #print(data_page_url)

    #Gets the data on the webpage
    page_source = driver.page_source

    #Names the file the station's name
    file_path = save_directory + "\\" + StationLocation + ".txt"

    #Saves the data to the directory that this file is
    with open(file_path, "w", encoding="utf-8") as file:
        file.write(page_source)
