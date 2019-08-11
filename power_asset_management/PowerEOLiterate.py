#!/usr/bin/env python3

# Import packages

import pandas as pd 
from datetime import datetime
from datetime import timedelta

# Convert csv's to datafranes, git path may need to be changed depending on your local setup

powerdf = pd.read_csv('~/git/sit/data/power.csv')
blessed = pd.read_csv('~/git/sit/data/blessed_tech.csv')

# Convert model names from dataframe to a list

dfValues = blessed['MODEL'].values

# Convert installation dates from strings to datetime 

powerdf['INSTALLATION_DATE'] = pd.to_datetime(powerdf['INSTALLATION_DATE'], format='%d/%m/%Y')

# Convert the lifetime of power equipment models from years to days as a timedelta object

blessed['LIFE'] = pd.to_timedelta(blessed['LIFE'], unit='Y')

# Extract relevant columns from dataframe pre-merge 

cutblessed = blessed[['MODEL','LIFE']]

# Create an empty list

power_list_appended = []

# Iterate through first dataframe with model list and create new dataframe where models match 

for i in dfValues:
	powerdf_matched = powerdf.loc[powerdf['MODEL'] == i]

    #  Append eached matched model dataframe created into a list

	power_list_appended.append(powerdf_matched)

    #  Concatenate list's into a single dataframe containing all matched models

	powerdf_matched_combined = pd.concat(power_list_appended)

	#  Merge dataframe on model name to add new column with model lifetimes

	
	powerdf_lifetime = pd.merge(cutblessed, powerdf_matched_combined,  how='outer', left_on=['MODEL'], right_on = ['MODEL'])

    #  Add the lifetime of models to their installation dates to produce a new column that gives the end of life for the specific model

	powerdf_lifetime['EOL'] = powerdf_lifetime['INSTALLATION_DATE'] + powerdf_lifetime['LIFE']

	#  Produce a boolean value in a new column to show if the end of life of the model is less than the current datetime

	powerdf_lifetime['REPLACEMENT'] = powerdf_lifetime.EOL < datetime.now()

	#  Output all true boolean values to a new dataframe with relevant columns

	powerdf_final = powerdf_lifetime.loc[powerdf_lifetime['REPLACEMENT'] == True, ['CODE','MANUFACTURER','MODEL','QUANTITY','INSTALLATION_DATE','EOL','REPLACEMENT']]
    
    #  Sort the dataframe alphabetically by site code

	powerdf_final_sorted = powerdf_final.sort_values(by=['CODE'], ascending=True)

    #  Print output of final datframe showing all equipment models last are end of lifetime
	
	with pd.option_context('display.max_rows', None, 'display.max_columns', None):
		print(powerdf_final_sorted)

	#  Convert final dataframe to a csv	

	powerdf_final_sorted.to_csv("EOL_power_assets.csv", index=False)

	
		


		


