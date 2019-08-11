#!/usr/bin/env python3

import pandas as pd 
from datetime import datetime
from datetime import timedelta

SunXtenderLife = timedelta(weeks=416)
Haze40Life = timedelta(weeks=416)

powerdf = pd.read_csv('~/git/sit/data/power.csv')

print('''\nPower asset options: SunXtender PVX1040T, Haze 40\n''')

Asset = input("Enter power asset: ")

# Add zero padding

#powerdf['INSTALLATION_DATE'] = pd.to_timedelta(powerdf['INSTALLATION_DATE'])

# Convert dates from strings to datetime 

powerdf['INSTALLATION_DATE'] = pd.to_datetime(powerdf['INSTALLATION_DATE'], format='%d/%m/%Y')

if Asset == "SunXtender PVX1040T":

	A = powerdf[powerdf['MODEL'].str.match('SunXtender', na=False)]
	A['LIFE'] = A['INSTALLATION_DATE'] + SunXtenderLife

	A['REPLACEMENT'] = A.LIFE < datetime.now()
	C = A.loc[A['REPLACEMENT'] == True, ['CODE','MANUFACTURER','MODEL','QUANTITY','INSTALLATION_DATE','LIFE','REPLACEMENT']]
	print(C)

if Asset == "Haze 40":

	Haze = powerdf[powerdf['MANUFACTURER'].str.match('Haze', na=False)]
	HazeModel = Haze[Haze['MODEL'].str.match('40', na=False)]
	HazeModel['LIFE'] = HazeModel['INSTALLATION_DATE'] + Haze40Life

	HazeModel['REPLACEMENT'] = HazeModel.LIFE < datetime.now()
	HazeEOL = HazeModel.loc[HazeModel['REPLACEMENT'] == True, ['CODE','MANUFACTURER','MODEL','QUANTITY','INSTALLATION_DATE','LIFE','REPLACEMENT']]
	print(HazeEOL)






# C = A[A['REPLACEMENT']]
# print(C)

# C = A.loc[A['REPLACEMENT'] == 'True', A[['CODE','TYPE']]
# print(C)p

# df_new = A.loc


# A['REPLACEMENT'] = np.where( A.LIFE <= datetime.now(), 'YES', 'NO') 
# print(A)

# EOL = A['LIFE'] < datetime.now()
# if EOL.all():
# 	print(A)


	





