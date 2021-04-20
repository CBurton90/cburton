import requests
import json
import pandas as pd
import numpy as np
from pandas.io.json import json_normalize
import matplotlib.pyplot as plt
import re

L = []

# retrieve all lines containing Trimble Alloy from a csv 
with open("/home/conradb/git/network/data/devices.csv",'r') as f:
    for line in f.readlines():
        # python can do regexes, but this is for s fixed string only
        if "Alloy" in line:
            idx1 = line.find('"')
            idx2 = line.find('"', idx1+1)
            field = line[idx1+1:idx2-1]
            L.append(field)

# create a list of only the first entries (site device locality name) delimited by comma
output = []
for string in L:
	a = string.split(',')[0]
	# print(a)
	output.append(a)
      

api_call = []
r = []
req = []

# build these into a list of api url calls
for i in output:
	url = "https://dapper-api.geonet.org.nz/data/fdmp?key="+i+"&fields=voltage&starttime=2021-01-01T00:00:00Z&endtime=2021-04-01T00:00:00Z"
	api_call.append(url)
	
print(api_call)

# loop through and request each api url and convert to Json, store all in a list
for n in api_call:
	try:
		print(n)
		print(requests.get(n))
		r = requests.get(n).json()
		req.append(r)
	except:
		pass

print(req[0])

# loop through list of requests and convert Json to a dataframe and clean 
for x in range(len(req)):
	df = pd.json_normalize(req[x]['results'])
	print(df)
	df = pd.json_normalize(data=req[x]['results'], record_path='records', meta=['domain', 'key', 'field'])
	df['timestamp'] = pd.to_datetime(df['timestamp'], unit='s')
	df["value"] = pd.to_numeric(df["value"])
	print(df.head())
	print(df.tail())

	# create a new dataframe for each that groups the percentage of values above and below 13500mV
	df_threshold = (
		df
		# fill nan values
		.pipe(lambda x: x.assign(
			time_delta=x["timestamp"].diff().fillna(pd.Timedelta(seconds=0))))
        # convert time to numeric (i.e., hours)
		.pipe(lambda x: x.assign(time_delta=x["time_delta"].dt.total_seconds() / (60 * 60)))
        # flag the values that are above the threshold
		.pipe(lambda x: x.assign(is_above=x["value"] > 13500))
         # drop unneeded columns
		.drop(["value", "domain", "key", "field"], axis="columns")
        # pivot around threshold flag and sum time
		.groupby("is_above").sum()
        # create percentage column
		.pipe(lambda x: x.assign(percent=x["time_delta"] / x["time_delta"].sum() * 100)
			)
		)
	df_threshold.reset_index()
	print(df_threshold)

	# choose only false values that are under 13500mV
	df_false = df_threshold.loc[False]
	print(df_false)


	# plot dataframes where more than 95% of values are under 13500mV
	if df_false.loc['percent'] > 95:
			ax = df.set_index("timestamp")["value"].plot.line(figsize=(10, 4))
			ax.axhline(13500, c="C1", ls="--")
			props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
			ax.text(0.01, 1.1, 'Percent below 13.5V (%s)' %(round(df_threshold.iloc[0,1], 3)), transform=ax.transAxes, fontsize=10, verticalalignment='top', bbox=props)
			plt.suptitle(df.iloc[0,3])
			# plt.show()
			plt.savefig(df.iloc[0,3]+'_voltage.png', format='PNG', dpi=400)
			plt.close()
    