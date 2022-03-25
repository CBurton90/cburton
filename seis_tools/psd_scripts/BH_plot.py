#!/usr/bin/env python3

import datetime
import os
from glob import glob

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.patheffects as pe
import numpy as np
import pandas as pd
from itertools import cycle
from cycler import cycler

from obspy import UTCDateTime, read
from obspy.clients.fdsn import Client
from obspy.signal import PPSD
from obspy.signal.spectral_estimation import get_nlnm, get_nhnm


start = UTCDateTime("2020-08-01")
end = UTCDateTime("2020-08-31")

stations = ['URZ']

path = "/home/conradb/git/cburton/seis_tools/psd_scripts/URZ_ppsd/"

datelist = pd.date_range(start.datetime, end.datetime, freq="D")

print(datelist)


ppsds1 = {}

# loop through stations (array reversed for colour visibility on plots i.e. darker blues/purples
# are assigned to WAZ/WEL at the top away from the MLNM) and load each daily ppsd and combine to
# produce yearly (2020) station ppsd's in a dictionary

for station in stations[::-1]:
	for day in datelist:
	    datestr = day.strftime("%Y-%m-%d")
	    print(datestr)
	    fn_pattern = "{}_*.npz".format(datestr)
	    print(fn_pattern)
	    for fn in glob(path+fn_pattern):
	    	print(fn)
	    	mseedid = fn.replace(".npz", "").split("_")[-1]
	    	print(mseedid)
	    	if mseedid not in ppsds1:
	    		ppsds1[mseedid] = PPSD.load_npz(fn, allow_pickle=True)
	    	else:
	            ppsds1[mseedid].add_npz(fn, allow_pickle=True)


start2 = UTCDateTime("2021-08-01")
end2 = UTCDateTime("2021-08-31")

datelist2 = pd.date_range(start2.datetime, end2.datetime, freq="D")

ppsds2 = {}

for station in stations[::-1]:
	for day in datelist2:
	    datestr = day.strftime("%Y-%m-%d")
	    print(datestr)
	    fn_pattern = "{}_*.npz".format(datestr)
	    print(fn_pattern)
	    for fn in glob(path+fn_pattern):
	    	print(fn)
	    	mseedid = fn.replace(".npz", "").split("_")[-1]
	    	print(mseedid)
	    	if mseedid not in ppsds2:
	    		ppsds2[mseedid] = PPSD.load_npz(fn, allow_pickle=True)
	    	else:
	            ppsds2[mseedid].add_npz(fn, allow_pickle=True)

# plotting parameters

num = len(stations)
styles = ['solid', 'dashed', 'dashdot', 'dotted']
num_styles = len(styles)

fig = plt.figure(figsize=(18,9))
ax = fig.add_subplot(111)

# trying a few different colormaps

# cm = plt.get_cmap('gist_rainbow')
# colors = [cm(1.*i/num) for i in range(num)]
# ax.set_prop_cycle(cycler('color', colors))

ax.set_prop_cycle('color',plt.cm.jet(np.linspace(0,1,num)))

# cycle through dictionary to grab individual station percentiles and modes 
# (modes also plotted) then combine in a larger array/list

for i, (mseedid, ppsd) in enumerate(ppsds1.items()):
	# ppsd.plot()
	p5 = ppsd.get_percentile(percentile=5)
	p95 = ppsd.get_percentile(percentile=95)
	mode = ppsd.get_mode()
	ax.plot(mode[0],mode[1], c='r', label=mseedid+' Aug 2020 (Guralp 3TB) Mode')


for i, (mseedid, ppsd) in enumerate(ppsds2.items()):
	# ppsd.plot()
	p5 = ppsd.get_percentile(percentile=5)
	p95 = ppsd.get_percentile(percentile=95)
	mode = ppsd.get_mode()
	lines = ax.plot(mode[0],mode[1], label=mseedid+' Aug 2021 (NMX T120-BH1) Mode')
	# lines[0].set_color(cm(i//num_styles*float(num_styles)/num))
	# lines[0].set_color(cm(i//num))
	lines[0].set_linestyle(styles[i%num_styles])

# get Petersons NHNM and NLNM for comparison

per, nlnm = get_nlnm()
per, nhnm = get_nhnm()

# print(lnm)
# plt.plot(p5[0],lnm, label='lowest combined NI NZNSN 5th percentiles')
# plt.plot(p5[0],hnm, label='highest combined NI NZNSN 95th percentiles')
# ax.plot(p5[0],mlnm, c='k', linestyle='dashed', label='NI MLNM')
ax.plot(per, nlnm, color='darkgrey', linewidth=2, linestyle='dashed')
ax.plot(per, nhnm, color='darkgrey', linewidth=2, linestyle='dashed', label='NLNM/NHNM')
plt.title("URZ Sensor Comparison", fontsize=15)
plt.xscale('log')
plt.xlim(0.02,200)
plt.xlabel('Period (s)')
plt.ylim(-200,-60)
plt.ylabel('dB[m^2/s^4/Hz]')
# plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.legend()
# plt.show()
plt.savefig('/home/conradb/git/cburton/seis_tools/psd_scripts/URZ_ppsd/URZ_comp.png', dpi=400, format='png')