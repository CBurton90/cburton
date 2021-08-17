#!/usr/bin/env python3

from obspy import UTCDateTime
from obspy import Stream
from obspy.clients.fdsn import Client
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from matplotlib.lines import Line2D
from scipy import signal
from scipy.signal import find_peaks
import sys

def get_and_remove_response(station, channel, location, output, t1, duration=100):
    # client = Client("http://service.iris.edu")
    client = Client("http://service.geonet.org.nz")
    st = client.get_waveforms(
        network="NZ", station=station, location=location,
        channel=channel, starttime=t1, endtime=t1 + duration)
    print(st)
    st.merge(fill_value='interpolate')
    print(st)
    tr = Stream()
    for n in range(len(st)):
    	
    	# print(st)
    	inv = client.get_stations(
    		network=st[n].stats.network, station=st[n].stats.station, 
    		location=st[n].stats.location, channel=st[n].stats.channel, 
    		level="response", startbefore=t1, endafter=t1 + duration)
    		# pre_filt = (0.005, 0.006, 30.0, 35.0)
    	st[n].remove_response(output=output, pre_filt=False, plot=False, 
    		water_level=60, inventory=inv)
    	tr += st[n]

    return tr

# pre_filt = [0.01,1,25,30]
tr = get_and_remove_response(station="GLKZ", channel="HHZ", location="10", output="VEL", t1=UTCDateTime(2021, 3, 4, 19, 28, 45))
tr += get_and_remove_response(station="RIZ", channel="HHZ", location="10", output="VEL", t1=UTCDateTime(2021, 3, 4, 19, 28, 45))
tr += get_and_remove_response(station="RIZ", channel="HNZ", location="20", output="ACC", t1=UTCDateTime(2021, 3, 4, 19, 28, 45))



fig, axs = plt.subplots(3, sharex=True, figsize=(20,15))
fig.suptitle('Mw8.1 GLKZ and RIZ vertical component comparison (PGV/PGA)', fontsize=16)



for i in range(len(tr)):
	peakspos, _ = find_peaks(tr[i].data, distance=tr[i].stats.npts)
	peaksneg, _ = find_peaks(-tr[i].data, distance=tr[i].stats.npts)
	np.diff(peakspos)
	np.diff(peaksneg)
	axs[i].plot(tr[i].times(), tr[i].data)
	axs[i].plot(peakspos/tr[i].stats.sampling_rate, tr[i][peakspos], "o")
	axs[i].plot(peaksneg/tr[i].stats.sampling_rate, tr[i][peaksneg], "o")
	axs[i].hlines(tr[i][peakspos], 0, peakspos/tr[i].stats.sampling_rate, linestyles='dashed', color='k')
	axs[i].hlines(tr[i][peaksneg], 0, peaksneg/tr[i].stats.sampling_rate, linestyles='dashed', color='r')
	line = Line2D([0,1],[0,1],linestyle='--', color='k')
	line2 = Line2D([0,1],[0,1],linestyle='--', color='r')
	axs[i].legend([line,line2], [tr[i][peakspos],tr[i][peaksneg]], loc=3, bbox_to_anchor=(0.01, 0.06), fontsize=10)
	props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
	axs[i].text(0.01, 0.85, tr[i].id, transform=axs[i].transAxes, fontsize=10,
        verticalalignment='top', bbox=props)
	axs[i].set_xlim((0,100))
	axs[2].set_xlabel('Seconds since 19:28:45 UTC')
	axs[0].set_ylabel('Velocity (m/s)', fontsize=10)
	axs[1].set_ylabel('Velocity (m/s)', fontsize=10)
	axs[2].set_ylabel('Acceleration (m/s^2)', fontsize=10)
	print(peakspos/tr[i].stats.sampling_rate, tr[i][peakspos])
	


fig.align_labels()
fig.tight_layout()
# plt.show()
plt.savefig('Raoul_Mw8.1_Z_comparison.png', format='PNG', dpi=400)
	



# tr.plot(equal_scale=False)
