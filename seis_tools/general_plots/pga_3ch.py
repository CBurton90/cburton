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
    client = Client("http://service.iris.edu")
    st = client.get_waveforms(
        network="IU", station=station, location=location,
        channel=channel, starttime=t1, endtime=t1 + duration)
    tr = Stream()
    for n in range(len(st)):
        st.merge()[n]
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
tr = get_and_remove_response(station="RAO", channel="H**", location="*", output="ACC", t1=UTCDateTime(2021, 3, 4, 17, 41, 0))
# tr += get_and_remove_response(station="RIZ", channel="H**", location="10", output="ACC", t1=UTCDateTime(2021, 3, 4, 17, 41, 30))



fig, axs = plt.subplots(3, sharex=True, figsize=(15,15))
fig.suptitle('RIZ HN* PGA', fontsize=14)



for i in range(len(tr)):
	peaks, _ = find_peaks(tr[i].data, distance=tr[i].stats.npts)
	np.diff(peaks)
	axs[i].plot(tr[i].times(), tr[i].data)
	axs[i].plot(peaks/tr[i].stats.sampling_rate, tr[i][peaks], "o")
	axs[i].hlines(tr[i][peaks], 0, peaks/tr[i].stats.sampling_rate, linestyles='dashed', color='k')
	line = Line2D([0,1],[0,1],linestyle='--', color='k')
	axs[i].legend([line], tr[i][peaks], loc=3)
	props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
	axs[i].text(0.05, 0.85, tr[i].id, transform=axs[i].transAxes, fontsize=10,
        verticalalignment='top', bbox=props)
	axs[i].set_xlim((0,20))
	axs[2].set_xlabel('Seconds since 17:41:30 UTC')
	axs[i].set_ylabel('Acceleration (m/s^2)')
	print(peaks/tr[i].stats.sampling_rate, tr[i][peaks])
	



plt.show()
# plt.savefig('RaoulSM_Mw7.4.png', format='PNG', dpi=400)
	



# tr.plot(equal_scale=False)
