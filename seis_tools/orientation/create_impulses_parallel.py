#!/usr/bin/env python3

from __future__ import print_function
from obspy import UTCDateTime
from obspy import Stream
from obspy import Trace
from obspy import read, read_inventory
from obspy.clients.fdsn import Client
from obspy.signal.cross_correlation import correlate
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal

def get_and_remove_response(station, channel, location, output, t1, duration):
    client = Client("http://service-nrt.geonet.org.nz")
    st = client.get_waveforms(
        network="NZ", station=station, location=location,
        channel=channel, starttime=t1, endtime=t1 + duration)
    tr = Stream()
    st.merge(fill_value='interpolate')
    print(st)
    for n in range(len(st)):
        
        inv = client.get_stations(
            network=st[n].stats.network, station=st[n].stats.station, 
            location=st[n].stats.location, channel=st[n].stats.channel, 
            level="response", startbefore=t1, endafter=t1 + duration)
        # st[n].detrend('linear')
        # st[n].taper(max_percentage=0.05, type='cosine')
        st[n].filter('bandpass', freqmin=0.1, freqmax=0.2, zerophase=True)
        st[n].remove_response(output=output, pre_filt=False, plot=False,
                       water_level=60, inventory=inv)
        tr += st[n]

    return tr

def equal_and_split(st, id, stacks):
	st = st.select(id=id)
	st = np.delete(st, 0, 1)
	rec_st = np.hsplit(st[0].data, stacks)

	return rec_st

def get_windows(n, Mt, olap):
    # Split a signal of length n into olap% overlapping windows each containing Mt terms 
    ists = []
    ieds = []

    ist = 0
    while 1:
        ied = ist + Mt
        if ied > n:
            break
        ists.append(ist)
        ieds.append(ied)
        ist += int(Mt * (1 - olap/100))
    return ists, ieds

def check_length(stream):
    """
    Forces all traces to have same number of samples.

    Traces must be one day long.
    :type stream:`~obspy.core.stream.Stream` object. 
    :param stream: Stream containing one or more day-long trace 
    :return: Stream of similar length traces 

    """
    pts = 24*3600*stream[0].stats.sampling_rate
    npts = []
    for trace in stream:
        npts.append(trace.stats.npts)
    npts = np.array(npts)
    if len(npts) == 0:
        return stream	
    index = np.where(npts != pts)
    index = list(index[0])[::-1]

    # remove short traces
    for trace in index:
        stream.pop(trace)

    return stream


t1 = UTCDateTime(2021, 6, 24, 6, 30, 0)
# duration = 2*48.2*60*60
# time = 4*24*60*60

duration = 4.3*24*60*60
time = 4*24*60*60

stime = UTCDateTime("2021-06-24T07:00:00")

rst = get_and_remove_response(station='RPZ', channel='HH*', location='10', output='VEL', t1=t1, duration=duration)
rst = rst.trim(starttime=stime, endtime=stime + time)

sst = get_and_remove_response(station='WVZ', channel='HH*', location='10', output='VEL', t1=t1, duration=duration)
sst = sst.trim(starttime=stime, endtime=stime + time)

sst.sort()
rst.sort()

print(rst)
print(sst)

sou_st_e = equal_and_split(st=sst, id="NZ.WVZ.10.HHE", stacks=time/1800)
sou_st_n = equal_and_split(st=sst, id="NZ.WVZ.10.HHN", stacks=time/1800)
sou_st_z = equal_and_split(st=sst, id="NZ.WVZ.10.HHZ", stacks=time/1800)

rec_st_2 = equal_and_split(st=rst, id="NZ.RPZ.10.HH2", stacks=time/1800)
rec_st_1 = equal_and_split(st=rst, id="NZ.RPZ.10.HH1", stacks=time/1800)
rec_st_z = equal_and_split(st=rst, id="NZ.RPZ.10.HHZ", stacks=time/1800)



for n in range(len(sou_st_z)):

    signal.detrend(sou_st_z[n],type='linear')
    window = signal.tukey(len(sou_st_z[n]),alpha=0.05)
    sou_st_z[n] * window
    signal.detrend(sou_st_n[n],type='linear')
    window = signal.tukey(len(sou_st_n[n]),alpha=0.05)
    sou_st_n[n] * window
    signal.detrend(sou_st_e[n],type='linear')
    window = signal.tukey(len(sou_st_e[n]),alpha=0.05)
    sou_st_e[n] * window

print("%s break")
print(sou_st_z[0])

for x in range(len(rec_st_z)):
	
	signal.detrend(rec_st_z[x],type='linear')
	window = signal.tukey(len(rec_st_z[x]),alpha=0.05)
	rec_st_z[x] * window
	signal.detrend(rec_st_1[x],type='linear')
	window = signal.tukey(len(rec_st_1[x]),alpha=0.05)
	rec_st_1[x] * window
	signal.detrend(rec_st_2[x],type='linear')
	window = signal.tukey(len(rec_st_2[x]),alpha=0.05)
	rec_st_2[x] * window

print(rec_st_z[0])
print(range(len(rec_st_z)))

z_ccf_windows = []
n_ccf_windows = []
e_ccf_windows = []
for i in range(len(rec_st_z)):
    z_corr = correlate(rec_st_z[i],sou_st_z[i],150*100)
    n_corr = correlate(rec_st_1[i],sou_st_z[i],150*100)
    e_corr = correlate(rec_st_2[i],sou_st_z[i],150*100)
    print("%s zcorr break")
    print(z_corr)
    z_ccf_windows.append(z_corr)
    n_ccf_windows.append(n_corr)
    e_ccf_windows.append(e_corr)
    # for a in range(len(z_corr)):
    #     z_stacked = np.add(z_corr[a],z_corr[a+1]) 
    #     z_stacked / len(rec_st_z)
    #     # Czz.append(z_stacked)
print("%s z stacked break")
print(z_ccf_windows)



zz = np.array(z_ccf_windows, dtype=np.float32)
nz = np.array(n_ccf_windows, dtype=np.float32)
ez = np.array(e_ccf_windows, dtype=np.float32)
z_stacked = zz.sum(axis=0) / len(zz)
n_stacked = nz.sum(axis=0) / len(nz)
e_stacked = ez.sum(axis=0) / len(ez)


print(len(zz))
print(range(len(zz)))

plt.plot(z_stacked)
plt.show()

print(len(z_stacked))
print(range(len(z_stacked)))




# Convert to NumPy character array
data = np.array(z_stacked, dtype='|S1')

# Fill header attributes
stats = {'network': 'NZ', 'station': 'RPZ', 'location': '10',
         'channel': 'HHZ', 'npts': len(z_stacked), 'sampling_rate': 100,
         'mseed': {'dataquality': 'D'}}
stats2 = {'network': 'NZ', 'station': 'RPZ', 'location': '10',
         'channel': 'HH1', 'npts': len(z_stacked), 'sampling_rate': 100,
         'mseed': {'dataquality': 'D'}}
stats3 = {'network': 'NZ', 'station': 'RPZ', 'location': '10',
         'channel': 'HH2', 'npts': len(z_stacked), 'sampling_rate': 100,
         'mseed': {'dataquality': 'D'}}
# set current time
stats['starttime'] = UTCDateTime("2021-06-24T07:00:00")
stats2['starttime'] = UTCDateTime("2021-06-24T07:00:00")
stats3['starttime'] = UTCDateTime("2021-06-24T07:00:00")
st = Stream([Trace(data=z_stacked, header=stats)])
st2 = Stream([Trace(data=n_stacked, header=stats2)])
st3 = Stream([Trace(data=e_stacked, header=stats3)])
# write as ASCII file (encoding=0)
st.write("RPZ_WVZ_ZZ.mseed", format='MSEED', encoding="FLOAT32", reclen=512)
st2.write("RPZ_WVZ_1Z.mseed", format='MSEED', encoding="FLOAT32", reclen=512)
st3.write("RPZ_WVZ_2Z.mseed", format='MSEED', encoding="FLOAT32", reclen=512)







