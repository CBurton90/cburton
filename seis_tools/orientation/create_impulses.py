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
    client = Client("http://service.geonet.org.nz")
    st = client.get_waveforms(
        network="NZ", station=station, location=location,
        channel=channel, starttime=t1, endtime=t1 + duration)
    tr = Stream()
    st.merge()
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


t1 = UTCDateTime(2021, 5, 31, 23, 58, 59)
duration = 2*48.2*60*60

st = get_and_remove_response(station='URZ', channel='HH*', location='10', output='VEL', t1=t1, duration=duration)
stime = UTCDateTime("2021-06-01T00:00:00")
st = st.trim(starttime=stime, endtime=stime + 2*48*60*60)


rt = get_and_remove_response(station='PUZ', channel='HH*', location='10', output='VEL', t1=t1, duration=duration)
rt = rt.trim(starttime=stime, endtime=stime + 2*48*60*60)



# st_local = read("/home/conradb/seis_tools/orientation/NIC_T120_BH-1/NIC_centaur-6_7978_20210514_000000.seed")
# st_local = check_length(st_local)

# for n in range(len(st_local)):
# 	# st_local[n].detrend('linear')
# 	# st_local[n].taper(max_percentage=0.05, type='cosine')
# 	st_local[n].filter('bandpass', freqmin=0.1, freqmax=0.2, zerophase=True)

# tr_z = st_local.select(id="NZ.NIC.10.HHZ")
# inv_z = read_inventory("/home/conradb/seis_tools/orientation/NIC_T120_BH-1/NZ.NIC.10.xml/NZ.NIC.10.HHZ.xml")
# print(inv_z)
# tr_z.remove_response(inventory=inv_z, pre_filt=False, output='VEL', water_level=60, plot=False)

# tr_1 = st_local.select(id="NZ.NIC.10.HH1")
# inv_1 = read_inventory("/home/conradb/seis_tools/orientation/NIC_T120_BH-1/NZ.NIC.10.xml/NZ.NIC.10.HH1.xml")
# tr_1.remove_response(inventory=inv_1, pre_filt=False, output='VEL', water_level=60, plot=False)

# tr_2 = st_local.select(id="NZ.NIC.10.HH2")
# inv_2 = read_inventory("/home/conradb/seis_tools/orientation/NIC_T120_BH-1/NZ.NIC.10.xml/NZ.NIC.10.HH2.xml")
# tr_2.remove_response(inventory=inv_2, pre_filt=False, output='VEL', water_level=60, plot=False)

# rt = Stream()

# rt += tr_2 + tr_1 + tr_z



st.sort()
rt.sort()

# stime = UTCDateTime("2021-05-14T00:00:00")

 

st_trim = st
rt_trim = rt

print(st_trim)
print(rt_trim)

# for i in range(len(st)):
#     st_trim += np.delete(st[i], 0)
# print(len(st_trim))



# stime2 = UTCDateTime("2021-05-11T14:52:30")
# st_trim2 = st.trim(stime, stime + 450)

# stime3 = UTCDateTime("2021-05-11T15:00:00")
# st_trim3 = st.trim(stime, stime + 450)

# stime4 = UTCDateTime("2021-05-11T15:07:30")
# st_trim4 = st.trim(stime, stime + 450)


# print(st_trim)




# z_corr = correlate(st_trim[3].data,st_trim[2].data,240*100)
# print(len(z_corr))

# z_corr2 = correlate(st_trim2[3].data,st_trim2[2].data,240*100)
# z_corr3 = correlate(st_trim3[3].data,st_trim3[2].data,240*100)
# z_corr4 = correlate(st_trim4[3].data,st_trim4[2].data,240*100)

# z_corr_stacked = (z_corr + z_corr2 + z_corr3 + z_corr4) / 4

# plt.plot(z_corr_stacked)
# plt.show()


# n = 100
# x = np.arange(n)







# for n in range(len(st)):
# 	for x in range(len(rt)):
# 		ists, ieds = get_windows(len(st[0]), Mt=900*100, olap=0) # windows of length 20 and 50% overlap
# 		for ist, ied in zip(ists, ieds):
# 			plt.plot(st[n][ist:ied])
# 			plt.show()


# a = []

# a = np.hsplit(rt_trim[n].data, 2)

# print(a)

 # a = []

st_1 = st_trim.select(id="NZ.URZ.10.HH1")
st_1 = np.delete(st_1, 0, 1)
st_2 = st_trim.select(id="NZ.URZ.10.HH2")
st_2 = np.delete(st_2, 0, 1)
st_z = st_trim.select(id="NZ.URZ.10.HHZ")
st_z = np.delete(st_z, 0, 1)

rt_n = rt_trim.select(id="NZ.PUZ.10.HHN")
rt_n = np.delete(rt_n, 0, 1)
rt_e = rt_trim.select(id="NZ.PUZ.10.HHE")
rt_e = np.delete(rt_e, 0, 1)
rt_z = rt_trim.select(id="NZ.PUZ.10.HHZ")
rt_z = np.delete(rt_z, 0, 1)

# rt_trim = np.delete(rt_trim, 0, 1)
# st_trim = np.delete(st_trim, 0, 1)

rt_split_e = np.hsplit(rt_e[0].data, 2*96)
rt_split_n = np.hsplit(rt_n[0].data, 2*96)
rt_split_z = np.hsplit(rt_z[0].data, 2*96)

st_split_2 = np.hsplit(st_2[0].data, 2*96)
st_split_1 = np.hsplit(st_1[0].data, 2*96)
st_split_z = np.hsplit(st_z[0].data, 2*96)

print(range(len(rt_split_z)))

for n in range(len(rt_split_z)):

    signal.detrend(rt_split_z[n],type='linear')
    window = signal.tukey(len(rt_split_z[n]),alpha=0.05)
    rt_split_z[n] * window
    signal.detrend(rt_split_n[n],type='linear')
    window = signal.tukey(len(rt_split_n[n]),alpha=0.05)
    rt_split_n[n] * window
    signal.detrend(rt_split_e[n],type='linear')
    window = signal.tukey(len(rt_split_e[n]),alpha=0.05)
    rt_split_e[n] * window

print("%s break")
print(rt_split_z[0])

for x in range(len(st_split_z)):

    signal.detrend(st_split_z[x],type='linear')
    window = signal.tukey(len(st_split_z[x]),alpha=0.05)
    st_split_z[x] * window
    signal.detrend(st_split_1[x],type='linear')
    window = signal.tukey(len(st_split_1[x]),alpha=0.05)
    st_split_1[x] * window
    signal.detrend(st_split_2[x],type='linear')
    window = signal.tukey(len(st_split_2[x]),alpha=0.05)
    st_split_2[x] * window

print(st_split_z[0])
print(range(len(st_split_z)))

z_ccf_windows = []
n_ccf_windows = []
e_ccf_windows = []
for i in range(len(st_split_z)):
    z_corr = correlate(st_split_z[i],rt_split_z[i],120*100)
    n_corr = correlate(st_split_1[i],rt_split_z[i],120*100)
    e_corr = correlate(st_split_2[i],rt_split_z[i],120*100)
    print("%s zcorr break")
    print(z_corr)
    z_ccf_windows.append(z_corr)
    n_ccf_windows.append(n_corr)
    e_ccf_windows.append(e_corr)
    # for a in range(len(z_corr)):
    #     z_stacked = np.add(z_corr[a],z_corr[a+1]) 
    #     z_stacked / len(st_split_z)
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
# data = np.array(z_stacked, dtype='|S1')

# Fill header attributes
stats = {'network': 'NZ', 'station': 'URZ', 'location': '10',
         'channel': 'HHZ', 'npts': len(z_stacked), 'sampling_rate': 100,
         'mseed': {'dataquality': 'D'}}
stats2 = {'network': 'NZ', 'station': 'URZ', 'location': '10',
         'channel': 'HH1', 'npts': len(z_stacked), 'sampling_rate': 100,
         'mseed': {'dataquality': 'D'}}
stats3 = {'network': 'NZ', 'station': 'URZ', 'location': '10',
         'channel': 'HH2', 'npts': len(z_stacked), 'sampling_rate': 100,
         'mseed': {'dataquality': 'D'}}
# set current time
stats['starttime'] = UTCDateTime("2021-06-01T00:00:00")
stats2['starttime'] = UTCDateTime("2021-06-01T00:00:00")
stats3['starttime'] = UTCDateTime("2021-06-01T00:00:00")
st = Stream([Trace(data=z_stacked, header=stats)])
st2 = Stream([Trace(data=n_stacked, header=stats2)])
st3 = Stream([Trace(data=e_stacked, header=stats3)])
# write as ASCII file (encoding=0)
st.write("URZ_PUZ_ZZ.mseed", format='MSEED', encoding="FLOAT32", reclen=512)
st2.write("URZ_PUZ_1Z.mseed", format='MSEED', encoding="FLOAT32", reclen=512)
st3.write("URZ_PUZ_2Z.mseed", format='MSEED', encoding="FLOAT32", reclen=512)







