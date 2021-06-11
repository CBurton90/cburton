#!/usr/bin/env python3

from obspy import UTCDateTime
from obspy import Stream
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


t1 = UTCDateTime(2021, 5, 13, 23, 59, 59)
duration = 1.1*24*60*60

st = get_and_remove_response(station='WEL', channel='HH*', location='10', output='VEL', t1=t1, duration=duration)
stime = UTCDateTime("2021-05-14T00:00:00")
st = st.trim(starttime=stime, endtime=stime + 24*60*60)

st = check_length(st)

st_local = read("/home/conradb/seis_tools/orientation/NIC_T120_BH-1/NIC_centaur-6_7978_20210514_000000.seed")
st_local = check_length(st_local)

for n in range(len(st_local)):
	# st_local[n].detrend('linear')
	# st_local[n].taper(max_percentage=0.05, type='cosine')
	st_local[n].filter('bandpass', freqmin=0.1, freqmax=0.2, zerophase=True)

tr_z = st_local.select(id="NZ.NIC.10.HHZ")
inv_z = read_inventory("/home/conradb/seis_tools/orientation/NIC_T120_BH-1/NZ.NIC.10.xml/NZ.NIC.10.HHZ.xml")
print(inv_z)
tr_z.remove_response(inventory=inv_z, pre_filt=False, output='VEL', water_level=60, plot=False)

tr_1 = st_local.select(id="NZ.NIC.10.HH1")
inv_1 = read_inventory("/home/conradb/seis_tools/orientation/NIC_T120_BH-1/NZ.NIC.10.xml/NZ.NIC.10.HH1.xml")
tr_1.remove_response(inventory=inv_1, pre_filt=False, output='VEL', water_level=60, plot=False)

tr_2 = st_local.select(id="NZ.NIC.10.HH2")
inv_2 = read_inventory("/home/conradb/seis_tools/orientation/NIC_T120_BH-1/NZ.NIC.10.xml/NZ.NIC.10.HH2.xml")
tr_2.remove_response(inventory=inv_2, pre_filt=False, output='VEL', water_level=60, plot=False)

rt = Stream()

rt += tr_2 + tr_1 + tr_z

st.sort()
# rt.sort()

# stime = UTCDateTime("2021-05-14T00:00:00")
st_trim = st
rt_trim = rt

print(st_trim)
print(rt_trim)


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

for n in range(len(rt_trim)):


	rt_split_2 = np.hsplit(rt_trim[0].data, 144)
	rt_split_1 = np.hsplit(rt_trim[1].data, 144)
	rt_split_z = np.hsplit(rt_trim[2].data, 144)
	print(rt_split_z)
	signal.detrend(rt_split_z,type='linear')
	# window = signal.windows.tukey(len(rt_split_z),alpha=0.05)
	# rt_split_z[n] * window
	# plt.plot(rt_split[2])
	# plt.plot(rt_split[1])
	# plt.plot(rt_split[2])
	print("%s break")
print(rt_split_z[0])

for x in range(len(st_trim)):

	st_split_e = np.hsplit(st_trim[0].data, 144)
	st_split_n = np.hsplit(st_trim[1].data, 144)
	st_split_z = np.hsplit(st_trim[2].data, 144)
	signal.detrend(st_split_z,type='linear')
	# window = signal.windows.tukey(len(st_split[x]),alpha=0.05)
	# st_split[x] * window
	# plt.plot(st_split[2])
	# plt.plot(rt_split[1])
	# plt.plot(rt_split[2])
	# print(st_split[2])

print(st_split_z[0])	

for i in range(len(st_split_z)):
	z_corr = correlate(st_split_z[i],rt_split_z[i],120*100)
	print(z_corr)
	z_stacked = z_corr
	z_stacked += z_corr[i]
print(z_stacked)

plt.plot(z_stacked)
plt.show()





