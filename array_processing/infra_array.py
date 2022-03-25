# functions defined by Stephen Arrowsmith (sarrowsmith@smu.edu)

import numpy as np
from obspy.signal.array_analysis import *
from obspy.core import AttribDict
from obspy import read
from obspy.clients.fdsn import Client
from obspy import UTCDateTime
from pyproj import Geod
from numpy.linalg import inv
import matplotlib.pyplot as plt
import warnings, utm
import utm
from sklearn.neighbors import NearestNeighbors
warnings.filterwarnings("ignore")

def get_array_coords(st, inv, ref_station):
    '''
    Returns the array coordinates for an array, in km with respect to the reference array provided
    
    Inputs:
    st - ObsPy Stream object containing array data
    ref_station - A String containing the name of the reference station
    
    Outputs:
    X - [Nx2] NumPy array of array coordinates in km
    stnm - [Nx1] list of element names
    
    Stephen Arrowsmith (sarrowsmith@smu.edu) modified by C.Burton
    '''
    
    X = np.zeros((len(st), 2))
    stnm = []
    for i in range(0, len(st)):
        E, N, _, _ = utm.from_latlon(inv[i].stations[0].latitude, inv[i].stations[0].longitude)
        X[i,0] = E; X[i,1] = N
        stnm.append(st[i].stats.station)

    # Adjusting to the reference station, and converting to km:
    ref_station_ix = np.where(np.array(stnm) == ref_station)[0][0]    # index of reference station
    X[:,0] = (X[:,0] - X[ref_station_ix,0])
    X[:,1] = (X[:,1] - X[ref_station_ix,1])
    X = X/1000.
    
    return X, stnm

def gc_backzimuth(st, evlo, evla):
    '''
    Computes the Great-Circle backazimuth and epicentral distance for an array in st given a known event location
    
    Inputs:
    st - ObsPy Stream object containing array data
    evlo - A float containing the event longitude
    evla - A float containing the event latitude
    
    Outputs:
    baz - Backazimuth (degrees from North)
    dist - Great-circle distance (km)
    '''
    
    g = Geod(ellps='sphere')
    
    a12,a21,dist = g.inv(evlo,evla,st[0].stations[0].longitude,st[0].stations[0].latitude); dist = dist/1000.
    
    return a21%360., dist

client = Client("http://service.geonet.org.nz")
start_t = UTCDateTime("2022-01-15T04:10:00.000")
end_t = start_t + 60*60*6

# start_t = UTCDateTime("2022-01-15T04:00:00.000")
# end_t = start_t + 60*60

#InfraBSU microphone LC31
# stations = [('NTVZ'),('PREZ'),('KHEZ'),('WIZ'),('COVZ'),('WSRZ')]
# stations = [('NTVZ'),('PREZ'),('WIZ')]

# # Setra270 LC30
stations = ['COVZ', 'IVVZ', 'FWVZ', 'WHVZ', 'TRVZ', 'TOVZ', 'TMVZ', 'WTVZ', 'OTVZ', 'NGZ', 'ETVZ', 'MAVZ', 'WNVZ', 'NTVZ', 'SNVZ', 'NOVZ', 'KRVZ']

st = Stream()

for station in stations:
	try:
		tr = client.get_waveforms("NZ", station, '30', "HDF", start_t, end_t)
		tr.merge(fill_value='interpolate')
		tr.taper(type='cosine', max_percentage=0.05, max_length=60)
		tr.filter('bandpass', freqmin=0.1, freqmax=5)
		print(tr)
		st += tr
	except Exception as e:
		print(station)
		print(e)
		# try:
		# 	tr = client.get_waveforms("NZ", station, '32', "HDF", start_t, end_t)
		# 	# tr.taper(type='cosine', max_percentage=0.05, max_length=60)
		# 	# tr.filter('bandpass', freqmin=0.1, freqmax=5)
		# 	print(tr)
		# 	st += tr
		# except Exception as e:
		# 	print(e)
		# 	try:
		# 		tr = client.get_waveforms("NZ", station, '33', "HDF", start_t, end_t)
		# 		# tr.taper(type='cosine', max_percentage=0.05, max_length=60)
		# 		# tr.filter('bandpass', freqmin=0.1, freqmax=5)
		# 		print(tr)
		# 		st += tr
		# 	except Exception as e:
		# 		print(e)

# st2 = Stream()

# for station in stations:
# 	try:
# 		tr = client.get_waveforms("NZ", station, '1*', "HHZ", start_t, end_t)
# 		tr.merge(fill_value='interpolate')
# 		tr.taper(type='cosine', max_percentage=0.05, max_length=60)
# 		# tr.filter('bandpass', freqmin=1/2, freqmax=15)
# 		print(tr)
# 		st2 += tr
# 	except Exception as e:
# 		print(station)
# 		print(e)

# st += st2
st.sort()
print(st)


inv = []
for n in range(len(st)):
	inv += client.get_stations(
		network=st[n].stats.network, station=st[n].stats.station, 
        location=st[n].stats.location, channel=st[n].stats.channel, 
        level="channel", startbefore=start_t, endafter=end_t)

# print(inv[8].stations[0].latitude)

for i in range(0, len(st)):
	print(range(0, len(st)))
	print(range(len(st)))
	E, N, _, _ = utm.from_latlon(inv[i].stations[0].latitude, inv[i].stations[0].longitude)
	print(E, N)


# client = Client('IRIS')
# tr = client.get_waveforms("AU","I07H*","*","BDF", start_t, end_t)
#tr.taper(type='cosine', max_percentage=0.05, max_length=60)
#tr.filter('bandpass', freqmin=0.1, freqmax=5)

st.plot(equal_scale=False)

X, stn = get_array_coords(st, inv, 'COVZ')

plt.plot(X[:,0], X[:,1], 'o')
for i in range(0, len(stn)):
    plt.text(X[i,0], X[i,1], stn[i])
plt.xlabel('Distance (km)')
plt.ylabel('Distance (km)')
plt.gca().set_aspect('equal')
plt.show()

# Defining a start time to compute an FK plot for:
time_start = st[0].stats.starttime+6000
time_end = st[0].stats.starttime+12000

baz, dist = gc_backzimuth(inv, -175.38, -20.57)
print('GC backazimuth =', baz)



