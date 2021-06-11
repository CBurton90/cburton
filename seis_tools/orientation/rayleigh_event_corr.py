#!/usr/bin/env python3

from obspy import UTCDateTime
from obspy import Stream
from obspy.clients.fdsn import Client
from obspy import read, read_inventory
from obspy.signal.cross_correlation import correlate
from obspy.geodetics.base import gps2dist_azimuth
import numpy as np
from scipy.signal import hilbert
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap


def get_and_remove_response(station, channel, location, output, t1, duration):
    client = Client("http://service.geonet.org.nz")
    st = client.get_waveforms(
        network="NZ", station=station, location=location,
        channel=channel, starttime=t1, endtime=t1 + duration)
    tr = Stream()
    # st.merge()
    print(st)
    for n in range(len(st)):
        
        inv = client.get_stations(
            network=st[n].stats.network, station=st[n].stats.station, 
            location=st[n].stats.location, channel=st[n].stats.channel, 
            level="response", startbefore=t1, endafter=t1 + duration)
        pre_filt = (0.01, 0.02, 0.04, 0.08)
        st[n].detrend('linear')
        st[n].taper(max_percentage=0.1, type='cosine')
        st[n].remove_response(output=output, pre_filt=pre_filt, plot=False,
                       water_level=60, inventory=inv)
        tr += st[n]
    
    return tr

t1 = UTCDateTime("2021-05-31T04:13:00")
duration = 60*6

st = get_and_remove_response(station="URZ", channel="HH*", location="10", output="VEL", t1=t1, duration=duration)
st.trim(starttime=t1+10, endtime=t1+duration-10)

# print(st)
# st.plot()

URZ_lat = -38.25925
URZ_long = 177.11089
event_lat = -45.30
event_long = 166.81


back_azi = gps2dist_azimuth(URZ_lat, URZ_long, event_lat, event_long)
event_back_azi = back_azi[1]




# rotating CLOCKWISE, correlating with ZZ90:
maxSrz= []
#coherences=[]
thetas = np.linspace(0,2*np.pi,360)
for i_,theta in enumerate(thetas):
    Rad_rot =  np.cos(theta)*st.select(id="NZ.URZ.10.HH1") - np.sin(theta)*st.select(id="NZ.URZ.10.HH2")
    Trv_rot = np.sin(theta)*st.select(id="NZ.URZ.10.HH1") + np.cos(theta)*st.select(id="NZ.URZ.10.HH2")
    Hil_rad = np.imag(hilbert(Rad_rot))
    Z = st.select(id="NZ.URZ.10.HHZ")
    Hil_Z = np.imag(hilbert(Z))
    #numpy
    Ncorr = np.correlate(Rad_rot[0].data, Hil_Z[0].data)
    Szz = np.correlate(Hil_Z[0], Hil_Z[0])
    #obspy
    # Ncorr = correlate(Rad_rot[0].data, Hil_Z[0].data, 0)
    # Szz = correlate(Hil_Z[0], Hil_Z[0], 0)
    #trying reverse corr for order and hil-tf
    # Ncorr = correlate(Z[0].data, Hil_rad[0].data, 0)
    # Ncorr = correlate(Hil_Z[0].data, Rad_rot[0].data, 0)
    maxSrz.append(max(Ncorr)/max(Szz))
    
    
# print(len(maxSrz))
maxmaxSrz= max(maxSrz)
# print(maxmaxSrz)
rotangle_index = maxSrz.index(maxmaxSrz)
# print(rotangle_index)
rotangle = thetas[rotangle_index]
radians = rotangle
degrees = int(rotangle*180/np.pi)
corr_coeff = maxmaxSrz
print(radians)
print(degrees)
print(corr_coeff)

correction = ((event_back_azi - degrees) + 360) % 360
print(correction)



# # PLOTTING



Rotated =  np.cos(rotangle)*st.select(id="NZ.URZ.10.HH1") - np.sin(rotangle)*st.select(id="NZ.URZ.10.HH2")
Final = np.imag(hilbert(Rotated))

tr_1 = st.select(id="NZ.URZ.10.HH1")
t = tr_1[0].times()




fig = plt.figure(figsize=(25,20))
fig.suptitle('T120-BH1 orientation at URZ with event 2021p405872 (GeoNet)', fontsize=18)


# Plot the basemap of event and station

plt.subplots_adjust(wspace= 0.3, hspace= 0.2)

ax1 = fig.add_subplot(2, 2, 1)
m = Basemap(llcrnrlon=165.,llcrnrlat=-50.,urcrnrlon=180.,urcrnrlat=-30.,\
            rsphere=(6378137.00,6356752.3142),\
            resolution='l',projection='merc',\
            lat_0=-40.,lon_0=170.,lat_ts=None)
# nylat, nylon are lat/lon of New York
nylat = URZ_lat; nylon = URZ_long
# lonlat, lonlon are lat/lon of London.
lonlat = event_lat; lonlon = event_long
# draw great circle route between NY and London
m.drawgreatcircle(nylon,nylat,lonlon,lonlat,linewidth=2,color='b')
m.drawcoastlines()
m.fillcontinents()
# draw parallels
m.drawparallels(np.arange(-60,-20,5),labels=[1,1,0,1])
# draw meridians
m.drawmeridians(np.arange(-180,180,5),labels=[1,1,0,1])
ax1.set_title('5.5 MLv to URZ (2021-05-31T04:09:01)', pad=10)


# Plot the polar plot showing event back azi and sensor orientation

ax2 = fig.add_subplot(2, 2, 3, projection='polar')

r = [1, 1,]
theta = [np.radians(correction), np.radians(event_back_azi)]
labels = ['HH1','Event']


# fig = plt.figure(figsize=(5, 5))
# ax2 = fig.subplot(111, projection = 'polar')



for i in range(len(r)):
    ax2.plot([0, theta[i]], [0, r[i]], linewidth=2, label=labels[i])
    


ax2.set_yticklabels([])
ax2.set_theta_zero_location('N')

# Go clockwise
ax2.set_theta_direction(-1)
# # Start from the top
# ax2.set_theta_offset(np.pi/2)

# Connect two points with a curve

for curve in [[[0, correction], [0.4, 0.4]]]:
    curve[0] = np.deg2rad(curve[0])
    x = np.linspace( curve[0][0], curve[0][1], 500)
    y = interp1d( curve[0], curve[1])( x)
    ax2.plot(x, y, linewidth=2, label='Orientation'+ ' ' + str(np.round(correction,0)))

for curve in [[[correction, event_back_azi], [0.6, 0.6]]]:
    curve[0] = np.deg2rad(curve[0])
    x = np.linspace( curve[0][0], curve[0][1], 500)
    y = interp1d( curve[0], curve[1])( x)
    ax2.plot(x, y, linewidth=2, label='BAZ measured'+ ' = ' + str(np.round(degrees,0)))

for curve in [[[0, event_back_azi], [0.8, 0.8]]]:
    curve[0] = np.deg2rad(curve[0])
    x = np.linspace( curve[0][0], curve[0][1], 500)
    y = interp1d( curve[0], curve[1])( x)
    ax2.plot(x, y, linewidth=2, label='BAZ expected'+ ' = ' + str(np.round(event_back_azi,0)))

angle = np.deg2rad(337.5)
ax2.legend(loc='upper left', bbox_to_anchor=(.5 + np.cos(angle)/2, .5 + np.sin(angle)/2))


# Plot the HH1 waveform, rotated and Hilbert transformed HH1, HHZ, and CCF


props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax3 = fig.add_subplot(4,2,2)
ax3.plot(tr_1[0].times(), tr_1[0].data)
ax3.text(0.01, 0.95, tr_1[0].id, transform=ax3.transAxes, fontsize=8,
        verticalalignment='top', bbox=props)
ax3.set_ylabel('m/s')
ax4 = fig.add_subplot(4,2,4)
ax4.plot(Z[0].times(), Hil_Z[0].data)
ax4.text(0.01, 0.95, Z[0].id + ' Hilbert-transformed', transform=ax4.transAxes, fontsize=8,
        verticalalignment='top', bbox=props)
ax4.set_ylabel('m/s')
ax5 = fig.add_subplot(4,2,6)
ax5.plot(t, Rotated[0])
ax5.text(0.01, 0.95, tr_1[0].id + ' rotated', transform=ax5.transAxes, fontsize=8,
        verticalalignment='top', bbox=props)
ax5.set_ylabel('m/s')
ax5.set_xlabel('Seconds since 2021-05-31T04:13:00')
ax6 = fig.add_subplot(4,2,8)
ax6.plot(maxSrz)
ax6.plot(degrees, corr_coeff, "o")
ax6.vlines(degrees, min(maxSrz)-0.5, corr_coeff, linestyles='dashed', color='k', label='R = '+str(np.round(corr_coeff,2)))
ax6.set_ylim(min(maxSrz)-0.3,max(maxSrz)+0.8)
ax6.text(0.01, 0.95, 'Cross Correlation Function', transform=ax6.transAxes, fontsize=8,
        verticalalignment='top', bbox=props)
ax6.text(degrees, corr_coeff + 0.1, str(degrees) + ' deg (clockwise)')
ax6.legend(loc='best')
ax6.set_xlabel('Degrees')
ax6.set_ylabel('Normalized Corr Coeff')




# fig.tight_layout()
# plt.show()
plt.savefig('T120-BH1_orient_URZ_2021p405872.png', format='PNG', dpi=400)
    



























    # find the angle with the maximum correlation:	
   #  maxmaxSrz= max(maxSrz)
  	# rotangle_index = maxSrz.index(maxmaxSrz)
   # 	rotangle = thetas[rotangle_index]




# t1 = UTCDateTime(2021, 4, 1, 19, 0, 0)
# duration = 100







# def rotatehorizontal(stream, angle1, angle2):
#     """
#     A function to rotate the horizontal components of a seismometer from 
#     radial and transverse into E and North components.
#     """
#     if stream[0].stats.channel in set(['HHE', 'HHN']):
#         stream.sort(['channel'], reverse=True)
#         angle1, angle2 = angle2, angle1
#         print(stream)
#         print(angle1)
#         print(anngle2)
#     theta_r1 = math.radians(angle1)
#     theta_r2 = math.radians(angle2)
#     swapSecond = False
#     if (angle2 >= 180. and angle2 <= 360.) or angle2 == 0.:
#         swapSecond = True 
#     # if the components are swaped swap the matrix
#     if theta_r1 > theta_r2 and swapSecond:
#         if debugRot:
#             print('Swap the components: ' + str((360. - angle1) - angle2)
#         stream.sort(['channel'], reverse=True)
#         theta_r1, theta_r2 = theta_r2, theta_r1
#         print(stream)
#     # create new trace objects with same info as previous
#     rotatedN = stream[0].copy()
#     rotatedE = stream[1].copy()
#     # assign rotated data
#     rotatedN.data = stream[0].data*math.cos(-theta_r1) +\
#         stream[1].data*math.sin(-theta_r1)
#     rotatedE.data = -stream[1].data*math.cos(-theta_r2-math.pi/2.) +\
#         stream[0].data*math.sin(-theta_r2-math.pi/2.)
#     rotatedN.stats.channel = 'LHN'
#     rotatedE.stats.channel = 'LHE'
#     # return new streams object with rotated traces
#     streamsR = Stream(traces=[rotatedN, rotatedE])
#     return streamsR        