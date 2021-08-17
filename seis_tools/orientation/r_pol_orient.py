#!/usr/bin/env python3

import argparse
import cartopy.feature
import cartopy.crs as ccrs
from obspy import UTCDateTime
from obspy import Stream
from obspy.clients.fdsn import Client
from obspy import read, read_inventory
from obspy.signal.cross_correlation import correlate
from obspy.geodetics.base import gps2dist_azimuth
from obspy.taup.taup_geo import calc_dist_azi
from obspy.taup import TauPyModel
import numpy as np
from scipy.signal import hilbert
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap


def get_and_remove_response(station, channel, location, output, t1, duration):
    client = Client("http://service-nrt.geonet.org.nz")
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



parser = argparse.ArgumentParser()
parser.add_argument('--time', type=str, required=True)
parser.add_argument('--ev_lat', type=float, required=True)
parser.add_argument('--ev_long', type=float, required=True)
parser.add_argument('--ev_code', type=str, required=True)
parser.add_argument('--Mww', type=str, required=True)
parser.add_argument('--site', type=str, required=True)
parser.add_argument('--site_lat', type=float, required=True)
parser.add_argument('--site_long', type=float, required=True)
parser.add_argument('--location', type=str, required=True)
parser.add_argument('--rad_comp', type=str, required=True)
parser.add_argument('--trv_comp', type=str, required=True)
# parser.add_argument('--site', type=str, required=True)
args = parser.parse_args()

t1 = args.time
event_lat = args.ev_lat
event_long = args.ev_long
ev_code = args.ev_code
Mww = args.Mww
site = args.site
location = args.location
radial = args.rad_comp
transverse = args.trv_comp
URZ_lat = args.site_lat 
URZ_long = args.site_long


# epi_dist = calc_dist_azi(event_lat, event_long, URZ_lat, URZ_long, 6371, 0)
# print(epi_dist)
# model = TauPyModel(model="iasp91")
# arrivals = model.get_travel_times(source_depth_in_km=44, distance_in_degree=epi_dist[0])
# print(arrivals)

back_azi = gps2dist_azimuth(URZ_lat, URZ_long, event_lat, event_long)
event_back_azi = back_azi[1]
great_circ_dist = back_azi[0]
trav_tm = great_circ_dist / 4000

t1 = UTCDateTime(t1) + trav_tm
duration = 60*10



st = get_and_remove_response(station=site, channel="HH*", location=location, output="VEL", t1=t1-25, duration=duration)
st.trim(starttime=t1-20, endtime=t1+duration-1)

# print(st)
# st.plot()


# rotating ANTI-CLOCKWISE, correlating with ZZ90:
maxSrz= []
#coherences=[]
thetas = np.linspace(0,2*np.pi,360)
for i_,theta in enumerate(thetas):
    Rad_rot =  np.cos(theta)*st.select(id="NZ."+site+"."+location+"."+radial) - np.sin(theta)*st.select(id="NZ."+site+"."+location+"."+transverse)
    Trv_rot = np.sin(theta)*st.select(id="NZ."+site+"."+location+"."+radial) + np.cos(theta)*st.select(id="NZ."+site+"."+location+"."+transverse)
    Hil_rad = np.imag(hilbert(Rad_rot))
    Z = st.select(id="NZ."+site+"."+location+".HHZ")
    Hil_Z = np.imag(hilbert(Z))
    
    # numpy
    
    # hilbert tform of vertical corr with radial
    Ncorr = np.correlate(Hil_Z[0].data,Rad_rot[0].data)
    Szz = np.correlate(Hil_Z[0], Hil_Z[0])
    
    # # hilbert tfrom of radial corr with vertical
    # Ncorr = np.correlate(Z[0].data,Hil_rad[0].data)
    # Szz = np.correlate(Z[0].data, Z[0].data)


    # Snn = np.correlate(Hil_rad[0].data, Hil_rad[0].data)
    
    # Snn = np.correlate(Hil_rad[0],Hil_rad[0])
    #obspy
    # Ncorr = correlate(Rad_rot[0].data, Hil_Z[0].data, 0)
    # Szz = correlate(Hil_Z[0], Hil_Z[0], 0)
    
    #trying reverse corr for order and hil-tf
    # Ncorr = correlate(Z[0].data, Hil_rad[0].data, 0)
    # Szz = correlate(Hil_Z[0], Hil_Z[0], 0)
    # Ncorr = correlate(Hil_Z[0].data, Rad_rot[0].data, 0)

    maxSrz.append(max(Ncorr)/max(Szz))
    # maxSrz.append(max(Ncorr)/np.sqrt((max(Szz)*max(Snn))))
    # maxSrz.append(max(Ncorr))

    
    
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

#clockwise-hh1-rot
# correction = ((event_back_azi - degrees) + 360) % 360
#anti-clockwise-hh1-rot
correction = ((event_back_azi + degrees) + 360) % 360
print(correction)



# # PLOTTING



Rotated =  np.cos(rotangle)*st.select(id="NZ."+site+"."+location+"."+radial) - np.sin(rotangle)*st.select(id="NZ."+site+"."+location+"."+transverse)
Rotated_Hil = np.imag(hilbert(Rotated))

tr_1 = st.select(id="NZ."+site+"."+location+"."+radial)
t = tr_1[0].times()




fig = plt.figure(figsize=(25,20))
fig.suptitle('Broadband seismometer orientation at '+site+' with event '+ev_code, fontsize=18)


# Plot the basemap of event and station

plt.subplots_adjust(wspace= 0.3, hspace= 0.2)

ax1 = fig.add_subplot(2, 2, 1,projection=ccrs.PlateCarree(central_longitude=180))

ax1.set_extent([110, 320, -60, 80], crs=ccrs.PlateCarree())

# add some features to make the map a little more polished
ax1.add_feature(cartopy.feature.LAND)
ax1.add_feature(cartopy.feature.OCEAN)
ax1.coastlines('50m')
ax1.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--', crs=ccrs.PlateCarree())
ax1.plot([URZ_long, event_long], [URZ_lat, event_lat], color='blue',  transform=ccrs.Geodetic())
ax1.set_title(Mww+' Mww to '+ site +' ('+args.time+')', pad=10)

# # Background map
# m=Basemap(llcrnrlon=170, llcrnrlat=-50, urcrnrlon=-60, urcrnrlat=70, projection='merc')
# m.drawmapboundary(fill_color='#A6CAE0', linewidth=0)
# m.fillcontinents(color='grey', alpha=0.7, lake_color='grey')
# m.drawcoastlines(linewidth=0.1, color="white")
# # nylat, nylon are lat/lon of New York
# nylat = URZ_lat; nylon = URZ_long
# # lonlat, lonlon are lat/lon of London.
# lonlat = event_lat; lonlon = event_long
# # draw great circle route between NY and London
# m.drawgreatcircle(nylon,nylat,lonlon,lonlat,linewidth=2,color='b')


# m = Basemap(llcrnrlon=-180.,llcrnrlat=-50.,urcrnrlon=180.,urcrnrlat=70.,\
#             rsphere=(6378137.00,6356752.3142),\
#             resolution='l',projection='merc',\
#             lat_0=0.,lon_0=-150.,lat_ts=None)

# # nylat, nylon are lat/lon of New York
# nylat = URZ_lat; nylon = URZ_long
# # lonlat, lonlon are lat/lon of London.
# lonlat = event_lat; lonlon = event_long
# # draw great circle route between NY and London
# m.drawgreatcircle(nylon,nylat,lonlon,lonlat,linewidth=2,color='b')
# m.drawcoastlines()
# m.fillcontinents()
# # draw parallels
# m.drawparallels(np.arange(-50,70,20),labels=[1,1,0,1])
# # draw meridians
# m.drawmeridians(np.arange(-180,180,30),labels=[1,1,0,1])
# ax1.set_title(Mww+' Mww to '+ site +' ('+args.time+')', pad=10)







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
    ax2.plot(x, y, linewidth=2, label='Orientation'+ ' = ' + str(np.around(correction,0)))


#anticlockwise

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
ax4.text(0.01, 0.95, Z[0].id + ' Hilbert-transform', transform=ax4.transAxes, fontsize=8,
        verticalalignment='top', bbox=props)
ax4.set_ylabel('m/s')
ax5 = fig.add_subplot(4,2,6)
ax5.plot(t, Rotated[0])
ax5.text(0.01, 0.95, tr_1[0].id + ' rotated', transform=ax5.transAxes, fontsize=8,
        verticalalignment='top', bbox=props)
# ax5.set_xlabel('Seconds since 2021-06-06T23:21:10')
ax5.set_ylabel('m/s')
ax6 = fig.add_subplot(4,2,8)
ax6.plot(maxSrz)
ax6.plot(degrees, corr_coeff, "o")
ax6.vlines(degrees, min(maxSrz)-0.5, corr_coeff, linestyles='dashed', color='k', label='R = '+str(np.round(corr_coeff,2)))
ax6.set_ylim(min(maxSrz)-0.3,max(maxSrz)+0.8)
ax6.text(0.01, 0.95, 'Cross Correlation Function', transform=ax6.transAxes, fontsize=8,
        verticalalignment='top', bbox=props)
ax6.text(degrees, corr_coeff + 0.1, str(degrees) + ' deg (anti-clockwise)')
ax6.legend(loc='best')
ax6.set_xlabel('Degrees')
ax6.set_ylabel('Normalized Corr Coeff')




# fig.tight_layout()
plt.show()
# plt.savefig('Borehole_orient_'+site+'_'+ev_code+'_'+str(correction)+'.png', format='PNG', dpi=400)
    



























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