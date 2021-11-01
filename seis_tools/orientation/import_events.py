#!/usr/bin/env python3

import argparse
import cartopy.feature
import cartopy.crs as ccrs
import pandas as pd
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
from scipy.stats import circmean,circstd
import scipy.stats as sts
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

def get_and_remove_response(station, channel, location, output, t1, duration):
    client = Client("http://service.geonet.org.nz")
    st = client.get_waveforms(
        network="NZ", station=station, location=location,
        channel=channel, starttime=t1, endtime=t1 + duration)
    tr = Stream()
    # st.merge()
    # print(st)
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



# Param's

site = 'RIZ'
location = '10'
URZ_lat = -29.2449
URZ_long = -177.9289
radial = 'HHN'
transverse = 'HHE'
# t1 = UTCDateTime("2021-05-28T00:00:00")
t1 = UTCDateTime("2020-01-01T00:00:00")
t2 = UTCDateTime("2020-12-31T00:00:00")

# Get a catalog of events we want

client = Client("USGS")
cat = client.get_events(starttime=t1, endtime=t2, minmagnitude=6, maxdepth=50, latitude=URZ_lat, longitude=URZ_long, minradius=10, maxradius=180)
# print(cat.__str__(print_all=True))
# cat.plot()

# Loop through origin times and coords of events

angles = []
df = pd.DataFrame(columns = ['USGS eventID','Origin Time','Magnitude','Czr*','Orientation'])

# Plotting

fig = plt.figure(figsize=(25,20))
fig.suptitle(site+' Seismometer Orientation via Rayleigh-Wave Polarization', fontsize=18)

ax1 = fig.add_subplot(2, 2, 1,projection=ccrs.PlateCarree(central_longitude=180))
ax1.set_extent([60, 320, -80, 80], crs=ccrs.PlateCarree())
ax1.add_feature(cartopy.feature.LAND)
ax1.add_feature(cartopy.feature.OCEAN)
ax1.coastlines('50m')
ax1.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--', crs=ccrs.PlateCarree())
ax1.set_title('Selected Teleseisms (>Mw 6.0, >10deg, <50km depth)', pad=10)

ax2 = fig.add_subplot(2, 2, 3, projection='polar')
ax2.set_yticklabels([])
ax2.set_theta_zero_location('N')
# Go clockwise
ax2.set_theta_direction(-1)

ax3 = fig.add_subplot(1,2,2)

for n in range(len(cat)):
	back_azi = gps2dist_azimuth(URZ_lat, URZ_long, cat[n].origins[0].latitude, cat[n].origins[0].longitude)
	event_back_azi = back_azi[1]
	great_circ_dist = back_azi[0]
	trav_tm = great_circ_dist / 4000

	t1 = cat[n].origins[0].time + trav_tm
	duration = 60*10
	print(t1)

	try:
		st = get_and_remove_response(station=site, channel="HH*", location=location, output="VEL", t1=t1-25, duration=duration)
		st.trim(starttime=t1-20, endtime=t1+duration-1)
		# st.plot()

		# rotating ANTI-CLOCKWISE, correlating with ZZ90:
		maxSrz= []
		norm_Srz = []
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

		    Snn = np.correlate(Rad_rot[0].data, Rad_rot[0].data)
		    
		    # Snn = np.correlate(Hil_rad[0],Hil_rad[0])
		    #obspy
		    # Ncorr = correlate(Rad_rot[0].data, Hil_Z[0].data, 0)
		    # Szz = correlate(Hil_Z[0], Hil_Z[0], 0)
		    
		    #trying reverse corr for order and hil-tf
		    # Ncorr = correlate(Z[0].data, Hil_rad[0].data, 0)
		    # Szz = correlate(Hil_Z[0], Hil_Z[0], 0)
		    # Ncorr = correlate(Hil_Z[0].data, Rad_rot[0].data, 0)

	        # unbounded
		    maxSrz.append(max(Ncorr)/max(Szz))
		    # normalised
		    norm_Srz.append(max(Ncorr)/np.sqrt((max(Szz)*max(Snn))))
		    # maxSrz.append(max(Ncorr))

		    
		    
		# print(len(maxSrz))
		maxmaxSrz= max(maxSrz)
		max_norm = max(norm_Srz)
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

		if  max_norm > 0.80:
			print(max_norm)
			correction_rad = np.deg2rad(correction)
			angles.append(correction_rad)
			ax1.plot([URZ_long, cat[n].origins[0].longitude], [URZ_lat, cat[n].origins[0].latitude], color='blue',  transform=ccrs.Geodetic())
			ax2.plot([0, correction_rad], [0, 1], linewidth=2)
			df.loc[n] = [str(cat[n].resource_id)[-25:-15], cat[n].origins[0].time, cat[n].magnitudes[0].mag, maxmaxSrz, correction]
	
	except Exception as e:
		print('No data available from '+site+' '+str(t1))

		
		

print(angles)
print(np.unwrap(angles))


circmeanangle = int(np.round(circmean(angles)*180/np.pi))
circmeanstd = int(np.round(circstd(angles)*180/np.pi))
circmeanangle2 = int(np.round(circmean(np.unwrap(angles))*180/np.pi))
circmeanstd2 = int(np.round(circstd(np.unwrap(angles))*180/np.pi))
ax2.plot([0, circmean(angles)], [0, 1], linewidth=3, linestyle='dotted', color='k', label='Circular mean'+ ' = ' + str(np.round(circmeanangle,0)))
angle = np.deg2rad(202.5)
ax2.legend(loc='upper right', bbox_to_anchor=(.5 + np.cos(angle)/2, .5 + np.sin(angle)/2), fontsize=9)
print(circmeanangle,circmeanstd)
print(circmeanangle2, circmeanstd2)

print(df)

std_err_2 = circstd(angles) / np.sqrt(len(angles))
c_95_2 = sts.t.interval(0.95, len(angles)-1, loc=circmean(angles), scale=std_err_2)
print('conf int rad', c_95_2)
print('conf int deg', np.rad2deg(c_95_2))


std_err = circstd(np.unwrap(angles)) / np.sqrt(len(angles))
c_95 = sts.t.interval(0.95, len(angles)-1, loc=circmean(np.unwrap(angles)), scale=std_err)
print('unwrapped conf int rad', c_95)
print('unwrapped conf int deg', np.rad2deg(c_95))

c_95_edited = []

# rebuild confidence interval if there are negative radians from numpy's unwrap by adding 2pi

for n in range(len(c_95)):
	if c_95[n] < 0:
		a = c_95[n] + 2*np.pi
		print('conf int under zero (+2pi)', a)
		c_95_edited.append(a)
	else:
		print('conf int over zero', c_95[n])
		c_95_edited.append(c_95[n])

print('conf int deg (0-360)', np.rad2deg(c_95_edited))
lower_C95 = np.rad2deg(c_95_edited[0])
print(lower_C95)
upper_C95 = np.rad2deg(c_95_edited[1])
print(upper_C95)

# df2 = df[(df['Orientation'] >= int(lower_C95)) & (df['Orientation'] <= int(upper_C95))]
# print(df2)


ax3.axis('off')
ax3.axis('tight')
table = ax3.table(cellText=df.values, colLabels=df.columns, loc='center')
table.scale(1.32,1.45)
table.auto_set_font_size(False)
table.set_fontsize(8)
# plt.show()
plt.savefig(str(site)+'_orient_'+str(t1)+'.png', format='PNG', dpi=400)