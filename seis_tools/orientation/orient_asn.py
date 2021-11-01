#!/usr/bin/env python3

import argparse
from obspy.core.inventory.inventory import read_inventory
from obspy import read
from obspy.geodetics.base import gps2dist_azimuth
import os
import numpy as np
from scipy.signal import hilbert
from scipy.stats import circmean,circstd
from scipy.interpolate import interp1d
import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

def fold_ccfs(zz,nz,ez):
	st = read(zz)
	st += read(nz)
	st += read(ez)

	# From one trace, pick off some constants we'll need:
	starttime = st[0].stats.starttime
	delta  = st[0].stats.delta
	npts2 = st[0].stats.npts
	npts = int((npts2-1)/2)

	#Improve ccf SNR by folding
	st_folded= st.copy()
	for i_,tr in enumerate(st):
		causal =  st[i_].data[npts:-1]
		acausal = st[i_].data[npts:0:-1]
		st_folded[i_].data[0:npts] =  (causal+acausal)/2


	#Keep only the first half of the stream, as that is where all the action is.
	st_folded.trim(starttime=starttime, endtime=starttime+npts*delta)
	st_folded.taper(max_percentage=0.05, type='cosine')

	# Order ccf streams
	st_folded.sort() 
	ZZ = st_folded.select(id='NZ.'+site_A+'.10.HHZ')
	RZ = st_folded.select(id='NZ.'+site_A+'.10.HH1')
	TZ = st_folded.select(id='NZ.'+site_A+'.10.HH2')

	return ZZ,RZ,TZ

	

parser = argparse.ArgumentParser()

parser.add_argument('--site_ref', type=str, required=True)
parser.add_argument('--ref_lat', type=float, required=True)
parser.add_argument('--ref_long', type=float, required=True)

parser.add_argument('--site_source', type=str, required=True)
parser.add_argument('--source_lat', type=float, required=True)
parser.add_argument('--source_long', type=float, required=True)

args = parser.parse_args()

site_A = args.site_ref
lat_A = args.ref_lat
long_A = args.ref_long
site_B = args.site_source
lat_B = args.source_lat
long_B = args.source_long

path_zz = '/home/conradb/git/cburton/seis_tools/orientation/final_impulses/'+site_A+'_'+site_B+'_ZZ.mseed'
path_nz = '/home/conradb/git/cburton/seis_tools/orientation/final_impulses/'+site_A+'_'+site_B+'_1Z.mseed'
path_ez = '/home/conradb/git/cburton/seis_tools/orientation/final_impulses/'+site_A+'_'+site_B+'_2Z.mseed'


ZZ,RZ,TZ = fold_ccfs(zz=path_zz,nz=path_nz,ez=path_ez)

print(ZZ)
print(RZ)
print(TZ)

# fig, axs = plt.subplots(3)
# axs[0].plot(ZZ[0].data)
# axs[1].plot(RZ[0].data)
# axs[2].plot(TZ[0].data)

# st2 = read('URZ_TLZ_ZZ.mseed')
# plt.plot(st2[0].data)
# plt.show()


# Comparing the unknown radial comp with a 90 degree shifted ZZ in 1 degree steps


ZZ_hil = np.imag(hilbert(ZZ[0].data))
Szz = np.correlate(ZZ_hil,ZZ_hil)

# rotating anti-clockwise and correlating with ZZ_hil:

maxSrz= []

#coherences=[]
# build values from 0 to 2pi in 1/360 increments
thetas = np.linspace(0,2*np.pi,360)
for i_,theta in enumerate(thetas):
	RZ_rot =  np.cos(theta)*RZ[0].data - np.sin(theta)*TZ[0].data
	TZ_rot = np.sin(theta)*RZ[0].data + np.cos(theta)*TZ[0].data
	Srz = np.correlate(RZ_rot, ZZ_hil)
	Szz = np.correlate(ZZ_hil, ZZ_hil)
	maxSrz.append(max(Srz)/max(Szz))

# find the angle with the maximum correlation:
maxmaxSrz= max(maxSrz)
# index of list with maximum correlation
rotangle_index = maxSrz.index(maxmaxSrz)
# value in radians
rotangle = thetas[rotangle_index]
# value in degrees
degrees = int(rotangle*180/np.pi)
# correlation coefficient value
corr_coeff = maxmaxSrz
print(degrees)
print(corr_coeff)


back_azi = gps2dist_azimuth(lat_A, long_A, lat_B, long_B)
print(back_azi)
source_azi = back_azi[1]
print(source_azi)

#anti-clockwise-hh1-rot
correction = ((source_azi + degrees) + 360) % 360
print(correction)


# plotting

fig = plt.figure(figsize=(25,20))
fig.suptitle('Ambient Seismic Noise Orientation', fontsize=18)


# Plot the basemap of event and station

plt.subplots_adjust(wspace= 0.3, hspace= 0.2)

ax1 = fig.add_subplot(2, 2, 1)

m = Basemap(llcrnrlon=long_A-5,llcrnrlat=lat_A-2.5,urcrnrlon=long_A+5.,urcrnrlat=lat_A+2.5,\
            rsphere=(6378137.00,6356752.3142),\
            resolution='l',projection='merc',\
            lat_0=lat_A,lon_0=long_A,lat_ts=None)

# draw great circle route between refrence and source sites
m.drawgreatcircle(long_A,lat_A,long_B,lat_B,linewidth=2,color='b')
m.drawcoastlines()
m.fillcontinents()
# draw parallels
m.drawparallels(np.arange(-90,0,1),labels=[1,1,0,1])
# draw meridians
m.drawmeridians(np.arange(-180,180,2.5),labels=[1,1,0,1])
ax1.set_title(site_A+' to '+site_B, pad=10)




# Plot the polar plot showing event back azi and sensor orientation

ax2 = fig.add_subplot(2, 2, 3, projection='polar')

r = [1, 1,]
theta = [np.radians(correction), np.radians(source_azi)]
labels = ['HH1','Station Source']


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

for curve in [[[source_azi,correction], [0.6, 0.6]]]:
	# ax2.set_theta_direction(1)
	curve[0] = np.deg2rad(curve[0])
	x = np.linspace( curve[0][0], curve[0][1], 500)
	y = interp1d( curve[0], curve[1])( x)
	ax2.plot(x, y, linewidth=2, label='BAZ measured'+ ' = ' + str(np.round(degrees,0)))

for curve in [[[0, source_azi], [0.8, 0.8]]]:
    curve[0] = np.deg2rad(curve[0])
    x = np.linspace( curve[0][0], curve[0][1], 500)
    y = interp1d( curve[0], curve[1])( x)
    ax2.plot(x, y, linewidth=2, label='BAZ expected'+ ' = ' + str(np.round(source_azi,0)))

angle = np.deg2rad(337.5)
ax2.legend(loc='upper left', bbox_to_anchor=(.5 + np.cos(angle)/2, .5 + np.sin(angle)/2))


# Plot the 1Z impluse, RZ impulse, Hilbert transformed ZZ impulse and CCF


props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax3 = fig.add_subplot(4,2,2)
ax3.plot(RZ[0].data)
ax3.text(0.01, 0.95, 'G1z', transform=ax3.transAxes, fontsize=8,
        verticalalignment='top', bbox=props)
ax4 = fig.add_subplot(4,2,4)
ax4.plot(ZZ_hil)
ax4.text(0.01, 0.95,'H[Gzz]', transform=ax4.transAxes, fontsize=8,
        verticalalignment='top', bbox=props)
ax5 = fig.add_subplot(4,2,6)
ax5.plot(np.cos(rotangle)*RZ[0].data - np.sin(rotangle)*TZ[0].data)
ax5.text(0.01, 0.95,'Grz', transform=ax5.transAxes, fontsize=8,
        verticalalignment='top', bbox=props)


ax6 = fig.add_subplot(4,2,8)
ax6.plot(maxSrz)
ax6.plot(degrees, corr_coeff, "o")
ax6.vlines(degrees, min(maxSrz)-0.5, corr_coeff, linestyles='dashed', color='k', label='R = '+str(np.round(corr_coeff,2)))
ax6.set_ylim(min(maxSrz)-0.3,max(maxSrz)+0.8)
ax6.text(0.01, 0.95, 'Normalized Cross Correlation Function', transform=ax6.transAxes, fontsize=8,
        verticalalignment='top', bbox=props)
ax6.text(degrees, corr_coeff + 0.1, str(degrees) + ' deg (anti-clockwise)')
ax6.legend(loc='best')
ax6.set_xlabel('Degrees')
ax6.set_ylabel('Srz')




# fig.tight_layout()
plt.show()
# plt.savefig('asn_orient_'+str(site_A)+'_'+str(site_B)+'_'+str(correction)+'.png', format='PNG', dpi=400)

