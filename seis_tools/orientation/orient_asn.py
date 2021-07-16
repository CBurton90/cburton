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
	ZZ = st_folded.select(id="NZ.URZ.10.HHZ")
	RZ = st_folded.select(id="NZ.URZ.10.HH1")
	TZ = st_folded.select(id="NZ.URZ.10.HH2")

	return ZZ,RZ,TZ

	
ZZ,RZ,TZ = fold_ccfs(zz='URZ_TLZ_ZZ.mseed',nz='URZ_TLZ_1Z.mseed',ez='URZ_TLZ_2Z.mseed')

print(ZZ)

# fig, axs = plt.subplots(3)
# axs[0].plot(ZZ[0].data)
# axs[1].plot(RZ[0].data)
# axs[2].plot(TZ[0].data)

# st2 = read('URZ_TLZ_ZZ.mseed')
# plt.plot(st2[0].data)
# plt.show()


# we are going to compare the RT with a 90 degree
                # shifted ZZ


ZZ_hil = np.imag(hilbert(ZZ[0].data))
Szz = np.correlate(ZZ_hil,ZZ_hil)

# rotating ANTI-CLOCKWISE, correlating with ZZ_hil:
maxSrz= []
#coherences=[]
thetas = np.linspace(0,2*np.pi,360)
for i_,theta in enumerate(thetas):
	RZ_rot =  np.cos(theta)*RZ[0].data - np.sin(theta)*TZ[0].data
	TZ_rot = np.sin(theta)*RZ[0].data + np.cos(theta)*TZ[0].data
	Srz = np.correlate(RZ_rot, ZZ_hil)
	# Szz = np.correlate(ZZ_hil, ZZ_hil)
	maxSrz.append(max(Srz)/max(Szz))

# find the angle with the maximum correlation:
maxmaxSrz= max(maxSrz)
rotangle_index = maxSrz.index(maxmaxSrz)
rotangle = thetas[rotangle_index]
degrees = int(rotangle*180/np.pi)
corr_coeff = maxmaxSrz
print(degrees)
print(corr_coeff)


URZ_lat = -38.2592
URZ_long = 177.1109
TLZ_lat = -38.3294
TLZ_long = 175.538

back_azi = gps2dist_azimuth(URZ_lat, URZ_long, TLZ_lat, TLZ_long)
print(back_azi)
source_azi = back_azi[1]

print(source_azi)

#anti-clockwise-hh1-rot
correction = ((source_azi + degrees) + 360) % 360
print(correction)

# plt.plot(maxSrz)
# plt.show()














fig = plt.figure(figsize=(25,20))
fig.suptitle('Ambient Seismic Noise Orientation', fontsize=18)


# Plot the basemap of event and station

plt.subplots_adjust(wspace= 0.3, hspace= 0.2)

ax1 = fig.add_subplot(2, 2, 1)

m = Basemap(llcrnrlon=170.,llcrnrlat=-40.,urcrnrlon=180.,urcrnrlat=-35.,\
            rsphere=(6378137.00,6356752.3142),\
            resolution='l',projection='merc',\
            lat_0=-37.5,lon_0=175.,lat_ts=None)

# nylat, nylon are lat/lon of New York
nylat = URZ_lat; nylon = URZ_long
# lonlat, lonlon are lat/lon of London.
lonlat = TLZ_lat; lonlon = TLZ_long
# draw great circle route between NY and London
m.drawgreatcircle(nylon,nylat,lonlon,lonlat,linewidth=2,color='b')
m.drawcoastlines()
m.fillcontinents()
# draw parallels
m.drawparallels(np.arange(-40,-35,1),labels=[1,1,0,1])
# draw meridians
m.drawmeridians(np.arange(-180,180,2.5),labels=[1,1,0,1])
ax1.set_title('URZ to TLZ', pad=10)




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
#plt.savefig('T120-BH1_orient_URZ_TLZ.png', format='PNG', dpi=400)

