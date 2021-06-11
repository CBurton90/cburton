#!/usr/bin/env python3

from obspy import UTCDateTime
from obspy import Stream
from obspy.clients.fdsn import Client
from obspy import read, read_inventory
from obspy.signal.cross_correlation import correlate
import numpy as np
import matplotlib.pyplot as plt


def get_and_remove_response(stream, inventory, output, t1, duration):
	st = stream
	tr = Stream()
	st.merge()
	# print(st)
	for n in range(len(st)):
		st[n].trim(starttime=t1, endtime=duration)
		st[n].detrend('linear')
		st[n].taper(max_percentage=0.05, type='cosine')
		inv = inventory
		pre_filt = (0.05, 0.1, 0.2, 0.4)
		st[n].remove_response(output=output, pre_filt=pre_filt, plot=False,
			water_level=60, inventory=inv)
		tr += st[n]

	return tr

t1 = UTCDateTime("2021-05-26T14:20:00")

st = read("/home/conradb/seis_tools/orientation/URZ_data/Borehole/2021.146.URZ.10-HH*.NZ.D.mseed")
inv = read_inventory("/home/conradb/seis_tools/orientation/URZ_data/Borehole/NZ.URZ.10.HH*.xml")
st = get_and_remove_response(stream=st, inventory=inv, output='VEL', t1=t1, duration=t1 + 60*4)


st_2 = read("/home/conradb/seis_tools/orientation/URZ_data/Surface_sensor/SS01_centaur-6_8003_20210526_011856.seed")
st_2 = st_2.select(id="NZ.SS01.10.HH*")
inv2 = read_inventory("/home/conradb/seis_tools/orientation/URZ_data/Surface_sensor/NZ.SS01.10.HH*.xml")
st_2 = get_and_remove_response(stream=st_2, inventory=inv2, output='VEL', t1=t1, duration=t1 + 60*4)

print(st)
print(st_2)

# rotating CLOCKWISE, correlating with ZZ90:
maxSrz= []
#coherences=[]
thetas = np.linspace(0,2*np.pi,360)
for i_,theta in enumerate(thetas):
	Rad_rot =  np.cos(theta)*st.select(id="NZ.URZ.10.HH1") - np.sin(theta)*st.select(id="NZ.URZ.10.HH2")
	Trv_rot = np.sin(theta)*st.select(id="NZ.URZ.10.HH1") + np.cos(theta)*st.select(id="NZ.URZ.10.HH2")
	N_surface = st_2.select(id="NZ.SS01.10.HHN")
	# Ncorr = np.correlate(Rad_rot[0].data, N_surface[0].data)
	Ncorr = correlate(Rad_rot[0].data, N_surface[0].data, 0)
	# Snn = np.correlate(N_surface[0].data, N_surface[0].data)
	# maxSrz.append(max(Ncorr)/max(Snn))
	maxSrz.append(max(Ncorr))
# print(len(maxSrz))
maxmaxSrz= max(maxSrz)
rotangle_index = maxSrz.index(maxmaxSrz)
# print(rotangle_index)
rotangle = thetas[rotangle_index]
radians = rotangle
degrees = int(rotangle*180/np.pi)
corr_coeff = maxmaxSrz







# PLOTTING



Final =  np.cos(rotangle)*st.select(id="NZ.URZ.10.HH1") - np.sin(rotangle)*st.select(id="NZ.URZ.10.HH2")

fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, sharex=False, figsize=(15,15))
fig.suptitle('T120-BH1 Orientation URZ', fontsize=18)

tr_1 = st.select(id="NZ.URZ.10.HH1")
t = tr_1[0].times()


props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax1.plot(tr_1[0].times(), tr_1[0].data)
ax1.text(0.01, 0.95, tr_1[0].id, transform=ax1.transAxes, fontsize=8,
        verticalalignment='top', bbox=props)
ax2.plot(N_surface[0].times(), N_surface[0].data)
ax2.text(0.01, 0.95, N_surface[0].id, transform=ax2.transAxes, fontsize=8,
        verticalalignment='top', bbox=props)
ax3.plot(t, Final[0])
ax3.text(0.01, 0.95, tr_1[0].id + ' rotated', transform=ax3.transAxes, fontsize=8,
        verticalalignment='top', bbox=props)
ax4.plot(maxSrz)
ax4.plot(degrees, corr_coeff, "o")
ax4.vlines(degrees, min(maxSrz)-0.5, corr_coeff, linestyles='dashed', color='k')
ax4.set_ylim(min(maxSrz)-0.2,max(maxSrz)+0.4)
ax4.text(0.01, 0.95, 'Cross Correlation Function', transform=ax4.transAxes, fontsize=8,
        verticalalignment='top', bbox=props)
ax4.text(degrees, corr_coeff + 0.1, str(degrees) + ' deg (clockwise)')


fig.tight_layout()
plt.show()
# plt.savefig('T120-BH1_orient_NIC.png', format='PNG', dpi=400)
    



























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