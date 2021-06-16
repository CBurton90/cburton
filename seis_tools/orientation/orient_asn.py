from obspy.core.inventory.inventory import read_inventory
from obspy import read
import os
import numpy as np
from scipy.signal import hilbert
from scipy.stats import circmean,circstd
import matplotlib.pyplot as plt

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
	# st_folded.taper(max_percentage=0.1, type='cosine')

	# Order ccf streams
	st_folded.sort() 
	ZZ = st_folded.select(id="NZ.URZ.10.HHZ")
	RZ = st_folded.select(id="NZ.URZ.10.HH1")
	TZ = st_folded.select(id="NZ.URZ.10.HH2")

	return ZZ,RZ,TZ

	




ZZ,RZ,TZ = fold_ccfs(zz='URZ_PUZ_ZZ.mseed',nz='URZ_PUZ_1Z.mseed',ez='URZ_PUZ_2Z.mseed')

print(ZZ)

# fig, axs = plt.subplots(3)
# axs[0].plot(ZZ[0].data)
# axs[1].plot(RZ[0].data)
# axs[2].plot(TZ[0].data)

# st2 = read('URZ_PUZ_ZZ.mseed')
# plt.plot(st2[0].data)
# plt.show()


# we are going to compare the RT with a 90 degree
                # shifted ZZ


ZZ_hil = np.imag(hilbert(ZZ[0].data))
Szz = np.correlate(ZZ_hil,ZZ_hil)

# rotating CLOCKWISE, correlating with ZZ_hil:
maxSrz= []
#coherences=[]
thetas = np.linspace(0,2*np.pi,360)
for i_,theta in enumerate(thetas):
	RZ_rot =  np.cos(theta)*RZ[0].data - np.sin(theta)*TZ[0].data
	TZ_rot = np.sin(theta)*RZ[0].data + np.cos(theta)*TZ[0].data
	Srz = np.correlate(RZ_rot, ZZ_hil)
	Szz = np.correlate(ZZ_hil, ZZ_hil)
	maxSrz.append(max(Srz)/max(Szz))

plt.plot(maxSrz)
plt.show()