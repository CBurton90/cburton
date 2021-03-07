#!/usr/bin/env python3
from obspy.signal.spectral_estimation import get_nlnm, get_nhnm
from obspy.clients.fdsn import Client
from obspy.core import UTCDateTime, Stream
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from scipy import signal
import sys
# mpl.rc('font',family='serif')
# mpl.rc('font',serif='Times')
# mpl.rc('text', usetex=True)
# mpl.rc('font',size=18)

def rotate(data1, data2, data3):
    # create new trace objects with same info as previous
    rotatedZ = data1.copy()
    rotatedN = data2.copy()
    rotatedE = data3.copy()
    rotatedE.data = (1/np.sqrt(6))*(-data1.data*2 + data2.data + data3.data)
    rotatedN.data = (1/np.sqrt(6))*(np.sqrt(3)*data2.data - np.sqrt(3)*data3.data)
    rotatedZ.data = (1/np.sqrt(6))*(data1.data + data2.data + data3.data)
    st2 = Stream()
    st2=Stream(traces=[rotatedZ, rotatedN, rotatedE])
    return st2

client = Client("http://service.geonet.org.nz" , timeout=20)
#client = Client()
stime = UTCDateTime('2021-008T00:00:00')
etime = UTCDateTime('2021-010T00:00:00')
st = client.get_waveforms("NZ", "RACZ", "10",
                           "EHZ", stime, etime, attach_response = True)
st += client.get_waveforms("NZ", "HBAZ", "10",
                           "EHZ", stime, etime, attach_response = True)
st += client.get_waveforms("NZ", "ETAZ", "10",
                           "EH*", stime, etime, attach_response = True)
# st += client.get_waveforms("NZ", "WMSZ", "10",
#                            "EHZ", stime, etime, attach_response = False)
# st += client.get_waveforms("NZ", "WDSZ", "10",
#                            "EHZ", stime, etime, attach_response = False)


st.merge()
st.sort()
print(st)
per, nlnm = get_nlnm()
per, nhnm = get_nhnm()
NFFT=4096
fig = plt.figure(1,figsize=(10,12))
# plt.subplot(2,1,1)
for idx, tr in enumerate(st.select(channel="E**")):
    label= (tr.id).replace('.',' ')

    f, p = signal.welch(tr.data, fs = tr.stats.sampling_rate,
                     nperseg=NFFT, noverlap=256)
    amp, f= tr.stats.response.get_evalresp_response(tr.stats.delta, NFFT,
                                                    output='ACC')
    #(m/s^2)^2/Hz
    p /= np.abs(amp)**2
    p=10.*np.log10(p)
    plt.semilogx(f,p, label=label)
# plt.semilogx(1./per, nlnm, color='k', linewidth=2)
# plt.semilogx(1./per, nhnm, color='k', linewidth=2, label='NLNM/NHNM')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Amplitude (dB rel. 1 $(m/s^2)^2/Hz$)')
plt.xlim((1./10, 50.))
plt.vlines(4.5, -200, 250, colors='k', linestyles='dashed', label="4.5Hz")
plt.vlines(2, -200, 250, colors='c', linestyles='dashed', label="2Hz")
plt.legend(ncol=2, bbox_to_anchor=(0,1.02, 1,.12), loc='lower left', fontsize=10)
plt.ylim((-200, 250))

#Rotate to uvw
# st2 = rotate(st[1], st[2], st[0])

# plt.subplot(2,1,2)
# for idx, tr in enumerate(st2.select(channel="LH*")):
#     label='STS-2 ' + (tr.id).replace('.',' ')
#     label = label.replace('LHZ','LHU')
#     label = label.replace('LHN', 'LHV')
#     label = label.replace('LHE', 'LHW')
#     f, p = signal.welch(tr.data, fs = tr.stats.sampling_rate,
#                      nperseg=NFFT, noverlap=256)
#     amp, f= tr.stats.response.get_evalresp_response(tr.stats.delta, NFFT,
#                                                     output='ACC')
#     # (m/s^2)^2/Hz
#     p /= np.abs(amp)**2
#     p=10.*np.log10(p)
#     plt.semilogx(f,p, label=label)
# plt.semilogx(1./per, nlnm, color='k', linewidth=2)
# plt.semilogx(1./per, nhnm, color='k', linewidth=2, label='NLNM/NHNM')
# plt.xlabel('Frequency (Hz)')
# plt.ylabel('Amplitude (dB rel. 1 $(m/s^2)^2/Hz$)')
# plt.xlim((1./(1000.), 1./2.))

# plt.legend(ncol=2, loc=9, fontsize=18)
# plt.ylim((-200,-80))


plt.show()
# plt.savefig('Obsidian_compare_TOZ.png', format='PNG', dpi=400)