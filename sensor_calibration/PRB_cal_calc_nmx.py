import matplotlib


#from mtspec import *
from pylab import *
from numpy.fft import rfft

#from obspy.imaging import spectrogram
from obspy import read
from obspy.core import UTCDateTime, Stream
#from obspy.signal.spectral_estimation import psd
from obspy.signal.invsim import paz_to_freq_resp
from scipy import signal
#from obspy.signal import pazToFreqResp
#from obspy.signal.psd import welch_taper
#from obspy.signal.psd import welch_window

def cmg3t(sps, nfft, scale=1.0, norm = None, error = 5.0):

    poles = [(-0.00589 + 0.00589j) * (2 * pi), (-0.00589 - 0.00589j) * (2 * pi), (-73.2 + 37.6j) * ( 2 * pi), (-73.2 - 37.6j) * (2 * pi)]
    zeros = [0j, 0j, (146.5 + 0j) * (2 * pi)]

    (R, f) = paz_to_freq_resp(poles, zeros, scale, 1.0 / sps, nfft, freq=True);

    if norm is not None:
        offset = int(norm * float(nfft) / sps)
        R = R / abs(R[offset])

    Rmin = R * (100.0 - error) / 100.0
    Rmax = R * (100.0 + error) / 100.0

    return R, Rmin, Rmax, f

def t120qa(sps, nfft, scale=1.0, norm = None, error = 5.0):

    poles = [(-0.036614 + 0.037059j), (-0.036614 - 0.037059j), (-32.55 + 0j), (-142 + 0j), (-364 + 404j), (-364-404j), (-1260 + 0j), (-4900 + 5200j), (-4900 - 5200j), (-7100 + 1700j), (-7100 - 1700j)]
    zeros = [0j, 0j, (-31.63 + 0j), (-160 + 0j), (-350 + 0j), (-3177 + 0j)]

    (R, f) = paz_to_freq_resp(poles, zeros, scale, 1.0 / sps, nfft, freq=True);

    if norm is not None:
        offset = int(norm * float(nfft) / sps)
        R = R / abs(R[offset])

    Rmin = R * (100.0 - error) / 100.0
    Rmax = R * (100.0 + error) / 100.0

    return R, Rmin, Rmax, f


def ms_psd(infile, nfft, norm = None):
    print(infile)
    st = infile
    print(st)

    #Pxx, f = psd(st[0].data.astype(np.float64), NFFT = nfft, Fs = st[0].stats.sampling_rate, window = window_hanning, detrend = detrend_mean, noverlap=0)
    f, Pxx = signal.welch(st[0].data.astype(np.float64), fs = st[0].stats.sampling_rate, nperseg=nfft, noverlap=256)

#    if norm is not None:
#        offset[0] = int(norm[0] * float(nfft) / st[0].stats.sampling_rate)
#        offset[1] = int(norm[1] * float(nfft) / st[0].stats.sampling_rate)

    return Pxx, f, st[0].stats.sampling_rate

def ms_spec(infile, reffile, nfft, norm = None):
    Pxx, f, sps = ms_psd(infile, nfft)
    Rxx, f, sps = ms_psd(reffile, nfft)

    Syy = []
    Szz = []
    for i in range(0, len(f)):
        Syy.append(abs(Pxx[i]) * pow(2.0 * pi * f[i], 2.0));
#        Syy.append(abs(Pxx[i]) * pow(2.0 * pi * f[i], 2.0));
#        Syy.append(abs(Pxx[i]) * pow(2.0 * pi * f[i], 2.0));

        Szz.append(sqrt((Syy[i] / Rxx[i])))

        #Szz.append(abs((Pxx[i] * pow(2.0 * pi * f[i], 2.0) / Rxx[i])))
        #Szz.append(sqrt(abs((Pxx[i] * pow(2.0 * pi * f[i], 2.0)/ Rxx[i]))))
        #Szz.append(abs(Syy[i] / Rxx[i]))
#
    if norm is not None:
        offset = int(norm * float(nfft) / sps)
        Szz = Szz / abs(Szz[offset])

    #return Szz, f
    return Szz, f
    #return Rxx, f

sps = 100.0
#nfft = 32768
nfft = 32768 * 4
#sps = [100.0, 100.0]
##nfft = [16384, 32768]
#nfft = [16384, 32768]


starttime = UTCDateTime("2021-11-08T23:08:00")
endtime = UTCDateTime("2021-11-09T02:08:00")

# read in the calibration outputs for the 3ch's and the sine calibration signal

path = '/home/conradb/git/cburton/sensor_calibration/2021_Calibs/'

st = read(path+'RPZ_centaur-6_7981_20211108_230752.seed')
cal = read(path+'20211108_RPZ_PRB_CAL.seed')
st.trim(starttime=starttime, endtime=endtime)
cal.trim(starttime=starttime, endtime=endtime)
print(st)
print(cal)

ref = cal
hhz = st.select(component="Z")
hh1 = st.select(component="1")
hh2 = st.select(component="2")

fig = plt.figure(figsize=(10, 5))
ax = fig.add_subplot(111)
ax.set_title('RPZ Nominal Response Comparison with 5% Bounds Shown')
ax.grid(True)
ax.set_ylabel('Spectral Ratio')
ax.set_xlabel('Hz')
ax.set_ylim([0.1,2.0])
ax.set_xlim([1.0 / 400.0,1.0])

Rxx_hhz, f = ms_spec(hhz, ref, nfft, norm = 1.0)
Rxx_hh1, f = ms_spec(hh1, ref, nfft, norm = 1.0)
Rxx_hh2, f = ms_spec(hh2, ref, nfft, norm = 1.0)
#ax.loglog(f[1:-1], Rxx_hhz[1:-1], 'r', f[1:-1], Rxx_hh1[1:-1], 'g', f[1:-1], Rxx_hh2[1:-1], 'b')
ax.plot(f[1:-1], Rxx_hh2[1:-1], 'b', label='U')
ax.plot(f[1:-1], Rxx_hh1[1:-1], 'g', label='V')
ax.plot(f[1:-1], Rxx_hhz[1:-1], 'r', label='W')


ax.legend()


R, Rmin, Rmax, f = t120qa(sps, nfft, norm = 1.0)
print(R, Rmin, Rmax, f)

ax.loglog(f[1:-1],abs(R[1:-1]), color="black")
l, = ax.loglog(f[1:-1],abs(Rmin[1:-1]), '--', color="black")
l.set_dashes([1,1])
l, = ax.loglog(f[1:-1],abs(Rmax[1:-1]), '--', color="black")
l.set_dashes([1,1])

##leg = plt.gca().get_legend()
##ltext  = leg.get_texts()
##llines = leg.get_lines()
##plt.setp(ltext, fontsize='small')
##plt.setp(llines, linewidth=1.5)
#
##ax = fig.add_subplot(212)
##ax.grid(True)
##ax.set_ylabel('Spectral Ratio')
##ax.set_xlabel('Hz')
###ax.set_ylim([0.5,2])
#
#Rxx_lhz, f = ms_spec(hhz[1], ref[1], nfft[1], norm = 0.1)
#Rxx_lh1, f = ms_spec(hh1[1], ref[1], nfft[1], norm = 0.1)
#Rxx_lh2, f = ms_spec(hh2[1], ref[1], nfft[1], norm = 0.1)
#ax.loglog(f[1:-1], Rxx_lhz[1:-1], 'r', f[1:-1], Rxx_lh1[1:-1], 'g', f[1:-1], Rxx_lh2[1:-1], 'b')
#plt.legend(('HHZ', 'HH1', 'HH2'), 'lower right', shadow=False, fancybox=False)
#
#R, Rmin, Rmax, f = cmg3t(sps[1], nfft[1], norm = 0.1)
##
#ax.loglog(f[1:-1],abs(R[1:-1]), color="black")
#l, = ax.loglog(f[1:-1],abs(Rmin[1:-1]), '--', color="black")
#l.set_dashes([1,1])
#l, = ax.loglog(f[1:-1],abs(Rmax[1:-1]), '--', color="black")
#l.set_dashes([1,1])
#
##leg = plt.gca().get_legend()
##ltext  = leg.get_texts()
##llines = leg.get_lines()
##plt.setp(ltext, fontsize='small')
##plt.setp(llines, linewidth=1.5)

#plt.show()

# savefig('rpz-results.pdf', format='pdf')
plt.savefig('RPZ_PRB_cal_20211108.png', format='PNG', dpi=400)