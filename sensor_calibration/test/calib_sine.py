#!/usr/bin/env python

from obspy import read
from obspy.core import UTCDateTime
from obspy.signal import freqattributes
from obspy.signal import filter
import matplotlib.pyplot as plt
from scipy.fft import fft, rfft
from scipy import fftpack
from scipy import signal
import numpy as np


def sens_ptp_3TB(cal_in,cal_out,Vpp,bits_in,bits_out,cal_const,ex_sens):
	cal_in_dt = signal.detrend(cal_in[0], type='linear')
	cal_counts_ptp = np.ptp(cal_in_dt, axis=0)
	cal_volts_ptp = cal_counts_ptp * (Vpp/(2**bits_in))
	equiv_vel_z = cal_volts_ptp / (2 * np.pi * 1 * 51000 * cal_const)
	equiv_vol_z = equiv_vel_z * ex_sens
	equiv_counts_z = equiv_vol_z / (Vpp/(2**bits_out))
	cal_out_dt = signal.detrend(cal_out[0], type='linear')
	chan_counts_ptp = np.ptp(cal_out_dt, axis=0)
	chan_volts_ptp = chan_counts_ptp * (Vpp/(2**bits_out))
	sens = chan_volts_ptp / equiv_vel_z
	ratio = chan_counts_ptp / equiv_counts_z
	return sens, ratio

def sens_rms_3TB(cal_in,cal_out,Vpp,bits_in,bits_out,cal_const,ex_sens):
	ref_dt = signal.detrend(cal_in[0], type='linear')
	ref_sqrd = np.square(ref_dt)
	cal_counts_rms = np.sqrt(np.mean(ref_sqrd, dtype='int64'))
	cal_volts_rms = cal_counts_rms * (Vpp/(2**bits_in)) 
	equiv_vel_z = cal_volts_rms / (2 * np.pi * 1 * 51000 * cal_const)
	equiv_vol_z = equiv_vel_z * ex_sens
	equiv_counts_z = equiv_vol_z / (Vpp/(2**bits_out))
	chan_dt = signal.detrend(cal_out[0], type='linear')
	chan_sqrd = np.square(chan_dt)
	chan_counts_rms = np.sqrt(np.mean(chan_sqrd, dtype='int64'))
	chan_volts_rms = chan_counts_rms * (Vpp/(2**bits_out))
	sens = chan_volts_rms / equiv_vel_z
	ratio = chan_counts_rms / equiv_counts_z
	return sens, ratio

	
start = UTCDateTime("2020-11-10T02:23:00.000")
end = UTCDateTime("2020-11-10T02:25:00.000")


ref_all = read("/home/conradb/calib/RPZ/RPZ_HH2_SINE_calin.mseed", format="MSEED")
hhz_all = read("/home/conradb/calib/RPZ/RPZ_HHZ_SINE_calout.mseed", format="MSEED")
hh1_all = read("/home/conradb/calib/RPZ/RPZ_HH1_SINE_calout.mseed", format="MSEED")



ref = ref_all.trim(starttime=start, endtime=end)
hhz = hhz_all.trim(starttime=start, endtime=end)
hh1 = hh1_all.trim(starttime=start, endtime=end)

# HHZ params

Vpp = 40
bits_in = 26 
bits_out = 26
cal_const = 0.02705
ex_sens = 2 * 4974

# HHZ via PTP method

hhz_ptp_sens, ratio_ptp_hhz = sens_ptp_3TB(ref,hhz,Vpp,bits_in,bits_out,cal_const,ex_sens)
print("(PTP) HHZ sensitivity in V/m/s is %s" % hhz_ptp_sens)
print("counts ratio (measured / expected) is %s" % ratio_ptp_hhz)

# HHZ via RMS method

hhz_rms_sens, ratio_rms_hhz = sens_rms_3TB(ref,hhz,Vpp,bits_in,bits_out,cal_const,ex_sens)
print("(RMS) HHZ sensitivity in V/m/s is %s" %hhz_rms_sens)
print("counts ratio (measured / expected) is %s" % ratio_rms_hhz)

# HHZ params

Vpp = 40
bits_in = 26 
bits_out = 26
cal_const = 0.02872
ex_sens = 2 * 4973

# HH1 via PTP method

hh1_ptp_sens, ratio_ptp_hh1 = sens_ptp_3TB(ref,hh1,Vpp,bits_in,bits_out,cal_const,ex_sens)
print("(PTP) HH1 sensitivity in V/m/s is %s" % hh1_ptp_sens)
print("counts ratio (measured / expected) is %s" % ratio_ptp_hh1)

# HH1 via RMS method

hh1_rms_sens, ratio_rms_hh1 = sens_rms_3TB(ref,hhz,Vpp,bits_in,bits_out,cal_const,ex_sens)
print("(RMS) HH1 sensitivity in V/m/s is %s" %hh1_rms_sens)
print("counts ratio (measured / expected) is %s" % ratio_rms_hh1)


## Guralp Instructions 

# hhz_volts_rms = hhz_counts_ptp * (40/(2**26))
# sens_z2 = hhz_volts_rms / equiv_vel_z
# print(sens_z2)

## Ncalib

calib_hhz = 1000 / (2 * np.pi * hhz_ptp_sens * 1.6777216)
print(calib_hhz)
diff = (calib_hhz / 0.078856) 
print(1-diff)

## Ncalib try 2

calib_hh1 = 1000 / (hh1_ptp_sens * 1.6777216 * 2 * np.pi)
print(calib_hh1)
diff2 = (calib_hh1 / 0.078856) 
print(1-diff2)

sens_z_nom = 1000 / (0.078856 * 1.6777216 * 2 * np.pi)
print(sens_z_nom)

# c =  hhz_counts_ptp / equiv_vel_z
# d = c / 1670000




