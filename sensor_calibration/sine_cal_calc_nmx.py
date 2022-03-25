#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score
from obspy import UTCDateTime
from obspy.clients.fdsn import Client as FDSN_Client
from obspy import read

# FDSN

# client = FDSN_Client("GEONET")
# stime = UTCDateTime("2020-11-10T02:22:00.000")
# hhz = client.get_waveforms("NZ", "RPZ","10", "HHZ", stime, stime + 10)
# hh1 = client.get_waveforms("NZ", "RPZ","10", "HH1", stime, stime + 10)
# hh2 = client.get_waveforms("NZ", "RPZ","10", "HH2", stime, stime + 10)

# print(hh2[0].stats.npts)

# cal_sig = hh2[0].data * (40/(2**26))
# hhz_v = hhz[0].data * (40/(2**26))
# hh1_v = hh1[0].data * (40/(2**26))

# Local

freq = 1
Vpp = 40
bits= 24

sensor_sensitivty = 1202.5 # V/m/s
cal_coil_const = -0.01023 # (m/s^2)/V
dig_sens_V = 1 / (Vpp/(2**bits)) # counts/V
#dig_pre_amp = dig_sens_V / 1e6 # counts/uV
dig_pre_amp = 0.4
dig_gain = 1
dig_sens_uV = dig_pre_amp * dig_gain #count/uV
nom_sys_sens = (sensor_sensitivty * dig_sens_uV) / 1e3 # count/nm/s
NCALIB = 1 / (nom_sys_sens* 2 * np.pi * freq)
print(NCALIB)

y_comp = '1'
x_comp = '2'
 


print(dig_sens_uV)



# start and end of the sine calibration trimmed to 1 minute (out of 5)

starttime = UTCDateTime("2021-11-09T23:20:30")
endtime = UTCDateTime("2021-11-09T23:21:30")

# read in the calibration outputs for the 3ch's and the sine calibration signal

path = '/home/conradb/git/cburton/sensor_calibration/2021_Calibs/'

st = read(path+'RPZ_centaur-6_7981_20211109_231757.seed')
cal = read(str(path)+'20211109_RPZ_sine_CAL.seed')
st.trim(starttime=starttime, endtime=endtime)
cal.trim(starttime=starttime, endtime=endtime)
print(st)
print(cal)

# select the components in individual streams

cal_sig = cal.select(component="A")[0].data
hhz = st.select(component="Z")[0].data
hh1 = st.select(component=y_comp)[0].data
hh2 = st.select(component=x_comp)[0].data

# print(cal_sig.max())
# print(hhz.max())
# print(hh1.max())
# print(hh2.max())

Xc = np.linspace(0,(len(cal_sig)-1)/100,len(cal_sig))

# Define the sine function

def calc_sine(x,a,b,c,d):
    return a * np.sin(b* ( x + c)) + d

# An inital guess of amplitude, omega, phase, and offset for each calibration signal and each channel

hc = [cal_sig.max()/2.0, 2.*np.pi, 0.0, 0.0]
hz = [hhz.max()/2.0, 2.*np.pi, 0.0, 0.0]
h1 = [hh1.max()/2.0, 2.*np.pi, 0.0, 0.0]
h2 = [hh2.max()/2.0, 2.*np.pi, 0.0, 0.0]

# Find the optimal parameters and covariance

copt,ccov = curve_fit(calc_sine,Xc,cal_sig,p0=hc)
print('cal peak amp in counts is ', copt[0])
opt2,cov2 = curve_fit(calc_sine,Xc,hh2,p0=h2)
opt1,cov1 = curve_fit(calc_sine,Xc,hh1,p0=h1)
zopt,zcov = curve_fit(calc_sine,Xc,hhz,p0=hz)
print('UVW peak amp in counts is ', opt2[0], opt1[0], zopt[0])
U_rms = opt2[0] / np.sqrt(2)
V_rms = opt1[0] / np.sqrt(2)
W_rms = zopt[0] / np.sqrt(2)
print('UVW amp rms in counts is ', U_rms, V_rms, W_rms)

# convert the calibrtion peak amplitude to volts

#cal_amp = copt[0] * (Vpp/(2**bits))
cal_amp = 1
print('cal peak amp in volts is', cal_amp)

# calibration peak amplitude and rms in micro m per s (um/s)

cal_ums = abs((cal_amp / (2 * np.pi * freq * (1/cal_coil_const)))) * 1e6
cal_rms = cal_ums / np.sqrt(2)

# measured system sensitivty in counts/nm/s

U_meas = U_rms / (cal_rms * 1000)
V_meas = V_rms / (cal_rms * 1000)
W_meas = W_rms / (cal_rms * 1000)

# UVW calibration factor (CALIB) in nm/count

U_CALIB = 1 / (U_meas * 2 * np.pi * freq)
V_CALIB = 1 / (V_meas * 2 * np.pi * freq)
W_CALIB = 1 / (W_meas * 2 * np.pi * freq)
print('UVW CALIB', U_CALIB, V_CALIB, W_CALIB)

per_nom_U = (U_CALIB - NCALIB) / NCALIB *100
per_nom_V = (V_CALIB - NCALIB) / NCALIB *100
per_nom_W = (W_CALIB - NCALIB) / NCALIB *100
print('UVW percent from Nominal', per_nom_U, per_nom_V, per_nom_W)

# Accuracy of our fitted function

Accuracy_cal =r2_score(cal_sig,calc_sine(Xc,*copt))
print ("r2_score", Accuracy_cal)

# XYZ sys sens in counts/nm/s

X_meas = np.sqrt(1/6) * ((2 * (np.sqrt(1/6)*2)*U_meas) + (-1 * (np.sqrt(1/6)*-1)*V_meas) + (-1 * (np.sqrt(1/6)*-1)*W_meas))
Y_meas = np.sqrt(1/6) * ((0 * (np.sqrt(1/6)*0)*U_meas) + (np.sqrt(3) * (np.sqrt(1/6)*np.sqrt(3))*V_meas) + (-1*np.sqrt(3) * (np.sqrt(1/6)*-1*np.sqrt(3))*W_meas))
Z_meas = np.sqrt(1/6) * ((np.sqrt(2) * (np.sqrt(1/6)*np.sqrt(2))*U_meas) + (np.sqrt(2) * (np.sqrt(1/6)*np.sqrt(2))*V_meas) + (np.sqrt(2) * (np.sqrt(1/6)*np.sqrt(2))*W_meas))

# XYZ calibration factor (CALIB) in nm/count

X_CALIB = 1 / (X_meas * 2 * np.pi * freq) #HHE/HH2
Y_CALIB = 1 / (Y_meas * 2 * np.pi * freq) #HHN/HH1
Z_CALIB = 1 / (Z_meas * 2 * np.pi * freq) #HHZ
print('XYZ CALIB', X_CALIB, Y_CALIB, Z_CALIB)

per_nom_X = (X_CALIB - NCALIB) / NCALIB *100
per_nom_Y = (Y_CALIB - NCALIB) / NCALIB *100
per_nom_Z = (Z_CALIB - NCALIB) / NCALIB *100
print('XYZ percent from Nominal', per_nom_X, per_nom_Y, per_nom_Z)

senssssss = Z_meas * 1000 / dig_sens_uV
print(senssssss)

# Plot the fit of our optimal parameters

line1 = Line2D(range(10), range(10), marker='o',color="darkblue")
line2 = Line2D(range(10), range(10), marker='o',color="firebrick")

amp_cal = round(cal_amp, 3)
sens_z = round(Z_meas, 3)
sens_1 = round(Y_meas, 3)
sens_2 = round(X_meas, 3)

line_label = "Amplitude: " + str(amp_cal) + "V"
z_label = "Sensitivity: " + str(sens_z) + "count/nm/s" 
label_1 = "Sensitivity: " + str(sens_1) + "count/nm/s"
label_2 = "Sensitivity: " + str(sens_2) + "count/nm/s"

fig, ax = plt.subplots(4,1,figsize=(20,10))
fig.suptitle('Calibration Sine Regression TEST (10/11/2020)', fontsize=25, fontweight='bold')

# plot calibration signal

ax[0].set_title(cal[0].stats.station + "_" + cal[0].stats.location + "_" + cal[0].stats.channel + "(Calibration Signal)", fontsize=17)
ax[0].set_ylabel('Volts', fontsize=17)
ax[0].scatter(Xc,cal[0].data,color="darkblue")
ax[0].plot(Xc,calc_sine(Xc,*copt),c="firebrick")
ax[0].legend((line1,line2),('Data','Sine Function'),numpoints=1, loc=2, bbox_to_anchor=(1.01,1))
ax[0].text(1.01, 0.6, line_label, transform=ax[0].transAxes)

# plot vertical comp

ax[1].set_title(st.select(component="Z")[0].stats.station + "_" + st.select(component="Z")[0].stats.location + "_" + st.select(component="Z")[0].stats.channel, fontsize=17)
ax[1].set_ylabel('Volts', fontsize=17)
ax[1].scatter(Xc,hhz,color="darkblue")
ax[1].plot(Xc,calc_sine(Xc,*zopt),c="firebrick")
ax[1].text(1.01, 0.6, z_label, transform=ax[1].transAxes)

# plot hh1/hhn comp

ax[2].set_title(st.select(component=y_comp)[0].stats.station + "_" + st.select(component=y_comp)[0].stats.location + "_" + st.select(component=y_comp)[0].stats.channel, fontsize=17)
ax[2].set_ylabel('Volts', fontsize=17)
ax[2].set_xlabel('Time (s)', fontsize=17)
ax[2].scatter(Xc,hh1,color="darkblue")
ax[2].plot(Xc,calc_sine(Xc,*opt1),c="firebrick")
ax[2].text(1.01, 0.6, label_1, transform=ax[2].transAxes)

# plot hh2/hhe comp

ax[3].set_title(st.select(component=x_comp)[0].stats.station + "_" + st.select(component=x_comp)[0].stats.location + "_" + st.select(component=x_comp)[0].stats.channel, fontsize=17)
ax[3].set_ylabel('Volts', fontsize=17)
ax[3].set_xlabel('Time (s)', fontsize=17)
ax[3].scatter(Xc,hh2,color="darkblue")
ax[3].plot(Xc,calc_sine(Xc,*opt2),c="firebrick")
ax[3].text(1.01, 0.6, label_2, transform=ax[3].transAxes)



plt.show()
# plt.savefig('RPZ_sine_results.png', dpi=300)