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

client = FDSN_Client("GEONET")
stime = UTCDateTime("2020-11-10T02:22:00.000")
hhz = client.get_waveforms("NZ", "RPZ","10", "HHZ", stime, stime + 10)
hh1 = client.get_waveforms("NZ", "RPZ","10", "HH1", stime, stime + 10)
hh2 = client.get_waveforms("NZ", "RPZ","10", "HH2", stime, stime + 10)

print(hh2[0].stats.npts)

cal_sig = hh2[0].data * (40/(2**26))
hhz_v = hhz[0].data * (40/(2**26))
hh1_v = hh1[0].data * (40/(2**26))

Xc = np.linspace(0,(len(cal_sig)-1)/100,len(cal_sig))


# Define the sine function

def calc_sine(x,a,b,c,d):
    return a * np.sin(b* ( x + c)) + d

# An inital guess of amplitude, omega, phase, and offset for each calibration signal and each channel

hc = [cal_sig.max()/2.0, 2.*np.pi, 0.0, 0.0]
hz = [hhz_v.max()/2.0, 2.*np.pi, 0.0, 0.0]
h1 = [hh1_v.max()/2.0, 2.*np.pi, 0.0, 0.0]

# Find the optimal parameters and covariance

copt,ccov = curve_fit(calc_sine,Xc,cal_sig,p0=hc)
print(copt[0])
zopt,zcov = curve_fit(calc_sine,Xc,hhz_v,p0=hz)
print(zopt[0])
opt1,cov1 = curve_fit(calc_sine,Xc,hh1_v,p0=h1)
print(opt1[0])

# Acccuracy of our fitted function

Accuracy_cal =r2_score(cal_sig,calc_sine(Xc,*copt))
print (Accuracy_cal)

# Calculate the sensitivities of each component

# The equivalent input velocity is derived from the voltage amplitude of the calibration
# signal, resistor in  

rz =  (copt[0] * 2.0 * np.pi * 1.0 * 51000.0 * 0.02705 / zopt[0])
print(rz)

tyty = zopt[0] / (copt[0] / (2.0 * np.pi * 1.0 * 51000.0 * 0.02705))
print(tyty)

r1 = (copt[0] * 2.0 * np.pi * 1.0 * 51000.0 * 0.02872 / opt1[0])
print(r1)

calib_hhz = 1000 / (2 * np.pi * rz * 1.6777216)
print(calib_hhz)
calib_hh1 = 1000 / (2 * np.pi * r1 * 1.6777216)
print(calib_hh1)



# Plot the fit of our optimal parameters



line1 = Line2D(range(10), range(10), marker='o',color="darkblue")
line2 = Line2D(range(10), range(10), marker='o',color="firebrick")

amp_cal = round(copt[0], 3)
sens_z = round(rz, 2)
sens_1 = round(r1, 2)

line_label = "Amplitude: " + str(amp_cal) + "V"
z_label = "Sensitivity: " + str(sens_z) + "V/m/s" 
label_1 = "Sensitivity: " + str(sens_1) + "V/m/s"

fig, ax = plt.subplots(3,1,figsize=(20,10))
fig.suptitle('Calibration Sine Regression RPZ (10/11/2020)', fontsize=25, fontweight='bold')
ax[0].set_title(hh2[0].stats.station + "_" + hh2[0].stats.location + "_" + hh2[0].stats.channel + "(Calibration Signal)", fontsize=17)
ax[0].set_ylabel('Volts', fontsize=17)
ax[0].scatter(Xc,cal_sig,color="darkblue")
ax[0].plot(Xc,calc_sine(Xc,*copt),c="firebrick")
ax[0].legend((line1,line2),('Data','Sine Function'),numpoints=1, loc=2, bbox_to_anchor=(1.01,1))
ax[0].text(1.01, 0.6, line_label, transform=ax[0].transAxes)
ax[1].set_title(hhz[0].stats.station + "_" + hhz[0].stats.location + "_" + hhz[0].stats.channel, fontsize=17)
ax[1].set_ylabel('Volts', fontsize=17)
ax[1].scatter(Xc,hhz_v,color="darkblue")
ax[1].plot(Xc,calc_sine(Xc,*zopt),c="firebrick")
ax[1].text(1.01, 0.6, z_label, transform=ax[1].transAxes)
ax[2].set_title(hh1[0].stats.station + "_" + hh1[0].stats.location + "_" + hh1[0].stats.channel, fontsize=17)
ax[2].set_ylabel('Volts', fontsize=17)
ax[2].set_xlabel('Time (s)', fontsize=17)
ax[2].scatter(Xc,hh1_v,color="darkblue")
ax[2].plot(Xc,calc_sine(Xc,*opt1),c="firebrick")
ax[2].text(1.01, 0.6, label_1, transform=ax[2].transAxes)


# plt.show()
# plt.savefig('RPZ_sine_results.png', dpi=300)