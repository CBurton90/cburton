#!/usr/bin/env python3

from __future__ import print_function
from obspy import UTCDateTime
from obspy import Stream
from obspy import Trace
from obspy import read, read_inventory
from obspy.clients.fdsn import Client
from obspy.signal.cross_correlation import correlate
# from multiprocessing import Pool, cpu_count, get_context, set_start_method
import multiprocessing
from multiprocessing import cpu_count
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
from timeit import default_timer as timer



def get_and_remove_response(station, channel, location, output, t1, duration):
    client = Client("http://service.geonet.org.nz")
    st = client.get_waveforms(
        network="NZ", station=station, location=location,
        channel=channel, starttime=t1, endtime=t1 + duration)
    tr = Stream()
    st.merge(fill_value='interpolate')
    print(st)
    for n in range(len(st)):
        
        inv = client.get_stations(
            network=st[n].stats.network, station=st[n].stats.station, 
            location=st[n].stats.location, channel=st[n].stats.channel, 
            level="response", startbefore=t1, endafter=t1 + duration)
        # st[n].detrend('linear')
        # st[n].taper(max_percentage=0.05, type='cosine')
        st[n].filter('bandpass', freqmin=0.1, freqmax=0.2, zerophase=True)
        st[n].remove_response(output=output, pre_filt=False, plot=False,
                       water_level=60, inventory=inv)
        tr += st[n]

    return tr

def equal_and_split(st, id, stacks):
	st = st.select(id=id)
	st = np.delete(st, 0, 1)
	rec_st = np.hsplit(st[0].data, stacks)

	return rec_st

def multi_run_wrapper(args):
    return get_and_remove_response(*args)


def get_windows(n, Mt, olap):
    # Split a signal of length n into olap% overlapping windows each containing Mt terms 
    ists = []
    ieds = []

    ist = 0
    while 1:
        ied = ist + Mt
        if ied > n:
            break
        ists.append(ist)
        ieds.append(ied)
        ist += int(Mt * (1 - olap/100))
    return ists, ieds

def check_length(stream):
    """
    Forces all traces to have same number of samples.

    Traces must be one day long.
    :type stream:`~obspy.core.stream.Stream` object. 
    :param stream: Stream containing one or more day-long trace 
    :return: Stream of similar length traces 

    """
    pts = 24*3600*stream[0].stats.sampling_rate
    npts = []
    for trace in stream:
        npts.append(trace.stats.npts)
    npts = np.array(npts)
    if len(npts) == 0:
        return stream	
    index = np.where(npts != pts)
    index = list(index[0])[::-1]

    # remove short traces
    for trace in index:
        stream.pop(trace)

    return stream




start = timer()

t1 = UTCDateTime(2018, 6, 4, 6, 30, 0)

duration = 4*24*60*60
time = 23.5*24*60*60

# duration = 60*60
# time = 3*60*60

stime = UTCDateTime("2018-06-04T06:45:00")

shift = 120

args = [('URZ', 'HH*', '10', 'VEL', t1, duration),
    ('KNZ', 'HH*', '10', 'VEL', t1, duration),
    ('URZ', 'HH*', '10', 'VEL', t1+4*24*60*60, duration),
    ('KNZ', 'HH*', '10', 'VEL', t1+4*24*60*60, duration),
    ('URZ', 'HH*', '10', 'VEL', t1+8*24*60*60, duration),
    ('KNZ', 'HH*', '10', 'VEL', t1+8*24*60*60, duration),
    ('URZ', 'HH*', '10', 'VEL', t1+12*24*60*60, duration),
    ('KNZ', 'HH*', '10', 'VEL', t1+12*24*60*60, duration),
    ('URZ', 'HH*', '10', 'VEL', t1+16*24*60*60, duration),
    ('KNZ', 'HH*', '10', 'VEL', t1+16*24*60*60, duration),
    ('URZ', 'HH*', '10', 'VEL', t1+20*24*60*60, duration),
    ('KNZ', 'HH*', '10', 'VEL', t1+20*24*60*60, duration)]





def main():
    # pool = multiprocessing.Pool(processes=2, maxtasksperchild=1)
    pool = multiprocessing.get_context('spawn').Pool(processes=12)
    # pool = multiprocessing.Pool(processes=4)
    # rst,sst,rst2,sst2,rst3,sst3,rst4,sst4,rst5,sst5,rst6,sst6,rst7,sst7,rst8,sst8,rst9,sst9,rst10,sst10 = pool.starmap(
    #     get_and_remove_response, args, chunksize=1)
    rst,sst,rst2,sst2,rst3,sst3,rst4,sst4,rst5,sst5,rst6,sst6 = pool.starmap(
        get_and_remove_response, args)
    pool.close()
    pool.join()


    rst += rst2 +rst3 + rst4 + rst5 + rst6 
    sst += sst2 +sst3 + sst4 + sst5 + sst6 


    rst.merge(fill_value='interpolate')
    sst.merge(fill_value='interpolate')


    rst = rst.trim(starttime=stime, endtime=stime + time)

    sst = sst.trim(starttime=stime, endtime=stime + time)   

    station = rst[0].stats.station
    source = sst[0].stats.station
    loc = rst[0].stats.location

    sst.sort()
    rst.sort()

    print(rst)
    print(sst)

    sou_st_e = equal_and_split(st=sst, id="NZ."+source+"."+loc+".HHE", stacks=time/1800)
    sou_st_n = equal_and_split(st=sst, id="NZ."+source+"."+loc+".HHN", stacks=time/1800)
    sou_st_z = equal_and_split(st=sst, id="NZ."+source+"."+loc+".HHZ", stacks=time/1800)

    rec_st_2 = equal_and_split(st=rst, id="NZ."+station+"."+loc+".HH2", stacks=time/1800)
    rec_st_1 = equal_and_split(st=rst, id="NZ."+station+"."+loc+".HH1", stacks=time/1800)
    rec_st_z = equal_and_split(st=rst, id="NZ."+station+"."+loc+".HHZ", stacks=time/1800)



    for n in range(len(sou_st_z)):

        signal.detrend(sou_st_z[n],type='linear')
        window = signal.tukey(len(sou_st_z[n]),alpha=0.05)
        sou_st_z[n] * window
        signal.detrend(sou_st_n[n],type='linear')
        window = signal.tukey(len(sou_st_n[n]),alpha=0.05)
        sou_st_n[n] * window
        signal.detrend(sou_st_e[n],type='linear')
        window = signal.tukey(len(sou_st_e[n]),alpha=0.05)
        sou_st_e[n] * window

    print("%s break")
    print(sou_st_z[0])

    for x in range(len(rec_st_z)):
        signal.detrend(rec_st_z[x],type='linear')
        window = signal.tukey(len(rec_st_z[x]),alpha=0.05)
        rec_st_z[x] * window
        signal.detrend(rec_st_1[x],type='linear')
        window = signal.tukey(len(rec_st_1[x]),alpha=0.05)
        rec_st_1[x] * window
        signal.detrend(rec_st_2[x],type='linear')
        window = signal.tukey(len(rec_st_2[x]),alpha=0.05)
        rec_st_2[x] * window

    print(rec_st_z[0])
    print(range(len(rec_st_z)))

    z_ccf_windows = []
    n_ccf_windows = []
    e_ccf_windows = []
    for i in range(len(rec_st_z)):
        z_corr = correlate(rec_st_z[i],sou_st_z[i],shift*100)
        n_corr = correlate(rec_st_1[i],sou_st_z[i],shift*100)
        e_corr = correlate(rec_st_2[i],sou_st_z[i],shift*100)
        print("%s zcorr break")
        print(z_corr)
        z_ccf_windows.append(z_corr)
        n_ccf_windows.append(n_corr)
        e_ccf_windows.append(e_corr)
        # for a in range(len(z_corr)):
        #     z_stacked = np.add(z_corr[a],z_corr[a+1]) 
        #     z_stacked / len(rec_st_z)
        #     # Czz.append(z_stacked)
    print("%s z stacked break")
    print(z_ccf_windows)



    zz = np.array(z_ccf_windows, dtype=np.float32)
    nz = np.array(n_ccf_windows, dtype=np.float32)
    ez = np.array(e_ccf_windows, dtype=np.float32)
    z_stacked = zz.sum(axis=0) / len(zz)
    n_stacked = nz.sum(axis=0) / len(nz)
    e_stacked = ez.sum(axis=0) / len(ez)


    print(len(zz))
    print(range(len(zz)))

    end = timer()
    print(end - start)

    # plt.plot(z_stacked)
    # plt.show()

    print(len(z_stacked))
    print(range(len(z_stacked)))

    # Convert to NumPy character array
    data = np.array(z_stacked, dtype='|S1')

    # Fill header attributes
    stats = {'network': 'NZ', 'station': station, 'location': loc,
             'channel': 'HHZ', 'npts': len(z_stacked), 'sampling_rate': 100,
             'mseed': {'dataquality': 'D'}}
    stats2 = {'network': 'NZ', 'station': station, 'location': loc,
             'channel': 'HH1', 'npts': len(z_stacked), 'sampling_rate': 100,
             'mseed': {'dataquality': 'D'}}
    stats3 = {'network': 'NZ', 'station': station, 'location': loc,
             'channel': 'HH2', 'npts': len(z_stacked), 'sampling_rate': 100,
             'mseed': {'dataquality': 'D'}}
    # set current time
    stats['starttime'] = stime
    stats2['starttime'] = stime
    stats3['starttime'] = stime
    st = Stream([Trace(data=z_stacked, header=stats)])
    st2 = Stream([Trace(data=n_stacked, header=stats2)])
    st3 = Stream([Trace(data=e_stacked, header=stats3)])
    # write as ASCII file (encoding=0)
    st.write(station+"_"+source+"_ZZ.mseed", format='MSEED', encoding="FLOAT32", reclen=512)
    st2.write(station+"_"+source+"_1Z.mseed", format='MSEED', encoding="FLOAT32", reclen=512)
    st3.write(station+"_"+source+"_2Z.mseed", format='MSEED', encoding="FLOAT32", reclen=512)

if __name__ == "__main__":

    multiprocessing.set_start_method("spawn")

    main()