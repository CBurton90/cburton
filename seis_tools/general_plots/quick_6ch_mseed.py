#!/usr/bin/env python3

from obspy import UTCDateTime
from obspy import Stream
from obspy.clients.fdsn import Client

def get_and_remove_response(station, channel, location, output, t1, duration=60):
    client = Client("http://service-nrt.geonet.org.nz")
    st = client.get_waveforms(
        network="NZ", station=station, location=location,
        channel=channel, starttime=t1, endtime=t1 + duration)
    tr = Stream()
    for n in range(len(st)):
        st.merge()[n]
        inv = client.get_stations(
            network=st[n].stats.network, station=st[n].stats.station, 
            location=st[n].stats.location, channel=st[n].stats.channel, 
            level="response", startbefore=t1, endafter=t1 + duration)
        # pre_filt = (0.005, 0.006, 30.0, 35.0)
        st[n].remove_response(output=output, pre_filt=False, plot=False,
                       water_level=60, inventory=inv)
        tr += st[n]
    
    return tr

# pre_filt = [0.01,1,25,30]
tr = get_and_remove_response(station="RIZ", channel="HN*", location="20", output="ACC", t1=UTCDateTime(2021, 3, 4, 17, 41, 30))
# tr += get_and_remove_response(station="RIZ", channel="H**", location="10", output="ACC", t1=UTCDateTime(2021, 3, 4, 17, 41, 30))

tr.plot(equal_scale=False)
