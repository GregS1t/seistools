#!/usr/bin/env python
#
#
#	   @name	: rotation.py
#      @date    : 2022/05/11
#	   @author  : Greg Sainton (sainton@ipgp.fr)
#	   @purpose : module dedicated ease the rotation of the data
#
###############################################################################

import numpy as np
import numpy.ma as ma
from numpy.linalg import inv         # load function to invert matrices


from obspy import UTCDateTime, Stream

from modules.calibration import get_dip_azimuth

def check_4_rotation(chan_list,inv):
    """
    The following function is checking if the channels from the 
    input channels list is good to be used for rotation. 
    1> The number of channel in the list must be = 3
    3> The network, station dans location must be the same
    4> Channels must be ??U, ..V, ??W
    5> Starttime, Endtime, sampling and number of points must be the same
    ----
    INPUT:
        @chan_list: list of channels to be rotated
        @inv      : inventory object to get dip and azimuth for each channels
    ----
    OUTPUT:
        @rot_chan_list: list of rotated channels
        
    """
    if len(chan_list) == 3:
       # Check for validity of the data
        station_l = [cha.stats.station for cha in chan_list]
        network_l = [cha.stats.network for cha in chan_list]
        location_l = [cha.stats.location for cha in chan_list]
        channel_l = [cha.stats.channel for cha in chan_list]
        if len(set(network_l)) != 1 \
            or len(set(station_l)) !=1 or len(set(location_l)) !=1:
            print("Network, station or location are not the same for all channels.")
            return False
        else:
            axes = [orient[2] for orient in channel_l]
            axes.sort()
            if (axes == ['U', 'V', 'W']):
                stime_l = [cha.stats.starttime.timestamp for cha in chan_list]
                etime_l = [cha.stats.endtime.timestamp for cha in chan_list]
                sps_l = [cha.stats.sampling_rate for cha in chan_list]
                #print(stime_l, etime_l, sps_l)
                if len(set(stime_l)) != 1 \
                    or len(set(etime_l)) !=1 or len(set(sps_l)) !=1:
                    print("Invalid date range")
                    return False
                else:
                    return True
            else:
                print("It's not UVW channels")
                return False
    else:
        print("There is no exactly 3 channels...")
        return False
    
def get_unmasked_index_interval(trace, starttrim, endtrim):
    '''
    Function to find masked value at the beginning or at the end or 
    at the end of the data array.
    Other cases should be covered by "gap" cases
    ----
    INPUT:
        @trace: Trace object - input trace where to found masked data
        @starttrim: UTCDateTime to start trimming
        @endtrim: UTCDateTime to end trimming
    ----
    OUTPUT:
        @return: list with two values, start index and end index of 
        the unmasked data 
    
    '''
    temp_trace = trace.copy()
    temp_trace.trim(starttime=starttrim, endtime=endtrim, 
                    pad=False, nearest_sample=True, 
                    fill_value=None)
    
    
    #Look for masked value at the beginning
    for i in range(0, len(temp_trace.data)):
        if ma.is_masked(temp_trace.data[i])==False:
            break
    #Look for masked value at the end
    for j in range(1, len(temp_trace.data)):
        if ma.is_masked(temp_trace.data[len(temp_trace.data)-j])==False:
            break        
    return [i, len(temp_trace.data)-j], temp_trace

#==================================================================================
def get_unmasked_channels_interval(chan2trim, starttrim, endtrim):
    '''
    Function to return the common interval with unmasked data for 
    all of the channels
    in the channelTab list.
    ----
    INPUT:
        @chan2trim: list of Traces
        @starttrim: UTCDateTime to start trimming
        @endtrim: UTCDateTime to end trimming
    ----
    OUTPUT:
        @return: list with to UTCDateTime with the min and the max 
        common interval.

    '''
    try:
        import numpy as np
    except ImportError:
        pass
    
    # First part is to get common interval with unmasked values
    lower_bounds = []
    upper_bounds = []

    lower_bounds_time = []
    upper_bounds_time = []

    for i in range (0, len(chan2trim)):
        interv = []
        interv, new_chan = get_unmasked_index_interval(chan2trim[i],
                                                       starttrim, endtrim)
        lower_bounds.append(interv[0])
        upper_bounds.append(interv[1])
       
        #time vector
        t_temp = []
        
        starttime_raw = new_chan.stats.starttime
        npts_raw = new_chan.stats.npts
        sprate_raw = new_chan.stats.sampling_rate
        
        t_temp = np.arange(starttime_raw, starttime_raw + npts_raw/sprate_raw,
                           1/sprate_raw)
        t_tstamp = [t.timestamp for t in t_temp]
        
        #save the time vector for synchronization.
        
        lower_bounds_time.append(t_tstamp[interv[0]])
        upper_bounds_time.append(t_tstamp[interv[1]])
        
    #print(lower_bounds_time, upper_bounds_time)
    max_lower_bounds_time = max(lower_bounds_time)
    min_lower_bounds_time = min(upper_bounds_time)
    #print("Common interval in time without masked values: ", max_lower_bounds_time, min_lower_bounds_time)
    UTCxmin2trim = UTCDateTime(max_lower_bounds_time)
    UTCxmax2trim = UTCDateTime(min_lower_bounds_time)
    #print("In UTC: ", UTCxmin2trim, UTCxmax2trim)
    return [UTCxmin2trim, UTCxmax2trim]
        
          
def get_sync_channels_interval(channelTab, starttrim, endtrim):
    '''
    This function returns the common and synchronised interval 
    for a list of channels
    ----
    INPUT:
        @channelTab: list of Traces
        @starttrim: UTCDateTime to start trimming
        @endtrim: UTCDateTime to end trimming
    ----
    OUTPUT:
        @return: UTCDateTime interval of synchronized date of all channels
    '''
    try:
        import numpy as np
    except ImportError:
        pass
    
    
    UTCxmin2trim, UTCxmax2trim = get_unmasked_channels_interval(channelTab,
                                                                starttrim,
                                                                endtrim)
    #print("Unmasked range in UTC: ", UTCxmin2trim, UTCxmax2trim)
    
    #trim the data on this common interval without masked data
    trimmedTab = []
    t_prime = []
    
    for i in range (0, len(channelTab)):
        cp_channel = channelTab[i].copy()
        cp_channel.trim(starttime=UTCxmin2trim, endtime=UTCxmax2trim, 
                        pad=False, nearest_sample=True, fill_value=None)
        trimmedTab.append(cp_channel)
        #print(cp_channel.stats)
        starttime_raw = cp_channel.stats.starttime
        npts_raw = cp_channel.stats.npts
        sprate_raw = cp_channel.stats.sampling_rate
        
        t_temp = np.arange(starttime_raw, starttime_raw + npts_raw/sprate_raw,
                           1/sprate_raw)
        t_tstamp = [t.timestamp for t in t_temp]
        t_prime.append(t_tstamp)
    
    # look for the common synchronized interval
    
    temp_data = t_prime[0]
    #print(temp_data)
    np.set_printoptions(formatter={'float': '{: 0.3f}'.format})
    for i in range(0, len(t_prime)-1):
        #print(t_prime[i+1])
        temp_data = np.intersect1d(temp_data, t_prime[i+1])
        float_casted_data = temp_data.astype(float)
   
    if len(float_casted_data) >1:
        return [UTCDateTime(float_casted_data[0]),
                UTCDateTime(float_casted_data[-1])]
    else:
        return None

def sync_trim_channels(channelTab, starttrim, endtrim):
    """
    This function trim the channels in the synchronized mode.
    It works only if there is a common interval of time without
    any interpolation of the data. Otherwise, the function return None
    ----
    INPUT:
        @channelTab: list of Traces
        @starttrim: UTCDateTime to start trimming
        @endtrim: UTCDateTime to end trimming
    ----
    OUTPUT:
        @return: list of trimmed channels
    """
    if len(channelTab)>1:
        sync_interval = get_sync_channels_interval(channelTab, starttrim,
                                                   endtrim)
        if sync_interval is not None:
            trimmed_chan = []
            for i in range (0, len(channelTab)):
                cp_channel = channelTab[i].copy()
                cp_channel.trim(starttime=sync_interval[0],
                                endtime=sync_interval[1], \
                                pad=False, nearest_sample=True,
                                fill_value=None)
                trimmed_chan.append(cp_channel)
            return trimmed_chan
        else:
            print("Impossible to synchronize")
            return None
    else:
        trimmed_chan = []
        UTCxmin2trim, UTCxmax2trim = get_unmasked_channels_interval(channelTab,
                                                                    starttrim,
                                                                    endtrim)
        cp_channel = channelTab[0].copy()
        cp_channel.trim(starttime=UTCxmin2trim, endtime=UTCxmax2trim, \
                        pad=False, nearest_sample=True, fill_value=None)
        trimmed_chan.append(cp_channel)
        return trimmed_chan

def rotate_func_NEZ(cha_list, inventory):
    """
    Function which realize the rotation from UVW to ZNE.
    In this function, time series are supposed to be perfectly aligned
    with the same starttime, the same endtime and the same sampling rate.
    If the channels in cha_list are not aligned, one can try to use the function 
    'sync_trim_channels' which try to find the common part of thoses time series 
    between a given startime of trim and an endtime of trim
    
    ----
    INPUT
        @channel_list : list of 3 channels supposed to be U,V and W 
        @inv          : Inventory object generated by OBSPY from the dataless
    ----
    OUTPUT:
        @rotated_channel : list of 3 channels
    
    """

    channel2rotate = cha_list.copy()

    rot_param_l = []
    vec_data = []
    
    for cha in channel2rotate:
        rot_param = {}
        dip_azi = get_dip_azimuth(inventory,cha.stats.station, \
                cha.stats.location, \
                cha.stats.channel)
        rot_param["channel"] = cha.stats.channel
        rot_param["alpha"] = float(dip_azi.get("dip"))
        rot_param["beta"] = float(dip_azi.get("azimuth"))
        rot_param_l.append(rot_param)
        print(rot_param["channel"], rot_param["alpha"], rot_param["beta"])
        
        #vec_data.append(np.array(cha.data))
    # Fill the matrix to be rotated 
    #vec_data = [channel2rotate[0].data, channel2rotate[1].data, channel2rotate[2].data]
   
    
    for cha in channel2rotate:
        if cha.stats.channel[2] == "U":
            tr_N = cha.copy()
            tr_N.stats.channel = cha.stats.channel[:2]+"N"
            
        if cha.stats.channel[2] == "V":
            tr_E = cha.copy()
            tr_E.stats.channel = cha.stats.channel[:2]+"E"
                
        if cha.stats.channel[2] == "W":
            tr_Z = cha.copy()
            tr_Z.stats.channel = cha.stats.channel[:2]+"Z"
    
    
    vec_data = [tr_N.data,tr_E.data,tr_Z.data]
    
    #Define the rotation matrix XYZ to UVW
    A_NEZ_2_UVW = []
    # Fill the rotation matrix using input dip and azimuth from Inventory
    for i in range(0,len(rot_param_l)):
        A_NEZ_2_UVW.append([np.cos(np.deg2rad(rot_param_l[i].get("alpha")))*np.cos(np.deg2rad(rot_param_l[i].get("beta"))), \
                      np.cos(np.deg2rad(rot_param_l[i].get("alpha")))*np.sin(np.deg2rad(rot_param_l[i].get("beta"))), \
                      np.sin(np.deg2rad(rot_param_l[i].get("alpha")))])
    print("A_NEZ_2_UVW")
    print(A_NEZ_2_UVW)
    
    A_UVW_2_NEZ = inv(np.matrix(A_NEZ_2_UVW))
    print("A_UVW_2_NEZ")
    print(A_UVW_2_NEZ)
    rotate_vec_data = A_UVW_2_NEZ*vec_data
    
    
    tr_N.data = np.array(rotate_vec_data[0,:])[0,:]
    tr_E.data = np.array(rotate_vec_data[1,:])[0,:]
    tr_Z.data = np.array(rotate_vec_data[2,:])[0,:]
    
    tr_N._internal_add_processing_info("::ASPIC: Rotation from UVW to XYZ.")
    tr_E._internal_add_processing_info("::ASPIC: Rotation from UVW to XYZ.")
    tr_Z._internal_add_processing_info("::ASPIC: Rotation from UVW to XYZ.")
    
    rotated_stream = Stream([tr_N, tr_E, tr_Z])
    print(rotated_stream)
    rotated_stream.plot()
    return rotated_stream
