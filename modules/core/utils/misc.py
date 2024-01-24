#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 11 10:06:20 2019

@author: Greg
@email: sainton@ipgp.fr

@purpose: Place to put useful functions used multiple times in other scripts

- mkdir_p(pathdata, exportdir)
- create_tree_struct(outputdir= '', sol=None)

"""
import os
from obspy import UTCDateTime

def mkdir_p(pathdata, newdir):
    '''
    Function to create a directory if necessary
    This function is used to save exported files
    '''
    import os
    import errno
    fullpath = pathdata+"/"+newdir
    try:
        os.makedirs(fullpath)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(fullpath):
            # print("In fact, this directory already exists...")
            pass
        else:
            raise


# Define tree structure of directories
############################################
def create_tree_struct(outputdir='', sol=None):
    '''
    Function to create a basic folder structure according to a root dirtectory and
    a sol number
    @outputdir: root directory for the analysis.
    @sol: sol number

    '''
    if sol is not None:
        directory = "".join(["SOL",'{:04d}'.format(sol)])
    else:
        directory = ""

    if outputdir !='' and os.path.isdir(outputdir)==True:
            mkdir_p(outputdir, directory)
    else:
        here = os.getcwd()
        mkdir_p(here,outputdir+"/"+directory)

    savedir = outputdir+"/"+directory

    # Directory to save the raw data
    if savedir !='' and os.path.isdir(savedir)==True:
        mkdir_p(savedir, "Raw")
        raw_dir = "/".join([savedir,"Raw"])
    else:
        mkdir_p(outputdir, directory+"/Raw")
        raw_dir = "/".join([outputdir, directory, "Raw"])

    return directory, savedir, raw_dir

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
        @return: list with two values, start index and end index of the unmasked data 

    '''
    import numpy.ma as ma
    temp_trace = trace.copy()
    temp_trace.trim(starttime=starttrim, endtime=endtrim,
                    pad=False, nearest_sample=True,
                    fill_value=None)

    # Look for masked value at the beginning
    for i in range(0, len(temp_trace.data)):
        if not ma.is_masked(temp_trace.data[i]):
            break
    # Look for masked value at the end
    for j in range(1, len(temp_trace.data)):
        if not ma.is_masked(temp_trace.data[len(temp_trace.data)-j]):
            break
    return [i, len(temp_trace.data)-j], temp_trace

# =============================================================================

def get_unmasked_channels_interval(chan2trim, starttrim, endtrim):
    '''
    Function to return the common interval with unmasked data for all
    of the channels in the channelTab list.
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

        # time vector
        t_temp = []

        starttime_raw = new_chan.stats.starttime
        npts_raw = new_chan.stats.npts
        sprate_raw = new_chan.stats.sampling_rate
        
        t_temp = np.arange(starttime_raw,
                           starttime_raw + npts_raw/sprate_raw,
                           1/sprate_raw)

        t_tstamp = [t.timestamp for t in t_temp]

        # save the time vector for synchronization.

        lower_bounds_time.append(t_tstamp[interv[0]])
        upper_bounds_time.append(t_tstamp[interv[1]])

    # print(lower_bounds_time, upper_bounds_time)
    max_lower_bounds_time = max(lower_bounds_time)
    min_lower_bounds_time = min(upper_bounds_time)
    # print("Common interval in time without masked values: ", max_lower_bounds_time, min_lower_bounds_time)
    UTCxmin2trim = UTCDateTime(max_lower_bounds_time)
    UTCxmax2trim = UTCDateTime(min_lower_bounds_time)
    #print("In UTC: ", UTCxmin2trim, UTCxmax2trim)
    return [UTCxmin2trim, UTCxmax2trim]

def get_sync_channels_interval(channelTab, starttrim, endtrim):
    '''
    This function returns the common and synchronised interval for a list of channels
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

    UTCxmin2trim, UTCxmax2trim = get_unmasked_channels_interval(channelTab, starttrim, endtrim)
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
        
        t_temp = np.arange(starttime_raw, starttime_raw + npts_raw/sprate_raw, 1/sprate_raw)
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
        return [UTCDateTime(float_casted_data[0]), UTCDateTime(float_casted_data[-1])]
    else:
        return None


def sync_trim_channels(channelTab, starttrim, endtrim):
    """
    This function trim the channels in the synchronized mode.
    It works only if there is a common interval of time without any interpolation
    of the data. Otherwise, the function return None
    ----
    INPUT:
        @channelTab: list of Traces
        @starttrim: UTCDateTime to start trimming
        @endtrim: UTCDateTime to end trimming
    ----
    OUTPUT:
        @return: list of trimmed channels
    """
    try:
        import numpy as np
    except ImportError:
        pass
    
    if len(channelTab)>1:
        sync_interval = get_sync_channels_interval(channelTab, starttrim, endtrim)
        if sync_interval is not None:
            trimmed_chan = []
            for i in range (0, len(channelTab)):
                cp_channel = channelTab[i].copy()
                cp_channel.trim(starttime=sync_interval[0], endtime=sync_interval[1], \
                                pad=False, nearest_sample=True, fill_value=None)    
                trimmed_chan.append(cp_channel)
            return trimmed_chan
        else:
            #print("Impossible to synchronize")
            return None
    else:
        trimmed_chan = []
        UTCxmin2trim, UTCxmax2trim = get_unmasked_channels_interval(channelTab, starttrim, endtrim)
        cp_channel = channelTab[0].copy()
        cp_channel.trim(starttime=UTCxmin2trim, endtime=UTCxmax2trim, \
                        pad=False, nearest_sample=True, fill_value=None)
        trimmed_chan.append(cp_channel)
        return trimmed_chan


def get_lg_from_id(trace, df_properties):
    """
    Retrieve the longname fron the ID of the channel
    which mean the information given in the stats
    This function is linked to a config file saved in
    the configuration file
    ---
    INPUT:
        trace : obpsy Trace object
    ---
    OUTPUT:
        longname: str
    """
    try:
        locaidchan = trace.stats.location+"_"+trace.stats.channel
        #print(locaidchan)
        return df_properties.loc[(df_properties["locidchan"] == locaidchan)
                        & (df_properties["Network"] == trace.stats.network),
                        "longnamesps"].iloc[0]
    except:
        return None


def trace_pprint(trace):
    """
    Function to print the detail of the traces with more details than the
    usual print(traces)
    """

    net = trace.stats.network
    sta = trace.stats.station
    loc = trace.stats.location
    cha = trace.stats.channel
    lg = get_lg_from_id(trace)
    sps = trace.stats.sampling_rate
    npts = trace.stats.npts
    startt = trace.stats.starttime
    endt = trace.stats.endtime
    print(lg+" | "+net+"."+sta+"."+loc+"."+cha+" | "\
          +str(startt)+" | "+str(endt)+" | "\
          +str(sps)+" Hz | "+ str(npts)+" samp")
