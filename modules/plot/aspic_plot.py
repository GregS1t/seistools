#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 15 16:01:50 2019

@author: Grégory Sainton
@email:  sainton@ipgp.fr



"""


import sys
import os
import numpy as np
from itertools import cycle
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import pandas as pd
from obspy import UTCDateTime

sys.path.insert(0, 'modules/time/')
sys.path.insert(0, 'modules/core/utils/')
from time_utils import utc2lmst, utc2sol, lmst2utc, sol_span_in_utc,\
                            parse_location_channel_codes, dayplot_set_x_ticks

CHANNEL_PROPERTIES_FILE = "./configuration/channel.properties2.0.csv"
if os.path.exists(CHANNEL_PROPERTIES_FILE):
    pass
else:
    print("Config file not found")

DF_PROPERTIES  = pd.read_csv(CHANNEL_PROPERTIES_FILE, delimiter=";")

class plot_options_process():
    """
    Class to handle parameters of plot.
    There is one plot object per processing

    """
    share_x_axis = True
    share_y_axis = False
    utc_or_lmst = "LMST"
    show_gaps = True
    show_overlaps = False
    display = True
    sol_scale = True
    forced_starttime = None
    forced_endtime = None
    save = True
    ptype = "subplot"

    def __init__(self, suptitle="", xtitle="", ytitle="",
                 display=True, save=False):
        self.suptitle = suptitle
        self.xtitle = xtitle
        self.ytitle = ytitle
        self.display = display
        self.save = save

    def __str__(self):
        """
        Pretty printable components of the plots
        """
        param_plot = "Suptitle: " + self.suptitle + " / Plot type: "+self.ptype+"\n"
        param_plot += "Xtitle: " + self.xtitle + " / Ytitle: " + self.ytitle+" / Units: "+str(self.utc_or_lmst)+"\n"
        param_plot += "Show gaps: "+ str(self.show_gaps)+ " / Show overlaps: "+ str(self.show_overlaps)+"\n"
        param_plot += "Display plots: "+ str(self.display)+" / Save plots: "+str(self.save)

        return param_plot.format(self=self)

    def copy(self):
        temp = plot_options_process()
        temp.share_x_axis = self.share_x_axis
        temp.share_y_axis = self.share_y_axis
        temp.utc_or_lmst = self.utc_or_lmst
        temp.show_gaps = self.show_gaps
        temp.show_overlaps = self.show_overlaps
        temp.display = self.display
        temp.sol_scale = self.sol_scale
        temp.forced_starttime = self.forced_starttime
        temp.forced_endtime = self.forced_endtime
        temp.save = self.save
        temp.ptype = self.ptype
        temp.suptitle = self.suptitle
        temp.xtitle = self.xtitle
        temp.ytitle = self.ytitle
        return temp



def get_traces_id(stream):
        """
        Plain function to create an ID for each channels using concatenation 
        of the Network, Station, Location and Channel codes
        ie : XB_ELYSE_02_BHU
        
        """
        list_id = []
        for trace in stream:
                trace.stats.id = "{}_{}_{}_{}".format(trace.stats.network, trace.stats.station, 
                                                        trace.stats.location,trace.stats.channel)
                list_id.append(trace.stats.id)
                #print(sorted(list_id))
        set_id = set(sorted(list_id, key=str))
        return set_id

def plot_by_sols_slice(trace, mDate, yshift=0, title=None, xlabel=None, ylabel="Time (LMST)"):
    """
    Function to plot time series after slicing by sols
    @trace: ASPIC_Channel object
    @mDate: MarsConverter object to convert from UTC to LMST
    @yshift: number to shift the different time series for a better rendering
    @title: string = title of the plots
    @xlabel: string = xaxis label
    @ylabel: string = yaxis label

    """

    from matplotlib.ticker import AutoMinorLocator

    starttime_sol = str(mDate.get_utc_2_lmst(utc_date=trace.starttime))
    posT_min = starttime_sol.find('T')
    # Here are the sols.
    sol_min = int(starttime_sol[0:posT_min])

    # Slice the data
    sliced_channels = trace.slice_in_sols(mDate)
    fig, axs = plt.subplots(figsize=(10,6))
    import math
    for i in range(0,len(sliced_channels)):
        print(sliced_channels[i])
        if i==0:
            print(sliced_channels[i].stats.starttime)
            first_date= mDate.get_utc_2_lmst(utc_date=sliced_channels[i].stats.starttime, output="decimal")
            first_date= math.modf(first_date)[0]
            print("Frac part of the SOL: ", first_date)
            first_date = first_date*mDate.SECOND_IN_A_DAY*mDate.SOL_RATIO
            print("In seconds: ", first_date)
            #posT_min = first_date.find('T')
            #first_time = first_date[posT_min+1:-1]
            #first_time2list = first_time.split(":")
            #first_time2second = float(first_time2list[0])*3600+float(first_time2list[1])*60+float(first_time2list[2])+\
            #        float(first_time2list[3])/1000.
            #print(first_time2second)
            t = np.arange(first_date,first_date+sliced_channels[i].stats.npts/sliced_channels[i].stats.sampling_rate, \
                                      sliced_channels[i].stats.delta)

        if i>0:
            t = np.arange(0,sliced_channels[i].stats.npts/sliced_channels[i].stats.sampling_rate, \
                                          sliced_channels[i].stats.delta)
        t/=(3600*mDate.SOL_RATIO)
        print(i, t[0], t[-1], len(t), sliced_channels[i].stats.endtime - sliced_channels[i].stats.starttime)
        axs.plot(t,sliced_channels[i].data+i*yshift, linewidth=0.6, label="SOL{:04d}".format(sol_min+i))

    axs.xaxis.set_tick_params(labelsize=9)
    axs.set(xlabel=xlabel)
    axs.set(ylabel=ylabel)
    axs.set_xlim(0, 24)
    axs.set_title(title)
    axs.legend(framealpha=0.01, loc='best', borderaxespad=0.,
               fontsize="x-small", ncol=3)
    axs.grid(True)
    axs.xaxis.set_major_locator(MaxNLocator(12, min_n_ticks=12))
    axs.xaxis.set_minor_locator(AutoMinorLocator(4))
    fig.tight_layout()
    plt.show()


def plot_channels(channel_list, xunits="UTC", yunits="Intensity (in DU)",
                  typeplot="subplot", decim=1, title=None, longname=False):
    """
    Plain function to plot channels from the list channel_list.
    User can choose if it's in overplot mode or un subplot mode.
    
    Using decim variable, user can decimate data before plotting them. There is
    no decimation besides the plots. Further processing will keep the full 
    resolution of the data. It's just a convienient way to speed up the plots 
    
    INPUT
    ----------
            @channel_list : list of channels, ie Trace object
                    Contains the list of channels to plot
            @xunits : str - The default is "UTC".
            @yunits : str - The default is "Intensity (in DU)".
            @decim  : int - The default value is 1
            @title: str - title of the plot
            @longname: boolean - if label of data are replace by longname
                                    or not.

    Returns
    -------
    

    """
    
    #Import useful lib
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter, MaxNLocator
    import sys

    import numpy as np
    from itertools import cycle

    if xunits == "LMST":
        from MarsConverter import MarsConverter
        mDate = MarsConverter()
    set_id = get_traces_id(channel_list)
    
    channel_list_decim = channel_list.copy()
    if decim != 1:

        if int(decim) <=16:
            for tr in channel_list_decim:
                tr.decimate(decim, no_filter=False, strict_length=False)
        else:
            sys.exit("Decimation factor must be lower than 16")

    if typeplot in ["subplot", "SP", "S"]:
        #plot_number = len(channel_list)
        plot_number = len(set_id)
        # ----------------------------------------------------------
        # Case where the plots are not on the same figure
        #
        figs, axs = plt.subplots(plot_number, 1, sharex=True,
                                 sharey=True, figsize=(10, 3))
        figs.subplots_adjust(hspace = .01)
        cycol = cycle('grcmk')  # Cycle of color to change for each time series
        index = 0
        for tr_id in set_id:
            color= next(cycol)
            for chan in channel_list_decim:
                l_id = tr_id.split("_")
                if (l_id[0] == chan.stats.network) and (l_id[1] == chan.stats.station) \
                    and (l_id[2] == chan.stats.location) and (l_id[3] == chan.stats.channel):
                    start_time = chan.stats.starttime
                    npts           = chan.stats.npts
                    sps                   = chan.stats.sampling_rate
                    delta           = chan.stats.delta
                    #print(start_time, npts , sps , delta)
                    #---------------------
                    # Define time vector
                    
                    t = np.arange(start_time, start_time + npts/sps, delta)
                    dte = [UTCDateTime(dat).timestamp for dat in t]
                    
                    #print(chan.stats.channel, len(t), len(chan.data[0:len(t)]))
                    #---------------------
                    #Plot the time series
                    label = chan.stats.location + chan.stats.channel 
                    axs[index].plot(dte, chan.data[0:len(t)], color = color, \
                                                    linestyle='-', marker='None', \
                                                    linewidth=1, label=label)
                    axs[index].legend()
                    #---------------------
                    # X-axis options
                    axs[index].grid(True)
                    axs[index].xaxis.set_tick_params(rotation=30, labelsize=6, \
                            length=6, width=2)
                            
                    axs[index].xaxis.set_major_locator(MaxNLocator(12))
                    axs[index].xaxis.label.set_size(6)
                    if xunits == "UTC":
                            axs[index].xaxis.set_major_formatter(FuncFormatter(lambda tstamp, \
                                                    pos : str(UTCDateTime(tstamp))[:-8]))
                            axs[index].set(xlabel="Date (UTC)")
                    elif xunits == "LMST":
                            axs[index].xaxis.set_major_formatter(FuncFormatter(lambda tstamp, \
                                    pos : mDate.get_utc_2_lmst(UTCDateTime(tstamp), \
                                                            output="date")[:-7]))
    
                            axs[index].set(xlabel="Date (LMST)")
    
                    # Y-axis options
                    axs[index].xaxis.label.set_size(8)
                    axs[index].set(ylabel=yunits)
            index+=1
        if title is not None:
            plt.title(str(title))
        plt.show()
    elif typeplot in ["overplot", "OP", "O"]:
        figs = plt.figure(figsize=(10, 3))
        axs = figs.add_subplot(111)
        axs.titlesize = 10
        cycol = cycle('grcmk')  # Cycle of color to change for each time series

        for index, chan in enumerate(channel_list_decim):
            start_time = chan.stats.starttime
            npts = chan.stats.npts
            sps = chan.stats.sampling_rate
            delta = chan.stats.delta
            # ---------------------
            # Define time vector
            t = np.arange(start_time, start_time + npts/sps, delta)
            dte = [UTCDateTime(dat).timestamp for dat in t]
            # ---------------------
            # Plot the time series
            label = chan.stats.location + chan.stats.channel 
            axs.plot(dte, chan.data[0:len(t)], color=next(cycol),
                     label=label,
                     linestyle='-', marker='None',
                     linewidth=1)
        axs.legend()
        #---------------------
        # X-axis options
        axs.xaxis.set_major_locator(MaxNLocator(12))
        axs.xaxis.label.set_size(6)
        if xunits == "UTC":
            axs.xaxis.set_major_formatter(FuncFormatter(lambda tstamp, \
                                    pos : str(UTCDateTime(tstamp))[:-8]))
            axs.set(xlabel="Date (UTC)")
        elif xunits == "LMST":
            axs.xaxis.set_major_formatter(FuncFormatter(lambda tstamp, \
                    pos : mDate.get_utc_2_lmst(UTCDateTime(tstamp), \
                                            output="date")[:-7]))

            axs.set(xlabel="Date (LMST)")

                # Y-axis options
        axs.xaxis.label.set_size(8)
        axs.set(ylabel=yunits)
        axs.grid(True)
        axs.xaxis.set_tick_params(rotation=30, labelsize=6, \
                        length=6, width=2)
        del(channel_list_decim)
        if title is not None:
            plt.title(str(title))
        plt.show()
    return figs, axs 



def plot_channels_v2(channel_list, xunits="UTC",
                     yunits="Intensity (in DU)",
                     typeplot="subplot", decim=1, 
                    title=None, longname=False):
    """
    Plain function to plot channels from the list channel_list.
    User can choose if it's in overplot mode or un subplot mode.

    Using decim variable, user can decimate data before plotting them. There
    is no decimation besides the plots. Further processing will keep the full
    resolution of the data. It's just a convienient way to speed up the plots

    INPUT
    ----------
        @channel_list : list of channels, ie Trace object
            Contains the list of channels to plot
        @xunits : str - The default is "UTC".
        @yunits : str - The default is "Intensity (in DU)".
        @decim  : int - The default value is 1

    Returns
    -------

    """

    set_id = get_traces_id(channel_list)

    channel_list_decim = channel_list.copy()
    if decim != 1:
        if int(decim) <= 16:
            for tr in channel_list_decim:
                tr.decimate(decim, no_filter=False,
                            strict_length=False)
        else:
            sys.exit("Decimation factor must be lower than 16")

    if typeplot.upper() in ["SUBPLOT", "SP", "S"]:

        plot_number = len(set_id)
        # ----------------------------------------------------------
        # Case where the plots are not on the same figure
        #
        figs, axs = plt.subplots(plot_number, 1, sharex=True,
                                 sharey=True, figsize=(10,6))
        figs.subplots_adjust(hspace = .01)
        cycol = cycle('grcmk')             # Cycle of color to change for each time series
        index = 0
        for tr_id in set_id:
            color= next(cycol)
            for chan in channel_list_decim:
                l_id = tr_id.split("_")
                if (l_id[0] == chan.stats.network) and (l_id[1] == chan.stats.station) \
                    and (l_id[2] == chan.stats.location) and (l_id[3] == chan.stats.channel):
                    start_time = chan.stats.starttime
                    npts = chan.stats.npts
                    sps           = chan.stats.sampling_rate
                    delta       = chan.stats.delta
                    # print(start_time, npts , sps , delta)
                    #---------------------
                    # Define time vector

                    t = np.arange(start_time, start_time + npts/sps, delta)
                    dte = [UTCDateTime(dat).timestamp for dat in t]
                    # print(chan.stats.channel, len(t), len(chan.data[0:len(t)]))
                    # ---------------------
                    # Plot the time series
                    #label_plot = get_lg_from_id(chan, DF_PROPERTIES) 
                    label_plot = chan.stats.location + "." + chan.stats.channel 
                    axs[index].plot(dte, chan.data[0:len(t)], color = color, \
                                    linestyle='-', marker='None', label=label_plot,\
                                    linewidth=1)
                    #axs[index].legend([tr_id]) 
                    #---------------------
                    # X-axis options
                    axs[index].grid(True)
                    # Y-axis options
                    axs[index].xaxis.label.set_size(8)
                    axs[index].set(ylabel=yunits)
                    axs[index].legend(bbox_to_anchor=(1.02, 1),
                                      loc='upper left', borderaxespad=0.)
                    dayplot_set_x_ticks(start_time, start_time + npts/sps,
                                        ax=axs[index], print_sol=False)
            index+=1
            
        ax_lmst = axs[index-1].twiny()
        ax_lmst.xaxis.set_ticks_position("bottom")
        ax_lmst.xaxis.set_label_position("bottom")
        ax_lmst.spines["bottom"].set_position(("axes", -0.23))
        ax_lmst.set_frame_on(True)
        ax_lmst.patch.set_visible(False)
        for sp in ax_lmst.spines.values():
            sp.set_visible(False)
        ax_lmst.spines["bottom"].set_visible(True)

        sol_start_time, sol = utc2lmst(start_time)
        sol_end_time, _ = utc2lmst(start_time + npts/sps)
        dayplot_set_x_ticks(sol_start_time, sol_end_time,
                                  ax=ax_lmst, print_sol=True)
        if title is not None:
            plt.suptitle(str(title))

        plt.show()
    elif typeplot.upper() in ["OVERPLOT", "OP", "O"]:
        figs = plt.figure(figsize=(10, 3))
        axs = figs.add_subplot(111)
        axs.titlesize = 10
        cycol = cycle('grcmk') # Cycle of color to change for each time series
    
        for index, chan in enumerate(channel_list_decim):
            start_time = chan.stats.starttime
            npts       = chan.stats.npts
            sps        = chan.stats.sampling_rate
            delta      = chan.stats.delta
            #---------------------
            # Define time vector
            t = np.arange(start_time, start_time + npts/sps, delta)
            dte = [UTCDateTime(dat).timestamp for dat in t]
            #---------------------
            #Plot the time series
            label_plot = chan.stats.location + "." + chan.stats.channel 
            axs.plot(dte, chan.data[0:len(t)], color= next(cycol), \
                            label = label_plot, \
                            linestyle='-', marker='None', \
                            linewidth=1)

        #---------------------
        # X-axis options
        axs.legend(bbox_to_anchor=(1.02, 1),
                        loc='upper left', borderaxespad=0.)
        dayplot_set_x_ticks(start_time, start_time + npts/sps,
                            ax=axs, print_sol=False)

            # Y-axis options
        axs.grid(True)
        ax_lmst = axs.twiny()
        ax_lmst.xaxis.set_ticks_position("bottom")
        ax_lmst.xaxis.set_label_position("bottom")
        ax_lmst.spines["bottom"].set_position(("axes", -0.23))
        ax_lmst.set_frame_on(True)
        ax_lmst.patch.set_visible(False)
        for sp in ax_lmst.spines.values():
            sp.set_visible(False)
        ax_lmst.spines["bottom"].set_visible(True)

        sol_start_time, sol = utc2lmst(start_time)
        sol_end_time, _ = utc2lmst(start_time + npts/sps)
        dayplot_set_x_ticks(sol_start_time, sol_end_time,
                            ax=ax_lmst, print_sol=True)
        
        del(channel_list_decim)
        if title is not None:
            plt.suptitle(str(title))
        plt.show()
    return figs, axs 