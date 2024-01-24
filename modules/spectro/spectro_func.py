#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 11 16:26:23 2019

@author  : greg
@email   : sainton`at`ipgp.fr
@purpose : plot spectrogram for a given miniseed file


This file contains various function to plot spectrograms. 

    - estimate_specgram_trace : estimate the spectrograms of one Trace.
    - estimate_specgram_stream: estimate the spectrograms of each Traces of a Stream.
    - plot_spectrogram        : to plot a single trace spectrogram.
    - plot_spectrogram_multi  : to overplot spectrograms. Useful for data with gaps.
    - spectro_from_mseed      : to plot spectrogram from mseed file.
    - spectro_from_trace      : to plot spectrogram from Trace object.

Those functions need obspy library to be installed to manage mseed files

"""

import os
import sys

from itertools import cycle

import math

import numpy as np
from obspy import UTCDateTime, read

from matplotlib.ticker import MaxNLocator, FuncFormatter 
from matplotlib.colors import Normalize
import matplotlib.pyplot as plt
from matplotlib import mlab              # lib used to estimate spectrogram

from time_utils import utc2lmst, utc2sol, lmst2utc, sol_span_in_utc,\
                            parse_location_channel_codes, dayplot_set_x_ticks


#==============================================================================
# Function to estimate the nearest power of 2 of a number

def _nearest_pow_2(x):
    """
    Function to estimate the nearest power of 2 of a number
    ----
    INPUT:
        @x: int - number for which you are looking the power of 2
    ----
    OUTPUT:
        either a or b 
    """

    a = math.pow(2, math.ceil(np.log2(x)))
    b = math.pow(2, math.floor(np.log2(x)))
    if abs(a - x) < abs(b - x):
        return a
    else:
        return b

# ======================================================================
# Function of estimate spectrograms for a Stream of data
#
def estimate_specgram_stream(stream, wlen=200, per_lap=0.9, mult=8):
    """
    Compute and plot a spectrogram of data in Stream Data are split into 
    NFFT length segments and the spectrum of each section is computed. 
    The windowing function window is applied to each segment, and the 
    amount of overlap of each segment is specified with noverlap.
    
    ----
    INPUT:
        @stream  : Trace object - Trace with the time serie
        @wlen   : int/float - Window length for fft in seconds
        @per_lap: float - Percentage of overlap of sliding window, 
                    ranging from 0 to 1.
        @mult   : float - Pad zeros to length mult * wlen. 
                  This will make the spectrogram smoother.
    ----
    OUTPUT:
        @spectrum : list of list - 2-D array, columns are the periodograms 
                                of successive segments.
        @freq     : list of list - 1-D array, frequencies corresponding 
                                to the rows in spectrum.
        @time     : list of list - 1-D array, the times corresponding to 
                                midpoints of segments 
                                (i.e the columns in spectrum).
        
    """
    sta = stream.split()
    
    specgram_list = [] 
    freq_list = [] 
    time_list = []
    
    for trace in sta:
        specgram, freq, time =  estimate_specgram_trace(trace, wlen=wlen, \
                                                  per_lap=per_lap, \
                                                  mult = mult)
        specgram_list.append(specgram)
        freq_list.append(freq)
        time_list.append(time)
    
    return specgram_list, freq_list, time_list


# =============================================================================
# Function of estimate spectrograms for a Trace of data
# 
def estimate_specgram_trace(trace, wlen=200, per_lap=0.9, \
                            mult=8.0, detrend='linear'):    
    """
    Compute and plot a spectrogram of data in x. Data are split into 
    NFFT length segments and the spectrum of each section is computed. 
    The windowing function window is applied to each segment, and the 
    amount of overlap of each segment is specified with noverlap.
    
    ----
    INPUT:
        @trace  : Trace object - Trace with the time serie
        @wlen   : int/float - Window length for fft in seconds
        @per_lap: float - Percentage of overlap of sliding window, 
                    ranging from 0 to 1.
        @mult   : float - Pad zeros to length mult * wlen. 
                  This will make the spectrogram smoother.
    ----
    OUTPUT:
        @spectrum : array_like - 2-D array, columns are the periodograms 
                                of successive segments.
        @freq     : array_like - 1-D array, frequencies corresponding 
                                to the rows in spectrum.
        @time     : array_like - 1-D array, the times corresponding to 
                                midpoints of segments 
                                (i.e the columns in spectrum).

    """

    samp_rate = trace.stats.sampling_rate
    data= trace.data
    
    # set wlen from samp_rate if not specified otherwise
    if not wlen:
        wlen = samp_rate / 100.

    npts = len(data)

    # nfft needs to be an integer, otherwise a deprecation will be raised
    # XXX add condition for too many windows => calculation takes for ever

    nfft = int(_nearest_pow_2(wlen * samp_rate))
    
    if nfft > npts:
        nfft = int(_nearest_pow_2(npts / 8.0))

    if mult is not None:
        mult = int(_nearest_pow_2(mult))
        mult = mult * nfft

    nlap = int(nfft * float(per_lap))

    # spectrogram estimation
    specgram, freq, time = mlab.specgram(data, Fs=samp_rate, NFFT=nfft, \
                                         pad_to=mult, noverlap=nlap, \
                                         detrend=detrend)

    return specgram, freq, time

# =============================================================================
# Plot spectrogram according to the result of the spectrogram 
#    estimated in the function estimate_specgram_trace
def plot_spectrogram(specgram, freq, time, trace, \
                     dbscale=True, clip=[0.0,1.0], log=True, \
                     fmin=0, fmax=1, vmin_param=None, \
                     vmax_param=None, xunits = "UTC", yunits="Counts"):
    """
    INPUT:
        @spectrum : array_like - 2-D array, columns are the periodograms 
                                of successive segments.
        @freq     : array_like - 1-D array, frequencies corresponding 
                                to the rows in spectrum.
        @time     : array_like - 1-D array, the times corresponding to 
                                midpoints of segments 
                                (i.e the columns in spectrum).
        @trace    : Trace object used to calculate the spectrogram
    ----
    OUTPUT:
        @figs : Figure object
        @axs  : Axes object 
    
    """

    #----------
    if dbscale:
        specgram = 10 * np.log10(specgram[1:, :])
    else:
        specgram = np.sqrt(specgram[1:, :])
    start_time = trace.stats.starttime
    freq = freq[1:]
    dt = trace.stats.delta
    npts = trace.stats.npts
    sps = trace.stats.sampling_rate

    vmin, vmax = clip
    if vmin < 0 or vmax > 1 or vmin >= vmax:
        msg = "Invalid parameters for clip option."
        raise ValueError(msg)
    
    if vmin_param is not None:
        vmin = vmin_param
        
    if vmax_param is not None:
        vmax = vmax_param
    
    norm = Normalize(vmin, vmax, clip=True)

    # calculate half bin width
    halfbin_time = (time[1] - time[0]) / 2.0
    halfbin_freq = (freq[1] - freq[0]) / 2.0

    fig = plt.figure(figsize=(10, 3))
    axs = fig.add_subplot(111)
    
    if log:
        # pcolor expects one bin more at the right end
        freq = np.concatenate((freq, [freq[-1] + 2 * halfbin_freq]))
        time = np.concatenate((time, [time[-1] + 2 * halfbin_time]))
        time_shift = np.add(time, UTCDateTime(trace.stats.starttime).timestamp)
        # center bin
        time -= halfbin_time
        freq -= halfbin_freq
        # Log scaling for frequency values (y-axis)
        axs.set_yscale('log')
        # Plot times
        axs.pcolormesh(time_shift, freq, specgram, norm=norm, cmap="jet")
    else:
        time_shift = np.add(time,UTCDateTime(trace.stats.starttime).timestamp)
        # this method is much much faster!
        specgram = np.flipud(specgram)
        # center bin
        extent = (time_shift[0] - halfbin_time, time_shift[-1] + halfbin_time,
                      freq[0] - halfbin_freq, freq[-1] + halfbin_freq)
        
        axs.imshow(specgram, interpolation="nearest", extent=extent, cmap="jet")

    
    # set correct way of axis, whitespace before and 
    # after with window length
    axs.axis('tight')
    axs.set_xlim(UTCDateTime(trace.stats.starttime).timestamp, \
                 UTCDateTime(trace.stats.endtime).timestamp)
    axs.grid(False)
    axs.xaxis.set_major_locator(MaxNLocator(12))
    if xunits == "LMST":
        from MarsConverter import MarsConverter
        mDate = MarsConverter()
        axs.xaxis.set_major_formatter(FuncFormatter(lambda tstamp, \
                                pos : mDate.get_utc_2_lmst(UTCDateTime(tstamp), \
                                output="date")[:-10]))
        axs.set(xlabel="Date (LMST)")
    else: 
        axs.xaxis.label.set_size(8)
        #axs.set(ylabel=yunits)
        axs.legend(bbox_to_anchor=(1.02, 1),
                          loc='upper left', borderaxespad=0.)
        dayplot_set_x_ticks(start_time, start_time + npts/sps,
                            ax=axs, print_sol=False)

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

    axs.set_ylabel('Frequency [Hz]')
    title = True
    if title:
        axs.set_title("Spectrogram of "+trace.stats.location+"/"+trace.stats.channel)

    #Set f limits
    if fmin is None: 
        fmin = 1/(npts*dt)
    if fmax is None:
        fmax = 1/(2*dt)

    axs.set_ylim(fmin, fmax)

    axs.set_ylabel('Frequency [Hz]')  
    axs.xaxis.label.set_size(8)
    axs.xaxis.set_tick_params(rotation=20, labelsize=8)  
    fig.autofmt_xdate()

    mappable = axs.collections[0]

    #Colorbar on the right of the spectrogram
    cbar = plt.colorbar(mappable=mappable, norm=norm)
    if dbscale:
        ylab = r"$($"+yunits+"$)^2/Hz$"+" [dB]"
    else:
        ylab = yunits+"$/\sqrt{Hz}$"

    cbar.set_label(ylab)
    plt.show()    
    

    return fig, axs


#========================================================================
# Function to plot spectrograms with gaps

def plot_spectrogram_multi(specgram_list, freq_list, time_list, trace, \
                        dbscale=True, clip=[0,10], log=True, fmin=0, fmax=1, \
                        vmin_param=None, vmax_param=None, \
                        xunits = "UTC", yunits = "Counts", plotTS = False):
    '''
    Function to plot spectrogram where signal contains some gaps.
    This function is overplotting the spectrogram of each part of the
    time serie.
    To make it works, one need first to run "estimate_specgram_stream"
    to calculate "specgram_list", "freq_list" and "time_list"
    eg: specgram_l, freq_l, time_l = \
                estimate_specgram_stream(rotated_stream, wlen=2000, per_lap=0.2, mult = 8)
    ----
    INPUT:
        @specgram_list: list of spectrogram computed with the function 
            estimate_specgram_stream(stream, wlen=200, per_lap=0.9, mult=8)
            The list contains as many elements as elements in the initial time serie
        @freq_list: list of vector with frequency
        @time_list: list of vector with time
        @channel:list of traces (after splitting)
        @dbscale: boolean to choose is the scale is in decibels or not
        @fmin: minimum frequency
        @fmax: maximum frequency
        @vmin_param: lower limit of the spectrogram (ie: -200db)
        @vmax_param: uppe limit of the spectrogram (ie: 0db)
        @xunits: (str) can take the value LMST or UTC
        @yunits= (str) can take any values. Used for labelling axes
    
    ----    
    OUTPUT:
        @figs: Figure object to be plotted
        @axs:  Axes to put on the figure.
        
    '''
    
    nb_plt = len(specgram_list)
    plottitle = "Spectrogram"
    starttime= trace[0].stats.starttime
    endtime =  trace[-1].stats.endtime
    dt = trace[0].stats.delta
    npts = trace[0].stats.npts
    
    plt.suptitle("Spectrogram", fontsize = 7)


    if plotTS in [False,"False"]:
        figs = plt.figure(figsize=(10, 6))
        for i in range(0, nb_plt):
            
            # Hack used to overlap axes containing all the spectrograms
            # The first one is a classical axe.
            # Others are completly overlapped with a transparent background
            # axes are completly shared.
            
            if i == 0:
                ax1 = figs.add_axes([0.2,0.2,0.8,0.7], zorder=0)
            else:    
                ax1 = figs.add_axes([0.2,0.2,0.8,0.7], zorder=0, \
                        sharex=ax1, sharey=ax1)
                #set a complete transparency of the background
                ax1.patch.set_alpha(0)
            
            specgram = specgram_list[i]
            freq = freq_list[i]
            time = time_list[i]

            if dbscale:
                specgram = 10 * np.log10(specgram[1:, :])
            else:
                specgram = np.sqrt(specgram[1:, :])

            freq = freq[1:]
            vmin, vmax = clip

            if vmin_param is not None:
                vmin = vmin_param

            if vmax_param is not None:
                vmax = vmax_param

            norm = Normalize(vmin, vmax, clip=True)

            # calculate half bin width
            halfbin_time = (time[1] - time[0]) / 2.0
            halfbin_freq = (freq[1] - freq[0]) / 2.0

            if log:
                # pcolor expects one bin more at the right end
                freq = np.concatenate((freq, [freq[-1] + 2 * halfbin_freq]))
                time = np.concatenate((time, [time[-1] + 2 * halfbin_time]))
                time_shift = np.add(time,UTCDateTime(trace[i].stats.starttime).timestamp)
    
                # center bin
                time -= halfbin_time
                freq -= halfbin_freq
                # Log scaling for frequency values (y-axis)
                ax1.set_yscale('log')
                # Plot times
                ax1.pcolormesh(time_shift, freq, specgram, norm=norm, cmap="jet")
            else:
                time_shift = np.add(time,UTCDateTime(trace[i].stats.starttime).timestamp)
                # this method is much much faster!
                specgram = np.flipud(specgram)
                # center bin
                extent = (time_shift[0] - halfbin_time, time_shift[-1] + halfbin_time,
                              freq[0] - halfbin_freq, freq[-1] + halfbin_freq)
    
                ax1.imshow(specgram, interpolation="nearest", extent=extent, cmap="jet")
    
            # set the xlimit with the largest time range defined by the 
            # starttime of the first element
            # and the endtime of the last element of the channel
            ax1.set_xlim(starttime.timestamp, endtime.timestamp)
    
            # command to set a given number of ticks
            ax1.xaxis.set_major_locator(MaxNLocator(12))
    
            # desativate the grid
            ax1.grid(False)

            if xunits == "LMST" :
                if xunits == "LMST": 
                    from MarsConverter import MarsConverter
                    mDate = MarsConverter()
                ax1.xaxis.set_major_formatter(FuncFormatter(lambda tstamp, \
                                        pos : mDate.get_utc_2_lmst(UTCDateTime(tstamp), \
                                        output="date")[:-10]))
                ax1.set(xlabel="Date (LMST)")
            
            else: 
                ax1.xaxis.set_major_formatter(FuncFormatter(lambda tstamp, \
                                            pos : str(UTCDateTime(tstamp))[:-8]))
                ax1.set(xlabel="Date (UTC)")
            
            ax1.set_ylabel('Frequency [Hz]')
    
            #Set frequency limits
            if fmin is None: 
                fmin = 1/(npts*dt)
            if fmax is None:
                fmax = 1/(2*dt)
            ax1.set_ylim(fmin, fmax)

            ax1.xaxis.label.set_size(8)
            ax1.xaxis.set_tick_params(rotation=20, labelsize=8)  

            mappable = ax1.collections[0]

            #Colorbar on the right of the spectrogram
            cbar = plt.colorbar(mappable=mappable, norm=norm)
            if dbscale:
                ylab = r"$($"+yunits+"$)^2/Hz$"+" [dB]"
            else:
                ylab = yunits+"$/\sqrt{Hz}$"
            
            cbar.set_label(ylab)
            
            plottitle = "Spectrogram of "+str(trace[0].stats.location)+"/"+str(trace[0].stats.channel)
            figs.suptitle(plottitle, fontsize=8)
            #axs = [ax1]    
            figs.autofmt_xdate()

    else: # Case where the time serie is added above the spectrogram
        print("Plot avec TS")
        from itertools import cycle
        from matplotlib import gridspec
        from mpl_toolkits.axes_grid1.inset_locator import inset_axes
        
        figs = plt.figure(figsize=(10, 5))
        figs.subplots_adjust(hspace=0.01)

        gs = gridspec.GridSpec(2, 2, width_ratios=[15, 1]) 
        ax1 = plt.subplot(gs[0])
        ax1.set(ylabel = yunits)
        ax2 = plt.subplot(gs[2], sharex = ax1)
        #ax2 = plt.subplot(gs[2])
        #ax3 = plt.subplot(gs[3])
        
        cycol = cycle('grcmk')             # Cycle of color to change for each time series
        channel_list_decim = trace.copy()
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
            ax1.plot(dte, chan.data[0:len(t)], color= next(cycol), \
                            label = chan.stats.channel, \
                            linestyle='-', marker='None', \
                            linewidth=1)
        
        ax1.legend() 
        ax1.tick_params(
            axis='x',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            bottom='off',      # ticks along the bottom edge are off
            top='off',         # ticks along the top edge are off
            labelbottom='off'  # labels along the bottom edge are off
            )
        for ax in plt.gcf().axes:
            try:
                ax.label_outer()
            except:
                pass
        #---------------------
        # X-axis options
        #ax1.xaxis.set_major_locator(MaxNLocator(12))
        #ax1.axes.get_xaxis().set_ticklabels([])
        #ax1.xaxis.label.set_size(6)
        #if xunits == "UTC":
        #    ax1.xaxis.set_major_formatter(FuncFormatter(lambda tstamp, \
        #                pos : str(UTCDateTime(tstamp))[:-8]))
        #    ax1.set(xlabel="Date (UTC)")
        #elif xunits == "LMST":
        #    ax1.xaxis.set_major_formatter(FuncFormatter(lambda tstamp, \
        #        pos : mDate.get_utc_2_lmst(UTCDateTime(tstamp), \
        #                    output="date")[:-7]))

            #ax1.set(xlabel="Date (LMST)")
    
            # Y-axis options
        #ax2.xaxis.label.set_size(8)
        ax1.set(ylabel=yunits)
        ax1.grid(True)
        
        #ax1.set_xticklabels([])
        #ax2.xaxis.set_tick_params(rotation=30, labelsize=6, \
        #        length=6, width=2)

        for i in range(0, nb_plt):
        
        # Hack used to overlap axes containing all the spectrograms
        # The first one is a classical axe.
        # Others are completly overlapped with a transparent background
        # axes are completly shared.
        
            #if i == 0:
            #    ax2 = figs.add_axes([0.2,0.1,0.7,0.5], zorder=0)
            #    ax1 = figs.add_axes([0.2,0.6,0.7,0.5], zorder=0)
            #else:    
            #    ax2 = figs.add_axes([0.2,0.2,0.7,0.5], zorder=0, \
            #            sharex=ax2, sharey=ax2)
            #    ax1 = figs.add_axes([0.2,0.6,0.7,0.5], zorder=0)
            #    #set a complete transparency of the background
            ax2.patch.set_alpha(0)
            
            specgram = specgram_list[i]
            freq = freq_list[i]
            time = time_list[i]
        
            if dbscale:
                specgram = 10 * np.log10(specgram[1:, :])
            else:
                specgram = np.sqrt(specgram[1:, :])
    
            freq = freq[1:]
    
            vmin, vmax = clip
            #if vmin < 0 or vmax > 1 or vmin >= vmax:
            #    msg = "Invalid parameters for clip option."
            #    raise ValueError(msg)
            
    
            if vmin_param is not None:
                vmin = vmin_param
    
            if vmax_param is not None:
                vmax = vmax_param
    
            norm = Normalize(vmin, vmax, clip=True)
    
            # calculate half bin width
            halfbin_time = (time[1] - time[0]) / 2.0
            halfbin_freq = (freq[1] - freq[0]) / 2.0
    
            if log:
                # pcolor expects one bin more at the right end
                freq = np.concatenate((freq, [freq[-1] + 2 * halfbin_freq]))
                time = np.concatenate((time, [time[-1] + 2 * halfbin_time]))
                time_shift = np.add(time,UTCDateTime(trace[i].stats.starttime).timestamp)
    
                # center bin
                time -= halfbin_time
                freq -= halfbin_freq
                # Log scaling for frequency values (y-axis)
                ax2.set_yscale('log')
                # Plot times
                ax2.pcolormesh(time_shift, freq, specgram, norm=norm, cmap="jet")
            else:
                time_shift = np.add(time,UTCDateTime(trace[i].stats.starttime).timestamp)
                # this method is much much faster!
                specgram = np.flipud(specgram)
                # center bin
                extent = (time_shift[0] - halfbin_time, time_shift[-1] + halfbin_time,
                              freq[0] - halfbin_freq, freq[-1] + halfbin_freq)
    
                ax2.imshow(specgram, interpolation="nearest", extent=extent, cmap="jet")
    
            # set the xlimit with the largest time range defined by the 
            # starttime of the first element
            # and the endtime of the last element of the channel
            ax2.set_xlim(starttime.timestamp, endtime.timestamp)
    
            # command to set a given number of ticks
            ax2.xaxis.set_major_locator(MaxNLocator(12))
    
            # desativate the grid
            ax2.grid(False)
    
            if xunits == "LMST" :
                if xunits == "LMST": 
                    from MarsConverter import MarsConverter
                    mDate = MarsConverter()
                ax2.xaxis.set_major_formatter(FuncFormatter(lambda tstamp, \
                                        pos : mDate.get_utc_2_lmst(UTCDateTime(tstamp), \
                                        output="date")[:-10]))
                ax2.set(xlabel="Date (LMST)")
            
            else: 
                ax2.xaxis.set_major_formatter(FuncFormatter(lambda tstamp, \
                                            pos : str(UTCDateTime(tstamp))[:-8]))
                ax2.set(xlabel="Date (UTC)")
            
            ax2.set_ylabel('Frequency [Hz]')

            #Set frequency limits
            if fmin is None: 
                fmin = 1/(npts*dt)
            if fmax is None:
                fmax = 1/(2*dt)
            ax2.set_ylim(fmin, fmax)

            ax2.xaxis.label.set_size(8)
            ax2.xaxis.set_tick_params(rotation=20, labelsize=8)  

        mappable = ax2.collections[0]

        axins = inset_axes(ax2,
                 width="2%",  # width = 5% of parent_bbox width
                 height="100%",  # height : 50%
                 loc='lower left',
                 bbox_to_anchor=(1.05, 0., 1, 1),
                 bbox_transform=ax2.transAxes,
                 borderpad=0)

        #Colorbar on the right of the spectrogram
        cbar = plt.colorbar(mappable= mappable, cax= axins, norm=norm,
                            orientation="vertical", pad=-1.0)
        if dbscale:
            ylab = r"$($"+yunits+"$)^2/Hz$"+" [dB]"
        else:
            ylab = yunits+"$/\sqrt{Hz}$"
        
        cbar.set_label(ylab)
            
        plottitle = "Spectrogram of "+str(trace[0].stats.location)+"/"+str(trace[0].stats.channel)
        figs.suptitle(plottitle, fontsize=8)
        #axs = [ax2]    
        figs.autofmt_xdate()
        plt.show()



#==============================================================================
# Plot a spectrogram for each traces of a mseed file
def spectro_from_mseed(filename):
    """
    Plot spectrogram for every trace of a mseed file.
    The function does care if the trace are different peaces of a same channel
    
    ----
    INPUT:
        filename: string - path to a mseed file.
    
    """
    
    stream = read(filename)
    for tr in stream:
        spectro_from_trace(tr)

#==============================================================================
# Plot a spectrogram of a single trace
def spectro_from_trace(trace, wlen=2000, \
                    per_lap=0.9, mult = 8, fmin=0.01, fmax=10,\
                    vmin_param=-200, vmax_param=-140, \
                    xunits = "UTC", yunits=r"$m.s^{-2}$"):
    """
    Function to plot the spectrogram if a single trace.
    ----
    INPUT:
        @trace: Trace object 
        @trace  : Trace object - Trace with the time serie
        @wlen   : int/float - Window length for fft in seconds
        @per_lap: float - Percentage of overlap of sliding window, 
                    ranging from 0 to 1.
        @mult   : float - Pad zeros to length mult * wlen. 
                  This will make the spectrogram smoother.
        @fmin: minimum frequency
        @fmax: maximum frequency
        @vmin_param: lower limit of the spectrogram (ie: -200db)
        @vmax_param: uppe limit of the spectrogram (ie: 0db)
        @xunits: (str) can ONLY take the value 'LMST' or 'UTC'
        @yunits= (str) can take any values. Used for labelling axes
    
    """

    import matplotlib.pyplot as plt
    specgram, freq, time= estimate_specgram_trace(trace, wlen=wlen, \
                                                  per_lap=per_lap, \
                                                  mult=mult, \
                                                  detrend='linear')
    
    figs, axs = plot_spectrogram(specgram, freq, time, \
                        trace, \
                        fmin=fmin, fmax=fmax,\
                        vmin_param=vmin_param, vmax_param=vmax_param, \
                        xunits = xunits, yunits=yunits)
    plt.show()

def psd_process_on_channel_list(channel_list, nperseg = None, nfft = None,
                                freq_choice = "PSD"):
    
    """
    Function to apply the Power spectral density (PSD) on each channels of a list. 
    If choice was ASD, the function return the square root of the signal
    
    In this function we are using the Welsh method to estimate PSD. 
    More information here : 
    https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.signal.welch.html
    
    ----
    INPUT:
    
    @channel_list: stream - input station 
    @str_freq_type: string - the type of frequency process (so far welsh)
    @str_options: string - string containing the 
                parameters and their value (eg: 'freq=5')
    @option_dict: dictionnary with all the options 
    @freq_choice: string - "PSD" or ASD
    
    ----
    OUTPUT:
        @freq_channels
        @density_channels
        @index_in_channelTab
    
    """
    import math
    from scipy import signal
    import numpy as np
    #from itertools import cycle
    
    freq_channels = []
    density_channels = []
    f = []
    out_y = []
    index_in_channelTab = []
    for i in range(0, len(channel_list)):
        cptrace = channel_list[i].copy()
        datatrace = cptrace.data
        
        freqs = cptrace.stats.sampling_rate
        nfft_custom = 2**math.ceil(math.log(len(datatrace))/math.log(2.))
        nperseg = int(nperseg) if nperseg != None else len(datatrace)
        nfft = int(nfft) if nfft != None  else nfft_custom

        # processing of the power spectral density
        f, out_y = signal.welch(datatrace, fs=freqs, window='hanning', 
                                nperseg=nperseg, 
                                #noverlap=noverlap, 
                                nfft=nfft)
        freq_channels.append(f)

        # Compute the amplitude spectral density
        if freq_choice == "ASD":
            out_y = np.sqrt(out_y)

        density_channels.append(out_y)
        index_in_channelTab.append(i)

    return freq_channels, density_channels, index_in_channelTab


def create_freq_plot(channelList,
                    freq_choice = "PSD",
                    typeplot = "overplot",
                    xunits="Hertz", yunits = "m.s^-2"):
    '''
    This function is creating loglog plots using matplotlib.
    
    @selected_traces_freq: list of channels to be plotted
    @selected_traces_density: parallel list of channelsTab containing 
    all the informations about the channels
    @channelTab_infos: informations about the channels
    @frequency_Options: Informations about the filter applied
    @ytitle: global title for y axis
    @plottitle: global title of the plot
    @return: figs and axs, elements of matplotlib
    '''
    #selected_traces_freq, selected_traces_density, index_in_channelTab

    freq_channels, density_channels, index_in_channelTab = \
                psd_process_on_channel_list(channelList,
                                            nperseg = 2000,
                                            freq_choice = "PSD")

    #cycle of colors to change color for each signal
    cycol = cycle('grcmk')
    xtitle = "Hertz"
 
    #fig initialization
    nb_plt = len(freq_channels)
    
    # Subplot choice
    if typeplot.lower() =="subplot" and nb_plt>1:
        figs, axs = plt.subplots(nb_plt, 1, sharex=True,
                                 sharey= True,
                                 figsize=(11,7))
        if freq_choice == "PSD":
            figs.suptitle("Power Spectral Density")
        elif freq_choice == "ASD":
            figs.suptitle("Amplitude Spectral Density")
        # Remove horizontal space between axes
        figs.subplots_adjust(hspace=0.2)
        figs.canvas.set_window_title('Frequency')
        for index in range(0, nb_plt):
            nameChannel = channelList[index].stats.location+"."+\
                channelList[index].stats.channel
            axs[index].loglog(freq_channels[index], density_channels[index],
                              color= next(cycol), label = nameChannel, 
                              linestyle='-',marker='None', linewidth=0.7)

            if freq_choice == "PSD":
                axs[index].set(ylabel=r"$($"+yunits+r"$)^2/Hz$")
            elif freq_choice == "ASD":
                axs[index].set(ylabel=yunits+r"$/\sqrt{Hz}$")
                

            axs[index].set(xlabel=xunits)
            axs[index].grid(True, which='both', axis='both')
            axs[index].tick_params(labelsize='medium', width=1)
            axs[index].legend()
        return figs, axs
    # Overplot choice
    else:
        figs = plt.figure(figsize=(10, 6))
        if freq_choice == "PSD":
            figs.suptitle("Power Spectral Density")
        elif freq_choice == "ASD":
            figs.suptitle("Amplitude Spectral Density")
            
        axs = figs.add_subplot(111)
        axs = figs.gca()
        
        for index in range(0, nb_plt):
            nameChannel = channelList[index].stats.location+"."+\
                channelList[index].stats.channel
            axs.loglog(freq_channels[index],
                       density_channels[index],
                       color= next(cycol),
                       label = nameChannel,
                       linestyle='-',
                       marker='None',
                       linewidth=0.7)
            
            if freq_choice == "PSD":
                axs.set(ylabel=r"$($"+yunits+r"$)^2/Hz$")
            elif freq_choice == "ASD":
                axs.set(ylabel=yunits+r"$/\sqrt{Hz}$")
            
            axs.grid(True, which='both', axis='both')
            
            axs.set(xlabel=xunits)
            axs.tick_params(labelsize='medium', width=2)
            axs.legend()
    plt.show()

#========================================================================
# Function to plot spectrograms with gaps

def plot_spectrogram_multi2(specgram_list, freq_list, time_list, channel,
                           dbscale=True, clip=[0,10], log=True, fmin=0, fmax=1,
                           vmin_param=None, vmax_param=None, xunits = "UTC",
                           yunits = "Counts"):
    '''
    Function to plot spectrogram where signal contains some gaps.
    This function is overplotting the spectrogram of each part of the
    time serie.
    To make it works, one need first to run "estimate_specgram_stream"
    to calculate "specgram_list", "freq_list" and "time_list"
    eg: specgram_l, freq_l, time_l = \
                estimate_specgram_stream(rotated_stream, wlen=2000, per_lap=0.2, mult = 8)
    ----
    INPUT:
        @specgram_list: list of spectrogram computed with the function 
            estimate_specgram_stream(stream, wlen=200, per_lap=0.9, mult=8)
            The list contains as many elements as elements in the initial time serie
        @freq_list: list of vector with frequency
        @time_list: list of vector with time
        @channel:list of traces (after splitting)
        @dbscale: boolean to choose is the scale is in decibels or not
        @fmin: minimum frequency
        @fmax: maximum frequency
        @vmin_param: lower limit of the spectrogram (ie: -200db)
        @vmax_param: uppe limit of the spectrogram (ie: 0db)
        @xunits: (str) can take the value LMST or UTC
    
    ----    
    OUTPUT:
        @figs: Figure object to be plotted
        @axs:  Axes to put on the figure.
        
    '''
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator, FuncFormatter 
    from matplotlib.colors import Normalize
    
    nb_plt = len(specgram_list)
    plottitle = "Spectrogram"
    figs = plt.figure(figsize=(10, 3))
    starttime= channel[0].stats.starttime
    endtime =  channel[-1].stats.endtime
    #import colormaps as cmaps 
    from matplotlib import cm as cmaps
    for i in range(0, nb_plt):
        
        # Hack used to overlap axes containing all the spectrograms
        # The first one is a classical axe.
        # Others are completly overlapped with a transparent background
        # axes are completly shared.
        

        
        if i == 0:
            ax1 = figs.add_axes([0.2,0.2,0.8,0.7], zorder=0)
        else:    
            ax1 = figs.add_axes([0.2,0.2,0.8,0.7], zorder=0, sharex=ax1, sharey=ax1)
            #set a complete transparency of the background
            ax1.patch.set_alpha(0)
        
        
        specgram = specgram_list[i]
        freq = freq_list[i]
        time = time_list[i]
    
        if dbscale:
            specgram = 10 * np.log10(specgram[1:, :])
        else:
            specgram = np.sqrt(specgram[1:, :])

        freq = freq[1:]

        vmin, vmax = clip
        if vmin < 0 or vmax > 1 or vmin >= vmax:
            msg = "Invalid parameters for clip option."
            raise ValueError(msg)
        #_range = float(specgram.max() - specgram.min())

        if vmin_param is not None:
            vmin = vmin_param

        if vmax_param is not None:
            vmax = vmax_param

        norm = Normalize(vmin, vmax, clip=True)

        # calculate half bin width
        halfbin_time = (time[1] - time[0]) / 2.0
        halfbin_freq = (freq[1] - freq[0]) / 2.0

        if log:
            # pcolor expects one bin more at the right end
            freq = np.concatenate((freq, [freq[-1] + 2 * halfbin_freq]))
            time = np.concatenate((time, [time[-1] + 2 * halfbin_time]))
            time_shift = np.add(time,UTCDateTime(channel[i].stats.starttime).timestamp)

            # center bin
            time -= halfbin_time
            freq -= halfbin_freq
            # Log scaling for frequency values (y-axis)
            ax1.set_yscale('log')
            # Plot times
            #ax1.pcolormesh(time_shift, freq, specgram, norm=norm, cmap="jet")
            ax1.pcolormesh(time_shift, freq, specgram, norm=norm, cmap=cmaps.viridis)
            #cmap=cmaps.viridis
        else:
            time_shift = np.add(time,UTCDateTime(channel[i].stats.starttime).timestamp)
            # this method is much much faster!
            specgram = np.flipud(specgram)
            # center bin
            extent = (time_shift[0] - halfbin_time, time_shift[-1] + halfbin_time,
                          freq[0] - halfbin_freq, freq[-1] + halfbin_freq)

            #ax1.imshow(specgram, interpolation="nearest", extent=extent, cmap="jet")
            ax1.imshow(specgram, interpolation="nearest", extent=extent, cmap=cmaps.viridis)

        # set the xlimit with the largest time range defined by the starttime of the first element
        # and the endtime of the last element of the channel
        ax1.set_xlim(starttime.timestamp, endtime.timestamp)

        # command to set a given number of ticks
        ax1.xaxis.set_major_locator(MaxNLocator(12))

        # desativate the grid
        ax1.grid(False)

        if xunits == "LMST" :
            if xunits == "LMST": 
                from MarsConverter import MarsConverter
                mDate = MarsConverter()
            ax1.xaxis.set_major_formatter(FuncFormatter(lambda tstamp, \
                                    pos : mDate.get_utc_2_lmst(UTCDateTime(tstamp), \
                                    output="date")[:-10]))
            ax1.set(xlabel="Date (LMST)")
        
        else: 
            ax1.xaxis.set_major_formatter(FuncFormatter(lambda tstamp, \
                                        pos : str(UTCDateTime(tstamp))[:-8]))
            ax1.set(xlabel="Date (UTC)")
        
        ax1.set_ylabel('Frequency [Hz]')

        #Set frequency limits
        if fmin is None: 
            fmin = 1/(npts*dt)
        if fmax is None:
            fmax = 1/(2*dt)
        ax1.set_ylim(fmin, fmax)

        ax1.xaxis.label.set_size(8)
        ax1.xaxis.set_tick_params(rotation=20, labelsize=8)  

        mappable = ax1.collections[0]

        #Colorbar on the right of the spectrogram
        cbar = plt.colorbar(mappable=mappable, norm=norm)
        if dbscale:
            ylab = r"$($"+yunits+"$)^2/Hz$"+" [dB]"
        else:
            ylab = yunits+"$/\sqrt{Hz}$"
        
        cbar.set_label(ylab)
        
        plottitle = "Spectrogram of "+str(channel[0].stats.location)+"/"+str(channel[0].stats.channel)
        figs.suptitle(plottitle, fontsize=8)
        axs = [ax1]    
        figs.autofmt_xdate()
    return figs, axs



#==============================================================================
# MAIN PROGRAM

if __name__ == "__main__":

    # Part added to manage MarsConverter lib
    try: 
        MARSCONVERTER = os.environ["MARSCONVERTER"]
        print("MARSCONVERTER: ", MARSCONVERTER)
    except KeyError:
        sys.exit("MARSCONVERTER env path is not defined. Please do it first.")
    else:
        sys.path.insert(0, MARSCONVERTER)
    
    from MarsConverter import MarsConverter
    
    print('Welcome in spectrogram program.')
    print("Now is ", UTCDateTime.now())
    mDate = MarsConverter()
    marsDateNow = mDate.get_utc_2_lmst()
    posT = marsDateNow.find('T')

    print("in LMST, now, it is ", marsDateNow)
    