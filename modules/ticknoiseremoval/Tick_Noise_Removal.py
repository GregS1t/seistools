#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 11:58:41 2020

@author: n.compaire

Modified : G.Sainton - copy function Remove_Tick_Noise_02_BH to 
					Remove_Tick_Noise_02_BH_2_trace to return a Trace object


"""

# import packages

import os
import sys
import copy
import numpy as np

import obspy
from obspy.core import UTCDateTime
from obspy import read, read_inventory
from obspy.core import Stream, Trace
from obspy.signal.cross_correlation import correlate




###############################################################################
# Tick-Noise Removal on VBB 20 sps (02.BH)
###############################################################################

TICK_NOISE_DATADIR = os.path.join(os.curdir,'modules/ticknoiseremoval/')

def detick_noise_02_BH_stream(stream_data):
	"""
	Function with remove the tick of into the InSight data @1Hz
	This function is calling the function with detick each single 
	Trace of the input Stream object.

	---
	INPUT:
		@stream_data: Stream object - 20sps data
	
	OUTPUT:
		@return: deticked 20Hz data 

	"""
	if len(stream_data)>0:
		deticked_stream = Stream()
		for chan in stream_data:
			dt_chan = Remove_Tick_Noise_02_BH_2_trace(chan)
			print(dt_chan)
			deticked_stream.append(dt_chan)
		return deticked_stream
	else:
		print("NONE")
		return None


def Remove_Tick_Noise_02_BH(tr):
	"""
	-- Tick-Noise Removal function --
	Only for 20sps VBB data (02.BH). Remove the estimated tick-noise from the data.
	
	Packages needed : obspy, numpy (as np), copy 
	
	INPUT Variables :
	tr : trace obspy
	Tick_Noise_Datadir : (.../.../) directory of the Tick-Noise files 

	OUTPUT Variables:
	Detick_tr : cleaned timeseries (numpy array)

	"""
	
	
	Compo = tr.stats.channel[2]
	print("Compo : {}".format(Compo))
	
	Tick_Noise = np.load(TICK_NOISE_DATADIR + 'Tick_Noise_' + Compo + '.npy')
	
	fe = tr.stats.sampling_rate
	
	npts_in_tr = tr.stats.npts
	
	nbr_of_s = npts_in_tr/fe
	
	if nbr_of_s % 1 > 0:
		new_nbr = int((np.floor(nbr_of_s) + 1)*fe)
		npts_to_add = new_nbr - npts_in_tr
	else:
		npts_to_add = 0
		new_nbr = copy.deepcopy(npts_in_tr)
	
	temp_data = np.concatenate([tr.data, np.zeros(npts_to_add)])
	
	Long_Tick = np.tile(Tick_Noise, int(new_nbr/fe))
	
	
	corr = correlate(temp_data, Long_Tick, int(fe))
	
	ind_max=np.where(corr==np.max(corr))[0][0] - int(fe)
	
	Roll_Long_Tick = np.roll(Long_Tick, ind_max)
	
	Detick_tr = temp_data - Roll_Long_Tick
	if npts_to_add!=0:
		Detick_tr = Detick_tr[0:-npts_to_add]
	print(tr)
	print("npts_in_tr: ", npts_in_tr)
	print("Detick_tr: ", len(Detick_tr))

	return Detick_tr


def Remove_Tick_Noise_02_BH_2_trace(tr):
	"""
	-- Tick-Noise Removal function --
	Only for 20sps VBB data (02.BH). Remove the estimated tick-noise from the data.
	
	Packages needed : obspy, numpy (as np), copy 
	
	INPUT Variables :
	tr : trace obspy
	Tick_Noise_Datadir : (.../.../) directory of the Tick-Noise files 

	OUTPUT Variables:
	Detick_tr : cleaned timeseries (numpy array)

	"""
	
	
	Compo = tr.stats.channel[2]
	#print("Compo : {}".format(Compo))
	
	Tick_Noise = np.load(TICK_NOISE_DATADIR + 'Tick_Noise_' + Compo + '.npy')
	
	fe = tr.stats.sampling_rate
	
	npts_in_tr = tr.stats.npts
	
	nbr_of_s = npts_in_tr/fe
	
	if nbr_of_s % 1 > 0:
		new_nbr = int((np.floor(nbr_of_s) + 1)*fe)
		npts_to_add = new_nbr - npts_in_tr
	else:
		npts_to_add = 0
		new_nbr = copy.deepcopy(npts_in_tr)
	
	temp_data = np.concatenate([tr.data, np.zeros(npts_to_add)])
	
	Long_Tick = np.tile(Tick_Noise, int(new_nbr/fe))
	
	
	corr = correlate(temp_data, Long_Tick, int(fe))
	
	ind_max=np.where(corr==np.max(corr))[0][0] - int(fe)
	
	Roll_Long_Tick = np.roll(Long_Tick, ind_max)
	
	Detick_tr = temp_data - Roll_Long_Tick
	if npts_to_add!=0:
		Detick_tr = Detick_tr[0:-npts_to_add]
	#print(tr)
	#print("npts_in_tr: ", npts_in_tr)
	#print("Detick_tr: ", len(Detick_tr))
	trace = Trace(Detick_tr)
	trace.stats = tr.stats.copy()
	
	return trace



###############################################################################
# Estimate and Remove the tick-noise in the others cases
###############################################################################

def clean_tic_noise_1s_Automat(ts, dt):
	"""
	This function removes the 1s signal coming from the temperature
	recording at 1sps.
	
	
	From a code by P. Lognonn�/B. Kenda (IPGP). Modified and Translated in python by N. Compaire.
	
	Packages needed : numpy (as np), copy
	
	INPUT variables: ts timeseries
					 dt time step (seconds)
	
	OUTPUT : 
	stack : estimated Tick-Noise
	ts_clean : cleaned timeseries
	
	"""
	p = 0.95
	T = 1 # 1 second
	
	fs=1/dt  # sampling frequency
	
	df_s = 1/T
	df_d = (1/dt)/len(ts)

	r=np.remainder(df_s/df_d,1)
	
	vardf = (1/T)/(1/dt)
	add_npts = int(np.round((1-r)/vardf))
	
	ts_2 = np.concatenate((ts,np.zeros(add_npts)))
	
	uf = copy.deepcopy(ts_2)
	
	# Stack the signal over T #
	N = len(uf)
	NP = int(T*fs)
	factor1 = int((N/NP)*p)
	factor2 = int(N/NP)
	istack = 0
	stack = np.zeros(NP)
	for i in range(factor1):
		uf_mean = uf[i*NP:(i+1)*NP] - np.mean(uf[i*NP:(i+1)*NP])
		for j in range(NP):
			stack[j] = stack[j] + uf_mean[j]
		istack += 1
	stack /= istack
	
	# Remove the stacked signal over moving windows
	
	ts_clean = np.zeros(len(ts_2))
	
	for i in range(factor2):
		for j in range(NP):
			ts_clean[i*(NP)+j] = ts_2[i*(NP)+j] - stack[j]
	
	if add_npts!=0:
		ts_clean = ts_clean[0:-add_npts]
	
	return stack, ts_clean

