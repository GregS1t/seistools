#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 14:47:55 2019

@author: greg
"""



def decomposeFIR(sps_start, sps_final, fir_list = [5,4,2]):
	"""
	Function which decompose the ratio between two sps
	@sps_start (float): sampling rate at the beginning
	@sps_stop (float) : sampling rate at the end
	
	If the ratio is not integer, function returns an empty list
	
	@return (list) : list of FIR coef
	"""
	fir_list.sort(reverse=True)
	fir_2_apply = []
	sps_ratio = sps_start / sps_final
	if sps_ratio == int(sps_ratio):
		print("Ratio entier: ")
	else:
		print("Cannot decompose ratio: ", sps_ratio)
		return []
		
	i = 0
	ratio = sps_ratio
	while ratio > 1 : 
		if ratio<fir_list[i]:
			i+=1
		else:
			ratio = ratio/fir_list[i]
			fir_2_apply.append(fir_list[i])

	return fir_2_apply	# return empty if decomposition is impossible


def read_filter_file(file, print_coef=False):
	"""
	Function to read filter file
	@file : input text file
	"""
	import os
	
	#print("File: ", file)
	Lines = None
	if os.path.exists (file) :
		txt = open(file)
		Lines = txt.readlines()
	else :
		print("EBox filter not found.")
		exit (-1)

	Filtre = []
	
	if Lines is not None:
		for line in Lines  :
			if (line.find ("B061F09 ") == 0 ) :
				fact = float (line[17:])
				Filtre.append (fact)
			elif isinstance(line,float)  :
				fact = float (line)
				Filtre.append (fact)

		if print_coef == True:
			print(Filtre)
		return Filtre
	return None


def read_fir_config_file(fir_configfile):
	'''
	Function which read the FIR config file and return a dict structure 
	with informations
	
	@fir_configfile: str - input filename
	@return: dict - information about filter (name, type, gain, pathfile)
	'''
	#open and read XML file
	#xml lib to manage xml files 
	from lxml import etree 
	try:
		tree = etree.parse(fir_configfile) 
	except etree.XMLSyntaxError:
		print("Invalid XML file.")
		return None   
	except IOError:
		print('File not found: ', fir_configfile)
		return None
	except:
		print('Error')
		
	#parse the structure to extract informations
	root = tree.getroot()
	config_list_dict = []
	for el in root: 
		dict_fir = {}
		for sub_el in el:
			dict_fir[str(sub_el.tag)] = sub_el.text
		config_list_dict.append(dict_fir)
	return config_list_dict


def decimate_FIR_scipy(trace, decim_fact, trim = False):
	"""
	Function to decimate data with the scipy filter. It's a 30 points filter
	----
	INPUT:
		@trace: Trace object 
		@decim_fact: (int) decimation factor
		@trim: (boolean) to tell the code to trim data after decimation
	OUTPUT:
		@return the decimated trace
	
	"""
	from scipy.signal import decimate
	cptrace      = trace.copy()
	cptrace.data = decimate(trace.data, decim_fact, ftype="fir")
	cptrace.stats.sampling_rate = trace.stats.sampling_rate / float(decim_fact)
	
	return cptrace



def FIR_filter(trace, filt2apply, decim_fact, decim = False, verbose = False):
	"""
	Function which apply a FIR filter on input data. 
	decim_fact must be provided. If no gain_fact provided, the default one has
	value equal to 1.
	Data are automatically trimmed to remove the beginnning and the end of 
	the signal since no filtering is made.
	----
	INPUT:
		@trace  (Trace)    : - input data  (Trace is an OBPSY object)
		@filt2apply (list) : - input FIR coeff
		@decim_fact (int)  : - decimation factor in case of simultaneous decim
		@gain_fact (int)   : - gain factor
		@decim (boolean)   : - boolean to tell if decimation is done or not
	----
	OUTPUT:
		@return  (Trace)   : - data corrected from FIR
	"""
	import numpy as np
	
	len_filt2apply = len(filt2apply)
	#make a copy of the trace not to overwrite the original
	filtered_trace = trace.copy()
	
	if len(filt2apply) !=0:
		half_len_filt = int((len_filt2apply-1)/2)
		filter_data = filtered_trace.copy()

		#apply the filter
		for i in range(half_len_filt, (len(filtered_trace.data)-\
											 half_len_filt)-1):
				filter_data.data[i] = np.dot(filt2apply,\
					trace.data[i-half_len_filt:i+half_len_filt+1])
		
		#trim the new signal on the edges
		starttime = filter_data.stats.starttime
		npts = filter_data.stats.npts
		sampling_rate =filter_data.stats.sampling_rate
		delta = filter_data.stats.delta
		
		t = np.arange(starttime, starttime + npts/sampling_rate, delta)
		#print("vector t:", t[0], t[-1])
		t_trim = t[half_len_filt:-half_len_filt]
		#print("len(t_trim)",len(t_trim))
		
		filtered_trace.data=filter_data[half_len_filt:-half_len_filt]
		filtered_trace.stats.starttime = t_trim[0]
		filtered_trace.stats.npts = len(t_trim)
		if verbose:
			print("be4 decim", filtered_trace.stats.starttime, filtered_trace.stats.endtime, \
					 filtered_trace.stats.sampling_rate,  filtered_trace.stats.npts)
		
		if decim: 
			#decimate the filtered trace
			#no filter applied.
			if verbose: 
				print("    - Decimation applied with factor {}".format(int(decim_fact)))
			filtered_trace_d = filtered_trace.copy()
			filtered_trace_d.decimate(int(decim_fact), no_filter=False, strict_length=False)
			if verbose:
				print("after decim:", filtered_trace_d.stats.starttime, filtered_trace_d.stats.endtime, \
						filtered_trace_d.stats.npts, filtered_trace_d.stats.sampling_rate)
			
		return filtered_trace_d 
	else:
		print("No FIR filter coefficients")
		return None

def FIR_filter_as(as_channel, filt2apply, decim_fact, decim = False, verbose = False):
	"""
	Function which apply a FIR filter on input data. 
	decim_fact must be provided. If no gain_fact provided, the default one has
	value equal to 1.
	Data are automatically trimmed to remove the beginnning and the end of 
	the signal since no filtering is made.
	----
	INPUT:
		@trace  (Trace)    : - input data  (Trace is an OBPSY object)
		@filt2apply (list) : - input FIR coeff
		@decim_fact (int)  : - decimation factor in case of simultaneous decim
		@gain_fact (int)   : - gain factor
		@decim (boolean)   : - boolean to tell if decimation is done or not
	----
	OUTPUT:
		@return  (Trace)   : - data corrected from FIR
	"""
	#import numpy as np
	
	if len(filt2apply) !=0:
		#make a copy of the trace not to overwrite the original
		as_filtered_trace = as_channel.copy()
		as_filtered_trace.traces = FIR_filter(as_filtered_trace.traces, filt2apply, \
									decim_fact, decim = decim, verbose = verbose)
		
		if as_filtered_trace.traces is not None:
			as_filtered_trace.npts = as_filtered_trace.traces.stats.npts
			as_filtered_trace.sps = as_filtered_trace.traces.stats.sampling_rate
			as_filtered_trace.starttime = as_filtered_trace.traces.stats.starttime
			as_filtered_trace.endtime = as_filtered_trace.traces.stats.endtime
			return as_filtered_trace 
		else:
			print("No FIR filter coefficients")
			return None
	else:
		print("No FIR filter coefficients")
		return None



def get_filemane(config_list_dict, factor, gain = 1):
	"""
	Function to recover the filename giving the factor and the gain
	
	INPUT:
	@config_list_dict : dict to reach the filename with gain and factor
	@factor           : factor of the FIR
	@gain             : gain (defaultly = 1)
	
	OUTPUT:
	@filename         : string with the filename
	
	"""
	filename = None
	for fird in config_list_dict:
		if fird.get("type") == str(factor) and fird.get("gain") == str(gain):
			filename = fird.get("filename")
	return filename

#==============================================================================
#    Main part

if __name__ == "__main__":
	print("File to manage filter functions...")
	
	fir_configfile = "../../configuration/fir_filters.xml"
	config_list_dict = read_fir_config_file(fir_configfile)
	print(config_list_dict)
	fir_list = []
	for fir in config_list_dict:
		if fir.get("gain") == str(1):
			#print("Only FIR with gain = 1 kept here...")
			print("FIR {} filter saved in {} (gain = {})".format(fir.get("type")\
		 , fir.get("filename"), fir.get("gain")))
			fir_list.append(int(fir.get("type")))
	fir_list.sort(reverse=True)
	
	filename = get_filemane(config_list_dict, 2)
	print(filename)
	filename = get_filemane(config_list_dict, 4)
	print(filename)
	print("fir_list to be transmitted in 'decomposeFIR' function: ", fir_list)
	
	
	
