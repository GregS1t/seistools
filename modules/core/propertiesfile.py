#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 15:22:39 2019

@author: greg
"""

#==============================================================================
# Return the longname from location and channel code extracted 
#    data are stored into channel.properties file
#

def get_long_name_from_prop_list(longname_tab, channel, locid):
	'''
	Function which return a longname (eg: VBB_1_SCI_VEL_HG) giving a channel 
	and a location.
	Informations are read from the properties file.
	
	INPUT:
		@longname_tab: Pandas Structure with all the properties
		@channel: channel
		@locId: location
		@sps : sample per second
	OUTPUT:
		@return: string with the long name.
	'''
	import datetime
	pdlgname = longname_tab.loc[(longname_tab["chan"] == channel) \
							 & (longname_tab["locid"] ==locid), 'longnamesps']
	if len(pdlgname)==0:
		print(str(datetime.datetime.now()) + " - "+ locid\
				+"_"+channel+" not found in the file channel.properties.")
		longname = locid+"_"+channel
	else:
		longname = pdlgname.apply(str).iloc[0]
	#print(longname)
	return longname

#==============================================================================
# Return units for a given channel
#
def get_unit_from_prop_list(channel_properties_list, channel, locid):
	'''
	Function which return units (eg: VBB_1_SCI_VEL_HG) giving a channel 
	and a location.
	Informations are read from the properties file.
	
	INPUT:
		@longname_tab: Pandas Structure with all the properties
		@channel: channel
		@locId: location
		@sps : sample per second
	
	OUTPUT:
		@return: string with the long name.
	'''
	import datetime
	units = channel_properties_list.loc[(channel_properties_list["chan"] == channel) \
						  & (channel_properties_list["locid"] == locid), 'Units']
	if len(units)==0:
		print(str(datetime.datetime.now()) + " - "+ locid\
				+"_"+channel+" not found in the file channel.properties.")
		yunits = "Counts"
	else:
		yunits = units.apply(str).iloc[0]
	return yunits


def get_chanid4sensor(channel_properties_list, sensor):
	'''
	Function to return location code, channel code and longname for a given 
	sensor
	----
	INPUT
		@channel_properties_list: list of properties for every known channels
		@sensor: sensor type = VBBPOS, VBBVel, VBBTemp...
	
	----
	OUTPUT
		@pdlgsensor : pandaframe list of channels with the same input sensor
	
	'''
	pdlgsensor = channel_properties_list.loc[\
							(channel_properties_list["sensor"] == sensor), \
							('longnamesps', 'locid', 'chan')]
	
	return pdlgsensor

