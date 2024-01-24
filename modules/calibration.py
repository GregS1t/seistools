#
#
#	   @name	: calibration.py
#	   @author  : Greg Sainton (sainton@ipgp.fr)
#	   @purpose : module dedicated to read an get informations from the dataless
#
#
# 2019/07/04 - GS - Update functions to add startime and and time
#
#
###############################################################################

#obspy libs

import obspy.signal
from obspy.signal.invsim import paz_to_freq_resp

#numpy lib
import numpy as np

# matplotlib lib to plot things
import matplotlib.pyplot as plt


#==============================================================================

def funct_ML(coeffList, x):
	"""
	MAC LAURIN function
	----
	INPUT
		@coeffList : List - coefficients of the polynomial
		@x: float - input value
	----
	OUTPUT:
		@p_x: float - value of the polynomial
	"""
	p_x = 0
	for i in range (0, len(coeffList)):
		p_x+=coeffList[i]*x**i
	return p_x

def get_coeff_poly_resp_stage(inv, station_name, loc_name, channel_name, 
							  starttime=None, endtime=None):
	"""
	Function to return polynomial coefficients from a given location/channel
	----
	INPUT:
		@inv (inventory)      : inventory structure made with dataless
		@chan_name (string)   : channel name (eg: "MHU")
		@station_name (string): station name (eg: "KAWA2")

	----
	OUPUT
		@return: list with all the coefficients

	"""
	try:
		response_chan_sta = inv.select(channel=channel_name, location=loc_name, 
									  station=station_name, 
									  starttime=starttime, 
									  endtime=endtime)[0][0][0]
	except:
		#print("No response found with {}/{}/{}".format(loc_name, channel_name, 
		#												station_name))
		return None
	else:
		#print(response_chan_sta)
		resp = response_chan_sta.response
		stages = resp.response_stages
		for stage in stages:
			if isinstance(stage, obspy.core.inventory.response.PolynomialResponseStage):
				return stage._coefficients


def get_counts_2_Volts(inv, sta, loc, chan, starttime=None, endtime=None):
	"""
	Function to extract Gain from Counts to Volt for a given (station, channel)
	couple from an Inventory object. Inventory object instanciated by ObsPy 
	while reading a dataless file
	----
	INPUT:
		@inv : Inventory object
		@sta : string - station name
		@chan : string - channel name
	----
	OUTPUT
		@return : float - value of the gain

	"""
	try:
		response_chan_sta = inv.select(channel=chan, location=loc, \
								 station=sta, starttime=starttime, \
								 endtime=endtime)[0][0][0]
	except IndexError:
		print("It seems that there is an error in you channel or station name.\
					Please check")
		return None
	else:
		resp = response_chan_sta.response
		for stage in resp.response_stages:

			if stage.input_units == "V" and stage.output_units == "COUNTS":
				gain_C2V = stage.stage_gain
		return gain_C2V


def get_FIR_gains(inv, sta, loc, chan, starttime=None, endtime=None):
	"""
	Function to get the gains from FIR stages.

	----
	INPUT
		@inv : Inventory object
		@sta : string - station name
		@loc : string - location code
		@chan : string - channel name
	OUTPUT
	----
		@return : list - FIR gains
	"""
	try:
		response_chan_sta = inv.select(channel=chan, location=loc, \
								 station=sta, starttime=starttime, \
								 endtime=endtime)[0][0][0]
	except IndexError:
		return None
	else:
		response_chan_sta = inv.select(channel=chan, location=loc, \
								 station=sta, starttime=starttime, \
								 endtime=endtime)[0][0][0]
		resp = response_chan_sta.response
		stages = resp.response_stages
		numerators_list = []

		for stage in stages:
			if isinstance(stage, obspy.core.inventory.response.FIRResponseStage):
				numerators_list.append(stage.stage_gain)
		return numerators_list


def get_FIR_numerators(inv, sta, chan, loc, dec_factor, starttime=None, endtime=None):
	"""
	Function to get the gains from FIR stages.

	----
	INPUT:
		@inv: Inventory object
		@loc: location code
		@sta: string - station name
		@chan: string - channel name
		@dec_factor: decimation factor
		@starttime: starttime of the signal 
		@endtime: endtime of the signal
	----
	OUPUT
		@return: list - FIR numerators
	"""
	try:
		response_chan_sta = inv.select(channel=chan, location=loc, \
								 station=sta, starttime=starttime, \
								 endtime=endtime)[0][0][0]
	except IndexError:
		print("It seems that there is an error in you channel or station name. Please check")
		return None
	else:
		response_chan_sta = inv.select(channel=chan, location=loc, \
								 station=sta, starttime=starttime, \
								 endtime=endtime)[0][0][0]
		resp = response_chan_sta.response
		stages = resp.response_stages
		numerators_list = []
		for stage in stages:
			#print(stage)
			#Look for FIR Stages
			if isinstance(stage, obspy.core.inventory.response.CoefficientsTypeResponseStage):
				# Filters with decimation factors 2 or 4 have gain = 1.
				if stage.decimation_factor == dec_factor:
					#print("stage.decimation_factor", stage.decimation_factor)
					numerators_list= stage._numerator
			#if isinstance(stage, obspy.core.inventory.response.FIRResponseStage):
			#	print(stage.__dict__)
		return numerators_list


def get_units(inv, sta, loc, chan, starttime=None, endtime=None):
	"""
	Function to retrieve the physical unit of a given channel
	@inv : Inventory object
	@sta : string - station name
	@chan : string - channel name
	@return : string - physical unit of the channel

	"""
	try:
		response_chan_sta = inv.select(channel=chan, location=loc, \
								 station=sta, starttime=starttime, \
								 endtime=endtime)[0][0][0]
	except IndexError:
		print("GET_UNITS: It seems that there is an error in you channel or station name. Please check")
		return None
	else:
		response_chan_sta = inv.select(channel=chan, location=loc, \
								 station=sta, starttime=starttime, \
								 endtime=endtime)[0][0][0]
		resp = response_chan_sta.response
		stages = resp.response_stages

		for stage in stages:

			if isinstance(stage, obspy.core.inventory.response.PolynomialResponseStage):
				#print("Polynomial response...")
				stage_dict={}
				stage_dict["type"] = "Poly"
				stage_dict["xunit"] = stage.output_units
				stage_dict["yunit"] = stage.input_units
				stage_dict["yunit_desc"] = stage.input_units_description

			if isinstance(stage, obspy.core.inventory.response.PolesZerosResponseStage):
				#print("Poles and Zeros response...")
				stage_dict={}
				stage_dict["type"] = "PAZ"
				stage_dict["input_units"] = stage.output_units
				stage_dict["output_units"] = stage.input_units
				stage_dict["output_units_desc"] = stage.input_units_description

		return stage_dict


def print_dip_azimuth(inv, sta, loc, chan, starttime=None, endtime=None):

	try:
		response_chan_sta = inv.select(channel=chan, location=loc, \
								 station=sta, starttime=starttime, \
								 endtime=endtime)[0][0][0]
	except IndexError:
		#print("PRINT_DIP_AZIMUTH: It seems that there is an error in you channel or station name. Please check !")
		return None
	else:
		#print("dip="+str(response_chan_sta.dip))
		#print("azimuth="+str(response_chan_sta.azimuth))
		return {'dip':str(response_chan_sta.dip), 'azimuth':str(response_chan_sta.azimuth)}


def get_dip_azimuth(inv, sta, loc, chan, starttime=None, endtime=None):
	"""
	Function to return DIP and AZIMUTH for a given channel
	----
	INPUT:
		@inv        : inventory object
		@sta        : (str) station code (eg. 'ELYSE')
		@location   : (str) location code (eg. '02')
		@chan       : (str) channel code ('MHU')
		@starttime  : (UTCDateTime) starttime of the time serie
		@endtime    : (UTCDateTime) end of the time serie
	----
	OUTPUT:
		dictionnary with two keys : "dip" and "azimuth"
	
	"""
	
	try:
		response_chan_sta = inv.select(channel=chan, location=loc, \
								 station=sta, starttime=starttime, \
								 endtime=endtime)[0][0][0]
	except IndexError:
		#print("PRINT_DIP_AZIMUTH: It seems that there is an error in you channel or station name. Please check !")
		return None
	else:
		#print("dip="+str(response_chan_sta.dip))
		#print("azimuth="+str(response_chan_sta.azimuth))
		return {'dip':str(response_chan_sta.dip), 'azimuth':str(response_chan_sta.azimuth)}


def print_stage(inv, sta, loc, chan, starttime = None, endtime = None):
	"""
	Function to display all the stage of a channel in the dataless.
	----
	INPUT
		@inv: Inventory object
		@sta: string - station name
		@loc: location code
		@chan: string - channel name
		@starttime: starttime of the signal
		@endtime: endtime of the signal
	"""

	print("===========================")
	print("Stages for", loc," - ", chan)
	print("__________________________")

	try:
		response_chan_sta = inv.select(channel=chan, location=loc, \
								 station=sta, starttime=starttime, \
								 endtime=endtime)[0][0][0]
	except IndexError:
		print("It seems that there is an error in you channel or station name. Please check")
	else:
		resp = response_chan_sta.response

		if isinstance(resp.instrument_sensitivity, obspy.core.inventory.response.InstrumentSensitivity):
				print("Sensitivity object: ")
				print("Input units: ", resp.instrument_sensitivity.input_units)
				print("Output units: ", resp.instrument_sensitivity.output_units)


		for stage in resp.response_stages:

			print("	 stage sequence number: ",stage.stage_sequence_number)
			print("	 +++++++++++++++++++++++++++++++++++++++++++++++++++++")
			print(stage)

			#if stage.input_units == "M/S**2" and stage.output_units == "V":
				#stage.stage_gain = -stage.stage_gain
			if isinstance(stage, obspy.core.inventory.response.PolesZerosResponseStage):
				poles = resp.get_paz()._poles
				zeros = resp.get_paz()._zeros
				print("poles: {} \n zeros: {}".format(poles, zeros))



def plot_maclaurin(coeffML, start, end, step, xtitle, ytitle, title):
	"""
	Function to plot the Mac Laurin function giving the input coefficients.
	The size of the list of coefficients set the order of the polynom.
	----
	INPUT:
		@coeffML: list of coefficients to be used for the Mac Laurin function
		@start: float - lower bound of x
		@end: float - upper bound of x
		@step: float - step
		@xtitle: string - title of the x axis
		@ytitle: string - title of the y axis
		@title: string - title of the plot

	"""

	volt = np.arange(start, end, step)
	temp = []
	for v in volt:
		t= funct_ML(coeffML, v)
		temp.append(t)

	fig, ax = plt.subplots()
	ax.plot(volt, temp)
	ax.set(xlabel=xtitle, ylabel=ytitle, title=title)
	ax.grid()


def plot_responses(inv, channel_list, output, min_freq, starttime = None, endtime = None):
	"""
	Function to plot the response of list of channels
	@inv: dataless file read by Obspy
	@channe_list: dict containing the name of the channel with the name of the station.
	@output: string - VEL, DISP or ACC
	@min_freq: minimum frequency
	@return : figs and axs with plot structure
	"""
	from itertools import cycle
	import matplotlib.pyplot as plt
	import numpy as np

	resp_list=[]
	for channel in channel_list:
		chan = channel.get("channel")
		sta  = channel.get("station")
		loc  = channel.get("location")
		lgname = channel.get("longname")
		detector = channel.get("detector")
		sps = channel.get("sps")
		gaintype = channel.get("gaintype")
		out = channel.get("output")
		try:
			response_chan_sta = inv.select(channel=chan, location=loc , \
						station=sta, starttime = None, endtime = None)[0][0][0]
		except IndexError:
			#with dataless_message:
			print("It seems that there is an error in you channel or station name. Please check")
		else:
			dict_resp ={}
			dict_resp["channel"]=chan
			dict_resp["station"]=sta
			dict_resp["location"]=loc
			dict_resp["longname"]=lgname
			dict_resp["detector"]=detector
			dict_resp["response"]=response_chan_sta
			dict_resp["gaintype"]=gaintype
			dict_resp["output"]= out
			dict_resp["sps"]= sps
			print(dict_resp)
			resp_list.append(dict_resp)

	#figs.subplots_adjust(hspace=0.2)
	cycol = cycle('bgrcmk')
	cpx_resp_list=[]
	freq_list=[]
	longname=[]
	nyquist_l = []
	for i in range(0, len(resp_list)):

		resp = resp_list[i].get("response").response
		sampling_rate = None
		min_freq = min_freq
		output = output

		# detect sampling rate from response stages
		if sampling_rate is None:
			for stage in resp.response_stages[::-1]:
				if (stage.decimation_input_sample_rate is not None and
						stage.decimation_factor is not None):
					sampling_rate = (stage.decimation_input_sample_rate /
									 stage.decimation_factor)
					break
			else:
				msg = ("Failed to autodetect sampling rate of channel from "
					   "response stages. Please manually specify parameter "
					   "`sampling_rate`")
				raise Exception(msg)
		if sampling_rate == 0:
			#with response_message:

			msg = "Can not plot response for channel "
			print (msg, resp_list[i].get("longname"))
		else:
			t_samp = 1.0 / sampling_rate
			nyquist = sampling_rate / 2.0
			nyquist_l.append(nyquist)
			nfft = sampling_rate / min_freq

			try:
				cpx_resp, freq = resp.get_evalresp_response(
					t_samp=t_samp, nfft=nfft, output=output)
			except NotImplementedError:
				#with response_message:
				print("No polynomial responses for ", resp_list[i].get("longname"))
			else:
				cpx_resp_list.append(cpx_resp)
				freq_list.append(freq)
				longname.append(resp_list[i].get("longname"))

	if len(cpx_resp_list)>1:
		figs, axs = plt.subplots(len(cpx_resp_list), 2, figsize=(11.69,8.27), sharex='col')
		for i in range(0, len(cpx_resp_list)):
			color = next(cycol)
			#Response

			axs[i,0].loglog(freq_list[i], abs(cpx_resp_list[i]), lw=1, \
							  color=color , label= output+"_"+chan)

			title = ''
			if resp_list[i].get("detector") is not None:
				title += resp_list[i].get("detector")+"/"
			if resp_list[i].get("output") is not None:
				title += resp_list[i].get("output")+"/"
			if resp_list[i].get("sps") is not None:
				title += str(resp_list[i].get("sps"))
			title += "("+output+")"

			axs[i,0].set_title(title, fontsize=6)
			axs[i,0].set_xlabel('Freq (Hz)', fontsize=6)
			axs[i,0].set_ylabel('Amplitude', fontsize=6)
			#axs[i,0].set_xlim(0, 0.8*nyquist_l[i])
			axs[i,0].grid(True,which="both",ls="-")
			#Phase
			phase = np.angle(cpx_resp_list[i], deg=False)
			axs[i,1].set_title('Phase for '+ 'Transfer function for ' + title, fontsize=6)
			axs[i,1].semilogx(freq_list[i], phase, lw=1, color=color, linewidth = 0.5)
			axs[i,1].grid(True,which="both",ls="-")
			axs[i,1].set_xlabel('Freq (Hz)', fontsize=6)
			axs[i,1].set_ylabel('Phase [Rad]', fontsize=6)
			axs[i,1].set_yticks([0, np.pi / 2, np.pi, 3 * np.pi / 2, 2 * np.pi])
			axs[i,1].set_yticklabels(['$0$', r'$\frac{\pi}{2}$', r'$\pi$', r'$\frac{3\pi}{2}$', r'$2\pi$'])

	else:
		figs, axs = plt.subplots(nrows=1, ncols=2, figsize=(11.69,8.27))
		#Response
		#lw = 3
		#lines = axs[0].loglog(freq_list[0], abs(cpx_resp_list[0]),lw=1, label= output+"_"+chan, linewidth = 0.5)

		title = ''
		if resp_list[0].get("detector") is not None:
			title += resp_list[0].get("detector")+"/"
		if resp_list[0].get("output") is not None:
			title += resp_list[0].get("output")+"/"
		if resp_list[0].get("sps") is not None:
			title += str(resp_list[0].get("sps"))
		title += "("+output+")"

		axs[0].set_title(title, fontsize=6)

		axs[0].set_xlabel('Freq (Hz)')
		axs[0].set_ylabel('Amplitude')
		axs[0].grid(True,which="both",ls="-")
		#axs[0].set_xlim(0, 0.8*nyquist_l[i])
		#Phase
		phase = np.angle(cpx_resp_list[0], deg=False)
		axs[1].set_title('Phase for '+ title, fontsize=6)
		axs[1].semilogx(freq, phase, lw=1)
		axs[1].grid(True,which="both",ls="-")
		axs[1].set_xlabel('Freq (Hz)')
		axs[1].set_ylabel('Phase [Rad]')
		axs[1].set_yticks([0, np.pi / 2, np.pi, 3 * np.pi / 2, 2 * np.pi])
		axs[1].set_yticklabels(['$0$', r'$\frac{\pi}{2}$', r'$\pi$', r'$\frac{3\pi}{2}$', r'$2\pi$'])
	return figs, axs

#
# POLES AND ZEROS FUNCTIONS
#
#################################################
def get_paz(inv, chan_name, loc_name, station_name, starttime=None, endtime=None):
	"""
	Function to return the poles and zeros of a channel
	@inv: inventory structure made with dataless
	@chan_name: string containing a channel name (eg: "MHU")
	@station_name: string containing a station name (eg: "KAWA2")

	@return: list with the poles and zeros
	"""

	try:
		response_chan_sta = inv.select(channel=chan_name, location=loc_name, \
				 station=station_name, starttime=None, endtime=None)[0][0][0]
	except:
		print("Error while reading the response")
	else:
		resp = response_chan_sta.response
		try:
			poles = resp.get_paz()._poles
			zeros = resp.get_paz()._zeros
		except Exception:
			print ("Channel {}: No poles nor zeros founds".format(chan_name))
			return None
		else:
			return [poles, zeros]

#
# Plot Pole and Zero function
######################################################
def plot_paz(inv, sta, loc, chan, starttime=None, endtime=None):
	try:
		response_chan_sta = inv.select(channel=chan, location=loc,\
						 station=sta, starttime=None, endtime=None)[0][0][0]
	except IndexError:
		print("It seems that there is an error in you channel or station name. Please check")
	else:
		resp = response_chan_sta.response

		try:
			poles = resp.get_paz().poles
			zeros = resp.get_paz().zeros
		except:
			print("No PolesZerosResponseStage found.")
		else:
			#print(poles)
			#print(zeros)
			scale_fac = 0.4

			h, f = paz_to_freq_resp(poles, zeros, scale_fac, 0.005, 16384, freq=True)

			plt.figure()
			plt.subplot(121)
			plt.loglog(f, abs(h))
			plt.xlabel('Frequency [Hz]')
			plt.ylabel('Amplitude')
			plt.grid(True)
			plt.subplot(122)
			phase = 2 * np.pi + np.unwrap(np.angle(h))
			plt.semilogx(f, phase)
			plt.xlabel('Frequency [Hz]')
			plt.ylabel('Phase [rad]')
			# ticks and tick labels at multiples of pi
			plt.yticks(
				[0, np.pi / 2, np.pi, 3 * np.pi / 2, 2 * np.pi],
				['$0$', r'$\frac{\pi}{2}$', r'$\pi$', r'$\frac{3\pi}{2}$', r'$2\pi$'])
			plt.ylim(-0.2, 2 * np.pi + 0.2)
			# title, centered above both subplots
			plt.suptitle('Frequency Response of ' + chan)
			# make more room in between subplots for the ylabel of right plot
			plt.subplots_adjust(wspace=0.3)
			plt.show()

def plot_calibration(inv, channel_name, loc_name, station_name, starttime=None, endtime=None):
	"""
	Function to plot all the calibration curve available among checked channels
	----
	INPUT:
		@inv: dataless read by ObsPy
		@channel_name: string - channel code
		@loc_name: string - location code
		@station_name: string - station code
		
	"""

	Vmin = -25.0
	Vmax = 25.0
	Vstep = 0.001
	stage_data = []

	response_chan_sta = inv.select(channel=channel_name, location = loc_name, \
				station=station_name, starttime=None, endtime=None)[0][0][0]
	resp = response_chan_sta.response
	stages = resp.response_stages

	for stage in stages:
		if isinstance(stage, obspy.core.inventory.response.PolynomialResponseStage):
			#print(stage)
			stage_dict={}
			stage_dict["coefML"] = stage._coefficients
			stage_dict["xunit"] = stage.output_units
			stage_dict["yunit"] = stage.input_units
			stage_dict["yunit_desc"] = stage.input_units_description
			stage_data.append(stage_dict)
	nb_plot = len(stage_data)
	if nb_plot>1:
		fig, ax = plt.subplots(nb_plot, 1, figsize=(10,7))
		fig.subplots_adjust(hspace=0.5)
		for i in range(0, nb_plot):
			volt = np.arange(Vmin, Vmax, Vstep)
			temp = []
			for v in volt:
				t= funct_ML(stage_data[i].get("coefML"), v)
				temp.append(t)
			ax[i].plot(volt, temp)
			ax[i].set_title("Response from {} to {}.".format(stage_data[i].get("xunit"), stage_data[i].get("yunit_desc"), fontsize=6))
			ax[i].set_ylabel(ylabel=stage_data[i].get("yunit"), fontsize=6)
			ax[i].set_xlabel(xlabel=stage_data[i].get("xunit"), fontsize=6)

			ax[i].grid()
	else:
		fig, ax = plt.subplots(figsize=(10,7))
		volt = np.arange(Vmin, Vmax, Vstep)
		temp = []
		for v in volt:
			t= funct_ML(stage_data[0].get("coefML"), v)
			temp.append(t)
		ax.plot(volt, temp)
		ax.set_title("Response from {} to {}.".format(stage_data[0].get("xunit"), stage_data[0].get("yunit_desc"), fontsize=6))
		ax.set_ylabel(ylabel=stage_data[0].get("yunit"), fontsize=6)
		ax.set_xlabel(xlabel=stage_data[0].get("xunit"), fontsize=6)
		ax.grid()
	return fig, ax

def patch2CalibRotated(chan):
	"""
	Patch add to be able to invert calibration and 
	rotation process. Since there is no transfert 
	function for channels after rotation, we decided 
	to apply the following approximation:
		TF(**N) -> TF(**U)
		TF(**E) -> TF(**V)
		TF(**Z) -> TF(**W)
		
	Parameters:
		@chan (string)  = input channel name
	
	Return:
		@return (string) = modified channel name

	@author  : sainton@ipgp.fr
	@date    : aug 29th, 19'
	@version : 1.0
	"""
	
	chandict = {"N":"U","E":"V","Z":"W"}
	if chan[-1] in ["N", "E", "Z"]:
		newchan = chan[:-1]+chandict.get(chan[-1])
		#print("WARNING : Approximation made on TF({}) -> TF({})".format(chan,newchan))
		return newchan
	else:
		return chan


def hack_VBB_Pos_stage_gain(inv, channel_properties_list, verbose= False):
	"""
	Since obspy does not correct properly the stage gain,
	this function is used to correct the stage gain 
	For input units M/S**2, M/S and M for every Pos channels
	
	-> stage_gain is divided by 4
	
	In Pos channels -> only M/S**2 stages
	In Vel channels -> only M/S	stages
	
	----
	INPUT: 
		@inv: inventory object (from read_inventory command)
	----
	OUPUT:
		@inv: inventory object (from read_inventory command)
	
	"""
	from modules.core.propertiesfile import get_chanid4sensor
	result = get_chanid4sensor(channel_properties_list, 'VBBPOS')

	for chan in inv._networks[0]._stations[0].channels:
		for i in range(0, len(result)): 
			if chan._location_code+"_"+chan._code == \
					result["locid"].iloc[i]+"_"+result["chan"].iloc[i]:
				if verbose:
					print(result["longnamesps"].iloc[i])
				for stage in chan.response.response_stages:
					if stage.input_units in ["M/S**2", "M/S", "M"]:
						if verbose:
							print(stage.input_units)
							print("Before: ", stage.stage_gain) 
						stage.stage_gain = stage.stage_gain/4.
						if verbose:
							print("After: ", stage.stage_gain) 
	return inv



def hack_VBB_Pos_FIR_gain_2_one(inv, channel_properties_list, verbose= False):
	"""
	Since obspy does not correct properly the stage gain,
	this function is used to correct the stage gain 
	For input units M/S**2, M/S and M for every Pos channels

	-> stage_gain is divided by 4
	
	In Pos channels -> only M/S**2 stages
	In Vel channels -> only M/S    stages
	
	----
	INPUT: 
		@inv: inventory object (from read_inventory command)
	----
	OUPUT
		@inv: inventory object (from read_inventory command)

	"""
	from modules.core.propertiesfile import get_chanid4sensor
	from obspy.core.inventory.response import FIRResponseStage

	result = get_chanid4sensor(channel_properties_list, 'VBBPOS')

	for chan in inv._networks[0]._stations[0].channels:
		for i in range(0, len(result)): 
			if chan._location_code+"_"+chan._code == \
					result["locid"].iloc[i]+"_"+result["chan"].iloc[i]:
				if verbose:
					print(result["longnamesps"].iloc[i])
				for stage in chan.response.response_stages:
					print(stage)
					if isinstance(stage, FIRResponseStage):
						if stage.stage_gain > 1:
							stage.stage_gain = 1.0
							print(stage.stage_gain)
	return inv


if __name__ == "__main__":
	print("Welcome in the calibration module.")