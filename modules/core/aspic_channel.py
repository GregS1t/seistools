#import sys

#import obspy.core
from obspy import Stream
import fnmatch
import sys

sys.path.insert(0, '../..')


from modules.calibration import get_counts_2_Volts, get_units, funct_ML, \
						get_coeff_poly_resp_stage, get_FIR_gains
from modules.plot.aspic_plot import plot_options_process
from modules.io.read_xml import read_XML_script
from modules.core.propertiesfile import get_long_name_from_prop_list, \
										get_unit_from_prop_list

from modules.ticknoiseremoval.Tick_Noise_Removal import Remove_Tick_Noise_02_BH_2_trace
#from MarsConverter import MarsConverter


# Define path and filename to reach the config file in Ampere
conf_dir       = './configuration/'
amp_calib_file = 'config_calib_channels_ampere.xml'
path_tickdir   = "./modules/ticknoiseremoval/"

#################################################################
# Function to seach a channel into a table

def search(locidchan, All_calib_list):
	"""
	@locidchan: string format : eg "08/UEA"
	@All_calib_list: list of dict eg:
		{'id': '08/LEA', 'location': '08', 'channel': 'LEA', 'expression': 0.001518230532679, 'unit': 'mA'}
		{'id': '08/UEA', 'location': '08', 'channel': 'UEA', 'expression': 0.001518230532679, 'unit': 'mA'}

	@return: dict: eg {'id': '08/UEA', 'location': '08', 'channel': 'UEA', 'expression': 0.001518230532679, 'unit': 'mA'}
	"""

	for chan in All_calib_list:
		if locidchan in chan.get("id"):
			return chan


class AspicChannel(Stream):

	from obspy import Trace

	network  = None
	station  = None
	location = None
	channel  = None
	longname = None # readable name completed with channel.properties
	sps      = None # number in Herz
	npts     = None # number of points in the time serie
	output   = None # Vel or Pos completed using the longname
	gaintype = None # HG of LG completed using the longname
	mode	 = None # Sci or Eng completed using the longname
	detector = None # VBB, SP... completed using the longname
	orientation = None #U, V, W   completed using the longname
	starttime = None  # Always in UTCDateTime
	endtime = None # Always in UTCDateTime
	stream = Stream()
	traces = Trace()
	processing  = None
	activities = []
	gap_list = []
	overlap_list= []
	xunits = None
	yunits = None
	calibrated = False
	log_process = {}       #Log of what is done on the channels
	plot_options = plot_options_process()


	#Builder
	def __init__(self, network=None, station= None, location= None, \
				  channel= None, traces = None, channelprop = None):
		self.network  = network
		self.station  = station
		self.location = location
		self.channel  = channel
		self.gap_list = []
		self.overlap_list = []
		self.log_process = {}
		self.calibrated = False
		self.plot_options = plot_options_process()
		if channelprop is not None:
			self.longname = get_long_name_from_prop_list(channelprop, self.channel, self.location)
			self.yunits   = get_unit_from_prop_list(channelprop, self.channel, self.location)
			self.xunits   = "Time (UTC)"
			self.sps = float(channelprop.loc[(channelprop["chan"] == self.channel) & (channelprop["locid"]==self.location), 'sps'].apply(str).iloc[0])
			self.orientation = channelprop.loc[(channelprop["chan"] == self.channel) & (channelprop["locid"]==self.location), 'orientation'].apply(str).iloc[0]
			# Define mode
			if "Sci" in self.longname and "Scientific" not in self.longname:
				self.mode = "Sci"
			elif "Eng" in self.longname:
				self.mode = "Eng"
			elif "-HKT" in self.longname:
				self.mode = "HKT"
			elif "-PXT" in self.longname:
				self.mode = "PXT"
			elif "-Temp" in self.longname or "_Temp" in self.longname:
				self.mode = "Temp"
			elif "SEIS-AC" in self.longname:
				self.mode = "SEIS-AC"
			elif "SEIS-DC" in self.longname:
				self.mode = "SEIS-DC"
			elif "Pressure_Outside" in self.longname:
				self.mode = "Pres"
			elif "Magne" in self.longname:
				self.mode = "Mag"
			else:
				self.mode = "Undefined"
			#Define detector
			if "VBB" in self.longname:
				if  "-HKT" in self.longname or "-PXT" in self.longname \
						or "-Temp" in self.longname \
						or "_Temp" in self.longname:
					self.detector = "ThermalSensor"
				else:
					self.detector = "VBB"
			elif "SP" in self.longname:
				self.detector = "SP"
			elif "SEIS-DC" in self.longname:
				self.detector = "Leveling"
			elif "SEIS-AC" in self.longname:
				self.detector = "Acquisition"
			elif "Scientific_Temp" in self.longname or \
							"Atmosphere_Temp" in self.longname:
				self.detector = "ThermalSensor"
			elif "Pressure_Outside" in self.longname:
				self.detector = "PressureSensor"
			elif "Magne" in self.longname:
				self.detector = "MagSensor"
			else:
				self.detector = "Other"

			#Define output 
			if "Vel" in self.longname:
				self.output = 'Vel'
			elif "Pos" in self.longname or "Position" in self.longname:
				self.output = 'Pos'
			elif "Temp" in self.longname or "Temperature" in self.longname \
					or "HKT" in self.longname or "PXT" in self.longname:
				self.output = 'Temp'
			elif "Pressure_Outside" in self.longname:
				self.output = 'Pres'
			elif "Magne" in self.longname:
				self.output = 'Mag'
				
			else:
				if self.detector == "SP":
					self.output = 'Vel'
				else:
					self.output = self.detector

			#Define gain
			if "HG" in self.longname:
				self.gaintype = "HG"
			elif "LG" in self.longname:
				self.gaintype = "LG"
			else:
				self.gaintype = "Undefined"

			#Define if calibrated or not
			if self.yunits not in ["Counts", "Count"]:
				self.calibrated = True


		if traces is not None:
			self.stream = traces
			cpstream = self.stream.copy()
			self.traces = cpstream.merge(method=0, fill_value=None)
			self.overlap_list, self.gap_list = \
						self.get_gap_overlap_time_in_channel(location, channel)
			self.starttime = self.traces[0].stats.starttime
			self.endtime   = self.traces[0].stats.endtime
			self.sps       = self.traces[0].stats.sampling_rate 
			self.npts      = self.traces[0].stats.npts


	def __str__(self):

		list2join = ["=====================================================\n"\
					"longname: {self.longname}\n", \
					"------------------------------------------------------\n",
					 "network: {self.network}\n", "station: {self.station}\n", \
					 "location: {self.location}\n", "channel: {self.channel}\n"
					 "sps: {self.sps}\n","detector: {self.detector}\n", \
					 "output: {self.output}\n", "gaintype: {self.gaintype}\n"
					 "orientation: {self.orientation}\n", "mode: {self.mode}\n"
					 "starttime: {self.starttime}\n", \
					 "endtime: {self.endtime}\n", \
					 "npts= {self.npts}\n", \
					 "xunits: {self.xunits} \n", "yunits: {self.yunits} \n",
					 "processing: {self.processing}\n", \
					 "gaps: {self.gap_list}\n", \
					 "overlaps: {self.overlap_list}\n", \
					 "SEIS activity: {self.activities}\n",
					 "log_process: {self.log_process}\n",\
					 "calibrated: {self.calibrated}\n", \
					 "plot options: {self.plot_options}"]

		return "".join(list2join).format(self=self)


	def match_process(self, process_list):
		self.processing == None
		for proc in process_list:
			
			if proc.detector == self.detector and proc.mode == self.mode \
					and proc.gaintype == self.gaintype and proc.output == self.output:
				self.processing = proc
				
				pass
		if self.processing == None:
			print("No processing found for {}: {}/{}/{}/{}".format(self.longname,\
					self.network, self.station, self.location, self.channel))

	def plot(self, *args, **kwargs):
		self.stream.plot(*args, **kwargs)


	def copy(self):
		temp = AspicChannel(self.network,self.station, self.location, \
							self.channel)
		temp.longname       = self.longname
		temp.sps            = self.sps
		temp.orientation    = self.orientation
		temp.mode           = self.mode
		temp.detector       = self.detector
		temp.output         = self. output
		temp.gaintype       = self.gaintype
		temp.starttime      = self.starttime
		temp.endtime        = self.endtime
		temp.stream         = self.stream.copy()
		temp.traces         = self.traces.copy()
		temp.processing     = self.processing
		temp.activities     = list(self.activities)
		temp.gap_list       = list(self.gap_list)
		temp.overlap_list   = list(self.overlap_list)
		temp.xunits         = self.xunits
		temp.yunits         = self.yunits
		temp.log_process    = self.log_process
		temp.calibrated     = self.calibrated
		temp.plot_options   = self.plot_options.copy()
		temp.npts           = self.npts

		return temp


	def tra2st(self):
		"""
		To convert Trace object to Stream object.
		In fact, it's Stream object are preferable to use if many traces, 
		especially for the calibration.
		Just to make sure that filters, calibration... 
		are applied on each part of the signal

		@return: no return, it's a setter of the class Aspic_channel'

		"""
		tr_copy = self.traces.copy()
		sta = tr_copy.split()
		#print(sta)
		self.stream = Stream(traces=sta)


	def st2tra(self, method = 0, fill_value = None, interpolation_samples=0):
		"""
		To convert Stream object to Trace

		@return:no return, it's a setter of the class Aspic_channel'
		"""
		cp_stream = self.stream.copy()
		cp_stream.merge(method = method, fill_value = fill_value, \
				interpolation_samples = interpolation_samples)
		self.traces = cp_stream[0]
		self.npts = cp_stream[0].stats.npts


	def set_plot_options(self, fixed_sol = None, mDate = None):

		plot_options = plot_options_process()
		
		for proc in self.processing.processes:
			if proc.name == "plot":
				if "ptype" in proc.param_dict:
					plot_options.ptype = proc.param_dict["ptype"]

				if "utc_or_lmst" in proc.param_dict:
					plot_options.utc_or_lmst = proc.param_dict["utc_or_lmst"]

				if "showgaps" in proc.param_dict:
					plot_options.show_gaps = bool(proc.param_dict["showgaps"])

				if "showoverlaps" in proc.param_dict:
					plot_options.show_overlaps = bool(proc.param_dict["showoverlaps"])

				if "showactivity" in proc.param_dict:
					plot_options.show_activity = bool(proc.param_dict["showactivity"])

				if "solscale" in proc.param_dict:
					plot_options.sol_scale = bool(proc.param_dict["solscale"])

				if plot_options.sol_scale:
					if fixed_sol is not None and mDate is not None:
						plot_options.forced_starttime = mDate.get_lmst_to_utc(lmst_date=fixed_sol)
						plot_options.forced_endtime = mDate.get_lmst_to_utc(lmst_date=fixed_sol+1)
				if "save" in proc.param_dict:
					plot_options.save = proc.param_dict["save"]
				if "display" in proc.param_dict:
					plot_options.display = proc.param_dict["display"]
				self.plot_options = plot_options.copy()


	def set_location(self, locid):
		self.location = locid
		self.stream[0].stats.location = locid


	def set_channel(self, chan):
		self.channel = chan
		self.stream[0].stats.channel = chan



	def get_output_units(self):
		"""
		Function to retrieve the output units for calibration. Interesting in
		case of several calibration for the same channel

		@return: List of units

		"""
		for proc in self.processing.processes:
			if proc.name == "calibration":
				if "output" in proc.param_dict:
					output = proc.param_dict["output"]
					output_list = output.split(",")
					new_outputlist = []
					for out in output_list:
						newout = out.replace(" ", "")
						new_outputlist.append(newout)
					#print(new_outputlist)
					return new_outputlist
				else:
					return None
			else:
				return None


	def get_gap_overlap_time_in_channel(self, loc, chan):
		'''
		This function checks if there are some overlap in the trace.
		@channels:
		'''
		#print("Inside Get overlap")
		station = self.stream

		#print("type(station)", type(station))

		#Initialize list of gap and overlap
		overlap_list = []
		gap_list = []

		idx_in_station = get_indexes_of_channel_in_station(station, loc, chan)
		#print("Index = ", idx_in_station)
		if len(idx_in_station) >= 1:

			#
			# OVERLAP
			##############################
			for i in range(0, len(idx_in_station)):
				for j in range(i+1, len(idx_in_station)):
					#print("Combination: ", i, j)
					# Situation 1
					#
					#   ----------------
					#			 --------------
					if station[idx_in_station[i]].stats.starttime < station[idx_in_station[j]].stats.starttime and \
						station[idx_in_station[j]].stats.starttime < station[idx_in_station[i]].stats.endtime:

						if station[idx_in_station[i]].stats.endtime > station[idx_in_station[j]].stats.starttime and \
							station[idx_in_station[i]].stats.endtime<station[idx_in_station[j]].stats.endtime:
							overlap_list.append((station[idx_in_station[j]].stats.starttime, station[idx_in_station[i]].stats.endtime))
					# Situation 2
					#
					#   --------------------------------
					#			 --------------

					if station[idx_in_station[i]].stats.starttime < station[idx_in_station[j]].stats.starttime and \
						station[idx_in_station[i]].stats.starttime < station[idx_in_station[j]].stats.endtime:
						if station[idx_in_station[j]].stats.endtime < station[idx_in_station[i]].stats.endtime:
							overlap_list.append((station[idx_in_station[j]].stats.starttime, station[idx_in_station[j]].stats.endtime))

					# Situation 3
					#
					#			  ----------------------
					#	------------------
					if station[idx_in_station[i]].stats.starttime > station[idx_in_station[j]].stats.starttime and \
						station[idx_in_station[i]].stats.starttime > station[idx_in_station[j]].stats.starttime:
						if station[idx_in_station[i]].stats.endtime > station[idx_in_station[j]].stats.endtime:
							overlap_list.append((station[idx_in_station[i]].stats.starttime, station[idx_in_station[j]].stats.endtime))

					# Situation 4
					#
					#		  ------------
					#  ----------------------------
					if station[idx_in_station[j]].stats.starttime < station[idx_in_station[i]].stats.starttime and \
						station[idx_in_station[i]].stats.endtime < station[idx_in_station[j]].stats.endtime:
						overlap_list.append((station[idx_in_station[i]].stats.starttime, station[idx_in_station[i]].stats.endtime))

			#
			# GAP
			#######################

			for i in range(0, len(idx_in_station)-1):
				if station[idx_in_station[i]].stats.endtime.timestamp < station[idx_in_station[i+1]].stats.starttime.timestamp:
					gap_list.append((station[idx_in_station[i]].stats.endtime, station[idx_in_station[i+1]].stats.starttime))
			return overlap_list, gap_list

		else:
			return None, None

	def set_gap_overlap_time_in_channel(self):
		self.overlap_list, self.gap_list = self.get_gap_overlap_time_in_channel(self.location, self.channel)

	def analyze_overlap(self, plot_overlap = True):
		"""
		Function used to analyze overlaps in a signal.
		It returns 3 subplots:
			- one with the common plotted signal,
			- one with each traces colored differently
			- one with just the position of the trace to see overlaps

		"""
		import numpy as np
		import sys
		import matplotlib.pyplot as plt
		#from matplotlib.dates import (MINUTELY, MONTHLY, YEARLY, DAILY, DateFormatter,
		#                              rrulewrapper, RRuleLocator, drange, AutoDateLocator)
		from matplotlib.ticker import FuncFormatter, MaxNLocator
		from obspy import UTCDateTime


		# channel selection
		if self.stream is not None:
			try:
				channel_selected = self.stream
			except TypeError:
				sys.exit("Channel not found. Please check your parameters.")
				channel_selected = None

		else:
			channel_selected = None

		# merge the different traces into a single trace
		# in case of overlap, the overlap zone remains EMPTY

		if channel_selected is not None and len(channel_selected)>0:
			# Stacking with no overlaping (leave a gap between data -> masked data)
			channel_selected_no_overlaping = Stream()
			for chan in channel_selected:
				channel_selected_no_overlaping+=chan

			# Create a single merged trace
			channel_selected_no_overlaping.merge(method=0, fill_value=None)
			channel_selected_no_overlaping = channel_selected_no_overlaping[0]

			# Create the time vector
			channel_selected_no_overlaping_time = np.arange(
					channel_selected_no_overlaping.stats.starttime,
					channel_selected_no_overlaping.stats.starttime + \
					channel_selected_no_overlaping.stats.npts/channel_selected_no_overlaping.stats.sampling_rate,
					channel_selected_no_overlaping.stats.delta)

		else:
			sys.exit("Oups, your channel is empty... No more process is possible.")

		if channel_selected is not None and len(channel_selected)>0:
			if self.overlap_list == [] or self.gap_list == []:
				overlap_list, gap_list = \
						self.get_gap_overlap_time_in_channel(self.location,
															 self.channel)

			else:
				overlap_list = self.overlap_list
				#gap_list = self.gap_list

			#if len(overlap_list) != 0:
				#print("List of overlaps: ")
			#	for ov in overlap_list:
			#		print(ov)
			#if len(gap_list) != 0:
			#	#print("\n ******* \n List of gaps:")
			#	for ga in gap_list:
			#		print(ga)

		if plot_overlap:
			from matplotlib import collections  as mc
			from itertools import cycle # cycle over colors on plots



			figs, axs = plt.subplots(3, 1, sharex='col', sharey=False, figsize=(9,7))
			plottitle = channel_selected_no_overlaping.stats.location+"-"\
				+channel_selected_no_overlaping.stats.channel
			figs.suptitle(plottitle, fontsize=8)


			# Subplot 1 : Show only the signal like in ASPIC
			axs[0].set_title("Signal as represented in ASPIC")
			axs[0].plot(channel_selected_no_overlaping_time,
						channel_selected_no_overlaping, color= 'red',
						linestyle='-', marker=None, linewidth=0.4,
						alpha=0.6, label = "Raw data merged")
			axs[0].grid()
			axs[0].tick_params(labelsize='small', width=1)
			axs[0].xaxis.set_tick_params(rotation=30, labelsize=8)
			axs[0].xaxis.set_major_locator(MaxNLocator(10))
			axs[0].xaxis.set_major_formatter(FuncFormatter(
					lambda tstamp, pos : (str(UTCDateTime(tstamp))[:-4])))
			axs[0].set(xlabel="Time (UTC)")
			axs[0].set(ylabel="Counts")
			axs[0].legend(framealpha=0.01, loc='best', borderaxespad=0.,
					fontsize="x-small", ncol = 3)

			# Subplot 2 : Show the different pieces of the signal tilted to make overlaps appearing
			cycol = cycle('bgrcmk')
			for i in range(0, len(channel_selected)):
				starttime = channel_selected[i].stats.starttime
				#print("channel ", i," - starttime : ",  starttime)
				delta = channel_selected[i].stats.delta
				npts = channel_selected[i].stats.npts
				sampling_rate = channel_selected[i].stats.sampling_rate
				t=np.arange(starttime, starttime+npts/sampling_rate, delta)
				dte = [UTCDateTime(dat).timestamp for dat in t]
				axs[1].plot(dte, np.add(channel_selected[i].data,
					50000+i*0), color=next(cycol),
					linestyle='-', marker=None, linewidth=0.4,alpha = 0.5,
					label=channel_selected[i].stats.channel+" trace n°"+str(i))
			axs[1].set_title("Signal with trace tilted in Y")
			#axs[1].legend(framealpha=0.01, loc='best', borderaxespad=0.,
			#			  fontsize="x-small", ncol = 3)
			axs[1].grid()
			axs[1].tick_params(labelsize='small', width=1)
			axs[1].xaxis.set_tick_params(rotation=30, labelsize=8)
			axs[1].xaxis.set_major_locator(MaxNLocator(12))
			axs[1].xaxis.set_major_formatter(
				FuncFormatter(
							lambda tstamp, pos : (str(UTCDateTime(tstamp))[:-2])))
			#mDate = MarsConverter()
			#axs[1].xaxis.set_major_formatter(FuncFormatter(lambda tstamp, \
			#		pos : mDate.get_utc_2_lmst(UTCDateTime(tstamp), \
			#					output="date")[:-7]))
			
			#axs[1].set(xlabel="Time (LMST)")
			axs[1].set(xlabel="Time (UTC)")
			axs[1].set(ylabel="Counts")

			#Subplot 3 : Show only segment with overlaps and gaps
			i=0

			lines_chan = []
			for chan in channel_selected:
				i+=.1
				lines_chan.append([(chan.stats.starttime, i), (chan.stats.endtime, i)])

			cycol = cycle('bgrcmk')
			for overlap in overlap_list:
				plt.axvspan(overlap[0].timestamp, overlap[1].timestamp,
							fc=next(cycol), alpha=0.5)

			lc = mc.LineCollection(lines_chan, colors=next(cycol), linewidths=2)
			axs[2].set_title("Simplified plot with positions of signals and overlaps")
			axs[2].add_collection(lc)
			axs[2].autoscale()
			axs[2].margins(0.1)
			axs[2].grid()
			axs[2].xaxis.set_tick_params(rotation=30, labelsize=8)
			axs[2].xaxis.set_major_locator(MaxNLocator(12))
			axs[2].xaxis.set_major_formatter(FuncFormatter(lambda tstamp, \
									pos : (str(UTCDateTime(tstamp))[:-2])))
			axs[2].set(xlabel="Time (UTC)")

			#axs[2].xaxis.set_major_formatter(FuncFormatter(lambda tstamp, \
			#		pos : mDate.get_utc_2_lmst(UTCDateTime(tstamp), \
			#					output="date")[:-7]))
			
			#axs[2].set(xlabel="Time (LMST)")


			plt.show()
			figs.autofmt_xdate()



	def apply_process(self, inv=None):
		"""
		
		Function used to analyse the type of process with their parameters
		The function check if there are some gaps in the data. 
		In this case, some calibration could be aborted: detrend, filtering, 
						rotation.
						

		Parameters
		----------
		process : STR 
					Contain the name of the process read in the process XML file
			
		kwargs : dict, optional
			dictionnary of all addictional parameters

		Returns
		-------
		Function apply on self object

		"""
		for proc in self.processing.processes:
			if (proc.name) in ["calibration", "detrend", "filter", "rotation", "detick"]:
				if proc.name == "detrend":
					if len(self.gap_list)==0 and len(self.overlap_list)==0:
						if proc.param_dict is not None:
							
							for keys in proc.param_dict:
								if keys in ["dspline", "order"]: 
									proc.param_dict[keys] = int(proc.param_dict[keys])
							kwargs = proc.param_dict
							#kwargs["plot"] = True
							self.traces.detrend(**kwargs)
							self.tra2st()
							self.log_process["detrend"] = kwargs
					else: 
						print("DETREND: Gaps in traces.")
				if proc.name == "filter":
					filttype = proc.param_dict["type"]
					if len(self.gap_list)==0 and len(self.overlap_list)==0:
						if proc.param_dict is not None:
							for keys in proc.param_dict:
								if keys in ["freq", "freqmin", "freqmax", "df"]: 
									proc.param_dict[keys] = float(proc.param_dict[keys])
								if keys in ["corners"]: 
									proc.param_dict[keys] = int(proc.param_dict[keys])
								if keys in ["zerophase"]:
									if proc.param_dict[keys] == "True" : proc.param_dict[keys]=True
									if proc.param_dict[keys] == "False" : proc.param_dict[keys]=False
							
							kwargs = proc.param_dict.copy()
							kwargs.pop("type")
							self.traces.filter(filttype, **kwargs)
							self.tra2st()
							self.log_process["filter"] = kwargs
					else: 
						print("FILTER: Gaps in traces.")
				if proc.name == "calibration":
					
					if proc.param_dict is not None and inv is not None:
						self.remove_response_with_process(inv)
						self.log_process["calibration"] = proc.param_dict
						

				if proc.name == "detick":
					if self.sps == 20:
						self.traces = Remove_Tick_Noise_02_BH_2_trace(self.traces, path_tickdir)
						self.tra2st()
						self.log_process["detick"] = "Applied @ 20Hz"
					
	#==========================================================================
	# Calibration function using the XML file containing steps of the process
	#
	def remove_response_with_process(self, inv, out_units=None, plotopt=False):
		"""
		Function able to take is self process parameters to perform calibration.
		It analyses the parameter in self.processes
		1> function is looking for a calibration stage
		2> get the ouput units
		3> call the remove response function above to apply the calibration
				according to the input parameters
		@inv : Inventory object
		@out_units : string - used to bypass the units in process especially
				if there are more than one (ie. Volts, ACC)
		"""

		volts_l = ["Volt", "Volts", "V", "volts","volt"]
		sci_l = ["Disp", "DISP", "Vel", "VEL", "Velocity", \
						"Acceleration", "Acc", "ACC"]
		amp_l = ["Ampere", "Ampère", "Amp", "A", "AMPERE", "mA"]
		
		if self.processing is not None:
			if len(self.processing.processes) > 0:
				#print("In remove_response_with_process")
				for proc in self.processing.processes:
					if proc.name == "calibration":
						if proc.param_dict["output"] in volts_l or out_units in volts_l:
							#print(proc.param_dict["output"])
								self.remove_response(inv, output="Volts", plotopt= False)

						elif proc.param_dict["output"] in sci_l or out_units in sci_l:
							#print(proc.param_dict["output"])

							if "zeromean" in proc.param_dict.keys():
								if proc.param_dict["zeromean"] == "False":
									zero_mean = False
								elif proc.param_dict["zeromean"] == "True":
									zero_mean = True
							else:
								zero_mean = True
							if "prefilt" in proc.param_dict.keys():
								
								prefilt = proc.param_dict["prefilt"]
								if len(prefilt)!=4: sys.exit("prefilt option is an array of 4 values")
								for i in range(0, len(prefilt)):
									prefilt[i] = float(prefilt[i])
							else:
								prefilt = None
							if "waterlevel" in proc.param_dict.keys():
								waterlevel = proc.param_dict["waterlevel"]
							else:
								waterlevel = None

							if out_units is not None:								
								self.remove_response(inv, output=out_units, \
													zero_mean = zero_mean,
													prefilt = prefilt, \
													waterlevel=waterlevel, \
													plotopt= plotopt)
							else:
								self.remove_response(inv, \
									output=proc.param_dict["output"],\
									zero_mean = zero_mean,
									prefilt = prefilt, waterlevel=waterlevel,\
									plotopt= plotopt)
						elif proc.param_dict["output"] in amp_l or out_units in amp_l:
							self.remove_response(inv, output="A", plotopt= False)
							#print(proc.param_dict["output"])
						else:
							self.remove_response(inv, output="Other", plotopt= False)
							#print(proc.param_dict["output"])
					#else:
					#	print("No calibration stage in your XML file, please check...")
		else:
			print("No processes defined for {}/{}/{}/{}".format(self.detector,\
						self.output, self.mode, self.gaintype))


	def remove_response(self, inv, output=None, zero_mean = False, \
					 prefilt = None, waterlevel=None, plotopt = False):
		"""
		Function to calibrate the signal in physical units.
		It's based on ObsPy commands : 'remove_response'
		For the last 4 parameters, read "trace.remove_response" in obspy 
		documentation.
		Make a copy of your channel, calibration data takes place of the previous
		data.
		----
		INPUT:
			@inv: dataless read by ObsPy ('read_inventory' command)
			@output: string - output units
			@zero_mean: boolean - only for output in (DISP, ACC, VEL) 
						-> subtract mean value
			@prefilt: list of 4 arguments - only for output in (DISP, ACC, VEL)
			@waterlevel: frequency filter level
		----
		OUTPUT:
			'self' is calibrated

		"""
		import numpy as np
		import os, sys
		result = None
		#print("Option of plot:", plotopt)
		#print("Unités: ", self.yunits)
		if self.yunits in ["Counts", "Kelvin", "Kelvins", "Count", "K"]:
			if output in ["Disp", "DISP", "Vel", "VEL", "Velocity", \
				 "Acceleration", "Acc", "ACC"]:
				#print("    - Calibration in {}".format(output))
				#print("    - zero_mean = {} ".format(zero_mean))
				#print("    - Option of plot:", plotopt)
				# Classical calibration
				# --------------------------
				FIR_coeff_gains = get_FIR_gains(inv, self.station, \
									self.location, self.channel, \
									starttime=self.starttime, \
									endtime=self.endtime)
				
#				#FIR correction
#				
#				if FIR_coeff_gains is not None:
#					if len(FIR_coeff_gains)!=0:
#						for coef in FIR_coeff_gains:
#							for tr in self.stream:
#								tr.data = np.divide(tr.data, coef)
#							self.st2tra()
#						#print("Stream object: ", self.stream)
#
#						print("    - FIR coefficients correction for {}: {}".format(self.longname, FIR_coeff_gains))
#					else:
#						print("ERROR: FIR coefficients correction for {}. No gain FIR correction".format(self.longname))


				#PREFILT OPTION
				if prefilt is not None:
					# In case of passband filter
					#print("prefilt: ", prefilt)
					#print("Make sure it's in the good format")
					#if prefilt[0]=="[" and prefilt[-1]=="]":
					#	tmp = prefilt[1:len(prefilt)-1].split(',')
					#	prefilt = [float(val) for val in tmp]

					if waterlevel is not None:
						#print("waterlevel :{}".format(waterlevel))
						#print("Case with prefilt and waterlevel provided...")
						waterlevel = float(waterlevel)
						#Water level provided
						try:
							
							result = self.stream.remove_response(inventory=inv, \
										output=output, plot=plotopt, \
										zero_mean=zero_mean,\
										water_level = waterlevel, \
										pre_filt = prefilt)

							self.st2tra()
						except NotImplementedError:
							print(self.longname+ ": Impossible to calibrate this channel.")
						except ValueError:
							print(self.longname+ ": No matching response information found.")
					else:
						#No water level provided
						try:
							result = self.stream.remove_response(inventory=inv, \
											output=output, plot=plotopt, \
											zero_mean=zero_mean,\
											pre_filt=prefilt)
							self.st2tra()
						except NotImplementedError:
							print(self.longname+ ": Impossible to calibrate this channel.")
						except ValueError:
							print(self.longname+ ": No matching response information found.")
					#else:
					#	# Invalid format of the prefilt
					#	print("Invalid format of the prefilt parameters.")
					#	result = None
				else:
					# No "prefilt" passband filter provided
					# If no prefilt is provided, add a default prefilt with the 
					# following values
					# f1 = 1/(duration *10)
					# f2 = 1.5/(duration *10)
					# f3 = 0.8*sps/2
					# f4 = sps/4
					#f1 = 1/(self.stream[0].stats.npts*self.stream[0].stats.sampling_rate)
					#f2 = 1.5/(self.stream[0].stats.npts*self.stream[0].stats.sampling_rate)
					#f3 = 0.8*self.stream[0].stats.sampling_rate/2.0
					#f4 = self.stream[0].stats.sampling_rate/4.0
					#prefiltdefault = [f1, f2, f3, f4]
					#print("default prefilt: ",prefiltdefault)
					
					if waterlevel is not None:
						#Water level provided
						try:
							#result = self.stream.remove_response(inventory=inv, output=output, plot=plotopt, zero_mean=zero_mean,
							#			water_level=waterlevel, pre_filt = prefiltdefault)
							result = self.stream.remove_response(inventory=inv, output=output, plot=plotopt, zero_mean=zero_mean,
										water_level=waterlevel)
							self.st2tra()
						except NotImplementedError:
							print(self.longname+ ": Impossible to calibrate this channel.")
						except ValueError:
							print(self.longname+ ": No matching response information found.")
					else:
						#No water level
						# If no waterlevel is provided one take 145dB as default 
						# value.
						waterlevel = float(145)
						#print("default waterlevel: ", waterlevel)
						try:
							#print("plotopt", plotopt)
							#print("No prefilt, no waterlevel")
							#result = self.stream.remove_response(inventory=inv,\
							#				output=output, plot=plotopt,\
							#				zero_mean=zero_mean, water_level=waterlevel, \
							#				pre_filt = prefiltdefault)
							result = self.stream.remove_response(inventory=inv,\
											output=output, plot=plotopt,\
											zero_mean=zero_mean)
							self.st2tra()
						except NotImplementedError:
							print(self.longname+ ": Impossible to calibrate this channel.")
						except ValueError:
							print(self.longname+ ": No matching response information found.")

				# Units change in metadata
				if result is not None:

					self.calibrated = True
					if output == "DISP":
						self.yunits = r"$m$"
					elif output == "VEL":
						self.yunits= r"$m.s^{-1}$"
					elif output == "ACC":
						self.yunits = r"$m.s^{-2}$"

			elif output in ["Volt", "Volts", "V", "volts","volt"]:
				#print("    - Calibration in volts")
				#print("output", output, self.station, self.location, self.channel, self.starttime, self.endtime)
				# To calibrate our data in Ampere
				# -------------------------------
				# This part contains to corrections from counts to mA.
				# Since no stage is defined in the dataless, the calibration is
				# done using a special config file

				gain_c2v = get_counts_2_Volts(inv, self.station, self.location, self.channel, starttime=self.starttime, endtime=self.endtime )

				#get FIR coeff.
				# get list of coefficients
				FIR_coeff = get_FIR_gains(inv, self.station, self.location, self.channel, starttime=self.starttime, endtime=self.endtime)

				# get coefficients from the polynomial coef C_i to estimate the function Unit = f(V)
				#poly_coef = get_coeff_poly_resp_stage(inv, self.station, self.location, self.channel)


				#FIR correction
				####
				FIRCorr = True
				if FIRCorr:
					if FIR_coeff is not None:
						if len(FIR_coeff)!=0:
							for coef in FIR_coeff:
								for tr in self.stream:
									tr.data = np.divide(tr.data, coef)
							self.st2tra()
							#print("    - FIR coefficients correction for {}: {}".format(self.longname, FIR_coeff))
							result = self.stream.copy()
					else:
						#print("ERROR: FIR coefficients correction for {}: Impossible, error in coefs.")
						result = None

				#Counts to gain conversion
				####
				if gain_c2v is not None:
					if gain_c2v !=0:
						for tr in self.stream:
							tr.data = np.divide(tr.data, gain_c2v)
							
							#--------------------------------------------------
							# WARNING : here is added a division by 4 to get a 
							# correct amplitude in Volts
							#if self.output == "Pos" and self.detector == "VBB":
							#	tr.data = np.divide(tr.data, 4)
							#print(tr)
						self.st2tra()
						#print("    - Gain Counts to Volts for {}: {} ".format(self.longname, gain_c2v))
						result = self.traces.copy()
					else:
						#print("ERROR: Counts to Volts conversion for {}: Impossible, error in gain.".format(self.longname))
						result = None

				if result is not None:
					self.calibrated = True
					self.yunits = "Volts"

			elif output in ["Ampere", "Ampère", "Amp", "A", "AMPERE", "mA"]:
				#print("    - Calibration in Ampère resquested")

				if os.path.isfile(conf_dir+amp_calib_file) == True:
					try:
						script_tree_ampere = read_XML_script(conf_dir+amp_calib_file) # Read is an Obspy function
						root = script_tree_ampere.getroot()
						All_calib_list = []
						calib_dict = None
						for calibtype in root:
							param_dict_tmp = {}
							param_dict = {}
							for param in calibtype:
								if param.tag == "channel":
									chantab = param.text.split(",")
									param_dict_tmp[param.tag] = chantab
								else:
									param_dict_tmp[param.tag] = param.text
							for cha in param_dict_tmp["channel"]:
								param_dict["id"] = param_dict_tmp["location"]+"/"+cha
								param_dict["location"] = param_dict_tmp["location"]
								param_dict["channel"] = cha
								param_dict["expression"] = float(param_dict_tmp["expression"])
								param_dict["unit"] = param_dict_tmp["unit"]
								#print(param_dict)
							All_calib_list.append(param_dict)
						calib_dict = search(str(self.location)+"/"+self.channel, All_calib_list)
						if calib_dict is not None:
							#print(calib_dict)
							for tr in self.stream:
								tr.data = np.multiply(tr.data, calib_dict.get("expression"))
							#print(tr)
							self.st2tra()
							#print("    - Gain Counts to {} for {}: {} ".format(calib_dict.get("unit"), self.longname, calib_dict.get("expression")))
							result = self.traces.copy()
							if result is not None:
								self.calibrated = True
								self.yunits = calib_dict.get("unit")
						else:
							#print("No information to calibrate "+self.longname+" in Ampère.")
							result = None
					except FileNotFoundError:
						result = None
						sys.exit("{}: No such file or directory.".format(conf_dir+amp_calib_file))
					except TypeError:
						sys.exit("{}: Invalid file.".format(conf_dir+amp_calib_file))
					except SystemExit:
						sys.exit("{}: Invalid file.".format(conf_dir+amp_calib_file))
				else:
					sys.exit("{}: No such file or directory.".format(conf_dir+amp_calib_file))

			else:
				if self.yunits in ["K", "Kelvin", "Kelvins"]:
					#print("En kelvin")
					for tr in self.stream:
						tr.data = np.subtract(tr.data, 273.15)
					self.st2tra()
					result = self.traces.copy()
					self.calibrated = True
					self.yunits = "Celsius"
				else:
					# Calibration in every other physical units
					#print("output", output, self.station, self.location, self.channel, self.starttime, self.endtime)
					# get the counts to volts gain
					gain_c2v = get_counts_2_Volts(inv, self.station, self.location, self.channel, starttime=self.starttime, endtime=self.endtime)
					#print("gain_c2v", gain_c2v)
					# get list of coefficients
					FIR_coeff = get_FIR_gains(inv, self.station, self.location, self.channel, starttime=self.starttime, endtime=self.endtime)
					#print("FIR_coeff", FIR_coeff)
	
					# get polynomial coeff
					poly_coef = []
					poly_coef = get_coeff_poly_resp_stage(inv,  self.station, self.location, self.channel, starttime=self.starttime, endtime=self.endtime)
					#print("poly_coef", poly_coef)
	
					# get physical units from the dataless
					physical_units = get_units(inv, self.station	, self.location, self.channel, starttime=self.starttime, endtime=self.endtime)
					#print("physical_units", physical_units)
					#FIR correction
	
					if poly_coef != None and poly_coef != []:
	
						if len(FIR_coeff)!=0 or FIR_coeff != None:
							for coef in FIR_coeff:
								for tr in self.stream:
									tr.data = np.divide(tr.data, coef)
	
							#print("    - FIR coefficients correction for {}: {}".format(self.longname, FIR_coeff))
							self.st2tra()
							result = self.traces.copy()
						else:
							print("ERROR: FIR coefficients correction for {}: Impossible, error in coefs.")
							result = self.traces.copy()
	
						#Count to gain conversion
						if gain_c2v !=0 or gain_c2v != None:
							for tr in self.stream:
								#print(tr.stats)
								tr.data = np.divide(tr.data, gain_c2v)
							self.st2tra()
							#print("    - Gain Counts to Volts for {}: {} ".format(self.longname, gain_c2v))
							result = self.traces.copy()
						else:
							print("ERROR: Counts to Volts conversion for {}: Impossible, error in gain.")
	
	
						#Volts to physical units conversion
						#print("Volts to physical units conversion")
						if len(poly_coef)!=0 or (poly_coef)!=None:
							if len(poly_coef) == 1 and poly_coef[0] == 1.0:
								print("{} is already in physical units. Update of the y axis title.".format(self.longname))
							else:
								#print("{}: {}".format(self.longname, poly_coef))
								for i in range (0, len(self.traces.data)):
									self.traces.data[i]=funct_ML(poly_coef, self.traces.data[i])
								#print("{} is converted in {}.".format(self.longname, physical_units.get("yunit")))
								self.tra2st()
							self.calibrated = True
							if physical_units.get("yunit") == "C":
								self.yunits = "Celsius"
							else:
								self.yunits = physical_units.get("yunit")
					else:
						print("Wrong choice of units with {}".format(self.longname))
			# add information in MSEED File Header.
			if result is not None:
				if self.yunits is not None:
					self.log_process["calibration"] = "output_units=" + self.yunits
				else:
					self.log_process["calibration"] = "output_units=" + "ERROR"

		else:
			if self.calibrated:
				print(self.longname + " is already calibrated -> Please choose data in counts (DU)")


	def slice_in_sols(self, mDate):
		"""
		Function to slice a trace into sols it can work on traces with masked data
		Input parameters is a Trace object.

		@mDate: MarsConverter object to convert UTC to LMST
		@return : Stream object containing the different traces after slicing

		"""
		from obspy import Stream
		import numpy as np
		#Get the sols interval of the data
		starttime_sol = str(mDate.get_utc_2_lmst(utc_date=self.starttime))
		endtime_sol = str(mDate.get_utc_2_lmst(utc_date=self.endtime))
		posT_min = starttime_sol.find('T')
		posT_max = endtime_sol.find('T')

		# Here are the sols.
		sol_min = int(starttime_sol[0:posT_min])
		sol_max = int(endtime_sol[0:posT_max])

		# Create a list of UTCDateTime for each sols
		# ie, for InSight:
		#     s=226 -> UTC Time = 2019-07-16T10:17:35.978789Z
		#     s=227 -> UTC Time = 2019-07-17T10:57:11.222956Z
		#
		utc_time_list = []
		for s in range(sol_min+1, sol_max+1):
			s2utc = mDate.get_lmst_to_utc(lmst_date = s)
			utc_time_list.append(s2utc)

		# Don't forget the residual parts at the beginning and at the end
		utc_time_list.insert(0,self.starttime)
		utc_time_list.append(self.endtime)

		sliced_channels = self.copy()
		# Slice the trace into traces using the previous list
		trace_list = []
		for i in range(0, len(utc_time_list)-1):
			t = self.traces.slice(utc_time_list[i], utc_time_list[i+1])
			trace_list.append(t)


		# Create Stream object to store the data
		sliced_channels = Stream(traces=trace_list)
		self.stream = sliced_channels
		print("Masked data in output of slice methode: ", \
				isinstance(self.stream, np.ma.masked_array))

	def export(self, exfmt="mseed", filename = None ):
		"""
		Function to export data to various formats
		@exfmt : output format : mseed, mat, csv, pdf, jpg, eps
		"""
		import numpy as np

		if exfmt in ["mseed","miniseed","mSEED"]:
			if filename is not None:
				if isinstance(self.traces, np.ma.masked_array) == False:
					self.traces.write(str(filename), format='MSEED')
				else:
					sys.exit("EXPORT : Data are masked")

			else:
				filename  = str(self.starttime)+'.'
				filename += self.network+'.'
				filename += self.station+'.'
				filename += self.location+'.'
				filename += self.channel+'.'
				filename += "mseed"
				filename = filename.replace(":", "_")
				if isinstance(self.stream, np.ma.masked_array):
					self.stream.write(str(filename), format='MSEED')
				else:
					sys.exit("Data are masked")
			print("Data exported in file: ", filename)



################################

def get_indexes_of_channel_in_station(chan_list, locID, chan):
	'''
	Function which return a list of indexes for a given location/channel
	combination in a Stream variable.
	We have to assume that a channel can be in the file more than one time.
	So we return a list !

	@chan_list: Stream() variable with all the channels
	@locID: location
	@chan: channel
	@return: index_list
	'''
	# import fnmatch
	#from obspy import Stream, Trace

	index_list = []
	#print(station)
	for i in range(0,chan_list.count()):
		if(fnmatch.fnmatch(chan_list[i].stats.channel, chan) and
				 (fnmatch.fnmatch(chan_list[i].stats.location, locID))):
			index_list.append(i)
	return index_list

#  Generate a list of ASPICChannel to calibrate
#  if a channel is appearing more than one time
#  Trace are grouped in the same ASPIC Channel_object
#######################################################

def create_channel_list_2_process(trace, channel_list, Output,
								  channel_properties_list,
								  seis_activity_list_sol= None,
								  processes_list= None,
								  mDate= None, fixed_sol = None, 
								  nomerge=False):
	"""
	@trace: list of traces objects (obspy)
	@channel_list: list of channels from the XML files of requested channels

	"""

	import sys
	List_Channels_to_process = []
	#logger_aspic.info("Reading data from {}".format(Output))
	#print("Inside create_channel_list_2_process..." )
	# Loop over the list of channels
	if channel_list is not None and len(channel_list) >0:
		for chan in channel_list:
			#print(chan.get("location"), chan.get("channel"))
			#look for channel indexes in the trace object
			#(ie. 02.BHU @ the index 0     -> Only one trace))
			#     03.BHW @ the indexes 4,5 -> 2 traces -> One gap or one ovlap)
			index_list = get_indexes_of_channel_in_station(trace,
								chan.get("location"), chan.get("channel"))
			#print(index_list, chan.get("location"), chan.get("channel"))
			#if len(index_list)>0:
			#	print(chan," - " ,index_list)

			# Stream object created (Obspy object containing traces)
			sta = Stream()

			# Check if there are more than one trace per channel
			if len(index_list)>0:
				if len(index_list)>1 and nomerge==False:
	
					t = trace[index_list[0]]
	
					# Create Aspic_channel object
					new_aspic_channel = AspicChannel(network = t.stats.network, \
							 station= t.stats.station, location = t.stats.location, \
							 channel = t.stats.channel, \
							 channelprop = channel_properties_list)
					starttime = trace[index_list[0]].stats.starttime
					endtime   = trace[index_list[0]].stats.endtime
					
					for i in range(0, len(index_list)):
						sta_i= Stream(traces=trace[index_list[i]])
						sta +=sta_i
	
					new_aspic_channel.stream = sta
					copy_sta = sta.copy()
					# Create a unique trace for plots
					copy_sta.merge(method=0, fill_value=None)
					new_aspic_channel.traces = copy_sta[0]
	
					starttime = new_aspic_channel.traces.stats.starttime
					endtime = new_aspic_channel.traces.stats.endtime
	
					#Fill ASPIC Object
					new_aspic_channel.starttime = starttime
					new_aspic_channel.endtime = endtime
					new_aspic_channel.longname = get_long_name_from_prop_list(
						channel_properties_list, t.stats.channel, t.stats.location)
					#print("new_aspic_channel.longname: ", new_aspic_channel.longname)
					new_aspic_channel.sps = t.stats.sampling_rate
					# Check for gaps between traces
					new_aspic_channel.overlap_list, new_aspic_channel.gap_list = \
						new_aspic_channel.get_gap_overlap_time_in_channel(new_aspic_channel.location, new_aspic_channel.channel)
					# Match process from the XML file to the channel
					#print(new_aspic_channel)
					if processes_list is not None:
						new_aspic_channel.match_process(processes_list)
						new_aspic_channel.set_plot_options(mDate = mDate, fixed_sol = fixed_sol)
						#print("In aspic_channel, l889: ", new_aspic_channel.plot_options)
					# Add the lander activity if any
					if seis_activity_list_sol is not None:
						if len(new_aspic_channel.activities) == 0:
							new_aspic_channel.activities = list(seis_activity_list_sol)
					List_Channels_to_process.append(new_aspic_channel)
					#print(new_aspic_channel)
				
				elif len(index_list) == 1 or nomerge == True:
					#print("Number of traces for the same loc/chan combination: {}".format(len(index_list)))
					for i in range(0, len(index_list)):
						t = trace[index_list[i]]
						new_aspic_channel = AspicChannel(t.stats.network, t.stats.station, 
											 t.stats.location, t.stats.channel, 
											 channelprop=channel_properties_list)
						#print(new_aspic_channel)
						starttime = t.stats.starttime
						endtime   = t.stats.endtime
						new_aspic_channel.stream = Stream(traces=trace[index_list[i]])
						new_aspic_channel.traces = trace[index_list[i]]
						new_aspic_channel.starttime = starttime
						new_aspic_channel.endtime = endtime
						new_aspic_channel.npts    = t.stats.npts
						new_aspic_channel.longname = get_long_name_from_prop_list(channel_properties_list, t.stats.channel, t.stats.location)
						new_aspic_channel.sps = t.stats.sampling_rate
						# Match process from the XML file to the channel
						if processes_list is not None:
							new_aspic_channel.match_process(processes_list)
							new_aspic_channel.set_plot_options(mDate = mDate, fixed_sol = fixed_sol)
						if seis_activity_list_sol is not None:
							if len(new_aspic_channel.activities) == 0:
								new_aspic_channel.activities = list(seis_activity_list_sol)
						List_Channels_to_process.append(new_aspic_channel)
					#print(new_aspic_channel)

		return List_Channels_to_process

	else:
		#logger_aspic.error("Trace empty - no data in the miniseed file.")
		sys.exit("Trace is empty - no data in the miniseed file.")
		return None


def extract_channels_list_from_xml(script_tree):
	"""
	Function which read the tree structure with the list of channels 
	to download
	
	INPUT:
		@script_tree : tree - structure which contains the list of channels


	OUTPUT:
		@return : list of channels
	"""
	#print("extract_channels_list_from_xml")
	root = script_tree.getroot()
	channel_List = []
	for el in root:
	  #SELECTED CHANNELS
		if el.tag == 'channel':
			channel_dict={}
			for subel in el:
				channel_dict[subel.tag] = subel.text
				#print(subel.tag, subel.text)
			channel_List.append(channel_dict)
	return channel_List


def generate_channels_list(xml_config_file):

	#
	# Read the XML file with the list of channels to process
	##############################################

	channel_list = []
	import os
	import sys
	from modules.io.read_xml import read_XML_script
	if os.path.isfile(xml_config_file) == True:
		try:
			script_tree_channels = read_XML_script(xml_config_file)
			#print(script_tree_channels)
			channel_list = extract_channels_list_from_xml(script_tree_channels)
			#print("{}: channels list read.".format(xml_config_file))
		except FileNotFoundError:
			sys.exit("{}: No such file or directory.".format(xml_config_file))
		finally:
			return channel_list

	else:
		sys.exit("{}: No such file or directory.".format(xml_config_file))


#  Functions related with the management of the list
#  of properties
###################################################
def read_csv_panda(textfiletoread, delimiter):

	'''
	Function used to read the properties of the channels 
	(propertieFile variable)
	
	INPUT:
	@textfiletoread: input properties file
	@delimiter: to define the delimiter
	
	OUTPUT:
	@return: A pandas structure containing all the properties
	'''
	import pandas as pd
	import os
	import sys
	try:
		if os.stat(textfiletoread).st_size > 0:
			fulldata = pd.read_csv(textfiletoread, delimiter=delimiter, \
									dtype=str)
			return fulldata
		else:
			print("The file is empty.")
			sys.exit()
	except IOError:
		print("Could not read the file ", textfiletoread)
		sys.exit()



#==============================================================================
# Function to generate signal in an interval of time at a given frequency
#
def straight_ln(starttime, endtime, sps, startsig = None, endsig = None, \
							forced_val= None , dtype= "int32"):
	
	"""
	Function to generate signal to fill the gaps into a seismic signal
	----
	INPUT:
		@starttime : UTCDateTime object with start time of the gap
		@endtime   : UTCDateTime object with end time of the gap
		@startsig  : int - value of the signal @ starttime
		@endsig    : int - value of the signal @ endtime
		@sps       : sampling rate of the signal
	----
	OUTPUT: 
		@total_sig : list of data corresponding to the equation of a line
						
	"""
	import numpy as np
	
	np_pts = int((endtime-starttime)*sps)

	# estimate the slope of the linear function
	if forced_val is None: 
		#print("forced_val: ", )
		#print("linear filling of the gap")
		if startsig is not None and endsig is not None:
			if endtime != starttime:
				coef = (endsig - startsig)/((endtime.timestamp - starttime.timestamp))
				coef/=sps
			else:
				sys.exit("starttime and endtime must be different to estimate \
							the slope.")
		else:
			coef = 1
	
		total_sig = [startsig+i*coef for i in range(1, np_pts+1)]
		#print("len of the vec: ", len(total_sig))
		total_sig_np = np.array(total_sig)
	else:
		#print("fill with value = ",forced_val)
		total_sig_np = np.full(np_pts, forced_val)
		
	# cast the result to dtype to be compatible with the other traces
	# it prevents from merge (function of Stream class) issues
	total_sig_int = total_sig_np.astype(dtype)
	return total_sig_int

#==============================================================================
#Function to generate white noise in an interval of time at a given frequency
#
def gen_wn(starttime, endtime, sps, startsig = None, endsig = None, mult = None
		, verbose = False):
	"""
	Function to generate white noise to fill the gaps into a seismic signal
	
	INPUT:
		@starttime : UTCDateTime object with start time of the gap
		@endtime   : UTCDateTime object with end time of the gap
		@startsig  : int - value of the signal @ starttime
		@endsig    : int - value of the signal @ endtime
		@sps       : sampling rate of the signal

	OUTPUT: 
		@total_sig_int: list of data corresponding to the equation of a line
						
	"""
	from random import gauss
	import numpy as np
	np_pts = int((endtime-starttime)*sps)
	#print("starttime: ",starttime, starttime.timestamp)
	#print("endtime: ", endtime, endtime.timestamp)
	#print("signal limits: ", startsig, endsig)


	# estimate the slope of the linear function
	if startsig is not None and endsig is not None:
		if endtime != starttime:
			coef = (endsig - startsig)/((endtime.timestamp - starttime.timestamp))
			coef/=sps
		else:
			sys.exit("starttime and endtime must be different to estimate \
						the slope.")
	else:
		coef = 1

	if mult is None :
		mult = 0

	# Create a random noise
	if mult is not None and mult!=0:
		rndsig = [mult*gauss(0.0, 1.0) for i in range(np_pts)]
	else:
		rndsig = [gauss(0.0, 1.0) for i in range(np_pts)]

	series = [startsig+i*coef for i in range(0, np_pts)]

	# add linear equation and noise
	total_sig = np.add(series, rndsig)

	# cast the result to int32 to be compatible with the other traces
	# it prevents from merge (function of Stream class) issues
	total_sig_int = total_sig.astype("int32")

	return total_sig_int



#==============================================================================
# Function to fill gaps into a given input Stream
#
def fill_gaps_with_wn(channels, forced_val = None, verbose = False):
	"""
	Test function to fill gaps in Stream objects with white noise
	
	INPUT:
	@channels : aspic_channel object
	
	OUPUT:
	@return : aspic_channel_object

	"""
	from obspy import Trace
	
	channels.stream = pop_fully_overlaped_trace(channels.stream)
	channels.set_gap_overlap_time_in_channel()
	channels.st2tra()


	for gap in channels.gap_list:
		print("gap: ", gap)
		for t in channels.stream:
			print(t.stats.starttime)
			if t.stats.endtime == gap[0]:
				startsig = t.data[-1]
				print(t.data[-1])
			if t.stats.starttime == gap[1]:
				print(t.data[0])
				endsig = t.data[0]

		print(">>>>>>>>>>>>>>>>>>>>>>>>>>>")
		print("Start time of the gap: {}".format(gap[0]))
		print("End time of the gap: {}".format(gap[1]))
		print("Duration of the gap: {} seconds".format(gap[1]-gap[0]))

		# Extracting millisecond informations
		mil0 = int(str(gap[0]._get_microsecond())[2:])
		mil1 = int(str(gap[1]._get_microsecond())[2:])
		delta = mil1-mil0

		# Adjust the gap to get the same millisecond digit at start and end.
		if delta<0:
			#print("delta = ", delta)
			l = list(gap)
			l[1] = l[1] - delta/1e6
			gap = tuple(l)
		else:
			#print("delta = ", 10000-delta)
			l = list(gap)
			l[1] = l[1] - (10000-delta)/1e6
			gap = tuple(l)

		starttime = gap[0]
		endtime = gap[1]+1/channels.sps

		# Create a white noise array to fill the gap
		if channels.sps is not None :
			# generation of the white noise
			my_series_wn = gen_wn(starttime, endtime, channels.sps, \
						 startsig = startsig, \
						 endsig = endsig, mult = 1, verbose = False)
			#print("len(my_series_wn)",len(my_series_wn))
			# convert white noise into a Trace object
			tr_wn = Trace()
			tr_wn.stats.station = channels.station
			tr_wn.stats.network = channels.network
			tr_wn.stats.location = channels.location
			tr_wn.stats.channel = channels.channel
			tr_wn.data = my_series_wn
			tr_wn.stats.starttime = starttime
			tr_wn.stats.npts = len(my_series_wn)
			tr_wn.stats.sampling_rate = channels.sps

			channels.stream.append(tr_wn)
			channels.stream.sort(keys = ['starttime'])

			channels_cp = channels.copy()


		else:
			sys.exit("sps is not defined for this channel")
	channels_cp.st2tra(method=1, fill_value="interpolate")
	return channels_cp

#==============================================================================
# Function to fill gaps into a given input Stream
#
def fill_gaps_with_line(channels, forced_val = None, verbose = False):
	"""
	Test function to fill gaps in Stream objects with white noise
	
	----
	INPUT:
		@channels: aspic_channel object
	----
	OUPUT:
		@return: aspic_channel_object

	"""
	from obspy import Trace 
	dtype = channels.stream[0].data.dtype

	channels_cp = channels.copy()

	channels_cp.stream = pop_fully_overlaped_trace(channels_cp.stream)
	channels_cp.set_gap_overlap_time_in_channel()
	channels_cp.st2tra()

	for gap in channels_cp.gap_list:
		#print("gap: ", gap)
		for t in channels_cp.stream:
			
			if t.stats.endtime == gap[0]:
				startsig = t.data[-1]
				
			if t.stats.starttime == gap[1]:
				endsig = t.data[0]

		# Extracting millisecond informations
		mil0 = int(str(gap[0]._get_microsecond())[2:])
		mil1 = int(str(gap[1]._get_microsecond())[2:])
		delta = mil1-mil0

		# Adjust the gap to get the same millisecond digit at start and end.
		if delta<0:
			#print("delta = ", delta)
			l = list(gap)
			l[1] = l[1] - delta/1e6
			gap = tuple(l)
		else:
			#print("delta = ", 10000-delta)
			l = list(gap)
			l[1] = l[1] - (10000-delta)/1e6
			gap = tuple(l)

		starttime = gap[0]
		endtime = gap[1] +(1/channels_cp.sps)
		npts = int((endtime-starttime)*channels_cp.sps)
		#print("npts", npts)
		
		if gap[0]+(npts*channels_cp.sps)<gap[1]:
			endtime = gap[0]+(npts+1)/channels_cp.sps
			
		#print("starttime, endtime: ", starttime, endtime)
		# Create a white noise array to fill the gap
		if channels_cp.sps is not None :
			# generation of the white noise
			my_series = straight_ln(starttime, endtime, channels_cp.sps, \
							startsig = startsig, endsig = endsig, \
							forced_val = forced_val, dtype = dtype)			
			# convert white noise into a Trace object
			#print("dtype of my_series: ", my_series.dtype)
			
			#cast_my_series = my_series.astype(dtype)
			cast_my_series = my_series.copy()
			#print("dtype of my_series after cast: ", cast_my_series.dtype)
			tr_wn = Trace()
			tr_wn.stats.station = channels_cp.station
			tr_wn.stats.network = channels_cp.network
			tr_wn.stats.location = channels_cp.location
			tr_wn.stats.channel = channels_cp.channel
			tr_wn.data = cast_my_series
			tr_wn.stats.starttime = starttime
			tr_wn.stats.npts = len(cast_my_series)
			tr_wn.stats.sampling_rate = channels_cp.sps
			channels_cp.stream.append(tr_wn)
			channels_cp.stream.sort(keys = ['starttime'])
		else:
			sys.exit("sps is not defined for this channel")
	channels_cp.st2tra(method=1,fill_value="interpolate")
	channels_cp.tra2st()

	return channels_cp

#==============================================================================
# merge the different traces into a single trace
# in case of overlap, the overlap zone remains EMPTY
def pop_fully_overlaped_trace(channel_selected):
	"""
	Function which look for fully overlaped signals and pop them up from the
	list to avoid to create some wrong gaps

	@channel_selected: list of channels
	@return: list of channels

	"""

	list2pop = []
	if channel_selected is not None and len(channel_selected)>0:

		for i in range(0, len(channel_selected)):
			for j in range(0, len(channel_selected)):
				if i!=j:
					if channel_selected[j].stats.starttime >= channel_selected[i].stats.starttime \
						and channel_selected[j].stats.endtime <= channel_selected[i].stats.endtime:
							print("tobe poped: ",channel_selected[j].stats.starttime, channel_selected[j].stats.endtime )
							list2pop.append(j)

	for idx in list2pop:
		channel_selected.pop(idx)

	return channel_selected


if __name__ == "__main__":
	print("Welcome in the aspicchannel class")