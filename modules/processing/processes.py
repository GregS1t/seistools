import os


#
# Class to define the list of processes attached to a detector
#
###################################################

class Processing():
	"""
	Class to handle the whole list of processes for each type of data
	A processing is defined by a bench of attributes
	- detector -> 'VBB', 'SP'...
	- mode	 -> 'Sci' or 'Eng' or None
	- output   -> 'Vel' or 'Pos'
	- gaintype -> 'HG' or 'LG'


	"""
	detector  = None
	mode	  = None
	output	= None
	gaintype  = None
	processes = []

	def __init__(self, detector=None, mode=None, output=None, gaintype=None):
		self.detector = detector
		self.mode = mode
		self.output = output
		self.gaintype = gaintype
		self.processes = []

	def __str__(self):
		"""
		Pretty print of the list of processes
		"""

		str2return  = " {self.detector} / {self.mode} / {self.output} / {self.gaintype} \n"
		for proc in self.processes:
			str2return+="- " + str(proc) + "\n"

		return str2return.format(self=self)

	def add_process(self, process):
		self.processes.append(process)

#
# Class to define a given process with its parameters
#
############################################

class Process: #define processing object
	"""
	Class to store all the informations about a single process found in the XML file.
	A process is defined by its name and a dictionnary of processes

	"""
	proctype = None
	param_dict = {}
	def __init__(self, process_name):
		self.name = process_name
		self.param_dict = {}

	def __str__(self):
		"""
		Pretty printable components of the process
		"""
		process_str =  self.name+": "
		for k, v in self.param_dict.items():
			process_str+=str(k)+"= "+str(v)+" "
		return process_str

	def add_option(self, key, value):
		self.param_dict[key] = value


def extract_process_from_xml(script_tree):

	"""
	Function which read the tree structure input with all the parameters of processing.
	This function defines here all the type of processing which can be taken into account.
	If the XML file contains unknown process it won't be used.

	@script_tree : tree - structure which contains all the paramaters for all the parameters
				of each process

	@return : The list of processing to apply for a given type of data

	"""
	root = script_tree.getroot()
	Process_number = 1
	All_process_list = []
	for processtype in root:
		proc = Processing()
		for process in processtype:
			#PROCESSING STEPS
			# TYPE OF DATA
			if process.tag == "datatype":
				for detail in process:
					if detail.tag == "detector":
						proc.detector = detail.text
					if detail.tag == "mode":
						proc.mode = detail.text
					if detail.tag == "output":
						proc.output = detail.text
					if detail.tag == "gaintype":
						proc.gaintype = detail.text
			# IF DETREND
			if process.tag == "detrend":
				#print(process.tag)
				pr = Process(process.tag)
				for options in process:
					if options.tag == "type":
						pr.proctype = options.text
						pr.add_option(options.tag, options.text)
					if options.tag == "parameters":
						for param in options:
							if param.text is not None:
								pr.add_option(param.tag, param.text)
				proc.add_process(pr)
				#print(pr)

			#
			# FILTER
			if process.tag == "filter":
				#print(process.tag)
				pr = Process(process.tag)
				for options in process:
					if options.tag == "type":
						pr.proctype = options.text
						pr.add_option(options.tag, options.text)
					if options.tag == "parameters":
						for param in options:
							if param.text is not None:
								pr.add_option(param.tag, param.text)
				proc.add_process(pr)
				#print(pr)
				
			if process.tag == "detick":
				pr = Process(process.tag)
				for options in process:
					if options.tag == "freq":
						pr.proctype = options.text
						pr.add_option(options.tag, options.text)
				proc.add_process(pr)
				#print(pr)

			# CALIBRATION
			if process.tag == "calibration":
				#print(process.tag)
				pr = Process(process.tag)
				for options in process:
					if options.tag == "parameters":
						for param in options:
							if param.tag == "output" :
								pr.add_option(param.tag, param.text)
							elif param.tag == "prefilt":
								prefilt = param.text
								#prefilt  = "[0.001,0.005,0.34,1]"
								prefilt = prefilt.replace("[", "")
								prefilt = prefilt.replace("]", "")
								prefilt = prefilt.split(",")
								pr.add_option(param.tag, prefilt)
							elif param.tag == "waterlevel":
								if param.text is not None or param.text!="":
									pr.add_option(param.tag, param.text)

							elif param.tag ==  "zeromean":
								if param.text is not None or param.text!="":
									pr.add_option(param.tag, param.text)
							else:
								if param.text is not None or param.text!="":
									pr.add_option(param.tag, param.text)
				proc.add_process(pr)
				#print(pr)
			# THRESHOLDS DECLARATION
			if process.tag == "thresholddef":
				pr = Process(process.tag)
				for options in process:
					if options.tag == "values":
						val = options.text
						thres = val.split(",")
						pr.add_option(options.tag, thres)
					else:
						if options.text is not None or options.text!="":
							pr.add_option(options.tag, options.text)
				proc.add_process(pr)
			#
			# ROTATION
			if process.tag == "rotation":
				pr = Process(process.tag)
				for options in process:
					if options.tag == "type" :
						pr.processtype = options.text
						pr.add_option(options.tag, options.text)

				proc.add_process(pr)
				#print(pr)

			if process.tag == "plot":
				#print(process.tag)
				pr = Process(process.tag)
				for options in process:
					if options.tag == "order":
						if options.text is not None or options.text!="":
							pr.add_option(options.tag, options.text)
					if options.tag == "ptype":
						if options.text is not None or options.text!="":
							pr.add_option(options.tag, options.text)
					if options.tag == "utc_or_lmst":
						if options.text is not None or options.text!="":
							pr.add_option(options.tag, options.text)
					if options.tag == "showoverlaps":
						if options.text is not None or options.text!="":
							pr.add_option(options.tag, options.text)
					if options.tag == "showgaps":
						if options.text is not None or options.text!="":
							pr.add_option(options.tag, options.text)
					if options.tag == "showactivity":
						if options.text is not None or options.text!="":
							pr.add_option(options.tag, options.text)
					if options.tag == "solscale":
						if options.text is not None or options.text!="":
							pr.add_option(options.tag, options.text)
				proc.add_process(pr)
				#print(pr)

			# PSD_WELSH
			##################
			if process.tag == "psdwelsh":
				pr = Process(process.tag)
				for options in process:
					if options.tag == "type":
						pr.processtype = options.text
						pr.add_option(options.tag, options.text)
					if options.tag == "parameters":
						for param in options:
							if param.text is not None:
								pr.add_option(param.tag, param.text)
								pass
				proc.add_process(pr)

			#
			# Spectrogram
			##################
			if process.tag == "spectrogram":
				pr = Process(process.tag)
				for options in process:
					if options.tag == "type":
						pr.processtype = options.text
						pr.add_option(options.tag, options.text)
					if options.tag == "parameters":
						for param in options:
							if param.text is not None:
								pr.add_option(param.tag, param.text)
								pass
				proc.add_process(pr)
		All_process_list.append(proc)
		Process_number+=1

	return All_process_list


def generate_processes_list(xml_process_file = None):

	from modules.io.read_xml import read_XML_script
	import sys
	#import lxml
	#
	# read the XML file to retrieve all the processings per type of channels
	###############################################
	if os.path.isfile(xml_process_file) == True:
		try:
			script_tree_processes = read_XML_script(xml_process_file)
			processes_list = extract_process_from_xml(script_tree_processes)
			#logger_aspic.info("{}: processes list read.".format(xml_process_file))
			#print("{}: processes list read.".format(xml_process_file))
		except FileNotFoundError:
			#logger_aspic.error("{}: No such file or directory.".format(xml_process_file))
			sys.exit("{}: No such file or directory.".format(xml_process_file))
		finally:
			return processes_list
	else:
		#logger_aspic.error("{}: No such file or directory.".format(xml_process_file))
		sys.exit("{}: No such file or directory.".format(xml_process_file))



#def calib_amp_param_from_xml(script_tree):
#	"""
#	Function to read the tree structure of the XML containing the calibration
#	parameters in Ampere.
#	This function is used for SEIS-DC* and SEIS-AC* channels where no calibra-
#	-tion stages are found in the dataless.
#	@script_tree : tree structure read from xml file with lxml library
#
#
#	"""
#	#for calib in root:
#	#	if calib.tag == "calibration":
#	#		for param in process:
#	#			if detail.tag == "detector":
#	#			if detail.tag == "mode":
#	#			if detail.tag == "output":
#	#			if detail.tag == "gaintype":

