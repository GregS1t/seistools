import os
import sys

def create_MSDS_postfile(channel_list, outputfile = "postfile.txt", 
						 station = None, outputdir= None, time_extension=None, 
						 starttime=None, endtime=None):
	"""
	Function to create a postfile used as input in a wget command to get data from MSDS server
	The postfile format is for example:

		XB ELYSE 02 MHU 2019-03-20T00:00:00.000000 2019-03-21T00:00:00.000000
		XB ELYSE 02 MHV 2019-03-20T00:00:00.000000 2019-03-21T00:00:00.000000
		XB ELYSE 67 SHU 2019-03-20T00:00:00.000000 2019-03-21T00:00:00.000000

	@channel_list : list dict with the channels to download
	@outputfile   : output postfile name -> if None, a default filename is 'postfile.txt'
	@station	  : station name in case you want to overwrite it
	@outputdir    : output directory name -> if None, a directory "exportASPIC" is created
	@timeextension: unique timestamp number to create a unique directory into 
					outputdir to make sure never overwrite data
	@starttime    : Datetime format : YYYY-MM-DDThh:mm:ss.ssssss
	@endtime      : Datetime format : YYYY-MM-DDThh:mm:ss.ssssss

"""
	import os
	if outputdir !='' and os.path.isdir(outputdir)==True:
		exportdir = outputdir
	else:
		directory = os.getcwd()
		exportdir = directory

	fullname = exportdir+"/"+outputfile
	#print("postfile", fullname)
	#print("from {} to {}".format(starttime, endtime))

	with open(fullname, 'w') as f:
		for chan in channel_list:
			if station is None:
				sta = chan["station"]
			else:
				sta = station
			f.write(chan["network"]+" "+sta+" "+chan["location"]+" "+chan["channel"]+" "+starttime+" "+endtime+"\n")
	f.close()
	return fullname


def get_data_from_portal(postfile, MSDSlogin, MSDSpasswd,  raw_dir_path= None, outputfile=None):
	"""
	Function used to get data from SEIS Data portail
	@postfile: text file with the channel id, the startdate and the enddate of the request
			ie:
				XB ELYSE 00 HHU 2019-04-06T17:38:02.052431Z 2019-04-07T18:17:37.296558Z
				XB ELYSE 00 BHU 2019-04-06T17:38:02.052431Z 2019-04-07T18:17:37.296558Z
				XB ELYSE 01 BHU 2019-04-06T17:38:02.052431Z 2019-04-07T18:17:37.296558Z
	@MSDSlogin: string containing email address
	@MSDSpasswd: string containing password
	@raw_dir_path: string containing the path to save the output file
	@outputfile: string containing the name of the output file.

	@return: exit code after processing
	"""

	MSDSAddress = "https://ws.seis-insight.eu/fdsnws/dataselect/1/query"

	#wgetline="wget -d --quiet --post-file="+"\""+postfile+"\""+" -O "+"\""+raw_dir_path+"/"+outputfile+"\""+" --http-user="+"\""+MSDSlogin+"\""+" --http-password="+"\""+MSDSpasswd+"\""+" https://ws.seis-insight.eu/fdsnws/dataselect/1/query"
	#wgetline="wget --no-check-certificate --quiet --post-file="+"""+postfile+"""+" -O "+"""+raw_dir_path+"/"+outputfile+"""+" --http-user="+MSDSlogin+" --http-password="+MSDSpasswd+" https://ws.seis-insight.eu/fdsnws/dataselect/1/query"
	
	
	if raw_dir_path is not None and outputfile is not None: 
		wgetline = "wget --no-check-certificate  --quiet -nv --post-file=\'{}\' -O \"{}/{}\" --http-user=\"{}\" --http-password=\"{}\" \"{}\"".format(\
                postfile, raw_dir_path, outputfile, MSDSlogin, MSDSpasswd, MSDSAddress)
	else:
		wgetline = "wget --no-check-certificate  -nv --quiet --post-file=\'{}\' --http-user=\"{}\" --http-password=\"{}\" \"{}\"".format(\
                postfile, raw_dir_path, outputfile, MSDSlogin, MSDSpasswd, MSDSAddress)
	
	try:
		wget_code = os.system(wgetline)	# Place where the request is executed.
	except:
		sys.exit("Download of the data did not worked.")

	return wget_code


def get_dataless_from_portal(MSDSlogin, MSDSpasswd, network="XB", station="ELYSE", 
							 raw_dir_path = None, outputfile = None, 
							 starttime = None, endtime = None, sol=""):
	"""
	Function to get dataless from SEIS Data Portal
	@MSDSlogin: string containing email address
	@MSDSpasswd: string containing password
	@raw_dir_path: string containing the path to save the output file
	@outputfile: string containing the name of the output file.
	@starttime: Starttime to extract the dataless
	@endtime: endtime to extract the dataless

	@return: exit code after processing
	"""
	#wgetline="wget -O '"+raw_dir_path+"/"+outputfile+"' --http-user="+MSDSlogin+" --http-password="+MSDSpasswd+\
	#	" 'https://ws.seis-insight.eu/fdsnws/station/1/query?network="+network+"&station="+station+\
	#		"&startTime="+starttime+"&endTime="+endtime+"&level=response&format=xml'"
	
	MSDSAddressDL = "https://ws.seis-insight.eu/fdsnws/station/1/query?"
	
	wgetline ="wget --no-check-certificate --quiet -O \"{}/{}\" --http-user=\"{}\" --http-password=\"{}\" \"{}network={}&station={}&startTime={}&endTime={}&level=response&format=xml\"".format(\
						raw_dir_path, outputfile,\
						MSDSlogin, MSDSpasswd, MSDSAddressDL,\
						network, station, starttime, endtime)
	try:
		print(wgetline)
		wget_code = os.system(wgetline)	# Place where the request is executed.

	except:
		#logger_aspic.error("Download of the dataless did not worked.")
		sys.exit("Download of the dataless did not worked.")
	return wget_code


def clean_dir(soldir, dataless=False, mseed=False, postfile=False):
	"""
	Function to clean the directory from old files. 
	This function keep the most recent file per type. 

	----
	INPUT : 
		@solstr: sol num  
		@dataless: (boolean) -> If True, directory is cleaned from old dataless files
		@mseed: (boolean) -> If True, directory is cleaned from old mseed files files
		@postfile: (boolean) -> If True, directory is cleaned from EVERY postfiles

	"""
	import os
	import sys
	import glob

	if dataless:
		#soldir = "{}/SOL{:04d}".format(workingdir, sol)
		if os.path.isdir(soldir)==False:
			print("{} - directory not found, please extract your data first".format(soldir))
		else:
			soldir_rawdata = soldir+"/"+"Raw"
			if os.path.isdir(soldir_rawdata)==False: 
				print("{} - directory not found, please extract your data first".format(soldir_rawdata))



		listdataless = glob.glob(soldir+"/"+"*.xml")
		#for file in listdataless: 
			#print(file)

		l_dataless_sci = []
		l_dataless_eng = []

		if listdataless is not None or len(listdataless)>0:
			for f in listdataless:
				#print("Filename: ", f)
				if "ELYSE" in f or "ELYS0" in f:
					l_dataless_sci.append(f)
				if "ELYHK" in f or "ELYH0" in f:
					l_dataless_eng.append(f)

		# Get most recent science dataless in this directory
		maxstamp = 0
		for f in l_dataless_sci:
			if f[-3:] == "xml":
				if int(f[-14:-4]) > maxstamp:
					maxstamp = int(f[-14:-4])
			elif f[-4:] == "seed":
				if int(f[-15:-5]) > maxstamp:
					maxstamp = int(f[-14:-4])
		for i,f in enumerate(l_dataless_sci):
			if str(maxstamp) in f:
				imax = i
				#print("most recent file is {} at index {}".format(f, imax))
				break

		for i in range(0,len(l_dataless_sci)):
			if i!=imax:
				#print(l_dataless_sci[i])
				os.remove(l_dataless_sci[i])

		#---------------------------------------------------------------------
		# Get most recent engineer dataless in this directory
		#
		maxstamp = 0
		for f in l_dataless_eng:
			if f[-3:] == "xml":
				if int(f[-14:-4]) > maxstamp:
					maxstamp = int(f[-14:-4])
			elif f[-4:] == "seed":
				if int(f[-15:-5]) > maxstamp:
					maxstamp = int(f[-14:-4])
		for i,f in enumerate(l_dataless_eng):
			if str(maxstamp) in f:
				imax = i
				#print("most recent file is {} at index {}".format(f, imax))
				break
		for i in range(0,len(l_dataless_eng)):
			if i!=imax:
				#print(l_dataless_eng[i])
				os.remove(l_dataless_eng[i])
	#----------------------------------------------------------
	#			 CLEAN DIR FROM OLD MSEED FILES
	#
	if mseed:
		if os.path.isdir(soldir)==False:
			print("{} - directory not found, please extract your data first".format(soldir))
		else:
			soldir_rawdata = soldir+"/"+"Raw"
			if os.path.isdir(soldir_rawdata)==False: 
				print("{} - directory not found, please extract your data first".format(soldir_rawdata))

		#==========================================================================
		# Check eng data in the directory Raw if exist
		#print("Looking for engineer data:")
		listfile = glob.glob(soldir_rawdata+"/"+"*.mseed")

		filetype_l = ["eng","sci"]
		for ftype in filetype_l: 
			if listfile is not None or len(listfile)>0:
				englist = glob.glob(soldir_rawdata+"/"+"*_{}_*.mseed".format(ftype))
				#for i, f_eng_file in enumerate(englist):
					#print(i, f_eng_file)
			else:
				sys.exit("No mseed file in this directory")
			#for each, get the most recent
			maxstamp = 0
			for f in englist:
				if int(f[-16:-6]) > maxstamp:
					maxstamp = int(f[-16:-6])

			f2open = None
			for i, f in enumerate(englist):
				if str(maxstamp) in f:
					f2open = f
					imax = i
					#print("most recent file is {} with type {} in place {}".format(f, ftype, imax))
					pass

			for i in range(0,len(englist)):
				if i!=imax:
					#print(englist[i])
					os.remove(englist[i])

	#----------------------------------------------------------
	#			 CLEAN DIR FROM ALL POSTFILES
	#
	if postfile:
		listpostfile = []
		listpostfile = glob.glob(soldir+"/"+"postfile_*.txt")
		if listpostfile is not None or len(listpostfile)>0:
			for file in listpostfile: 
				os.remove(file)
