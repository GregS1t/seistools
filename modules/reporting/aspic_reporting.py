#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 11 14:51:13 2019

@author: Grégory Sainton
@email : sainton@ipgp.fr

@purpose: module where one can fine several types of reporting



"""
from obspy import UTCDateTime
import sys

def create_daily_report(List_Channels_to_process, raw_data_plot, df, marsDate=None, tf=False, dataless=None, outfile = None):
	"""
	Function which plots all the data of channels from a list. 
	Transfer function can also be added to the report.
	----
	INPUT:
		@List_Channels_to_process: list - list of channels to process
		@raw_data_plot: plot_options_process object which contains the options of plot
		@df: pandas structures with all the details of the each channels
		@marsDate: MarsConverter object to convert date from UTC to LMST
		@tf: boolean to process transfert functions plots
		@dataless : inventory file to read the transfert functions
		@oudir : string with the output directory.
	----
	OUPUT:
	"""

	global mDate
	mDate = marsDate

	#Import useful lib
	import matplotlib.pyplot as plt
	from matplotlib.ticker import FuncFormatter, MaxNLocator
	from matplotlib import markers
	from matplotlib.backends.backend_pdf import PdfPages

	import numpy as np
	import math
	from itertools import cycle

	import sys
	sys.path.insert(0, './modules/')

	from modules.calibration import plot_responses

	#Open PDF file if option is asked (raw_data_plot.save = True)
	if outfile is not None:
		pdf_pages = PdfPages(outfile+'.pdf')
	else:
		import datetime
		tstamp = str(datetime.datetime.now().timestamp()).split('.')[0]

		pdf_pages = PdfPages("Reporting_"+tstamp+".pdf")

	
	
	#Sort the data by type of dectectors (SP, VBB, ... )
	ds = df.groupby(['detector'], sort= True)
	#print(ds.head())
	for name, group in ds:
		#print("Detector: ", name)
		#Sort the each group by output (Vel, Pos, Temp...)
		doutput = group.groupby(['output'], sort= True)
		for n,groupoutput in doutput:
			#print("Output: ", n)
			#Sort each group of output by sps
			dsps = groupoutput.groupby(['sps'], sort = True)
			for no,groupsps in dsps:
				#Sort each group with the same SPS with orientation
				#print('sps', no)
				dorient = groupsps.groupby(['orientation'], sort = True)
				plot_number = 0
				for nn, grouporient in dorient:
					plot_number+=1
				index = 0

				#Plot option to share x axis
				if raw_data_plot.share_x_axis:
					sharex = 'col'
					if raw_data_plot.sol_scale:
						min_date_plot = raw_data_plot.forced_starttime
						max_date_plot = raw_data_plot.forced_endtime


				else:
					sharex = False

				#Plot option to share x axis
				if raw_data_plot.share_y_axis:
					sharey = True
				else:
					sharey = False


				# Manage the case where there are more than one page per page.
				if plot_number>1:
					suptitle = ""
					figs, axs = plt.subplots(plot_number, 1, sharex=sharex, \
							  sharey=sharey, figsize=(11,7))
					figs.subplots_adjust(hspace = .01)
					cycol = cycle('grcmk')
					for nn, grouporient in dorient:
						cycol_gaps = cycle('by')
						#print("---> orientation", nn)
						#print("SUBPLOT N°", plot_number)
						for i in range(0, len(grouporient)):
							idx = grouporient.iloc[i]["index_in_list"]

							if List_Channels_to_process[idx].longname !='':
								nameChannel = List_Channels_to_process[idx].longname + " (" + List_Channels_to_process[idx].location + "/" +  \
								List_Channels_to_process[idx].channel +")"
							else:
								nameChannel = List_Channels_to_process[idx].location+ "/" + List_Channels_to_process[idx].channel


							start_time = List_Channels_to_process[idx].traces.stats.starttime
							t = np.arange(start_time, start_time + List_Channels_to_process[idx].traces.stats.npts/List_Channels_to_process[idx].traces.stats.sampling_rate, \
								 List_Channels_to_process[idx].traces.stats.delta)
							dte = [UTCDateTime(dat).timestamp for dat in t]
							#Plot the time series

							if len(dte) == len(List_Channels_to_process[idx].traces.data):
								axs[index].plot(dte, List_Channels_to_process[idx].traces.data, color= next(cycol), label = nameChannel, linestyle='-', marker='None', linewidth=1)
							else:
								delta = len(List_Channels_to_process[idx].traces.data) - len(dte)
								if delta > 0:
									temp = List_Channels_to_process[idx].traces.data[:-delta]
									axs[index].plot(dte, temp, color= next(cycol), label = nameChannel, linestyle='-', marker='None', linewidth=1)
								else:
									temp = dte[:delta]
									axs[index].plot(temp, List_Channels_to_process[idx].traces.data, color= next(cycol), label = nameChannel, linestyle='-', marker='None', linewidth=1)

							# X-axis options
							axs[index].xaxis.label.set_size(6)
							axs[index].xaxis.set_tick_params(rotation=30, labelsize=6, length=6, width=2)
							axs[index].xaxis.set_major_locator(MaxNLocator(12))

							if raw_data_plot.share_x_axis:
								axs[index].set_xlim(UTCDateTime(min_date_plot).timestamp, UTCDateTime(max_date_plot).timestamp)
							else:
								axs[index].set_xlim(dte[0], dte[-1])
							if raw_data_plot.utc_or_lmst == "LMST":
								axs[index].xaxis.set_major_formatter(FuncFormatter(lambda tstamp, pos : mDate.get_utc_2_lmst(UTCDateTime(tstamp), output="date")[:-7]))
								axs[index].set(xlabel="Date (LMST)")

							elif raw_data_plot.utc_or_lmst == "UTC":
								axs[index].xaxis.set_major_formatter(FuncFormatter(lambda tstamp, pos : str(mDate.format_tick_labels_UTC(UTCDateTime(tstamp))[:-7])[:-8]))
								axs[index].set(xlabel=raw_data_plot.xtitle)

							elif raw_data_plot.utc_or_lmst == "Both":
								axs[index].xaxis.set_major_formatter(FuncFormatter(lambda tstamp, pos : str(mDate.format_tick_labels_UTC(UTCDateTime(tstamp))[:-7])[:-8]))
								axs[index].set(xlabel=raw_data_plot.xtitle)

								if index == 0:
									ax2 = axs[0].twiny()
									ax2.xaxis.set_ticks_position('top') # set the position of the second x-axis to top
									ax2.xaxis.set_label_position('top') # set the position of the second x-axis to top
									ax2.spines['top'].set_position(('outward', 0))
									ax2.xaxis.set_tick_params(labelsize=6)

									ax2.set_xlim(axs[0].get_xlim())
									ax2.xaxis.set_major_formatter(FuncFormatter(lambda tstamp, pos : str(round(mDate.get_utc_2_lmst(UTCDateTime(tstamp), output="decimal"),2))))

							# Display Gaps
							######################
							color_gaps = next(cycol_gaps)
							if raw_data_plot.show_gaps == True:
								if List_Channels_to_process[idx].gap_list is not None or len(List_Channels_to_process[idx].gap_list)>0:
									#
									for gap in List_Channels_to_process[idx].gap_list:
										if (gap[0]!=0 and gap[1]!=0):
											x0_lower = gap[0]
											x0_lower = UTCDateTime(x0_lower)
											x1_upper = gap[1]
											x1_upper = UTCDateTime(x1_upper)

											axs[index].axvspan(x0_lower, x1_upper, fc=color_gaps, alpha=0.3)
											axs[index].plot(x0_lower, axs[index].get_ylim()[1], marker=markers.CARETLEFT, color = color_gaps)
											axs[index].plot(x1_upper, axs[index].get_ylim()[1], marker=markers.CARETRIGHT, color= color_gaps)
											axs[index].annotate(
												"Gap ({}/{})".format(List_Channels_to_process[idx].location, List_Channels_to_process[idx].channel),
												xy=(x1_upper, axs[index].get_ylim()[1]), xytext=(-10, -10),
												textcoords='offset points', ha='right', va='bottom',size=6, color = color_gaps)


							# Display overlaps
							######################
							if raw_data_plot.show_overlaps == True:
								if List_Channels_to_process[idx].overlap_list is not None or len(List_Channels_to_process[idx].overlap_list)>0:
									for overl in List_Channels_to_process[idx].overlap_list:
										if (overl[0]!=0 and overl[1]!=0):
											x0_lower = overl[0]
											x0_lower = UTCDateTime(x0_lower)
											x1_upper = overl[1]
											x1_upper = UTCDateTime(x1_upper)

											axs[index].axvspan(x0_lower, x1_upper, fc='red', alpha=0.3)
											axs[index].plot(x0_lower, axs[index].get_ylim()[1], marker=markers.CARETLEFT, color = 'red')
											axs[index].plot(x1_upper, axs[index].get_ylim()[1], marker=markers.CARETRIGHT, color= 'red')
											axs[index].annotate(
													"Overlap",
													xy=(x1_upper, axs[index].get_ylim()[1]), xytext=(-10, -10),
													textcoords='offset points', ha='right', va='bottom',size=6, color = "red")

							#
							# Display SEIS Activity
							########################
							if raw_data_plot.show_activity == True:
								if List_Channels_to_process[idx].activities is not None or len(List_Channels_to_process[idx].activities)>0:
									cycol2 = cycle('gcmk')
									for act in List_Channels_to_process[idx].activities:
										#SP Data
										if List_Channels_to_process[idx].detector == "SP":
											if "VBB" not in act.desc or "SP" in act.desc :
												col = next(cycol2)
												if "ON" in act.desc or "OFF" in act.desc:
													#print(act.starttime, act.endtime, act.desc)
													if "OFF" in act.desc:
														#color = "red"
														axs[index].plot(act.endtime, axs[index].get_ylim()[1], marker=markers.CARETLEFT, color = col)
													if "ON" in act.desc:
														#color = "green"
														axs[index].plot(act.endtime, axs[index].get_ylim()[1], marker=markers.CARETRIGHT, color= col)
													axs[index].vlines(act.endtime, 0, 1, transform=axs[index].get_xaxis_transform(), colors=col, alpha=0.3)
													axs[index].annotate(
															act.desc,
															xy=(act.endtime, axs[index].get_ylim()[1]), xytext=(-10, -10),
															textcoords='offset points', ha='right', va='top',size=5)

												else:
													col = next(cycol2)
													x0_lower = act.starttime
													x1_upper = act.endtime
													axs[index].axvspan(x0_lower, x1_upper, fc = col, alpha=0.5)
													axs[index].plot(act.starttime, axs[index].get_ylim()[1], marker=markers.CARETLEFT, color = col)
													axs[index].plot(act.endtime, axs[index].get_ylim()[1], marker=markers.CARETRIGHT, color = col)
													axs[index].annotate(
															act.desc,
															xy=(x1_upper, axs[index].get_ylim()[1]), xytext=(-10, -10),
															textcoords='offset points', ha='left', va='top',size=6, color=col)
										#VBB Data
										if List_Channels_to_process[idx].detector == "VBB":
											if "SP" not in act.desc or "VBB" in act.desc :
												col = next(cycol2)
												if "ON" in act.desc or "OFF" in act.desc:
													if "OFF" in act.desc:
														axs[index].plot(act.endtime, axs[index].get_ylim()[1], marker=markers.CARETLEFT, color = col)
													if "ON" in act.desc:
														axs[index].plot(act.endtime, axs[index].get_ylim()[1], marker=markers.CARETRIGHT, color= col)
													axs[index].vlines(act.endtime, 0, 1, transform=axs[index].get_xaxis_transform(), colors=col, alpha=0.3)
													axs[index].annotate(
															act.desc,
															xy=(act.endtime, axs[index].get_ylim()[1]), xytext=(-10, -10),
															textcoords='offset points', ha='right', va='top',size=5)

												else:
													col = next(cycol2)
													x0_lower = act.starttime
													x1_upper = act.endtime
													axs[index].axvspan(x0_lower, x1_upper, fc = col, alpha=0.5)
													axs[index].plot(act.starttime, axs[index].get_ylim()[1], marker=markers.CARETLEFT, color = col)
													axs[index].plot(act.endtime, axs[index].get_ylim()[1], marker=markers.CARETRIGHT, color = col)
													axs[index].annotate(
															act.desc,
															xy=(x1_upper, axs[index].get_ylim()[1]), xytext=(-10, -10),
															textcoords='offset points', ha='left', va='top',size=6, color=col)

										if List_Channels_to_process[idx].detector not in ["VBB", "SP"]:
											if "SP" not in act.desc and "VBB" not in act.desc :
												col = next(cycol2)
												if "ON" in act.desc or "OFF" in act.desc:
													if "OFF" in act.desc:
														axs[index].plot(act.endtime, axs[index].get_ylim()[1], marker=markers.CARETLEFT, color = col)
													if "ON" in act.desc:
														axs[index].plot(act.endtime, axs[index].get_ylim()[1], marker=markers.CARETRIGHT, color= col)
													axs[index].vlines(act.endtime, 0, 1, transform=axs[index].get_xaxis_transform(), colors=col, alpha=0.3)
													axs[index].annotate(
															act.desc,
															xy=(act.endtime, axs[index].get_ylim()[1]), xytext=(-10, -10),
															textcoords='offset points', ha='right', va='top',size=5)

												else:
													col = next(cycol2)
													x0_lower = act.starttime
													x1_upper = act.endtime
													axs[index].axvspan(x0_lower, x1_upper, fc = col, alpha=0.5)
													axs[index].plot(act.starttime, axs[index].get_ylim()[1], marker=markers.CARETLEFT, color = col)
													axs[index].plot(act.endtime, axs[index].get_ylim()[1], marker=markers.CARETRIGHT, color = col)
													axs[index].annotate(
															act.desc,
															xy=(x1_upper, axs[index].get_ylim()[1]), xytext=(-10, -10),
															textcoords='offset points', ha='left', va='top',size=6, color=col)


							#Y-axis options
							axs[index].yaxis.label.set_size(6)
							axs[index].yaxis.set_tick_params(labelsize=6)
							axs[index].set(ylabel=raw_data_plot.ytitle)

							# General plot options
							axs[index].legend(framealpha=0.01, loc='best', borderaxespad=0., fontsize="x-small")
							axs[index].grid(True)
							axs[index].tick_params(labelsize='x-small', width=2)
							figs.autofmt_xdate()
							suptitle += List_Channels_to_process[idx].longname+" "
							if len(suptitle)%100 > 90:
								suptitle += "\n"
						# Next subplot
						index +=1
					figs.suptitle(raw_data_plot.suptitle+"\n"+suptitle, fontsize=12)

				# Manage the case where there are only one plot
				else:
					suptitle = ""
					figs = plt.figure(figsize=(11, 7))
					axs = figs.add_subplot(111)
					axs.titlesize = 10

					cycol = cycle('bgrcmk')
					for nn, grouporient in dorient:

						for i in range(0, len(grouporient)):
							idx = grouporient.iloc[i]["index_in_list"]
							print("idx =", idx, List_Channels_to_process[idx].longname, List_Channels_to_process[idx].yunits )
							if List_Channels_to_process[idx].longname !='':
								nameChannel = List_Channels_to_process[idx].longname + " (" + List_Channels_to_process[idx].location + "/" +  \
								List_Channels_to_process[idx].channel +")"
							else:
								nameChannel = List_Channels_to_process[idx].location+ "/" + List_Channels_to_process[idx].channel

							start_time = List_Channels_to_process[idx].traces.stats.starttime
							t = np.arange(start_time, start_time + List_Channels_to_process[idx].traces.stats.npts/List_Channels_to_process[idx].traces.stats.sampling_rate, \
									 List_Channels_to_process[idx].traces.stats.delta)
							dte = [UTCDateTime(dat).timestamp for dat in t]
							#Plot the time series
							axs.plot(dte, List_Channels_to_process[idx].traces.data, color= next(cycol), label = nameChannel, linestyle='-', marker='None', linewidth=1)
							# X-axis options
							axs.xaxis.label.set_size(6)
							axs.xaxis.set_tick_params(rotation=30, labelsize=6, length=6, width=2)

							axs.xaxis.set_major_locator(MaxNLocator(12))

							if raw_data_plot.share_x_axis:
								axs.set_xlim(UTCDateTime(min_date_plot).timestamp, UTCDateTime(max_date_plot).timestamp)
							else:
								axs.set_xlim(dte[0], dte[-1])
							if raw_data_plot.utc_or_lmst == "LMST":
								axs.xaxis.set_major_formatter(FuncFormatter(lambda tstamp, pos : mDate.get_utc_2_lmst(UTCDateTime(tstamp), output="date")[:-7]))
								axs.set(xlabel="Date (LMST)")

							elif raw_data_plot.utc_or_lmst == "UTC":
								axs.xaxis.set_major_formatter(FuncFormatter(lambda tstamp, pos : str(mDate.format_tick_labels_UTC(UTCDateTime(tstamp))[:-7])[:-8]))
								axs.set(xlabel=raw_data_plot.xtitle)
							elif raw_data_plot.utc_or_lmst == "Both":
								axs.xaxis.set_major_formatter(FuncFormatter(lambda tstamp, pos : str(mDate.format_tick_labels_UTC(UTCDateTime(tstamp))[:-7])[:-8]))
								axs.set(xlabel=raw_data_plot.xtitle)

								if index == 0:
									ax2 = axs.twiny()
									ax2.xaxis.set_ticks_position('top') # set the position of the second x-axis to top
									ax2.xaxis.set_label_position('top') # set the position of the second x-axis to top
									ax2.spines['top'].set_position(('outward', 0))
									ax2.xaxis.set_tick_params(labelsize=6)
									ax2.set_xlim(axs.get_xlim())
									ax2.xaxis.set_major_formatter(FuncFormatter(lambda tstamp, pos : str(round(mDate.get_utc_2_lmst(UTCDateTime(tstamp), output="decimal"),2))))

							# Display Gaps
							if raw_data_plot.show_gaps == True:
								if List_Channels_to_process[idx].gap_list is not None or len(List_Channels_to_process[idx].gap_list)>0:
									cycol2 = cycle('bgrcmk')
									for gap in List_Channels_to_process[idx].gap_list:
										if (gap[0]!=0 and gap[1]!=0):
											x0_lower = gap[0]
											x0_lower = UTCDateTime(x0_lower)
											x1_upper = gap[1]
											x1_upper = UTCDateTime(x1_upper)
											axs.axvspan(x0_lower, x1_upper, fc='black', alpha=0.3)
											axs.plot(x0_lower, axs.get_ylim()[1], marker=markers.CARETLEFT, color = 'black')
											axs.plot(x1_upper, axs.get_ylim()[1], marker=markers.CARETRIGHT, color= 'black')
											axs.annotate(
													"Gap",
													xy=(x1_upper, axs.get_ylim()[1]), xytext=(-10, -10),
													textcoords='offset points', ha='right', va='bottom',size=6, color = "black")

							# Display overlaps
							if raw_data_plot.show_overlaps == True:
								if List_Channels_to_process[idx].overlap_list is not None or len(List_Channels_to_process[idx].overlap_list)>0:
									for overl in List_Channels_to_process[idx].overlap_list:
										if (overl[0]!=0 and overl[1]!=0):
											x0_lower = overl[0]
											x0_lower = UTCDateTime(x0_lower)
											x1_upper = overl[1]
											x1_upper = UTCDateTime(x1_upper)
											axs.axvspan(x0_lower, x1_upper, fc='red', alpha=0.3)
											axs.plot(x0_lower, axs.get_ylim()[1], marker=markers.CARETLEFT, color = 'red')
											axs.plot(x1_upper, axs.get_ylim()[1], marker=markers.CARETRIGHT, color= 'red')
											axs.annotate(
													"Overlap",
													xy=(x1_upper, axs.get_ylim()[1]), xytext=(-10, -10),
													textcoords='offset points', ha='right', va='bottom',size=6, color = "red")

							# Display SEIS Activity
							if raw_data_plot.show_activity == True:
								if List_Channels_to_process[idx].activities is not None or len(List_Channels_to_process[idx].activities)>0:
									cycol2 = cycle('gcmk')
									for act in List_Channels_to_process[idx].activities:
										if List_Channels_to_process[idx].detector == "SP":
											if "VBB" not in act.desc or "SP" in act.desc :
												col = next(cycol2)
												if "ON" in act.desc or "OFF" in act.desc:
													#print(act.starttime, act.endtime, act.desc)
													if "OFF" in act.desc:
														axs.plot(act.endtime, axs.get_ylim()[1], marker=markers.CARETLEFT, color = col)
													if "ON" in act.desc:
														axs.plot(act.endtime, axs.get_ylim()[1], marker=markers.CARETRIGHT, color= col)
													axs.vlines(act.endtime, 0, 1, transform=axs.get_xaxis_transform(), colors=col, alpha=0.3)
													axs.annotate(
															act.desc,
															xy=(act.endtime, axs.get_ylim()[1]), xytext=(-10, -10),
															textcoords='offset points', ha='right', va='top',size=5)

												else:
													col = next(cycol2)
													x0_lower = act.starttime
													x1_upper = act.endtime
													axs.axvspan(x0_lower, x1_upper, fc = col, alpha=0.5)
													axs.plot(act.starttime, axs.get_ylim()[1], marker=markers.CARETLEFT, color = col)
													axs.plot(act.endtime, axs.get_ylim()[1], marker=markers.CARETRIGHT, color = col)
													axs.annotate(
															act.desc,
															y=(x1_upper, axs.get_ylim()[1]), xytext=(-10, -10),
															textcoords='offset points', ha='left', va='top',size=6, color=col)

										if List_Channels_to_process[idx].detector == "VBB":
											if "SP" not in act.desc or "VBB" in act.desc :
												col = next(cycol2)
												if "ON" in act.desc or "OFF" in act.desc:
													if "OFF" in act.desc:
														axs.plot(act.endtime, axs.get_ylim()[1], marker=markers.CARETLEFT, color = col)
													if "ON" in act.desc:
														axs.plot(act.endtime, axs.get_ylim()[1], marker=markers.CARETRIGHT, color= col)
													axs.vlines(act.endtime, 0, 1, transform=axs.get_xaxis_transform(), colors=col, alpha=0.3)
													axs.annotate(
															act.desc,
															xy=(act.endtime, axs.get_ylim()[1]), xytext=(-10, -10),
															textcoords='offset points', ha='right', va='top',size=5)

												else:
													col = next(cycol2)
													x0_lower = act.starttime
													x1_upper = act.endtime
													axs.axvspan(x0_lower, x1_upper, fc = col, alpha=0.5)
													axs.plot(act.starttime, axs.get_ylim()[1], marker=markers.CARETLEFT, color = col)
													axs.plot(act.endtime, axs.get_ylim()[1], marker=markers.CARETRIGHT, color = col)
													axs.annotate(
															act.desc,
															xy=(x1_upper, axs.get_ylim()[1]), xytext=(-10, -10),
															textcoords='offset points', ha='left', va='top',size=6, color=col)

										if List_Channels_to_process[idx].detector not in ["VBB", "SP"]:
											if "SP" not in act.desc and "VBB" not in act.desc :
												col = next(cycol2)
												if "ON" in act.desc or "OFF" in act.desc:
													if "OFF" in act.desc:
														axs.plot(act.endtime, axs.get_ylim()[1], marker=markers.CARETLEFT, color = col)
													if "ON" in act.desc:
														axs.plot(act.endtime, axs.get_ylim()[1], marker=markers.CARETRIGHT, color= col)
													axs.vlines(act.endtime, 0, 1, transform=axs.get_xaxis_transform(), colors=col, alpha=0.3)
													axs.annotate(
															act.desc,
															xy=(act.endtime, axs.get_ylim()[1]), xytext=(-10, -10),
															textcoords='offset points', ha='right', va='top',size=5)

												else:
													col = next(cycol2)
													x0_lower = act.starttime
													x1_upper = act.endtime
													axs.axvspan(x0_lower, x1_upper, fc = col, alpha=0.5)
													axs.plot(act.starttime, axs.get_ylim()[1], marker=markers.CARETLEFT, color = col)
													axs.plot(act.endtime, axs.get_ylim()[1], marker=markers.CARETRIGHT, color = col)
													axs.annotate(
															act.desc,
															xy=(x1_upper, axs.get_ylim()[1]), xytext=(-10, -10),
															textcoords='offset points', ha='left', va='top',size=6, color=col)

							#Y-axis options
							axs.yaxis.label.set_size(6)
							axs.yaxis.set_tick_params(labelsize=6)
							axs.set(ylabel=raw_data_plot.ytitle)
							# General plot options
							axs.legend(framealpha=0.01, loc='best', borderaxespad=0., fontsize="x-small")
							axs.grid(True)
							axs.tick_params(labelsize='x-small', width=2)
							figs.autofmt_xdate()
							suptitle += List_Channels_to_process[idx].longname+" "
							if len(suptitle)> 100:
								suptitle += "\n"
				figs.suptitle(raw_data_plot.suptitle +"\n"+suptitle, fontsize=8)
				# Display / Save options
				raw_data_plot.display = False
				if raw_data_plot.display:
					plt.show()

				if raw_data_plot.save:
					pdf_pages.savefig(figs)


	#     Transfer function plots
	###########################################@
	if tf:
		if dataless is not None:
			try:
				from obspy import read_inventory
				inv_seed = read_inventory(dataless)
			except:
				print("Oups")
			finally:
				print("Transfer functions plots generation.")
				min_freq = 0.001
				ddet = df.groupby(['detector'], sort= True)
				channel_list = []
				#print("Number of pages: ", pdf_pages.get_pagecount())
				for ndet, groupdet in ddet:
					#print("******* DETECTOR : {} ********".format(ndet))
					#print(groupdet)
					doutput = groupdet.groupby(['output'], sort= True)
					for detector,groupoutput in doutput:
						channel_list = []
						dsps = groupoutput.groupby(['sps'], sort = True)
						for no,groupsps in dsps:
							#print("sps: ", no)
							#print(groupsps)
							dgain = groupsps.groupby(['gaintype'], sort = True)
							for ng,groupgain in dgain:
								chan_dict={}
								chan_dict["channel"] = groupgain["channel"].iloc[0]
								chan_dict["station"] = groupgain["station"].iloc[0]
								chan_dict["location"] = groupgain["location"].iloc[0]
								chan_dict["network"] = groupgain["network"].iloc[0]
								chan_dict["output"] = groupgain["output"].iloc[0]
								chan_dict["gaintype"] = groupgain["gaintype"].iloc[0]
								chan_dict["longname"] = groupgain["longname"].iloc[0]
								chan_dict["detector"] = groupgain["detector"].iloc[0]
								chan_dict["sps"] = groupgain["sps"].iloc[0]
								channel_list.append(chan_dict)
								 #print(len(channel_list))
						if ndet in ["VBB", "SP"] and detector not in ["Temp"]:
							if len(channel_list)>0 and channel_list is not None:
								for output in ["DISP", "VEL", "ACC"]:
									figs, axs = plot_responses(inv_seed, channel_list, output, min_freq, starttime = None, endtime = None)
									if raw_data_plot.save:
										pdf_pages.savefig(figs)


	if raw_data_plot.save:
		if pdf_pages is not None:
			print("Final number of pages :", pdf_pages.get_pagecount(), " saved in ", raw_data_plot.savefile+'.pdf')
			pdf_pages.close()


#==============================================================================
#    DAILY REPORT ENG
#

def daily_report_eng(List_Channels_to_process, df_channels, suptitle, \
					 display_plots= False, sol = None, marsDate=None, \
					 tf=False, dataless=None, outfile = None, plotcount = True):
	"""
	Function which plots all the channels from a list. Transfer function can 
	also be added to the report.
	-----
	INPUT:
		@List_Channels_to_process: list of ASPICchannels object 
			- list of channels to process
		@raw_data_plot: plot_options_process object whici contains 
			the options of plot
		@df: pandas structures with all the details of the each channels
		@marsDate: MarsConverter object to convert date from UTC to LMST
		@tf: boolean to process transfert functions plots
		@dataless : inventory file to read the transfert functions
		@oudir : string with the output directory.

	"""

	global mDate
	mDate = marsDate
	

	from itertools import cycle
	#Import useful lib
	import matplotlib.pyplot as plt
	from matplotlib.ticker import FuncFormatter, MaxNLocator
	from matplotlib import markers
	from matplotlib.backends.backend_pdf import PdfPages
	from matplotlib.gridspec import GridSpec


	import numpy as np

	#sys.path.insert(0, './modules/')


	if outfile is not None:
		pdf_pages = PdfPages(outfile+'.pdf')
	else:
		import datetime
		tstamp = str(datetime.datetime.now().timestamp()).split('.')[0]

		pdf_pages = PdfPages("Reporting_"+tstamp+".pdf")


	if sol is not None:
		sol_starttime_utc = mDate.get_lmst_to_utc(lmst_date = int(sol)) 
		sol_endtime_utc   = mDate.get_lmst_to_utc(lmst_date = int(sol+1))

	#Sort the data by type of dectectors (SP, VBB, ... )
	#ds = df_channels.groupby(['detector'], sort= True)
	ds = df_channels.groupby(['plt_order'], sort= True)

	for name, group in ds:
		print("<<<<<<<<<<<<<<<<<<<<<<<<<")
		print("Detector: ", name)
		#print("<<<<<<<<<<<<<<<<<<<<<<<<<")
		#Sort the each group by output (Vel, Pos, Temp...)
		doutput = group.groupby(['output'], sort= True)
		page_nb = 1
		for n,groupoutput in doutput:
			print("******************")
			print("\tOutput: {}".format(n))
			#Sort each group of output by sps
			dsps = groupoutput.groupby(['sps'], sort = True)
			#print("*******************")

			#We change page at each new sps

			for no,groupsps in dsps:
				#Sort each group with the same SPS with orientation
				print("===========================")
				print('\t\tsps', no, "at page", page_nb )
				#print("=========")
				dunits = groupsps.groupby(['yunits'], sort = False)

				nb_counts=0
				nb_calib=0
				#Sort each groups according to the output units and the orientation
				chan2plot_per_unit = []
				dict_indexes = {}
				#dict_indexes["Counts"] = 1    # default value in case no raw data
				raw_data_plot = None
				for n_u, groupunits in dunits:
					print("++++++++++++++++")
					print("\t\t\tYunits: ", n_u)
					#print("++++++++++++++++")
					dict_indexes[n_u] = []
					dorient = groupunits.groupby(['orientation'], sort = True)
					nb_grp_orient = 0
					chan2plot_per_orient = []

					for n_o, grouporient in dorient:
						print("---------------------")
						print("\t\t\t\tOrientations: ", n_o)
						#print("---------------------")

						nb_grp_orient += 1
						chan2plot_per_orient = []
						for i in range(0, len(grouporient)):
							idx = grouporient.iloc[i]["index_in_list"]
							print("idx =", idx, List_Channels_to_process[idx].longname, List_Channels_to_process[idx].yunits )
							if List_Channels_to_process[idx].yunits == "Counts":
								nb_counts +=1
							else:
								nb_calib+=1
							chan2plot_per_orient.append(idx)
						dict_indexes[n_u].append(chan2plot_per_orient)
						print("number of different orientations: ", nb_grp_orient)
						chan2plot_per_unit.append(chan2plot_per_orient)

					print("\t\t\tchannels with unit ",n_u , chan2plot_per_unit)

				print("\t\tChannels in counts: {} and channels calibrated: {}".format(nb_counts, nb_calib))
				# Loop to get the list of units
				data_units = []
				for k,v in dict_indexes.items():
					#print(k,v)
					data_units.append(k)
				print("List des unités disponibles: ", data_units)
				print(dict_indexes)
				#plt.tight_layout()
				if "Counts" in dict_indexes.keys():
					nb_subplt_counts = len(dict_indexes["Counts"])
				else: # case with no Counts data
					nb_subplt_counts = len(dict_indexes[data_units[0]])
					dict_indexes["Counts"] = dict_indexes[data_units[0]]
					data_units.append(data_units[0])
				#print("suptitle", suptitle)
				
				if raw_data_plot is None:
					#print("dict_indexes['Counts'][0]", dict_indexes["Counts"][0])
					if List_Channels_to_process[dict_indexes["Counts"][0][0]].plot_options is not None:
						#print("dict_indexes['Counts'][0][0]", dict_indexes["Counts"][0][0])
						raw_data_plot = List_Channels_to_process[dict_indexes["Counts"][0][0]].plot_options
						raw_data_plot.suptitle = suptitle
						#print("NEW <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<")
						#print(raw_data_plot)
						#print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
					else: 
						print(List_Channels_to_process[dict_indexes["Counts"][0][0]])
						sys.exit(List_Channels_to_process[dict_indexes["Counts"][0][0]], "has no plot_options")
						
				else:
					raw_data_plot.suptitle = suptitle
					#print("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<")
					#print(raw_data_plot)
					#print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
				
				
				# Loop over units 
				# Left plot is supposed to be always in Counts
				# That's why u_idx starts at 1 and not 0 (0 is supposed to be for 
				# Counts)
				for u_idx in range(1, len(data_units)):
					print("Unit: ", data_units[u_idx])
					fig = plt.figure(figsize=(11.69,7))
					fig.subplots_adjust(hspace = 0)
					
					gs = GridSpec(nb_subplt_counts, 2)

					if raw_data_plot.share_x_axis:
						#sharex = 'col'
						if raw_data_plot.sol_scale:
							min_date_plot = raw_data_plot.forced_starttime
							max_date_plot = raw_data_plot.forced_endtime
						else:
							min_date_plot = None
							max_date_plot = None
					cycol = cycle('grcmk')
					for i in range(nb_subplt_counts):
						ax = fig.add_subplot(gs[i,0])
						cycol_gaps = cycle('by')
						# Plot data in counts
						for idx in dict_indexes["Counts"][i]:
							print("idx =", idx, List_Channels_to_process[idx].longname, List_Channels_to_process[idx].yunits )
							if List_Channels_to_process[idx].longname !='':
								nameChannel = List_Channels_to_process[idx].longname + " ("\
								 + List_Channels_to_process[idx].location + "/" +  \
								List_Channels_to_process[idx].channel +")"
							else:
								nameChannel = List_Channels_to_process[idx].location+ "/" +\
								 List_Channels_to_process[idx].channel

							start_time = List_Channels_to_process[idx].traces.stats.starttime
							t = np.arange(start_time, start_time + \
								 List_Channels_to_process[idx].traces.stats.npts/List_Channels_to_process[idx].traces.stats.sampling_rate, \
										 List_Channels_to_process[idx].traces.stats.delta)
							dte = [UTCDateTime(dat).timestamp for dat in t]
							ax.plot(dte, List_Channels_to_process[idx].traces.data, \
										color= next(cycol),label = nameChannel, \
										linestyle='-', marker='None', linewidth=0.6)

							if i < nb_subplt_counts-1:
								ax.xaxis.set_major_formatter(plt.NullFormatter())
							else:
								if raw_data_plot.utc_or_lmst == "LMST" and mDate is not None:
									ax.xaxis.set_major_formatter(FuncFormatter(lambda tstamp, pos : mDate.get_utc_2_lmst(UTCDateTime(tstamp), output="date")[:-7]))
									ax.set(xlabel="Date (LMST)")

								elif raw_data_plot.utc_or_lmst == "UTC":
									ax.xaxis.set_major_formatter(FuncFormatter(lambda tstamp, pos : str(UTCDateTime(tstamp))[:-8]))
									ax.set(xlabel="Date (UTC)")
								ax.xaxis.label.set_size(6)
							#ax.xaxis.set_visible(False)
							ax.xaxis.set_tick_params(labelsize=6, length=3, width=1)
							ax.xaxis.set_major_locator(MaxNLocator(12))
							ax.legend(framealpha=0.01, loc='best', borderaxespad=0., fontsize="x-small")
							ax.grid(True)
							ax.tick_params(labelsize='x-small', width=1)
							#t = ax.yaxis.get_offset_text()
							#t.set_size(30)
							plt.setp(ax.get_xticklabels(), rotation=30, ha='right')

							#if raw_data_plot.share_x_axis:
							#	if min_date_plot is not None and max_date_plot is not None:
							#		ax.set_xlim(UTCDateTime(min_date_plot).timestamp, UTCDateTime(max_date_plot).timestamp)
							#	else:
							#		ax.set_xlim(dte[0], dte[-1])
							#else:
							#	ax.set_xlim(dte[0], dte[-1])
							ax.set_xlim(sol_starttime_utc.timestamp,sol_endtime_utc.timestamp)
							# Display Gaps
							######################
							color_gaps = next(cycol_gaps)
							if raw_data_plot.show_gaps == True:
								if List_Channels_to_process[idx].gap_list is not None or len(List_Channels_to_process[idx].gap_list)>0:
									#
									for gap in List_Channels_to_process[idx].gap_list:
										if (gap[0]!=0 and gap[1]!=0):
											x0_lower = gap[0]
											x0_lower = UTCDateTime(x0_lower)
											x1_upper = gap[1]
											x1_upper = UTCDateTime(x1_upper)

											ax.axvspan(x0_lower, x1_upper, fc=color_gaps, alpha=0.3)
											ax.plot(x0_lower, ax.get_ylim()[1], marker=markers.CARETLEFT, color = color_gaps)
											ax.plot(x1_upper, ax.get_ylim()[1], marker=markers.CARETRIGHT, color= color_gaps)
											ax.annotate(
												"Gap ({}/{})".format(List_Channels_to_process[idx].location, List_Channels_to_process[idx].channel),
												xy=(x1_upper, ax.get_ylim()[1]), xytext=(-10, -10),
												textcoords='offset points', ha='right', va='bottom',size=6, color = color_gaps)


							# Display overlaps
							######################
							if raw_data_plot.show_overlaps == True:
								if List_Channels_to_process[idx].overlap_list is not None or len(List_Channels_to_process[idx].overlap_list)>0:
									for overl in List_Channels_to_process[idx].overlap_list:
										if (overl[0]!=0 and overl[1]!=0):
											x0_lower = overl[0]
											x0_lower = UTCDateTime(x0_lower)
											x1_upper = overl[1]
											x1_upper = UTCDateTime(x1_upper)

											ax.axvspan(x0_lower, x1_upper, fc='red', alpha=0.3)
											ax.plot(x0_lower, ax.get_ylim()[1], marker=markers.CARETLEFT, color = 'red')
											ax.plot(x1_upper, ax.get_ylim()[1], marker=markers.CARETRIGHT, color= 'red')
											ax.annotate(
													"Overlap",
													xy=(x1_upper, ax.get_ylim()[1]), xytext=(-10, -10),
													textcoords='offset points', ha='right', va='bottom',size=6, color = "red")
					#---------------------
					# Calibrated plots
					if raw_data_plot.ptype == "overplot":
						ax = fig.add_subplot(gs[:,1])
						cycol = cycle('grcmk')
						for i in range(nb_subplt_counts):
							for idx in dict_indexes[data_units[u_idx]][i]:
								print("idx =", idx)
								if List_Channels_to_process[idx].longname !='':
									nameChannel = List_Channels_to_process[idx].longname + " ("\
									 + List_Channels_to_process[idx].location + "/" +  \
									List_Channels_to_process[idx].channel +")"
								else:
									nameChannel = List_Channels_to_process[idx].location+ "/" +\
									 List_Channels_to_process[idx].channel

								start_time = List_Channels_to_process[idx].traces.stats.starttime
								t = np.arange(start_time, start_time + List_Channels_to_process[idx].traces.stats.npts/List_Channels_to_process[idx].traces.stats.sampling_rate, \
											 List_Channels_to_process[idx].traces.stats.delta)
								dte = [UTCDateTime(dat).timestamp for dat in t]
								ax.plot(dte, List_Channels_to_process[idx].traces.data, color= next(cycol), label = nameChannel, linestyle='-', marker='None', linewidth=0.6)
						if raw_data_plot.utc_or_lmst == "LMST" and mDate is not None:
							ax.xaxis.set_major_formatter(FuncFormatter(lambda tstamp, pos : mDate.get_utc_2_lmst(UTCDateTime(tstamp), output="date")[:-7]))
							ax.xaxis.set_label("Date (LMST)")
							#ax.set(xlabel="Date (LMST)")
						elif raw_data_plot.utc_or_lmst == "UTC":
							ax.xaxis.set_major_formatter(FuncFormatter(lambda tstamp, pos : str(UTCDateTime(tstamp))[:-8]))
							ax.set(xlabel="Date (UTC)")

						ax.xaxis.label.set_size(6)
						ax.xaxis.set_tick_params(rotation=30, labelsize=6, length=3, width=1)
						ax.xaxis.set_major_locator(MaxNLocator(12))
						ax.legend(framealpha=0.01, loc='best', borderaxespad=0., fontsize="x-small")
						ax.grid(True)
						ax.tick_params(labelsize='x-small', width=1)
						t = ax.yaxis.get_offset_text()
						t.set_size(30)
						plt.setp(ax.get_xticklabels(), rotation=30, ha='right')
						# Display Gaps
						######################
						color_gaps = next(cycol_gaps)
						if raw_data_plot.show_gaps == True:
							if List_Channels_to_process[idx].gap_list is not None or len(List_Channels_to_process[idx].gap_list)>0:
								#
								for gap in List_Channels_to_process[idx].gap_list:
									if (gap[0]!=0 and gap[1]!=0):
										x0_lower = gap[0]
										x0_lower = UTCDateTime(x0_lower)
										x1_upper = gap[1]
										x1_upper = UTCDateTime(x1_upper)

										ax.axvspan(x0_lower, x1_upper, fc=color_gaps, alpha=0.3)
										ax.plot(x0_lower, ax.get_ylim()[1], marker=markers.CARETLEFT, color = color_gaps)
										ax.plot(x1_upper, ax.get_ylim()[1], marker=markers.CARETRIGHT, color= color_gaps)
										ax.annotate(
											"Gap ({}/{})".format(List_Channels_to_process[idx].location, List_Channels_to_process[idx].channel),
											xy=(x1_upper, ax.get_ylim()[1]), xytext=(-10, -10),
											textcoords='offset points', ha='right', va='bottom',size=6, color = color_gaps)
						# Display overlaps
						######################
						if raw_data_plot.show_overlaps == True:
							if List_Channels_to_process[idx].overlap_list is not None or len(List_Channels_to_process[idx].overlap_list)>0:
								for overl in List_Channels_to_process[idx].overlap_list:
									if (overl[0]!=0 and overl[1]!=0):
										x0_lower = overl[0]
										x0_lower = UTCDateTime(x0_lower)
										x1_upper = overl[1]
										x1_upper = UTCDateTime(x1_upper)

										ax.axvspan(x0_lower, x1_upper, fc='red', alpha=0.3)
										ax.plot(x0_lower, ax.get_ylim()[1], marker=markers.CARETLEFT, color = 'red')
										ax.plot(x1_upper, ax.get_ylim()[1], marker=markers.CARETRIGHT, color= 'red')
										ax.annotate(
												"Overlap",
												xy=(x1_upper, ax.get_ylim()[1]), xytext=(-10, -10),
												textcoords='offset points', ha='right', va='bottom',size=6, color = "red")


					else: #no overlaping plots on the right
						cycol = cycle('grcmk')
						print("No overlaping in the right")
						for i in range(nb_subplt_counts):
							ax = fig.add_subplot(gs[i,1])
							cycol_gaps = cycle('by')
							for idx in dict_indexes[data_units[u_idx]][i]:
								print("idx =", idx)
								if List_Channels_to_process[idx].longname !='':
									nameChannel = List_Channels_to_process[idx].longname + " ("\
									 + List_Channels_to_process[idx].location + "/" +  \
									List_Channels_to_process[idx].channel +")"
								else:
									nameChannel = List_Channels_to_process[idx].location+ "/" +\
									 List_Channels_to_process[idx].channel

								start_time = List_Channels_to_process[idx].traces.stats.starttime
								t = np.arange(start_time, start_time + List_Channels_to_process[idx].traces.stats.npts/List_Channels_to_process[idx].traces.stats.sampling_rate, \
											 List_Channels_to_process[idx].traces.stats.delta)
								dte = [UTCDateTime(dat).timestamp for dat in t]
								ax.plot(dte, List_Channels_to_process[idx].traces.data, \
										color= next(cycol),  label = nameChannel, \
										linestyle='-', marker='None', linewidth=0.6)

								if i < nb_subplt_counts-1:
									ax.xaxis.set_major_formatter(plt.NullFormatter())
								else:
									if raw_data_plot.utc_or_lmst == "LMST" and mDate is not None:
										ax.xaxis.set_major_formatter(FuncFormatter(lambda tstamp, pos : mDate.get_utc_2_lmst(UTCDateTime(tstamp), output="date")[:-7]))
										ax.set(xlabel="Date (LMST)")

									elif raw_data_plot.utc_or_lmst == "UTC":
										ax.xaxis.set_major_formatter(FuncFormatter(lambda tstamp, pos : str(UTCDateTime(tstamp))[:-8]))
										ax.set(xlabel="Date (UTC)")
									ax.xaxis.label.set_size(6)
								#ax.xaxis.set_visible(False)
								ax.xaxis.set_tick_params(labelsize=6, length=3, width=1)
								ax.xaxis.set_major_locator(MaxNLocator(12))
								ax.legend(framealpha=0.01, loc='best', borderaxespad=0., fontsize="x-small")
								ax.grid(True)
								ax.tick_params(labelsize='x-small', width=1)
								#t = ax.yaxis.get_offset_text()
								#t.set_size(30)
								plt.setp(ax.get_xticklabels(), rotation=30, ha='right')

#								if raw_data_plot.share_x_axis:
#									if 	min_date_plot is not None and max_date_plot is not None:
#										ax.set_xlim(UTCDateTime(min_date_plot).timestamp, UTCDateTime(max_date_plot).timestamp)
#									else:
#										ax.set_xlim(dte[0], dte[-1])
#								else:
#									ax.set_xlim(dte[0], dte[-1])
								ax.set_xlim(sol_starttime_utc.timestamp,sol_endtime_utc.timestamp)

								# Display Gaps
								######################
								color_gaps = next(cycol_gaps)
								if raw_data_plot.show_gaps == True:
									if List_Channels_to_process[idx].gap_list is not None or len(List_Channels_to_process[idx].gap_list)>0:
										#
										for gap in List_Channels_to_process[idx].gap_list:
											if (gap[0]!=0 and gap[1]!=0):
												x0_lower = gap[0]
												x0_lower = UTCDateTime(x0_lower)
												x1_upper = gap[1]
												x1_upper = UTCDateTime(x1_upper)

												ax.axvspan(x0_lower, x1_upper, fc=color_gaps, alpha=0.3)
												ax.plot(x0_lower, ax.get_ylim()[1], marker=markers.CARETLEFT, color = color_gaps)
												ax.plot(x1_upper, ax.get_ylim()[1], marker=markers.CARETRIGHT, color= color_gaps)
												ax.annotate(
													"Gap ({}/{})".format(List_Channels_to_process[idx].location, List_Channels_to_process[idx].channel),
													xy=(x1_upper, ax.get_ylim()[1]), xytext=(-10, -10),
													textcoords='offset points', ha='right', va='bottom',size=6, color = color_gaps)


								# Display overlaps
								######################
								if raw_data_plot.show_overlaps == True:
									if List_Channels_to_process[idx].overlap_list is not None or len(List_Channels_to_process[idx].overlap_list)>0:
										for overl in List_Channels_to_process[idx].overlap_list:
											if (overl[0]!=0 and overl[1]!=0):
												x0_lower = overl[0]
												x0_lower = UTCDateTime(x0_lower)
												x1_upper = overl[1]
												x1_upper = UTCDateTime(x1_upper)

												ax.axvspan(x0_lower, x1_upper, fc='red', alpha=0.3)
												ax.plot(x0_lower, ax.get_ylim()[1], marker=markers.CARETLEFT, color = 'red')
												ax.plot(x1_upper, ax.get_ylim()[1], marker=markers.CARETRIGHT, color= 'red')
												ax.annotate(
														"Overlap",
														xy=(x1_upper, ax.get_ylim()[1]), xytext=(-10, -10),
														textcoords='offset points', ha='right', va='bottom',size=6, color = "red")


					#fig.autofmt_xdate()
					fig.text(0.06, 0.5, data_units[0], ha='center', va='center', rotation='vertical')
					fig.text(0.5, 0.5, data_units[u_idx], ha='center', va='center', rotation='vertical')

					#				# To add the right panel
					#				ax_right = fig.add_subplot(gs[:, 1])
					#				ax_right.xaxis.set_label("LMST de fou")
					#
					fig.suptitle(raw_data_plot.suptitle)
					#plt.show()

					if display_plots:
						plt.show()

					if outfile is not None:
						pdf_pages.savefig(fig)

	if outfile is not None:
		if pdf_pages is not None:
			print("Final number of pages :", pdf_pages.get_pagecount(), " saved in ", outfile+'.pdf')
			pdf_pages.close()

#==============================================================================
#    WEEKLY REPORT ENG
#

def weekly_report(List_Channels_to_process, df_channels, suptitle, \
					 display_plots= False, sols = None, marsDate=None, \
					 tf=False, dataless=None, outfile = None, plotcount = True):
	"""
	Function which plots all the channels from a list. Transfer function can 
	also be added to the report.
	-----
	INPUT:
		@List_Channels_to_process: list of ASPICchannels object 
			- list of channels to process
		@raw_data_plot: plot_options_process object whici contains 
			the options of plot
		@df_channels: pandas structures with all the details of the each channels
		@marsDate: MarsConverter object to convert date from UTC to LMST
		@tf: boolean to process transfert functions plots
		@dataless : inventory file to read the transfert functions
		@oudir : string with the output directory.

	"""

	global mDate
	mDate = marsDate
	
	from itertools import cycle
	#Import useful lib
	import matplotlib.pyplot as plt
	from matplotlib.ticker import FuncFormatter, MaxNLocator
	from matplotlib import markers
	from matplotlib.backends.backend_pdf import PdfPages
	from matplotlib.gridspec import GridSpec

	import math

	import numpy as np

	#sys.path.insert(0, './modules/')


	if outfile is not None:
		pdf_pages = PdfPages(outfile+'.pdf')
	else:
		import datetime
		tstamp = str(datetime.datetime.now().timestamp()).split('.')[0]

		pdf_pages = PdfPages("Reporting_"+tstamp+".pdf")

	if sols is not None:
		sol_starttime_utc = mDate.get_lmst_to_utc(lmst_date = int(sols[0])) 
		sol_endtime_utc   = mDate.get_lmst_to_utc(lmst_date = int(sols[1]))

	#Sort the data by type of dectectors (SP, VBB, ... )
	#ds = df_channels.groupby(['detector'], sort= True)
	ds = df_channels.groupby(['plt_order'], sort= True)

	for name, group in ds:
		print("<<<<<<<<<<<<<<<<<<<<<<<<<")
		print("Detector: ", name)
		#print("<<<<<<<<<<<<<<<<<<<<<<<<<")
		#Sort the each group by output (Vel, Pos, Temp...)

		doutput = group.groupby(['output'], sort= True)
		page_nb = 1
		for n,groupoutput in doutput:
			print("******************")
			print("\tOutput: {}".format(n))
			#Sort each group of output by sps
			dsps = groupoutput.groupby(['sps'], sort = True)
			#print("*******************")

			#We change page at each new sps

			for no,groupsps in dsps:
				#Sort each group with the same SPS with orientation
				print("===========================")
				print('\t\tsps', no, "at page", page_nb )
				#print("=========")
				dunits = groupsps.groupby(['yunits'], sort = False)

				nb_counts=0
				nb_calib=0
				#Sort each groups according to the output units and the orientation
				chan2plot_per_unit = []
				dict_indexes = {}
				#dict_indexes["Counts"] = 1    # default value in case no raw data
				raw_data_plot = None
				nb_subplt = 0
				for n_u, groupunits in dunits:
					print("++++++++++++++++")
					print("\t\t\tYunits: ", n_u)
					#print("++++++++++++++++")
					dict_indexes[n_u] = []
					dorient = groupunits.groupby(['orientation'], sort = True)
					nb_grp_orient = 0
					chan2plot_per_orient = []

					for n_o, grouporient in dorient:
						print("---------------------")
						print("\t\t\t\tOrientations: ", n_o)
						#print("---------------------")

						nb_grp_orient += 1
						chan2plot_per_orient = []
						for i in range(0, len(grouporient)):
							idx = grouporient.iloc[i]["index_in_list"]
							print("idx =", idx, List_Channels_to_process[idx].longname, List_Channels_to_process[idx].yunits )
							if List_Channels_to_process[idx].yunits == "Counts":
								nb_counts +=1
							else:
								nb_calib+=1
							chan2plot_per_orient.append(idx)
						dict_indexes[n_u].append(chan2plot_per_orient)
						print("number of different orientations: ", nb_grp_orient)
						chan2plot_per_unit.append(chan2plot_per_orient)

					print("\t\t\tchannels with unit ",n_u , chan2plot_per_unit)

				print("\t\tChannels calibrated: {}".format(nb_calib))
				# Loop to get the list of units
				data_units = []
				for k,v in dict_indexes.items():
					#print(k,v)
					data_units.append(k)
				print("List des unités disponibles: ", data_units)
				print(dict_indexes)


				# Loop over units 
				for u_idx in data_units:
					print("Unit: ", u_idx, " number of channels for this unit: ", len(dict_indexes.get(u_idx)))
					if raw_data_plot is None:
						if List_Channels_to_process[dict_indexes[u_idx][0][0]].plot_options is not None:
							raw_data_plot = List_Channels_to_process[dict_indexes[u_idx][0][0]].plot_options
							raw_data_plot.suptitle = suptitle
						else: 
							print(List_Channels_to_process[dict_indexes[u_idx][0][0]])
							sys.exit(List_Channels_to_process[dict_indexes[u_idx][0][0]], "has no plot_options")
					else:
						raw_data_plot.suptitle = suptitle
					
					nb_subplt_counts = len(dict_indexes.get(u_idx))
					fig = plt.figure(figsize=(11.69,7))
					fig.subplots_adjust(hspace = 0)
					
					gs = GridSpec(nb_subplt_counts, 1)  #
					share_x_axis = True
					if share_x_axis:
						#sharex = 'col'
						solscale = True
						if solscale:
							min_date_plot = raw_data_plot.forced_starttime
							max_date_plot = raw_data_plot.forced_endtime
						else:
							min_date_plot = None
							max_date_plot = None

					if raw_data_plot.ptype == "overplot":
						ax = fig.add_subplot(gs[:,0])
						cycol = cycle('grcmk')
						for i in range(nb_subplt_counts):
							for idx in dict_indexes[u_idx][i]:
								print("idx =", idx)
								if List_Channels_to_process[idx].longname !='':
									nameChannel = List_Channels_to_process[idx].longname + " ("\
									 + List_Channels_to_process[idx].location + "/" +  \
									List_Channels_to_process[idx].channel +")"
								else:
									nameChannel = List_Channels_to_process[idx].location+ "/" +\
									 List_Channels_to_process[idx].channel

								start_time = List_Channels_to_process[idx].traces.stats.starttime
								t = np.arange(start_time, start_time + \
									List_Channels_to_process[idx].traces.stats.npts/List_Channels_to_process[idx].traces.stats.sampling_rate, \
											 List_Channels_to_process[idx].traces.stats.delta)
								dte = [UTCDateTime(dat).timestamp for dat in t]
								ax.plot(dte, List_Channels_to_process[idx].traces.data, color= next(cycol), label = nameChannel, linestyle='-', marker='None', linewidth=0.6)
						if raw_data_plot.utc_or_lmst == "LMST" and mDate is not None:
							ax.xaxis.set_major_formatter(FuncFormatter(lambda tstamp, pos : mDate.get_utc_2_lmst(UTCDateTime(tstamp), output="date")[:-7]))
							ax.xaxis.set_label("Date (LMST)")
							#ax.set(xlabel="Date (LMST)")
						elif raw_data_plot.utc_or_lmst == "UTC":
							ax.xaxis.set_major_formatter(FuncFormatter(lambda tstamp, pos : str(UTCDateTime(tstamp))[:-8]))
							ax.set(xlabel="Date (UTC)")

						ax.xaxis.label.set_size(6)
						ax.xaxis.set_tick_params(rotation=30, labelsize=6, length=3, width=1)
						ax.xaxis.set_major_locator(MaxNLocator(12))
						ax.legend(framealpha=0.01, loc='best', borderaxespad=0., fontsize="x-small")
						ax.grid(True)
						ax.tick_params(labelsize='x-small', width=1)
						ax.set_xlim(sol_starttime_utc.timestamp,sol_endtime_utc.timestamp)
						t = ax.yaxis.get_offset_text()
						t.set_size(30)
						plt.setp(ax.get_xticklabels(), rotation=30, ha='right')
						# Display Gaps
						######################
						cycol_gaps = cycle('by')
						color_gaps = next(cycol_gaps)

						if raw_data_plot.show_gaps == True:
							if List_Channels_to_process[idx].gap_list is not None or len(List_Channels_to_process[idx].gap_list)>0:
								#
								for gap in List_Channels_to_process[idx].gap_list:
									if (gap[0]!=0 and gap[1]!=0):
										x0_lower = gap[0]
										x0_lower = UTCDateTime(x0_lower)
										x1_upper = gap[1]
										x1_upper = UTCDateTime(x1_upper)

										ax.axvspan(x0_lower, x1_upper, fc=color_gaps, alpha=0.3)
										ax.plot(x0_lower, ax.get_ylim()[1], marker=markers.CARETLEFT, color = color_gaps)
										ax.plot(x1_upper, ax.get_ylim()[1], marker=markers.CARETRIGHT, color= color_gaps)
										ax.annotate(
											"Gap ({}/{})".format(List_Channels_to_process[idx].location, List_Channels_to_process[idx].channel),
											xy=(x1_upper, ax.get_ylim()[1]), xytext=(-10, -10),
											textcoords='offset points', ha='right', va='bottom',size=6, color = color_gaps)
						# Display overlaps
						######################
						if raw_data_plot.show_overlaps == True:
							if List_Channels_to_process[idx].overlap_list is not None or len(List_Channels_to_process[idx].overlap_list)>0:
								for overl in List_Channels_to_process[idx].overlap_list:
									if (overl[0]!=0 and overl[1]!=0):
										x0_lower = overl[0]
										x0_lower = UTCDateTime(x0_lower)
										x1_upper = overl[1]
										x1_upper = UTCDateTime(x1_upper)

										ax.axvspan(x0_lower, x1_upper, fc='red', alpha=0.3)
										ax.plot(x0_lower, ax.get_ylim()[1], marker=markers.CARETLEFT, color = 'red')
										ax.plot(x1_upper, ax.get_ylim()[1], marker=markers.CARETRIGHT, color= 'red')
										ax.annotate(
												"Overlap",
												xy=(x1_upper, ax.get_ylim()[1]), xytext=(-10, -10),
												textcoords='offset points', ha='right', va='bottom',size=6, color = "red")

						# Enhance sol
						######################
						#min_date_plot=min([chan.stats.starttime.timestamp for chan in channelBH])
						#max_date_plot=max([chan.stats.endtime.timestamp for chan in channelBH])
					
						#min_UTC = min([chan.stats.starttime for chan in channelBH])
						#max_UTC = max([chan.stats.endtime for chan in channelBH])
						
						first_SOL = mDate.get_utc_2_lmst(utc_date=min_date_plot, output="decimal")
						first_SOL = math.modf(first_SOL)[1]
					
						last_SOL = mDate.get_utc_2_lmst(utc_date=max_date_plot, output="decimal")
						last_SOL = math.modf(last_SOL)[1]
					
						list_utc_sol = []
						list_utc_sol.append(sol_starttime_utc)
						for sol in range(int(sols[0]), int(sols[1])):
							list_utc_sol.append(mDate.get_lmst_to_utc(lmst_date=sol).timestamp)
					
						list_utc_sol.append(sol_endtime_utc)
						print("from {} to {}".format(first_SOL, last_SOL))
						for i in range(0, len(list_utc_sol)-1, 2):
							ax.axvspan(list_utc_sol[i], list_utc_sol[i+1], fc='grey', alpha=0.02)
						

					else: #no overlaping plots on the right
						cycol = cycle('grcmk')
						print("No overlaping in the right")
						for i in range(nb_subplt_counts):
							ax = fig.add_subplot(gs[i,0])
							cycol_gaps = cycle('by')
							for idx in dict_indexes[u_idx][i]:
								print("idx =", idx)
								if List_Channels_to_process[idx].longname !='':
									nameChannel = List_Channels_to_process[idx].longname + " ("\
									 + List_Channels_to_process[idx].location + "/" +  \
									List_Channels_to_process[idx].channel +")"
								else:
									nameChannel = List_Channels_to_process[idx].location+ "/" +\
									 List_Channels_to_process[idx].channel

								start_time = List_Channels_to_process[idx].traces.stats.starttime
								end_time   = List_Channels_to_process[idx].traces.stats.endtime 
								t = np.arange(start_time, start_time + List_Channels_to_process[idx].traces.stats.npts/List_Channels_to_process[idx].traces.stats.sampling_rate, \
											 List_Channels_to_process[idx].traces.stats.delta)
								dte = [UTCDateTime(dat).timestamp for dat in t]
								ax.plot(dte, List_Channels_to_process[idx].traces.data, \
										color= next(cycol),  label = nameChannel, \
										linestyle='-', marker='None', linewidth=0.6)

								if i < nb_subplt_counts-1:
									ax.xaxis.set_major_formatter(plt.NullFormatter())
								else:
									if raw_data_plot.utc_or_lmst == "LMST" and mDate is not None:
										ax.xaxis.set_major_formatter(FuncFormatter(lambda tstamp, pos : mDate.get_utc_2_lmst(UTCDateTime(tstamp), output="date")[:-7]))
										ax.set(xlabel="Date (LMST)")

									elif raw_data_plot.utc_or_lmst == "UTC":
										ax.xaxis.set_major_formatter(FuncFormatter(lambda tstamp, pos : str(UTCDateTime(tstamp))[:-8]))
										ax.set(xlabel="Date (UTC)")
									ax.xaxis.label.set_size(6)
								#ax.xaxis.set_visible(False)
								ax.xaxis.set_tick_params(labelsize=6, length=3, width=1)
								ax.xaxis.set_major_locator(MaxNLocator(12))
								ax.legend(framealpha=0.01, loc='best', borderaxespad=0., fontsize="x-small")
								ax.grid(True)
								ax.tick_params(labelsize='x-small', width=1)
								#t = ax.yaxis.get_offset_text()
								#t.set_size(30)
								plt.setp(ax.get_xticklabels(), rotation=30, ha='right')

								if raw_data_plot.share_x_axis:
									if 	min_date_plot is not None and max_date_plot is not None:
										ax.set_xlim(UTCDateTime(sol_starttime_utc).timestamp, UTCDateTime(sol_endtime_utc).timestamp)
									else:
										ax.set_xlim(dte[0], dte[-1])
								else:
									ax.set_xlim(dte[0], dte[-1])
								#ax.set_xlim(sol_starttime_utc.timestamp,sol_endtime_utc.timestamp)

								# Display Gaps
								######################
								color_gaps = next(cycol_gaps)
								if raw_data_plot.show_gaps == True:
									if List_Channels_to_process[idx].gap_list is not None or len(List_Channels_to_process[idx].gap_list)>0:
										#
										for gap in List_Channels_to_process[idx].gap_list:
											if (gap[0]!=0 and gap[1]!=0):
												x0_lower = gap[0]
												x0_lower = UTCDateTime(x0_lower)
												x1_upper = gap[1]
												x1_upper = UTCDateTime(x1_upper)

												ax.axvspan(x0_lower, x1_upper, fc=color_gaps, alpha=0.3)
												ax.plot(x0_lower, ax.get_ylim()[1], marker=markers.CARETLEFT, color = color_gaps)
												ax.plot(x1_upper, ax.get_ylim()[1], marker=markers.CARETRIGHT, color= color_gaps)
												ax.annotate(
													"Gap ({}/{})".format(List_Channels_to_process[idx].location, List_Channels_to_process[idx].channel),
													xy=(x1_upper, ax.get_ylim()[1]), xytext=(-10, -10),
													textcoords='offset points', ha='right', va='bottom',size=6, color = color_gaps)


								# Display overlaps
								######################
								if raw_data_plot.show_overlaps == True:
									if List_Channels_to_process[idx].overlap_list is not None or len(List_Channels_to_process[idx].overlap_list)>0:
										for overl in List_Channels_to_process[idx].overlap_list:
											if (overl[0]!=0 and overl[1]!=0):
												x0_lower = overl[0]
												x0_lower = UTCDateTime(x0_lower)
												x1_upper = overl[1]
												x1_upper = UTCDateTime(x1_upper)

												ax.axvspan(x0_lower, x1_upper, fc='red', alpha=0.3)
												ax.plot(x0_lower, ax.get_ylim()[1], marker=markers.CARETLEFT, color = 'red')
												ax.plot(x1_upper, ax.get_ylim()[1], marker=markers.CARETRIGHT, color= 'red')
												ax.annotate(
														"Overlap",
														xy=(x1_upper, ax.get_ylim()[1]), xytext=(-10, -10),
														textcoords='offset points', ha='right', va='bottom',size=6, color = "red")
								
								# Enhance sol
								######################
								#min_date_plot=min([chan.stats.starttime.timestamp for chan in channelBH])
								#max_date_plot=max([chan.stats.endtime.timestamp for chan in channelBH])
							
								#min_UTC = min([chan.stats.starttime for chan in channelBH])
								#max_UTC = max([chan.stats.endtime for chan in channelBH])
								
								first_SOL = mDate.get_utc_2_lmst(utc_date=min_date_plot, output="decimal")
								first_SOL = math.modf(first_SOL)[1]
							
								last_SOL = mDate.get_utc_2_lmst(utc_date=max_date_plot, output="decimal")
								last_SOL = math.modf(last_SOL)[1]
							
								list_utc_sol = []
								list_utc_sol.append(sol_starttime_utc)
								for sol in range(int(sols[0]), int(sols[1])):
									list_utc_sol.append(mDate.get_lmst_to_utc(lmst_date=sol).timestamp)
							
								list_utc_sol.append(sol_endtime_utc)
								print("from {} to {}".format(first_SOL, last_SOL))
								for i in range(0, len(list_utc_sol)-1, 2):
									ax.axvspan(list_utc_sol[i], list_utc_sol[i+1], fc='grey', alpha=0.02)

					#fig.autofmt_xdate()
					fig.text(0.06, 0.5, u_idx, ha='center', va='center', rotation='vertical')
					#fig.text(0.5, 0.5, u_idx, ha='center', va='center', rotation='vertical')

					#				# To add the right panel
					#				ax_right = fig.add_subplot(gs[:, 1])
					#				ax_right.xaxis.set_label("LMST de fou")
					#
					fig.suptitle(raw_data_plot.suptitle)
					#plt.show()

					if display_plots:
						plt.show()

					if outfile is not None:
						pdf_pages.savefig(fig)

	if outfile is not None:
		if pdf_pages is not None:
			print("Final number of pages :", pdf_pages.get_pagecount(), " saved in ", outfile+'.pdf')
			pdf_pages.close()
