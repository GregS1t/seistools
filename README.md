# SEIS TOOLBOX

This repository is intended to host all the codes useful for the calibration and analysis of data from the NASA InSight mission and in particular the SEIS seismometer. 

All mission data are in SEED or MiniSEED format. 

`Feel free to add your own codes to share with the community.`


## Dependencies
All the programs developed are based on the use of the **[Obspy](https://docs.obspy.org/)** library.

Moreover, the codes allow calculations in Martian time. A library has also been developed to perform these conversions: **[MarsConverter](https://github.com/GregS1t/marstimeconverter)**

## Installation

First, clone this repository on your computer.
After installing the code you have to add the path of the directory into you PYTHONPATH

For example: 
```
export SEIS_TB=~/CODE/seistools/
export PYTHONPATH="$SEIS_TB/modules:$PYTHONPATH"
```

## Content
### NB_Data_Calibration.ipynb

Simple notebook to calibrate some channels (detick, detrend, filtering, calibration, rotation)
It also allows you to plot the time series, the spectrograms and ASD/PSD.

### NB_Deglitch_data_JS.ipynb

This notebook is based in the John Scholz algorithm `SEISGLITCH`.
It allows user to load the data from the SEIS DATA PORTAL and process the deglitching.
One need to have a SEIS Data Portal Account or use the Public portal to reach the mseed files.
Dataless is retrieved directly from the SEISGLITCH code.  
This part is connected to the SEISGLITCH code. 

