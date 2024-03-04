import xarray as xr
import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter
# Set font style to match latex document----------
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size':16})
# ------------------------------------------------
import cartopy.crs as ccrs
import datetime
from functions import *

rpath="/projects/NS9600K/astridbg/data/model/noresm_postprocessed/"
wpath="/projects/NS9600K/astridbg/master/figures/model/spatavg/"

# Default cases----------------
#case1 = "def_20210126"; case1nm = "CAM6"
case1 = "meyers92_20220210"; case1nm = "CAM5"
# Modified cases---------------
#case2 = "meyers92_20220210"; case2nm = "CAM5"
case2 = "andenes21_20220222"; case2nm = "Andenes 2021"
#------------------------------	
date1 = "2007-04-15_2010-03-15"
date2 = "2007-04-15_2010-03-15"

#------------------------------
# Two-dimensional fields
#------------------------------

#variables = ["SWCF","LWCF","NETCF","CLDTOT","CLDHGH","CLDMED","CLDLOW","TGCLDIWP","TGCLDLWP","TREFHT"]
variables = ["SWCFS","LWCFS"]

#------------------------------
# Areas to average over
#------------------------------

svalbard = [[9,28],[75,81]]
quttinirpaaq = [[-120+360,-60+360],[78,86]]
"""
#------------------------------
# Plotting whole period averages
#------------------------------
for var in variables:
	print(var)
	ds1 = xr.open_dataset(rpath+var+"_"+case1+"_"+date1+".nc")
	ds2 = xr.open_dataset(rpath+var+"_"+case2+"_"+date2+".nc")
	
	# Get start and end date of period
	date_start = str(ds1.time[0].values).split(" ")[0]
	date_end = str(ds1.time[-1].values).split(" ")[0]
	
	# Make time array
	datetimeindex = ds1.indexes["time"].to_datetimeindex()	

	# Get spatial average of cases over Arctic
	ds1_arct = ds1[var].sel(lat=slice(66.5,90)).mean("lat").mean("lon")
	ds2_arct = ds2[var].sel(lat=slice(66.5,90)).mean("lat").mean("lon")

	# Get spatial average of cases globally
	ds1_glob = ds1[var].mean("lat").mean("lon")
	ds2_glob = ds2[var].mean("lat").mean("lon")

        # Get spatial average of Svalbard
	ds1_sval = ds1[var].sel(lon=slice(svalbard[0][0],svalbard[0][1]),
                                 lat=slice(svalbard[1][0],svalbard[1][1])).mean("lat").mean("lon")
	ds2_sval = ds2[var].sel(lon=slice(svalbard[0][0],svalbard[0][1]),
                                 lat=slice(svalbard[1][0],svalbard[1][1])).mean("lat").mean("lon")

        # Get spatial average over Quttinirpaaq
	ds1_qutt = ds1[var].sel(lon=slice(quttinirpaaq[0][0],quttinirpaaq[0][1]),
                                 lat=slice(quttinirpaaq[1][0],quttinirpaaq[1][1])).mean("lat").mean("lon")
	ds2_qutt = ds2[var].sel(lon=slice(quttinirpaaq[0][0],quttinirpaaq[0][1]),
                                 lat=slice(quttinirpaaq[1][0],quttinirpaaq[1][1])).mean("lat").mean("lon")
	
	fig = plt.figure(1, figsize=[15,6],dpi=300)
	ax = plt.subplot(1,1,1)
	ax.set_title(ds1[var].long_name+" "+case2nm+"-"+case1nm+" "+date_start+"-"+date_end, fontsize=22)
	
	ax.plot(datetimeindex, ds1_arct, color="navy", linestyle="--", label=case1nm+" Arctic")
	ax.plot(datetimeindex, ds1_glob, color="orangered", linestyle="--", label=case1nm+" Global")
	ax.plot(datetimeindex, ds1_sval, color="darkgreen", linestyle="--", label=case1nm+" Svalbard")
	ax.plot(datetimeindex, ds1_qutt, color="darkmagenta", linestyle="--", label=case1nm+" Quttinirpaaq")
	ax.plot(datetimeindex, ds2_arct, color="cornflowerblue", label=case2nm+" Arctic")
	ax.plot(datetimeindex, ds2_glob, color="orange", label=case2nm+" Global")
	ax.plot(datetimeindex, ds2_sval, color="mediumseagreen", label=case2nm+" Svalbard")
	ax.plot(datetimeindex, ds2_qutt, color="magenta", label=case2nm+" Quttinirpaaq")
	
	ax.set_ylabel(ds1[var].units, fontsize=18)

	# Shrink current axis's height by 10% on the bottom
	box = ax.get_position()
	ax.set_position([box.x0, box.y0 + box.height * 0.1,
                 box.width, box.height * 0.9])

	# Put a legend below current axis
	ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.08), ncol=4)

	plt.grid(alpha=0.5)
	plt.savefig(wpath+var+"_"+case1+"_"+case2+".pdf",bbox_inches="tight")
	
	plt.clf()
"""

#------------------------------
# Plotting monthly mean averages
#------------------------------
for var in variables:
	print(var)
	ds1 = xr.open_dataset(rpath+var+"_"+case1+"_"+date1+".nc")
	ds2 = xr.open_dataset(rpath+var+"_"+case2+"_"+date2+".nc")
	
        # Get start and end date of period
	date_start = str(ds1.time[0].values).split(" ")[0]
	date_end = str(ds1.time[-1].values).split(" ")[0]

        # Get monthly mean
	ds1m = ds1.groupby("time.month").mean("time")
	ds2m = ds2.groupby("time.month").mean("time")
	
	# Make time array
	months = []
	for i in range(len(ds1m.month.values)):
		datetime_object = datetime.datetime.strptime(str(ds1m.month.values[i]), "%m")
		months.append(datetime_object.strftime("%b"))
	
        # Get spatial average over Arctic
	ds1_arct = computeWeightedMean(ds1m[var].sel(lat=slice(66.5,90)))
	ds2_arct = computeWeightedMean(ds2m[var].sel(lat=slice(66.5,90)))

        # Get spatial average of cases globally
	ds1_glob = computeWeightedMean(ds1m[var])
	ds2_glob = computeWeightedMean(ds2m[var])
	
	# Get spatial average over Svalbard
	ds1_sval = computeWeightedMean(ds1m[var].sel(lon=slice(svalbard[0][0],svalbard[0][1]),
				 lat=slice(svalbard[1][0],svalbard[1][1])))
	ds2_sval = computeWeightedMean(ds2m[var].sel(lon=slice(svalbard[0][0],svalbard[0][1]),
				 lat=slice(svalbard[1][0],svalbard[1][1])))	

	# Get spatial average over Quttinirpaaq
	ds1_qutt = computeWeightedMean(ds1m[var].sel(lon=slice(quttinirpaaq[0][0],quttinirpaaq[0][1]),
				 lat=slice(quttinirpaaq[1][0],quttinirpaaq[1][1])))
	ds2_qutt = computeWeightedMean(ds2m[var].sel(lon=slice(quttinirpaaq[0][0],quttinirpaaq[0][1]),
				 lat=slice(quttinirpaaq[1][0],quttinirpaaq[1][1])))


	fig = plt.figure(1, figsize=[15,6],dpi=300)
	ax = plt.subplot(1,1,1)
	ax.set_title(ds1[var].long_name+" "+case2nm+"-"+case1nm+" "+date_start+"-"+date_end, fontsize=22)

	ax.plot(months, ds1_arct, color="navy", linestyle="--", label=case1nm+" Arctic")
	ax.plot(months, ds1_glob, color="orangered", linestyle="--", label=case1nm+" Global")
	ax.plot(months, ds1_sval, color="darkgreen", linestyle="--", label=case1nm+" Svalbard")
	ax.plot(months, ds1_qutt, color="darkmagenta", linestyle="--", label=case1nm+" Quttinirpaaq")
	ax.plot(months, ds2_arct, color="cornflowerblue", label=case2nm+" Arctic")
	ax.plot(months, ds2_glob, color="orange", label=case2nm+" Global")
	ax.plot(months, ds2_sval, color="mediumseagreen", label=case2nm+" Svalbard")
	ax.plot(months, ds2_qutt, color="magenta", label=case2nm+" Quttinirpaaq")

	ax.set_ylabel(ds1[var].units, fontsize=18)

        # Shrink current axis's height by 10% on the bottom
	box = ax.get_position()
	ax.set_position([box.x0, box.y0 + box.height * 0.1,
                 box.width, box.height * 0.9])

        # Put a legend below current axis
	ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.08), ncol=4)
	
	plt.grid(alpha=0.5)
	plt.savefig(wpath+"monthlymean/"+var+"_"+case1+"_"+case2+".pdf",bbox_inches="tight")

	plt.clf()

