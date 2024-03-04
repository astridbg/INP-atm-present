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
import functions

rpath="/projects/NS9600K/astridbg/data/model/noresm_postprocessed/"
wpath="/projects/NS9600K/astridbg/master/figures/model/diff_height/"

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
# Three-dimensional fields
# for specific level
#------------------------------

level = "750"
#variables = ["AWNI", "FREQI","CLDICE"]
#variables = ["NIMEY","AWNI", "FREQI","CLDICE"]
variables = ["NIMEY"]

for var in variables:
    print(var)
    ds1 = xr.open_dataset(rpath+var+"_"+case1+"_"+date1+".nc")
    ds2 = xr.open_dataset(rpath+var+"_"+case2+"_"+date2+".nc")
        
    # Get start and end date of period
    date_start = str(ds1.time[0].values).split(" ")[0]
    date_end = str(ds1.time[-1].values).split(" ")[0]
    print(ds1["hyam"])
"""
        # Get the time average of cases over the whole period
	ds1m = ds1.mean("time")
	ds2m = ds2.mean("time")

	# Make horizontal averages:
	# - for the Arctic
	ds1_arct = functions.computeWeightedMean(ds1m[var].sel(lat=slice(66.5,90)))
	ds2_arct = functions.computeWeightedMean(ds2m[var].sel(lat=slice(66.5,90)))
	
	fig = plt.figure(1, figsize=[5,7],dpi=300)

	fig.suptitle(ds1[var].long_name+"\n"+date_start+r"$-$"+date_end+", Arctic", fontsize=20)
       	
	if var == "TH":
	   levels = ds1.ilev.values
	else:
	   levels = ds1.lev.values
	print(levels)
 
	plt.plot(ds1_arct, levels, label=case1nm)
	plt.plot(ds2_arct, levels, label=case2nm)
	
	plt.ylabel("hPa")
	plt.xlabel(ds1[var].units)
	plt.legend(loc="upper left")
	plt.grid(alpha=0.5)
	plt.gca().invert_yaxis()
	
	plt.savefig(wpath+var+"_height_"+case1+"_"+case2+".pdf",bbox_inches="tight")
	plt.clf()	
"""
