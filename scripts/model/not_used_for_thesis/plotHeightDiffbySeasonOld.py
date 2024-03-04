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
wpath="/projects/NS9600K/astridbg/master/figures/model/diff_byseason/"

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
# Areas to average over
#------------------------------

# area = [[min lon, max lon],[min lat, max lat]]
norsea = [[0,72],[67,76]] # Norwegian Sea
bfortsea = [[-155+360,-103+360],[74,86]] # Beaufort Sea
greenland = [[-57+360,-22+360],[70,80]] # Greenland


#------------------------------
# Three-dimensional fields
# for specific level
#------------------------------

#variables = ["AWNI", "FREQI","CLDICE"]
#variables = ["NIMEY","AWNI", "FREQI","CLDICE"]
variables = ["TH","CLOUD"]

for var in variables:
    print(var)
    ds1 = xr.open_dataset(rpath+var+"_"+case1+"_"+date1+".nc")
    ds2 = xr.open_dataset(rpath+var+"_"+case2+"_"+date2+".nc")
        
    # Get start and end date of period
    date_start = str(ds1.time[0].values).split(" ")[0]
    date_end = str(ds1.time[-1].values).split(" ")[0]

    # Group cases by season and mean over the period by season
    ds1m = ds1.groupby("time.season").mean("time")
    ds2m = ds2.groupby("time.season").mean("time")
	
    # Make horizontal averages:
    # - for the Arctic
    ds1_arct = functions.computeWeightedMean(ds1m[var].sel(lat=slice(66.5,90)))
    ds2_arct = functions.computeWeightedMean(ds2m[var].sel(lat=slice(66.5,90)))
    diff_arct = ds2_arct - ds1_arct
    # - for Norwegian Sea
    ds1_nsea = functions.computeWeightedMean(ds1m[var].sel(lon=slice(norsea[0][0],norsea[0][1]),
                                 lat=slice(norsea[1][0],norsea[1][1])))
    ds2_nsea = functions.computeWeightedMean(ds2m[var].sel(lon=slice(norsea[0][0],norsea[0][1]),
                                 lat=slice(norsea[1][0],norsea[1][1])))
    diff_nsea = ds2_nsea - ds1_nsea
    # - for Greenland
    ds1_green = functions.computeWeightedMean(ds1m[var].sel(lon=slice(greenland[0][0],greenland[0][1]),
                                 lat=slice(greenland[1][0],greenland[1][1])))
    ds2_green = functions.computeWeightedMean(ds2m[var].sel(lon=slice(greenland[0][0],greenland[0][1]),
                                 lat=slice(greenland[1][0],greenland[1][1])))
    diff_green = ds2_green - ds1_green
    # - for Beaufort Sea
    ds1_bsea = functions.computeWeightedMean(ds1m[var].sel(lon=slice(bfortsea[0][0],bfortsea[0][1]),
                                 lat=slice(bfortsea[1][0],bfortsea[1][1])))
    ds2_bsea = functions.computeWeightedMean(ds2m[var].sel(lon=slice(bfortsea[0][0],bfortsea[0][1]),
                                 lat=slice(bfortsea[1][0],bfortsea[1][1])))
    diff_bsea = ds2_bsea - ds1_bsea

	
    fig,axs = plt.subplots(nrows=2, ncols=2, figsize=[9,10],dpi=300, constrained_layout=True)
        
    # Set the projection to use for plotting
    ax1 = axs[0,0]
    ax2 = axs[0,1]
    ax3 = axs[1,0]
    ax4 = axs[1,1]

    fig.suptitle(ds1[var].long_name+"\n"+case2nm+"-"+case1nm+"\n"+date_start+r"$-$"+date_end, fontsize=20)

    if var == "TH":       	
        levels = ds1.ilev.values
        print(ds1.ilev)
    else:
        levels = ds1.lev.values
        
    for ax,season in zip([ax1, ax2, ax3, ax4], ["DJF", "MAM","JJA","SON"]):

        ax.plot(diff_arct.sel(season=season), levels, label="Arctic")
        ax.plot(diff_nsea.sel(season=season), levels, label="Norwegian Sea")
        ax.plot(diff_green.sel(season=season), levels, label="Greenland")
        ax.plot(diff_bsea.sel(season=season), levels, label="Beaufort Sea")
            
        ax.set_title(season)
        ax.set_ylabel("hPa")
        ax.set_xlabel(ds1[var].units)
        ax.legend(loc="upper left")
        ax.grid(alpha=0.5)
        ax.invert_yaxis()
    
    plt.savefig(wpath+var+"_heightdiff_"+case1+"_"+case2+".pdf",bbox_inches="tight")
    plt.clf()	

