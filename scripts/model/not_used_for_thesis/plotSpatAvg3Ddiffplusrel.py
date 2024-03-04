import xarray as xr
import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter
# Set font style to match latex document----------
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size':20})
# ------------------------------------------------
import cartopy.crs as ccrs
import datetime
from functions import *

rpath="/projects/NS9600K/astridbg/data/model/noresm_postprocessed/"
wpath="/projects/NS9600K/astridbg/master/figures/model/spatavg_reldiff/"

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

level = "850"
variables = ["NIMEY","AWNI","AWNICC","CLDICE"]
variables = ["RELHUM","Q"]

#------------------------------
# Areas to average over
#------------------------------

svalbard = [[9,28],[75,81]] # Svalbard
quttinirpaaq = [[-120+360,-60+360],[78,86]] # Quttinirpaaq
greenland = [[-57+360,-22+360],[70,80]] # Greenland
npole = [[0,360],[85,90]] # North Pole
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
	ax.set_title(ds1[var].long_name+" "+case2nm+"-"+case1nm+" "+date_start+r"$-$"+date_end, fontsize=22)
	
	ax.plot(datetimeindex, ds2_arct-ds1_arct, color="cornflowerblue", label="Arctic")
	ax.plot(datetimeindex, ds2_glob-ds1_glob, color="orange", label="Global")
	ax.plot(datetimeindex, ds2_sval-ds1_sval, color="mediumseagreen", label="Svalbard")
	ax.plot(datetimeindex, ds2_qutt-ds1_qutt, color="magenta", label="Quttinirpaaq")
	
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
    
    # Select level
    ds1 = ds1.sel(lev=level, method="nearest")
    ds2 = ds2.sel(lev=level, method="nearest")
    lev_name = str(int(ds1.lev.values))

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

    # Get spatial average over Greenland
    ds1_gren = computeWeightedMean(ds1m[var].sel(lon=slice(greenland[0][0], greenland[0][1]),
                               lat=slice(greenland[1][0],greenland[1][1])))
    ds2_gren = computeWeightedMean(ds2m[var].sel(lon=slice(greenland[0][0],greenland[0][1]),
                                 lat=slice(greenland[1][0],greenland[1][1])))
        
    # Get spatial average over North Pole
    ds1_npol = computeWeightedMean(ds1m[var].sel(lon=slice(npole[0][0],npole[0][1]),
                                 lat=slice(npole[1][0],npole[1][1])))
    ds2_npol = computeWeightedMean(ds2m[var].sel(lon=slice(npole[0][0],npole[0][1]),
                                 lat=slice(npole[1][0],npole[1][1])))
	
    fig,axs = plt.subplots(ncols=2,nrows=1, gridspec_kw={'width_ratios': [3, 1]}, figsize=[13,4],dpi=300,constrained_layout=True)
    ax = axs[0]
    ax2 = axs[1]
    fig.suptitle(ds1[var].long_name+", "+case2nm+"-"+case1nm+", "+lev_name+" hPa", fontsize=26)

    ax.plot(months, ds2_arct-ds1_arct, color="cornflowerblue", label="Arctic")
    ax.plot(months, ds2_sval-ds1_sval, color="mediumseagreen", label="Svalbard")
    ax.plot(months, ds2_qutt-ds1_qutt, color="magenta", label="Quttinirpaaq")
    ax.plot(months, ds2_gren-ds1_gren, color="crimson", label="Greenland")
    ax.plot(months, ds2_npol-ds1_npol, color="olive", label="North Pole")
    ax.set_ylabel(ds1[var].units)
    ax.tick_params(axis="x",labelsize=20)
    ax.grid(alpha=0.5)

    # Get average relative change
    
    rel_arct = ((ds2_arct-ds1_arct)/ds1_arct.where(ds1_arct!=0)).fillna(0)*100*np.sign(ds1_arct)
    rel_sval = ((ds2_sval-ds1_sval)/ds1_sval.where(ds1_sval!=0)).fillna(0)*100*np.sign(ds1_sval)
    rel_qutt = ((ds2_qutt-ds1_qutt)/ds1_qutt.where(ds1_qutt!=0)).fillna(0)*100*np.sign(ds1_qutt)
    rel_gren = ((ds2_gren-ds1_gren)/ds1_gren.where(ds1_gren!=0)).fillna(0)*100*np.sign(ds1_gren)
    rel_npol = ((ds2_npol-ds1_npol)/ds1_npol.where(ds1_npol!=0)).fillna(0)*100*np.sign(ds1_npol)
    
    rel_all = pd.DataFrame({"Arctic":rel_arct,"Svalbard":rel_sval,"Quttinirpaaq":rel_qutt,"Greenland":rel_gren,"North Pole":rel_npol})
    bplot=ax2.boxplot(rel_all,patch_artist=True,medianprops={"color":"black"})
    colors=["cornflowerblue","mediumseagreen","magenta","crimson","olive"]
    for patch, color in zip(bplot['boxes'], colors):
        patch.set_facecolor(color)
    ax2.set_ylabel("%")
    ax2.set_xticklabels([])
    
    # Shrink current axis's height by 15% on the bottom
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.15,
                     box.width, box.height * 0.85])
    box = ax2.get_position()
    ax2.set_position([box.x0, box.y0 + box.height * 0.15,
                     box.width, box.height * 0.85])

    # Put a legend below current axis
    ax.legend(loc='upper center', bbox_to_anchor=(0.75, -0.12), ncol=5, columnspacing=0.5, handlelength=1,handletextpad=0.4)
	
    plt.grid(alpha=0.5)
    plt.savefig(wpath+"monthlymean/"+var+"_avg_"+lev_name+"_"+case1+"_"+case2+".pdf",bbox_inches="tight")

    plt.clf()
