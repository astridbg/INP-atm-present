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

rpath="/projects/NS9600K/astridbg/data/model/noresm_postprocessed/"
wpath="/projects/NS9600K/astridbg/master/figures/model/spatial/"

# Case ----------------
#case = "def_20210126"; casenm = "CAM6"
#case = "meyers92_20220210"; casenm = "CAM5"
case = "andenes21_20220222"; casenm = "Andenes 2021"
#------------------------------	
date = "2007-04-15_2010-03-15"

#------------------------------
# Two-dimensional fields
#------------------------------

#variables = ["SWCF","LWCF","CLDTOT","CLDHGH","CLDMED","CLDLOW","TGCLDIWP","TGCLDLWP","TREFHT"]
variables = ["LWCFS","SWCFS","FLNS","FLNSC","FSNS","FSNSC"]
variables = ["PBLH"]
#------------------------------
# Shaping and plotting fields
#------------------------------
for var in variables:
	print(var)
	ds = xr.open_dataset(rpath+var+"_"+case+"_"+date+".nc")
	
	# Get start and end date of period
	date_start = str(ds.time[0].values).split(" ")[0]
	date_end = str(ds.time[-1].values).split(" ")[0]

	# Get time average of case over the whole period
	ds_mean = ds[var].mean("time")

	fig = plt.figure(1, figsize=[9,10],dpi=300)

	fig.suptitle(ds[var].long_name+"\n"+casenm+"_"+date_start+r"$-$"+date_end, fontsize=26)
	
	lev_extent = round(max(abs(np.min(ds_mean.values)), abs(np.max(ds_mean.values))),2)
	if lev_extent < 0.004:
	   lev_extent = 0.004
	levels = np.linspace(-lev_extent,lev_extent,25)
	
	# Set the projection to use for plotting
	ax = plt.subplot(1, 1, 1, projection=ccrs.Orthographic(0, 90))

	map = ds_mean.plot.pcolormesh(ax=ax, transform=ccrs.PlateCarree(), 
						cmap='coolwarm',levels=levels,
						add_colorbar=False)

	ax.coastlines()
	
	cb_ax = fig.add_axes([0.15, 0.07, 0.7, 0.04])

	cbar = plt.colorbar(map, cax=cb_ax, spacing = 'uniform', extend='both', orientation='horizontal', fraction=0.046, pad=0.06)
	cbar.ax.tick_params(labelsize=18)
	cbar.ax.set_xlabel(ds[var].units, fontsize=23)

	if lev_extent >= 4:
           cbar.ax.xaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}')) # No decimal places        
	elif 0.4 <= lev_extent < 4:
           cbar.ax.xaxis.set_major_formatter(StrMethodFormatter('{x:,.1f}')) # One decimal place
	elif 0.04 <= lev_extent < 0.4:
           cbar.ax.xaxis.set_major_formatter(StrMethodFormatter('{x:,.2f}')) # Two decimal places     
	elif 0.004 <= lev_extent < 0.04:
           cbar.ax.xaxis.set_major_formatter(StrMethodFormatter('{x:,.3f}')) # Three decimal places

	plt.savefig(wpath+var+"_"+case+".pdf",bbox_inches="tight")
	
	plt.clf()
