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
wpath="/projects/NS9600K/astridbg/master/figures/model/reldiff_all/"

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

variables = ["SWCF","LWCF","SWCFS","LWCFS","CLDTOT","CLDHGH","CLDMED","CLDLOW","TGCLDIWP","TGCLDLWP","TREFHT"]
variables = ["NETCFS"]
#------------------------------
# Shaping and plotting fields
#------------------------------
for var in variables:
    print(var)
    ds1 = xr.open_dataset(rpath+var+"_"+case1+"_"+date1+".nc")
    ds2 = xr.open_dataset(rpath+var+"_"+case2+"_"+date2+".nc")
	
    # Get start and end date of period
    date_start = str(ds1.time[0].values).split(" ")[0]
    date_end = str(ds1.time[-1].values).split(" ")[0]

    # Get difference between cases time averaged over the whole period
    ds1m = ds1[var].mean("time")
    ds2m = ds2[var].mean("time")
    diff = ds2m-ds1m
    reldiff = diff/ds1m.where(ds1m>0)*100
    
    lev_extent = round(max(abs(np.nanmin(reldiff.sel(lat=slice(66.5,90)).values)),
                           abs(np.nanmax(reldiff.sel(lat=slice(66.5,90)).values))))
    print(lev_extent)
    lev_extent = 7
    if lev_extent < 0.004:
        lev_extent = 0.004
    levels = np.linspace(-lev_extent,lev_extent,25)

    fig = plt.figure(1, figsize=[9,10],dpi=300)

    fig.suptitle(ds1[var].long_name+" "+case2nm+"-"+case1nm+"\n"+date_start+"-"+date_end, fontsize=26)
	
    # Set the projection to use for plotting
    ax = plt.subplot(1, 1, 1, projection=ccrs.Orthographic(0, 90))

    map = diff.plot.pcolormesh(ax=ax, transform=ccrs.PlateCarree(), 
						cmap='coolwarm',
                                                levels=levels,
						add_colorbar=False)

    ax.coastlines()
	
    cb_ax = fig.add_axes([0.15, 0.07, 0.7, 0.04])

    cbar = plt.colorbar(map, cax=cb_ax, spacing = 'uniform', extend='both', orientation='horizontal', fraction=0.046, pad=0.06)
    cbar.ax.tick_params(labelsize=18)
    cbar.ax.set_xlabel("%", fontsize=23)
    
    if lev_extent >= 4:
        cbar.ax.xaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}')) # No decimal places        
    elif 0.4 <= lev_extent < 4:
        cbar.ax.xaxis.set_major_formatter(StrMethodFormatter('{x:,.1f}')) # One decimal place
    elif 0.04 <= lev_extent < 0.4:
        cbar.ax.xaxis.set_major_formatter(StrMethodFormatter('{x:,.2f}')) # Two decimal places     
    elif 0.004 <= lev_extent < 0.04:
        cbar.ax.xaxis.set_major_formatter(StrMethodFormatter('{x:,.3f}')) # Three decimal places
    
    plt.savefig(wpath+var+"_"+case1+"_"+case2+".pdf",bbox_inches="tight")
	
    plt.clf()
