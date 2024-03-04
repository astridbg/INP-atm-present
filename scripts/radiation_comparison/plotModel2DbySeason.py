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

rpath="/home/astridbg/Documents/data/noresm_postprocessed/"
wpath="/home/astridbg/Documents/master/figures/radiation_comparison/"

# Default cases----------------
#case1 = "def_20210126"; case1nm = "CAM6"
case1 = "meyers92_20220210"; case1nm = "M92"
# Modified cases---------------
#case2 = "meyers92_20220210"; case2nm = "CAM5"
case2 = "andenes21_20220222"; case2nm = "A21"
#------------------------------	
date1 = "2007-04-15_2010-03-15"
date2 = "2007-04-15_2010-03-15"

#------------------------------
# Two-dimensional fields
#------------------------------

variables = ["SWCFS","LWCFS"]

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

    # Group cases by season and mean over the period by season
    ds1_seas = ds1.groupby("time.season").mean("time")  
    ds2_seas = ds2.groupby("time.season").mean("time")


    diff = ds2_seas[var]-ds1_seas[var]

    lev_extent = round(max(abs(np.min(diff.sel(lat=slice(66.5,90)).values)), 
                            abs(np.max(diff.sel(lat=slice(66.5,90)).values))),2)
    print(lev_extent)
    if lev_extent < 0.004:
        lev_extent = 0.004
    lev_extent = 18
    levels = np.linspace(-lev_extent,lev_extent,25)


    fig = plt.figure(1, figsize=[9,10],dpi=300)
    title = ds1[var].long_name+"\n"+case2nm+"-"+case1nm
    fig.suptitle(title, fontsize=22)
	
    # Set the projection to use for plotting
    ax1 = plt.subplot(2, 2, 1, projection=ccrs.Orthographic(0, 90))
    ax2 = plt.subplot(2, 2, 2, projection=ccrs.Orthographic(0, 90))
    ax3 = plt.subplot(2, 2, 3, projection=ccrs.Orthographic(0, 90))
    ax4 = plt.subplot(2, 2, 4, projection=ccrs.Orthographic(0, 90))
    #plt.subplots_adjust(top=0.85)

    for ax,season in zip([ax1, ax2, ax3, ax4], ["DJF", "MAM","JJA","SON"]):
    	
        functions.polarCentral_set_latlim([66.5,90], ax)
        diff = ds2_seas[var].sel(season=season) - ds1_seas[var].sel(season=season)
        map = diff.plot.pcolormesh(ax=ax, transform=ccrs.PlateCarree(), 
                                           cmap='coolwarm', levels=levels,
                                           add_colorbar=False)
        ax.set_title(season, fontsize=22)
        ax.coastlines()

	
    cb_ax = fig.add_axes([0.15, 0.07, 0.7, 0.04])

    cbar = plt.colorbar(map, cax=cb_ax, spacing = 'uniform', extend='both', orientation='horizontal', fraction=0.046, pad=0.06)
    cbar.ax.tick_params(labelsize=18)
    cbar.ax.set_xlabel(ds1[var].units, fontsize=18)
    
    # Customize number of decimal points for each colorbar extent
    if lev_extent >= 4:
        cbar.ax.xaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}')) # No decimal places        
    elif 0.4 <= lev_extent < 4:
        cbar.ax.xaxis.set_major_formatter(StrMethodFormatter('{x:,.1f}')) # One decimal place
    elif 0.04 <= lev_extent < 0.4:
        cbar.ax.xaxis.set_major_formatter(StrMethodFormatter('{x:,.2f}')) # Two decimal places     
    elif 0.004 <= lev_extent < 0.04:
        cbar.ax.xaxis.set_major_formatter(StrMethodFormatter('{x:,.3f}')) # Three decimal places
    
    plt.savefig(wpath+var+"_"+case1+"_"+case2+"_eqbar.pdf", bbox_inches='tight')
    plt.savefig(wpath+var+"_"+case1+"_"+case2+"_eqbar.png", bbox_inches='tight')

    plt.clf()