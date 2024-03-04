import xarray as xr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

rpath="/projects/NS9600K/astridbg/data/model/noresm_postprocessed/"
wpath="/projects/NS9600K/astridbg/master/figures/model/test/"

#case = "def_20210126"
#case = "meyers92_20220210" 				
case = "andenes21_20220222"
timelength = "2007-04-15_2010-03-15"

var = "FLDS"

ds = xr.open_dataset(rpath+var+"_"+case+"_"+timelength+".nc")

timepoint = 3
date = str(ds.time[timepoint].values).split(" ")[0]
#ds = ds.sel(lev="750", method="nearest").mean("time")

fig = plt.figure(1, figsize=[10,5])

ax = plt.subplot(1, 1, 1, projection=ccrs.PlateCarree())

ds[var].isel(time=timepoint).plot(ax=ax, transform=ccrs.PlateCarree(), cmap='viridis')
#ds[var].plot(ax=ax, transform=ccrs.PlateCarree(), cmap='viridis')

plt.title(var+" "+date+" "+case, fontsize=18)
#plt.title(var+" time mean "+case, fontsize=18)


ax.coastlines()

plt.savefig(wpath+var+"_"+case+"_"+date+"_test.png")
