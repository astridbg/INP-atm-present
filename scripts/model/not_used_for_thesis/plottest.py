import xarray as xr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

rpath="/projects/NS9600K/astridbg/data/model/noresm_postprocessed/"
wpath="/projects/NS9600K/astridbg/INP-atm-present/figures/model/test/"

#case1 = "meyers92_20220210"; case1nm="M92"
case1 = "M92_20240522"; case1nm = "M92"	

#case2 = "andenes21_20220222"; case2nm="A21"
case2 = "A21_20240522"; case2nm = "A21"
date1 = "2007-04-15_2010-03-15"
date2 = "2007-04-15_2010-03-15"

var = "TGCLDLWP"

ds1 = xr.open_dataset(rpath+var+"_"+case1+"_"+date1+".nc")
ds2 = xr.open_dataset(rpath+var+"_"+case2+"_"+date2+".nc")

fig = plt.figure(1, figsize=[10,5])

ax = plt.subplot(1, 1, 1, projection=ccrs.PlateCarree())

(ds2[var].mean("time")-ds1[var].mean("time")).plot(ax=ax, transform=ccrs.PlateCarree(), cmap='coolwarm')
plt.title(var+" "+case2nm+"-"+case1nm, fontsize=18)
#plt.title(var+" time mean "+case, fontsize=18)


ax.coastlines()

plt.savefig(wpath+var+"_"+case1+"_"+case2+"_test.png")
