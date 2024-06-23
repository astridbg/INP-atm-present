import xarray as xr
import pandas as pd
import numpy as np

rpath="/projects/NS9600K/astridbg/data/model/noresm_postprocessed/"
wpath="/projects/NS9600K/astridbg/data/model/noresm_postprocessed/"

# Case ------------------------
#case = "def_20210126"; casenm = "CAM6"
#case = "meyers92_20220210"; casenm = "M92"
case = "andenes21_20220222"; casenm = "A21"
#------------------------------	
date = "2007-04-15_2010-03-15"
#------------------------------
# Two-dimensional fields
#------------------------------

ds_ocn = xr.open_dataset(rpath+"OCNFRAC"+"_"+case+"_"+date+".nc")

open_sea = ds_ocn > 0.85
open_sea = open_sea.groupby("time.season").mean("time")
print(open_sea)
open_sea = open_sea >= 0.5
print(open_sea)
quit()
MIZ_found = False
PI_found = False
for n in range(len(ds.season)):
    for i in range(len(ds.lon)):
        print(i)
        MIZ_found = False
        PI_found = False
        for j in range(len(ds.lat.sel(lat=slice(66.5,90)))):
            print(j)
            print(ds["ICEFRAC"].isel(season=n, lon=i, lat=j).values)
            if ds["ICEFRAC"].isel(season=n, lon=i, lat=83+j).values >= 0.15 and MIZ_found==False:
                MIZ_EDGE[n, i] = ds.lat.sel(lat=slice(66.5,90))[j]
                MIZ_found = True
                print("Found MIZ!")
            if ds["ICEFRAC"].isel(season=n, lon=i, lat=83+j).values >= 0.80 and PI_found==False:
                PI_EDGE[n, i] = ds.lat.sel(lat=slice(66.5,90))[j]
                PI_found = True
                print("Found PI!")



ds_MIZ = xr.Dataset(
    data_vars=dict(
        ice_edge_lat=(["season", "lon"], MIZ_EDGE),
    ),
    coords=dict(
        lon=ds.lon,
        season=ds.season,
    ),
    attrs=dict(description="Marginal Ice Zone Edge (Latitude), 15 % Sea Ice Limit"),
)

ds_PI = xr.Dataset(
    data_vars=dict(
        ice_edge_lat=(["season", "lon"], PI_EDGE),
    ),
    coords=dict(
        lon=ds.lon,
        season=ds.season,
    ),
    attrs=dict(description="Pack Ice Edge (Latitude), 80 % Sea Ice Limit"),
)

ds_MIZ.to_netcdf(wpath+"MIZ_EDGE_"+case+"_byseason.nc")
ds_PI.to_netcdf(wpath+"PI_EDGE_"+case+"_byseason.nc")

"""
MIZ_EDGE = np.zeros((len(ds.time), len(ds.lon)))
PI_EDGE = np.zeros((len(ds.time), len(ds.lon)))

MIZ_found = False
PI_found = False
for n in range(len(ds.time)):
    for i in range(len(ds.lon)):
        MIZ_found = False
        PI_found = False
        for j in range(len(ds.lat.sel(lat=slice(66.5,90)))):
            print(j)
            print(ds["ICEFRAC"].isel(time=n, lon=i, lat=j).values)
            if ds["ICEFRAC"].isel(time=n, lon=i, lat=j).values >= 0.15 and MIZ_found==False:
                MIZ_EDGE[n, i] = ds.lat.sel(lat=slice(66.5,90))[j]
                MIZ_found = True
                print("Found MIZ!")
            if ds["ICEFRAC"].isel(time=n, lon=i, lat=j).values >= 0.80 and PI_found==False:
                PI_EDGE[n, i] = ds.lat.sel(lat=slice(66.5,90))[j]
                PI_found = True
                print("Found PI!")

ds_MIZ = xr.Dataset(
    data_vars=dict(
        MIZ_edge=(["time", "lon"], MIZ_EDGE),
    ),

    coords=dict(
        lon=ds.lon,
        time=ds.time,
    ),

    attrs=dict(description="Marginal Ice Zone Edge (Latitude), 15 % Sea Ice Limit"),

)

ds_PI = xr.Dataset(
    data_vars=dict(
        PI_edge=(["time", "lon"], PI_EDGE),
    ),

    coords=dict(
        lon=ds.lon,
        time=ds.time,
    ),

    attrs=dict(description="Pack Ice Edge (Latitude), 80 % Sea Ice Limit"),

)
ds_MIZ.to_netcdf(wpath+"MIZ_EDGE_"+case+"_"+date+".nc")
ds_PI.to_netcdf(wpath+"PI_EDGE_"+case+"_"+date+".nc")
"""