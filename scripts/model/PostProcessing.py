import xarray as xr
import pandas as pd
import numpy as np
import glob
import functions

rpath="/projects/NS9600K/astridbg/data/model/noresm_rawdata/cases/"
wpath="/projects/NS9600K/astridbg/data/model/noresm_postprocessed/"

#case = "def_20210126"
#case = "meyers92_20220210"
#case = "andenes21_20220222"
#case = "andenes21_20240322"
#case = "andenes21_20240322_biggsoff"
#case = "A21_20240522"
#case = "M92_20240522"
#case = "M92_4K_20240607"
#case = "A21_20240612"
#case = "M92_20240612"
#case = "M92_20241118"
#case = "A21_20241118"
#case = "Sze23_20250115"
#case = "A21_NoContactFreezing_20250117"
#case = "Sze23_NoContactFreezing_20250121"
#case = "NorESM2_3_slf_output_20250421"
case = "NorESM2_3_slf_output_RaFSIP_20250421"


casefolder="NF2000climo_f19_tn14_"+case
casefolder="N1850_f19_tn14_noresm2_3_slf_output_20250529"

all_files = glob.glob(rpath+casefolder+"/atm/hist/"+casefolder+".cam.h0.193*")
all_files.sort()
print("Files found")

ds = xr.open_mfdataset(all_files)
print("Dataset created")
#print([ds.data_vars.keys()])
print(ds["CLDTOT_ISCCP"])
print(ds["FISCCP1_COSP"])
quit()

#-----------------------------
# Postprocessing of model data
#-----------------------------

# Fix timestamp of model data
ds = functions.fix_cam_time(ds)

# Remove spinup months of data set
ds = ds.isel(time=slice(3,len(ds.time)))

#-----------------------------

print("Postprocessing completed")

#--------------------------------------------------------------------
# Store relevant variables intermediately to save time when plotting,
# change to desired units and create combined variables 
#--------------------------------------------------------------------
date = "2007-04-15_2010-03-15"
#date = "2007-04-15_2008-03-15"
#date = "2007-04-15_2007-12-15"
#date = "2007-04-15"


# For calipso-comparison
#variables = ["CLD_CAL","CLD_CAL_ICE", "CLD_CAL_LIQ", "CLD_CAL_NOTCS", "CLD_CAL_UN",
#             "CLDHGH_CAL", "CLDHGH_CAL_ICE", "CLDHGH_CAL_LIQ", "CLDHGH_CAL_UN",
#             "CLDLOW_CAL", "CLDLOW_CAL_ICE", "CLDLOW_CAL_LIQ", "CLDLOW_CAL_UN",
#             "CLDMED_CAL", "CLDMED_CAL_ICE", "CLDMED_CAL_LIQ", "CLDMED_CAL_UN",
#             "CLDTOT_CAL", "CLDTOT_CAL_ICE", "CLDTOT_CAL_LIQ", "CLDTOT_CAL_UN",
#             "CLDHGH_CAL", "CLDTOT_CALCS"]

variables = ["OCNFRAC","FREQI","NUMICE","LWCFS","SWCFS","NETCFS","CLOUD","CLDTOT", "TGCLDTWP","TREFHT", "PREC", "CLDLWEMS",
             "FSNT","FSNTC","FLNT", "FLNTC", "FSNS", "FLNS", 
             "CLDTOT_CAL", "CLDTOT_CAL_ICE", "CLDTOT_CAL_LIQ", "CLDTOT_CAL_UN",
             "CLDLOW_CAL", "CLDLOW_CAL_ICE", "CLDLOW_CAL_LIQ", "CLDLOW_CAL_UN",
             "CLDMED_CAL", "CLDMED_CAL_ICE", "CLDMED_CAL_LIQ", "CLDMED_CAL_UN",
             "CLDHGH_CAL", "CLDHGH_CAL_ICE", "CLDHGH_CAL_LIQ", "CLDHGH_CAL_UN",
             "NSACWIO",
             "TGCLDIWP","TGCLDLWP"] # These last two must be placed last, before TGCLTWP and CLDLEWMS


for var in variables:
    print("Started writing variable:")

    # Make combined data variables
    if var == "TGCLDTWP":
        ds = ds.assign(TGCLDTWP=ds["TGCLDIWP"]+ds["TGCLDLWP"])
        ds[var].values = ds[var].values*1e+3 # Change unit to grams per squared meter
        ds[var].attrs["units"] = "g/m$^2$"
        ds[var].attrs["long_name"] = "Total Grid-box Cloud Total Water Path"

    if var == "LWCFS":
        ds = ds.assign(LWCFS=ds["FLNSC"]-ds["FLNS"])
        ds[var].attrs["units"] = "W/m$^2$"
        ds[var].attrs["long_name"] = "Longwave cloud radiative effect at surface"

    if var == "SWCFS":
        ds = ds.assign(SWCFS=ds["FSNS"]-ds["FSNSC"])
        ds[var].attrs["units"] = "W/m$^2$"
        ds[var].attrs["long_name"] = "Shortwave cloud radiative effect at surface"
	
    if var == "AWNINONIMEY":
        ds = ds.assign(AWNINONIMEY=ds["AWNI"]-ds["NIMEY"])
        ds[var].values = ds[var].values*1e-3 # Change unit to number per litre
        ds[var].attrs["units"] = "1/L"
        ds[var].attrs["long_name"] = "Average cloud ice number concentration minus Meyer's contribution"
    
    if var == "AWNICC":
        ds = ds.assign(AWNICC=ds["AWNI"]/ds["FREQI"].where(ds["FREQI"]>0))
        ds[var] = ds[var].fillna(0)
        ds[var].values = ds[var].values*1e-3 # Change unit to number per litre
        ds[var].attrs["units"] = "1/L"
        ds[var].attrs["long_name"] = "Average cloud ice number concentration in cold clouds"

    if var == "AWNICLD":
        ds = ds.assign(AWNICLD=ds["AWNI"]/ds["CLOUD"].where(ds["CLOUD"]>0))
        ds[var] = ds[var].fillna(0)
        ds[var].values = ds[var].values*1e-3 # Change unit to number per litre
        ds[var].attrs["units"] = "1/L"
        ds[var].attrs["long_name"] = "Average cloud ice number concentration in clouds"
    
    if var == "NIMEYCC":
        ds = ds.assign(NIMEYCC=ds["NIMEY"]/ds["FREQI"].where(ds["FREQI"]>0))
        ds[var] = ds[var].fillna(0)
        ds[var].values = ds[var].values*1e-3 # Change unit to number per litre
        ds[var].attrs["units"] = "1/L"
        ds[var].attrs["long_name"] = "Activated Ice Number Concentation due to Meyers' parameterisation in cold clouds"

    if var == "NIMEYCLD":
        ds = ds.assign(NIMEYCLD=ds["NIMEY"]/ds["CLOUD"].where(ds["CLOUD"]>0))
        ds[var] = ds[var].fillna(0)
        ds[var].values = ds[var].values*1e-3 # Change unit to number per litre
        ds[var].attrs["units"] = "1/L"
        ds[var].attrs["long_name"] = "Activated Ice Number Concentation due to Meyers' parameterisation in clouds"
        
    if var == "NETCFS":
        ds = ds.assign(NETCFS=ds["FSNS"]-ds["FSNSC"]-ds["FLNS"]-(-ds["FLNSC"]))
        ds[var].attrs["units"] = "W/m$^2$"
        ds[var].attrs["long_name"] = "Net cloud radiative effect at surface"
    
    if var == "CLDLWEM":
        ds = ds.assign(CLDLWEM=1-np.exp(-0.13*ds["TGCLDLWP"]*1e+3))
        ds[var].attrs["units"] = "Emissivity"
        ds[var].attrs["long_name"] = "Cloud Longwave Emissivity"

    if var == "CLDLWEMS":
        ds = ds.assign(CLDLWEMS=1-np.exp(-0.158*ds["TGCLDLWP"]*1e+3))
        ds[var].attrs["units"] = "Emissivity"
        ds[var].attrs["long_name"] = "Cloud Longwave Emissivity"
    
    if var == "PREC":
        print(ds["PRECL"].attrs["units"])
        print(ds["PRECC"].attrs["units"])
        print("NB! Changing units to mm/month!")
        ds = ds.assign(PREC=(ds["PRECL"]+ds["PRECL"])*1e+3*60*60*24*30) # Sum precipitation rates
        ds[var].attrs["units"] = "mm/month" # Change unit to millimeter per month
        ds[var].attrs["long_name"] = "Total Precipitation"
	
    # Change to desired units
    if var == "NIMEY" or var == "AWNI":
        ds[var].values = ds[var].values*1e-3 # Change unit to number per litre and name
        ds[var].attrs["units"] = "1/L"

    if var == "T" or var == "TREFHT" or var=="TH":
        ds[var].values = ds[var].values - 273.15 # Change unit to degrees Celsius
        ds[var].attrs["units"] = r"$^{\circ}$C"

    if var == "CLDICE" or var == "CLDLIQ" or var == "Q" or var == "ICIMR" or var == "ICWMR": 
        ds[var].values = ds[var].values*1e+3 # Change unit to grams per kilograms
        ds[var].attrs["units"] = "g/kg"
	
    if var == "TGCLDIWP" or var == "TGCLDLWP":
        ds[var].values = ds[var].values*1e+3 # Change unit to grams per squared meter
        ds[var].attrs["units"] = "g/m$^2$"

    if var == "IWC":
        ds[var].values = ds[var].values*1e+3 # Change unit to grams per cubic meter
        ds[var].attrs["units"] = "g/m$^3$"
    
    if var == "PRECC" or var == "PRECL":
        ds[var].values = ds[var].values*1e+3*60*60*24*30 # Change unit to millimeter per month
        ds[var].attrs["units"] = "mm/month"

    # Change to more meaningful name

    if var == "TREFHT":
        ds[var].attrs["long_name"] = "Surface (2m) Temperature"
    if var == "NIMEY":
        ds[var].attrs["long_name"] = "Activated Ice Number Concentation due to Meyers' parameterisation"
    if var == "AWNI":
        ds[var].attrs["long_name"] = "Average cloud ice number concentration"
    if var == "ICIMR":
        ds[var].attrs["long_name"] = "In-cloud ice mixing ratio"
    if var == "ICWMR":
        ds[var].attrs["long_name"] = "In-cloud liquid water mixing ratio"
    
    print(ds[var].attrs["long_name"])
    print("Units: ", ds[var].attrs["units"])

    ds[var].to_netcdf(wpath+var+"_"+case+"_"+date+".nc")
