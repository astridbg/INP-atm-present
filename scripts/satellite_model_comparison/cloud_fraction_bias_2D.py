import xarray as xr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob
# Set font style to match latex document----------
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size':20})
# ------------------------------------------------
from functions import *

# Paths to data
calipso_path = "/projects/NS9600K/astridbg/data/observations/CALIPSO-GOCCP/"
noresm_path = "/projects/NS9600K/astridbg/data/model/noresm_postprocessed/"

# Paths to figures
fig_path = "/projects/NS9600K/astridbg/INP-atm-present/figures/satellite_model_comparison/"


#-------------------------------------
# Read data
#-------------------------------------
filenames = glob.glob(calipso_path +"20*/Map*")
filenames.sort()
calipso_ds = xr.open_mfdataset(filenames, compat='override',coords='all')

model_var = 'CLDHGH'
model_vars = [model_var+'_CAL',model_var+'_CAL_ICE',model_var+'_CAL_LIQ', model_var+'_CAL_UN']

case1 = "NorESM2_3_slf_output_20250421"; label1 = "NorESM2.3"
case2 = 'A21_NoContactFreezing_20250117'; label2 = 'A21'

ds1 = xr.open_mfdataset([noresm_path+var+'_'+case1+'_2007-04-15_2010-03-15.nc' for var in model_vars])
#M92 = xr.open_mfdataset([noresm_path+var+'_M92_20240612_2007-04-15_2010-03-15.nc' for var in model_vars])
ds2 = xr.open_mfdataset([noresm_path+var+'_'+case2+'_2007-04-15_2010-03-15.nc' for var in model_vars])

#-------------------------------------
# Compute Arctic average
#-------------------------------------

calipso_ds = calipso_ds.sel(time=slice('2007-04-15','2010-03-16'))
if model_var=='CLDTOT':
    calipso_var = 'cltcalipso'
if model_var=='CLDLOW':
    calipso_var = 'cllcalipso'
if model_var=='CLDMED':
    calipso_var = 'clmcalipso'
if model_var=='CLDHGH':
    calipso_var = 'clhcalipso'

if model_var=='CLDHGH':
    CALIOP_Aavg = computeWeightedMean(calipso_ds[calipso_var].sel(lat=slice(66.5,82))).groupby("time.month").mean("time")
else:
    CALIOP_Aavg = computeWeightedMean_CALIOP(calipso_ds[calipso_var].sel(latitude=slice(66.5,82))).groupby("time.month").mean("time")
CALIOP_ICE_Aavg = computeWeightedMean_CALIOP(calipso_ds[calipso_var+'_ice'].sel(latitude=slice(66.5,82))).groupby("time.month").mean("time")
CALIOP_LIQ_Aavg = computeWeightedMean_CALIOP(calipso_ds[calipso_var+'_liq'].sel(latitude=slice(66.5,82))).groupby("time.month").mean("time")
CALIOP_UN_Aavg = computeWeightedMean_CALIOP(calipso_ds[calipso_var+'_un'].sel(latitude=slice(66.5,82))).groupby("time.month").mean("time")

ds1_Aavg = computeWeightedMean(ds1.sel(lat=slice(66.5,82))).groupby("time.month").mean("time")
ds2_Aavg = computeWeightedMean(ds2.sel(lat=slice(66.5,82))).groupby("time.month").mean("time")

monthsn = np.arange(1,13,1)
months = ['Jan', 'Feb', 'Mar', 'Apr', 'May','Jun','Jul','Aug','Sep','Oct','Nov', 'Dec']

#-------------------------------------
# Plot CLOUD COVER BIAS
#-------------------------------------


plt.figure(figsize=[13,4])
plt.hlines(0, 1,12,linewidth=3.5,color='black')


plt.fill_between(CALIOP_Aavg.month[:6], (ds1_Aavg[model_var+'_CAL'].values*0.01-CALIOP_Aavg.values)[:6], alpha=0.3,label=label1, color='tab:blue')
plt.fill_between(CALIOP_Aavg.month[:6], (ds2_Aavg[model_var+'_CAL'].values*0.01-CALIOP_Aavg.values)[:6], alpha=0.3,label=label2, color='tab:orange')

plt.fill_between(CALIOP_Aavg.month[5:8], (ds2_Aavg[model_var+'_CAL'].values*0.01-CALIOP_Aavg.values)[5:8], alpha=0.3, color='tab:orange')
plt.fill_between(CALIOP_Aavg.month[5:8], (ds1_Aavg[model_var+'_CAL'].values*0.01-CALIOP_Aavg.values)[5:8], alpha=0.3,color='tab:blue')

plt.fill_between(CALIOP_Aavg.month[7:], (ds1_Aavg[model_var+'_CAL'].values*0.01-CALIOP_Aavg.values)[7:], alpha=0.3, color='tab:blue')
plt.fill_between(CALIOP_Aavg.month[7:], (ds2_Aavg[model_var+'_CAL'].values*0.01-CALIOP_Aavg.values)[7:], alpha=0.3, color='tab:orange')

plt.fill_between([1,2],0,0, alpha=0.3,label="Total cloud",color='black')
plt.plot(1,0,linestyle='--',color='grey',label='Ice Cloud')
plt.plot(1,0,label='Liquid Cloud',color='grey')

plt.plot(CALIOP_Aavg.month, ds2_Aavg[model_var+'_CAL_LIQ'].values*0.01-CALIOP_LIQ_Aavg.values, color='tab:orange')
plt.plot(CALIOP_Aavg.month, ds1_Aavg[model_var+'_CAL_LIQ'].values*0.01-CALIOP_LIQ_Aavg.values, color='tab:blue')
plt.plot(CALIOP_Aavg.month, ds2_Aavg[model_var+'_CAL_ICE'].values*0.01-CALIOP_ICE_Aavg.values, color='tab:orange', linestyle='--')
plt.plot(CALIOP_Aavg.month, ds1_Aavg[model_var+'_CAL_ICE'].values*0.01-CALIOP_ICE_Aavg.values, color='tab:blue', linestyle='--')
#plt.plot(CALIOP_Aavg.month, ds2_Aavg['CLDTOT_CAL_UN'].values*0.01-CALIOP_UN_Aavg[:,0].values, color='tab:orange', linestyle='-.')
#plt.plot(CALIOP_Aavg.month, ds1_Aavg['CLDTOT_CAL_UN'].values*0.01-CALIOP_UN_Aavg[:,0].values, color='tab:blue', linestyle='-.')

plt.xticks(monthsn,months,rotation=45, ha='right')
plt.ylabel('Bias (Model-GOCCP)')
plt.grid(alpha=0.5)

# Shrink current axis by 20%
box = plt.gca().get_position()
plt.gca().set_position([box.x0, box.y0, box.width * 0.8, box.height])

# Put a legend to the right of the current axis
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

#plt.savefig(fig_path+'pdf/'+model_var+'cc_phase_bias.pdf', bbox_inches="tight")
plt.savefig(fig_path+'png/'+model_var+'_phase_bias_'+case1+'_'+case2+'.png', bbox_inches="tight")
plt.clf()

