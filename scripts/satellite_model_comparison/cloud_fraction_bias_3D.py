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
filenames = glob.glob(calipso_path +"20*/3D_CloudFraction*")
filenames.sort()
calipso_ds = xr.open_mfdataset(filenames, compat='override',coords='all')

model_vars = ["CLD_CAL","CLD_CAL_ICE", "CLD_CAL_LIQ", "CLD_CAL_UN"]
A21 = xr.open_mfdataset([noresm_path+var+'_A21_20240612_2007-04-15_2010-03-15.nc' for var in model_vars])
M92 = xr.open_mfdataset([noresm_path+var+'_M92_20240612_2007-04-15_2010-03-15.nc' for var in model_vars])

#-------------------------------------
# Compute Arctic average
#-------------------------------------

calipso_ds = calipso_ds.sel(time=slice('2007-04-15','2010-03-16'))
CALIOP_CLD = computeWeightedMean_CALIOP(calipso_ds['clcalipso'].sel(latitude=slice(66.5,82)))
CALIOP_ICE = computeWeightedMean_CALIOP(calipso_ds['clcalipso_ice'].sel(latitude=slice(66.5,82)))
CALIOP_LIQ = computeWeightedMean_CALIOP(calipso_ds['clcalipso_liq'].sel(latitude=slice(66.5,82)))
CALIOP_UN = computeWeightedMean_CALIOP(calipso_ds['clcalipso_un'].sel(latitude=slice(66.5,82)))
A21_Aavg = computeWeightedMean(A21.sel(lat=slice(66.5,82)))
M92_Aavg = computeWeightedMean(M92.sel(lat=slice(66.5,82)))

#-------------------------------------
# Plot time averaged data bias
#-------------------------------------
# Makte time average
CALIOP_ICE_tavg = CALIOP_ICE.mean('time').values*100
CALIOP_LIQ_tavg = CALIOP_LIQ.mean('time').values*100
CALIOP_UN_tavg = CALIOP_UN.mean('time').values*100
CALIOP_CLD_tavg = CALIOP_CLD.mean('time').values*100
A21_Aavg_tavg = A21_Aavg.mean('time')
M92_Aavg_tavg = M92_Aavg.mean('time')


plt.figure()
plt.plot(A21_Aavg_tavg['CLD_CAL_ICE'].values-CALIOP_ICE_tavg, CALIOP_ICE.altitude, label='A21', color='tab:orange')
plt.plot(M92_Aavg_tavg['CLD_CAL_ICE'].values-CALIOP_ICE_tavg, CALIOP_ICE.altitude, label='M92', color='tab:blue', linestyle='--')
plt.ylabel('Altitude')
plt.xlabel(r'%')
plt.grid(alpha=0.5)
plt.legend()
plt.savefig(fig_path+'pdf/vertical_ice_bias.pdf', bbox_inches="tight")
plt.savefig(fig_path+'png/vertical_ice_bias.png', bbox_inches="tight")
plt.clf()

plt.plot(A21_Aavg_tavg['CLD_CAL_LIQ'].values-CALIOP_LIQ_tavg, CALIOP_LIQ.altitude, label='A21', color='tab:orange')
plt.plot(M92_Aavg_tavg['CLD_CAL_LIQ'].values-CALIOP_LIQ_tavg, CALIOP_LIQ.altitude, label='M92', color='tab:blue', linestyle='--')
plt.ylabel('Altitude')
plt.xlabel(r'%')
plt.grid(alpha=0.5)
plt.legend()
plt.savefig(fig_path+'pdf/vertical_liq_bias.pdf', bbox_inches="tight")
plt.savefig(fig_path+'png/vertical_liq_bias.png', bbox_inches="tight")
plt.clf()

plt.plot(A21_Aavg_tavg['CLD_CAL'].values-CALIOP_CLD_tavg, CALIOP_CLD.altitude, label='A21', color='tab:orange')
plt.plot(M92_Aavg_tavg['CLD_CAL'].values-CALIOP_CLD_tavg, CALIOP_CLD.altitude, label='M92', color='tab:blue', linestyle='--')
plt.ylabel('Altitude')
plt.xlabel(r'%')
plt.grid(alpha=0.5)
plt.legend()
plt.savefig(fig_path+'pdf/vertical_cld_bias.pdf', bbox_inches="tight")
plt.savefig(fig_path+'png/vertical_cld_bias.png', bbox_inches="tight")
plt.clf()

"""
plt.plot(A21_Aavg_tavg['CLD_CAL_UN'].values-CALIOP_UN_tavg, CALIOP_UN.altitude, label='A21', color='tab:orange')
plt.plot(M92_Aavg_tavg['CLD_CAL_UN'].values-CALIOP_UN_tavg, CALIOP_UN.altitude, label='M92', color='tab:blue', linestyle='--')
plt.ylabel('Altitude')
plt.xlabel(r'%')
plt.grid(alpha=0.5)
plt.savefig(fig_path+'pdf/vertical_un_bias.pdf', bbox_inches="tight")
plt.savefig(fig_path+'png/vertical_un_bias.png', bbox_inches="tight")
plt.clf()
"""