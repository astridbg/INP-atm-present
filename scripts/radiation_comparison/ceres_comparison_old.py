import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd
from functions import *
import glob

# Paths to data
ceres_path = "/home/astridbg/Documents/data/CERES_EBAF-TOA_Ed4.1_Subset_200901-201012.nc"
noresm_path = "/home/astridbg/Documents/data/noresm_postprocessed/"

# Paths to figures
fig_path = "/home/astridbg/Documents/master/figures/radiation_comparison/ceres_comparison/"

# Model variables to consider
model_vars = ['FSUTOA', 'FLUT', 'FLUTC', 'FLNT', 'FLNTC', 'FSNT', 'FSNTC']
ceres_data = xr.open_dataset(ceres_path)
A21 = xr.open_mfdataset([noresm_path+var+'_andenes21_20220222_2007-04-15_2010-03-15.nc' for var in model_vars])
M92 = xr.open_mfdataset([noresm_path+var+'_meyers92_20220210_2007-04-15_2010-03-15.nc' for var in model_vars])

ceres_Aavg = computeWeightedMean(ceres_data.sel(lat=slice(66.5,90)))
A21_Aavg = computeWeightedMean(A21.sel(lat=slice(66.5,90)))
M92_Aavg = computeWeightedMean(M92.sel(lat=slice(66.5,90)))



plt.figure()
ceres_Aavg['toa_sw_all_mon'].plot(label='sw out all')
ceres_Aavg['toa_lw_all_mon'].plot(label='lw out all')
ceres_Aavg['toa_net_all_mon'].plot(label='net in all')
ceres_Aavg['toa_sw_clr_c_mon'].plot(label='sw out clr')
ceres_Aavg['toa_lw_clr_c_mon'].plot(label='lw out clr')
ceres_Aavg['toa_net_clr_c_mon'].plot(label='net in clr')
ceres_Aavg['solar_mon'].plot(label='sw in')
plt.legend()
plt.savefig(fig_path+'ceres.png')


# Outgoing longwave allsky
plt.figure()
plt.suptitle("Outgoing longwave all-sky")
ceres_Aavg['toa_lw_all_mon'].plot(label='CERES')
plt.plot(A21.time, A21_Aavg['FLUT'], label='A21')
plt.plot(M92.time, M92_Aavg['FLUT'], label='M92')
plt.legend()
plt.savefig(fig_path+'lw_all.png')

#Net and out
plt.figure()
plt.suptitle("Outgoing longwave and net longwave all-sky")
ceres_Aavg['toa_lw_all_mon'].plot(label='CERES out')
plt.plot(A21.time, A21_Aavg['FLNT'], label='A21 net')
plt.plot(A21.time, A21_Aavg['FLUT'], label='A21 out')
plt.legend()
plt.savefig(fig_path+'lw_net_out_all.png')


# Outgoing longwave clearsky
plt.figure()
plt.suptitle("Outgoing longwave clear-sky")
ceres_Aavg['toa_lw_clr_c_mon'].plot(label='CERES')
plt.plot(A21.time, A21_Aavg['FLUTC'], label='A21')
plt.plot(M92.time, M92_Aavg['FLUTC'], label='M92')
plt.legend()
plt.savefig(fig_path+'lw_clr.png')

# Outgoing longwave clouds
plt.figure()
plt.suptitle("Outgoing longwave clouds (all-sky - clear-sky)")
(ceres_Aavg['toa_lw_all_mon'] - ceres_Aavg['toa_lw_clr_c_mon']).plot(label='CERES')
plt.plot(A21.time, (A21_Aavg['FLUT']-A21_Aavg['FLUTC']), label='A21')
plt.plot(M92.time, (M92_Aavg['FLUT'] - M92_Aavg['FLUTC']), label='M92')
plt.legend()
plt.savefig(fig_path+'lw_clouds.png')

# Outgoing longwave clouds
plt.figure()
plt.suptitle("Net longwave clouds (all-sky - clear-sky)")
(ceres_Aavg['toa_lw_all_mon'] - ceres_Aavg['toa_lw_clr_c_mon']).plot(label='CERES out')
plt.plot(A21.time, (A21_Aavg['FLNT']-A21_Aavg['FLNTC']), label='A21')
plt.plot(M92.time, (M92_Aavg['FLNT'] - M92_Aavg['FLNTC']), label='M92')
plt.legend()
plt.savefig(fig_path+'lw_net_clouds.png')


# Outgoing shortwave
plt.figure()
plt.suptitle("Outgoing shortwave all-sky")
ceres_Aavg['toa_sw_all_mon'].plot(label='CERES')
plt.plot(A21.time, A21_Aavg['FSUTOA'], label='A21')
plt.plot(M92.time, M92_Aavg['FSUTOA'], label='M92')
plt.legend()
plt.savefig(fig_path+'sw_all.png')


# Net solar flux all-sky
plt.figure()
plt.suptitle("Net shortwave all-sky")
(ceres_Aavg['solar_mon']-ceres_Aavg['toa_sw_all_mon']).plot(label='CERES in-out')
plt.plot(A21.time, A21_Aavg['FSNT'], label='A21')
plt.plot(M92.time, M92_Aavg['FSNT'], label='M92')
plt.legend()
plt.savefig(fig_path+'sw_net_all.png')


# Net solar flux clear-sky
plt.figure()
plt.suptitle("Net shortwave clear-sky")
(ceres_Aavg['solar_mon']-ceres_Aavg['toa_sw_clr_c_mon']).plot(label='CERES in-out')
plt.plot(A21.time, A21_Aavg['FSNTC'], label='A21')
plt.plot(M92.time, M92_Aavg['FSNTC'], label='M92')
plt.legend()
plt.savefig(fig_path+'sw_net_clr.png')

# Net solar flux clouds
plt.figure()
plt.suptitle("Net shortwave clouds")
(ceres_Aavg['toa_sw_clr_c_mon'] - ceres_Aavg['toa_sw_all_mon']).plot(label='CERES in-out')
plt.plot(A21.time, (A21_Aavg['FSNT']-A21_Aavg['FSNTC']), label='A21')
plt.plot(M92.time, (M92_Aavg['FSNT'] - M92_Aavg['FSNTC']), label='M92')
plt.legend()
plt.savefig(fig_path+'sw_net_clouds.png')

# Net all-sky
plt.figure()
plt.suptitle("Net flux all-sky")
ceres_Aavg['toa_net_all_mon'].plot(label='CERES')
plt.plot(A21.time, (A21_Aavg['FSNT']-A21_Aavg['FLNT']), label='A21')
plt.plot(M92.time, (M92_Aavg['FSNT']-M92_Aavg['FLNT']), label='M92')
plt.legend()
plt.savefig(fig_path+'net_all.png')


# Net clear-sky
plt.figure()
plt.suptitle("Net clear-sky")
ceres_Aavg['toa_net_clr_c_mon'].plot(label='CERES')
plt.plot(A21.time, (A21_Aavg['FSNTC']-A21_Aavg['FLNTC']), label='A21')
plt.plot(M92.time, (M92_Aavg['FSNTC']-M92_Aavg['FLNTC']), label='M92')
plt.legend()
plt.savefig(fig_path+'net_clr.png')

# Net clouds
plt.figure()
plt.suptitle("Net flux clouds (all-sky - clear-sky)")
(ceres_Aavg['toa_net_all_mon'] - ceres_Aavg['toa_net_clr_c_mon']).plot(label='CERES')
plt.plot(A21.time, ((A21_Aavg['FSNT']-A21_Aavg['FLNT']) - (A21_Aavg['FSNTC']-A21_Aavg['FLNTC'])), label='A21')
plt.plot(M92.time, ((M92_Aavg['FSNT']-M92_Aavg['FLNT'])- (M92_Aavg['FSNTC']-M92_Aavg['FLNTC'])), label='M92')
plt.legend()
plt.savefig(fig_path+'net_clouds.png')

plt.show()

