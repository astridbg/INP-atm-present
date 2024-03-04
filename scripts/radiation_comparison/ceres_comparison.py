import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd
from functions import *
import glob

# Paths to data
ceres_path = "/home/astridbg/Documents/data/CERES_EBAF_Ed4.2_Subset_200704-201003.nc"
noresm_path = "/home/astridbg/Documents/data/noresm_postprocessed/"

# Paths to figures
fig_path = "/home/astridbg/Documents/master/figures/radiation_comparison/ceres_comparison/"

# Model variables to consider
model_vars = ['FSUTOA', 'FLUT', 'FLUTC', 'FLNT', 'FLNTC', 'FSNT', 'FSNTC', 'SWCFS', 'LWCFS', 'FSNS', 'FLNS', 'FSNSC', 'FLNSC', 'CLDTOT']


# Ceres uncertainty (https://ceres.larc.nasa.gov/documents/DQ_summaries/CERES_EBAF_Ed4.1_DQS.pdf)
se_toa_sw_all, se_toa_sw_clr, se_toa_sw_cre = [2.5, 5.4, 5.9]
se_toa_lw_all, se_toa_lw_clr, se_toa_lw_cre = [2.5, 4.6, 4.5]
se_toa_net_all, se_toa_net_clr, se_toa_net_cre = [3.5, 7.1, 7.4]

se_sfc_lwd_all, se_sfc_lwd_clr, se_sfc_lwd_cre = [9, 8, 9]
se_sfc_lwu_all, se_sfc_lwu_clr, se_sfc_lwu_cre = [15, 15, 17]
se_sfc_lwn_all, se_sfc_lwn_clr, se_sfc_lwn_cre = [17, 17, 18]
se_sfc_swd_all, se_sfc_swd_clr, se_sfc_swd_cre = [14, 6, 14]
se_sfc_swu_all, se_sfc_swu_clr, se_sfc_swu_cre = [11, 11, 14]
se_sfc_swn_all, se_sfc_swn_clr, se_sfc_swn_cre = [13, 13, 16]
se_sfc_net_all, se_sfc_net_clr, se_sfc_net_cre = [20, 21, 26]

ceres_data = xr.open_dataset(ceres_path)
A21 = xr.open_mfdataset([noresm_path+var+'_andenes21_20220222_2007-04-15_2010-03-15.nc' for var in model_vars])
M92 = xr.open_mfdataset([noresm_path+var+'_meyers92_20220210_2007-04-15_2010-03-15.nc' for var in model_vars])

print(ceres_data.variables)

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
ceres_Aavg['toa_cre_lw_mon'].plot(label='lw cre')
ceres_Aavg['toa_cre_sw_mon'].plot(label='sw cre')
plt.legend()
plt.savefig(fig_path+'ceres_toa.png')

plt.figure()
ceres_Aavg['sfc_sw_down_all_mon'].plot(label='sw down')
ceres_Aavg['sfc_sw_up_all_mon'].plot(label='sw up')
ceres_Aavg['sfc_net_sw_all_mon'].plot(label='sw net')
ceres_Aavg['sfc_sw_down_clr_c_mon'].plot(label='sw down clr')
ceres_Aavg['sfc_sw_up_clr_c_mon'].plot(label='sw up clr')
ceres_Aavg['sfc_net_sw_clr_c_mon'].plot(label='sw net clr')
plt.legend()
plt.savefig(fig_path+'ceres_sfc_sw.png')

plt.figure()
ceres_Aavg['sfc_lw_down_all_mon'].plot(label='lw down')
ceres_Aavg['sfc_lw_up_all_mon'].plot(label='lw up')
ceres_Aavg['sfc_net_lw_all_mon'].plot(label='lw net')
ceres_Aavg['sfc_lw_down_clr_c_mon'].plot(label='lw down clr')
ceres_Aavg['sfc_lw_up_clr_c_mon'].plot(label='lw up clr')
ceres_Aavg['sfc_net_lw_clr_c_mon'].plot(label='lw net clr')
plt.legend()
plt.savefig(fig_path+'ceres_sfc_lw.png')

plt.figure()
ceres_Aavg['sfc_cre_net_sw_mon'].plot(label='sw net')
ceres_Aavg['sfc_cre_net_lw_mon'].plot(label='lw up')
ceres_Aavg['sfc_cre_net_tot_mon'].plot(label='tot net')
plt.legend()
plt.savefig(fig_path+'ceres_sfc_cre.png')


### TOA FLUXES

# Outgoing longwave allsky
plt.figure()
plt.suptitle("Outgoing longwave all-sky")
plt.plot(ceres_Aavg.time, ceres_Aavg['toa_lw_all_mon'], label='CERES', color='tab:blue')
plt.plot(ceres_Aavg.time, ceres_Aavg['toa_lw_all_mon']+se_toa_lw_all, color='tab:blue', linestyle='--')
plt.plot(ceres_Aavg.time, ceres_Aavg['toa_lw_all_mon']-se_toa_lw_all, color='tab:blue', linestyle='--')
plt.plot(A21.time, A21_Aavg['FLUT'], label='A21', color='tab:orange')
plt.plot(M92.time, M92_Aavg['FLUT'], label='M92', color='tab:green')
plt.legend()
plt.xticks(rotation=45, ha='right')
plt.ylabel(r'W/m$^2$')
plt.savefig(fig_path+'lw_all.png')

# Outgoing longwave clearsky
plt.figure()
plt.suptitle("Outgoing longwave clear-sky")
plt.plot(ceres_Aavg.time, ceres_Aavg['toa_lw_clr_c_mon'], label='CERES', color='tab:blue')
plt.plot(ceres_Aavg.time, ceres_Aavg['toa_lw_clr_c_mon']+se_toa_lw_clr, color='tab:blue', linestyle='--')
plt.plot(ceres_Aavg.time, ceres_Aavg['toa_lw_clr_c_mon']-se_toa_lw_clr, color='tab:blue', linestyle='--')
plt.plot(A21.time, A21_Aavg['FLUTC'], label='A21', color='tab:orange')
plt.plot(M92.time, M92_Aavg['FLUTC'], label='M92', color='tab:green')
plt.legend()
plt.xticks(rotation=45, ha='right')
plt.ylabel(r'W/m$^2$')
plt.savefig(fig_path+'lw_clr.png')

# Outgoing longwave clouds
plt.figure()
plt.suptitle("Outgoing longwave clouds (all-sky - clear-sky)")
plt.plot(ceres_Aavg.time, ceres_Aavg['toa_lw_all_mon'] - ceres_Aavg['toa_lw_clr_c_mon'], label='CERES', color='tab:blue')
plt.plot(ceres_Aavg.time, ceres_Aavg['toa_lw_all_mon'] - ceres_Aavg['toa_lw_clr_c_mon']+se_toa_lw_cre, color='tab:blue', linestyle='--')
plt.plot(ceres_Aavg.time, ceres_Aavg['toa_lw_all_mon'] - ceres_Aavg['toa_lw_clr_c_mon']-se_toa_lw_cre, color='tab:blue', linestyle='--')
plt.plot(A21.time, (A21_Aavg['FLUT']-A21_Aavg['FLUTC']), label='A21', color='tab:orange')
plt.plot(M92.time, (M92_Aavg['FLUT'] - M92_Aavg['FLUTC']), label='M92', color='tab:green')
plt.legend()
plt.xticks(rotation=45, ha='right')
plt.ylabel(r'W/m$^2$')
plt.savefig(fig_path+'lw_clouds.png')

# Outgoing longwave clouds
plt.figure()
plt.suptitle("Net longwave clouds (all-sky - clear-sky)")
plt.plot(ceres_Aavg.time, ceres_Aavg['toa_lw_all_mon'] - ceres_Aavg['toa_lw_clr_c_mon'], label='CERES out', color='tab:blue')
plt.plot(ceres_Aavg.time, ceres_Aavg['toa_lw_all_mon'] - ceres_Aavg['toa_lw_clr_c_mon']+se_toa_lw_cre, color='tab:blue', linestyle='--')
plt.plot(ceres_Aavg.time, ceres_Aavg['toa_lw_all_mon'] - ceres_Aavg['toa_lw_clr_c_mon']-se_toa_lw_cre, color='tab:blue', linestyle='--')
plt.plot(A21.time, (A21_Aavg['FLNT']-A21_Aavg['FLNTC']), label='A21', color='tab:orange')
plt.plot(M92.time, (M92_Aavg['FLNT'] - M92_Aavg['FLNTC']), label='M92', color='tab:green')
plt.legend()
plt.xticks(rotation=45, ha='right')
plt.ylabel(r'W/m$^2$')
plt.savefig(fig_path+'lw_net_clouds.png')


# Outgoing shortwave
plt.figure()
plt.suptitle("Outgoing shortwave all-sky")
plt.plot(ceres_Aavg.time, ceres_Aavg['toa_sw_all_mon'], label='CERES', color='tab:blue')
plt.plot(ceres_Aavg.time, ceres_Aavg['toa_sw_all_mon']+se_toa_sw_all, color='tab:blue', linestyle='--')
plt.plot(ceres_Aavg.time, ceres_Aavg['toa_sw_all_mon']-se_toa_sw_all, color='tab:blue', linestyle='--')
plt.plot(A21.time, A21_Aavg['FSUTOA'], label='A21', color='tab:orange')
plt.plot(M92.time, M92_Aavg['FSUTOA'], label='M92', color='tab:green')
plt.legend()
plt.xticks(rotation=45, ha='right')
plt.ylabel(r'W/m$^2$')
plt.savefig(fig_path+'sw_all.png')


# Net solar flux all-sky
plt.figure()
plt.suptitle("Net shortwave all-sky")
plt.plot(ceres_Aavg.time, ceres_Aavg['solar_mon']- ceres_Aavg['toa_sw_all_mon'], label='CERES in-out', color='tab:blue')
plt.plot(ceres_Aavg.time, ceres_Aavg['solar_mon']- ceres_Aavg['toa_sw_all_mon']+se_toa_sw_all, color='tab:blue', linestyle='--')
plt.plot(ceres_Aavg.time, ceres_Aavg['solar_mon']- ceres_Aavg['toa_sw_all_mon']-se_toa_sw_all, color='tab:blue', linestyle='--')
plt.plot(A21.time, A21_Aavg['FSNT'], label='A21', color='tab:orange')
plt.plot(M92.time, M92_Aavg['FSNT'], label='M92', color='tab:green')
plt.legend()
plt.xticks(rotation=45, ha='right')
plt.ylabel(r'W/m$^2$')
plt.savefig(fig_path+'sw_net_all.png')


# Net solar flux clear-sky
plt.figure()
plt.suptitle("Net shortwave clear-sky")
plt.plot(ceres_Aavg.time, ceres_Aavg['solar_mon']- ceres_Aavg['toa_sw_clr_c_mon'], label='CERES in-out', color='tab:blue')
plt.plot(ceres_Aavg.time, ceres_Aavg['solar_mon']- ceres_Aavg['toa_sw_clr_c_mon']+se_toa_sw_clr, color='tab:blue', linestyle='--')
plt.plot(ceres_Aavg.time, ceres_Aavg['solar_mon']- ceres_Aavg['toa_sw_clr_c_mon']-se_toa_sw_clr, color='tab:blue', linestyle='--')
plt.plot(A21.time, A21_Aavg['FSNTC'], label='A21', color='tab:orange')
plt.plot(M92.time, M92_Aavg['FSNTC'], label='M92', color='tab:green')
plt.legend()
plt.xticks(rotation=45, ha='right')
plt.ylabel(r'W/m$^2$')
plt.savefig(fig_path+'sw_net_clr.png')

# Net solar flux clouds
plt.figure()
plt.suptitle("Net shortwave clouds (clear-sky - all-sky)")
plt.plot(ceres_Aavg.time, ceres_Aavg['toa_sw_clr_c_mon']- ceres_Aavg['toa_sw_all_mon'], label='CERES', color='tab:blue')
plt.plot(ceres_Aavg.time, ceres_Aavg['toa_sw_clr_c_mon']- ceres_Aavg['toa_sw_all_mon']+se_toa_sw_cre, color='tab:blue', linestyle='--')
plt.plot(ceres_Aavg.time, ceres_Aavg['toa_sw_clr_c_mon']- ceres_Aavg['toa_sw_all_mon']-se_toa_sw_cre, color='tab:blue', linestyle='--')
plt.plot(A21.time, (A21_Aavg['FSNT']-A21_Aavg['FSNTC']), label='A21', color='tab:orange')
plt.plot(M92.time, (M92_Aavg['FSNT'] - M92_Aavg['FSNTC']), label='M92', color='tab:green')
plt.legend()
plt.xticks(rotation=45, ha='right')
plt.ylabel(r'W/m$^2$')
plt.savefig(fig_path+'sw_net_clouds.png')

# Net all-sky
plt.figure()
plt.suptitle("Net flux all-sky")
plt.plot(ceres_Aavg.time, ceres_Aavg['toa_net_all_mon'], label='CERES', color='tab:blue')
plt.plot(ceres_Aavg.time, ceres_Aavg['toa_net_all_mon']+se_toa_net_all, color='tab:blue', linestyle='--')
plt.plot(ceres_Aavg.time, ceres_Aavg['toa_net_all_mon']-se_toa_net_all, color='tab:blue', linestyle='--')
plt.plot(A21.time, (A21_Aavg['FSNT']-A21_Aavg['FLNT']), label='A21', color='tab:orange')
plt.plot(M92.time, (M92_Aavg['FSNT']-M92_Aavg['FLNT']), label='M92', color='tab:green')
plt.legend()
plt.xticks(rotation=45, ha='right')
plt.ylabel(r'W/m$^2$')
plt.savefig(fig_path+'net_all.png')


# Net clear-sky
plt.figure()
plt.suptitle("Net clear-sky")
plt.plot(ceres_Aavg.time, ceres_Aavg['toa_net_clr_c_mon'], label='CERES', color='tab:blue')
plt.plot(ceres_Aavg.time, ceres_Aavg['toa_net_clr_c_mon']+se_toa_net_clr, color='tab:blue', linestyle='--')
plt.plot(ceres_Aavg.time, ceres_Aavg['toa_net_clr_c_mon']-se_toa_net_clr, color='tab:blue', linestyle='--')
plt.plot(A21.time, (A21_Aavg['FSNTC']-A21_Aavg['FLNTC']), label='A21', color='tab:orange')
plt.plot(M92.time, (M92_Aavg['FSNTC']-M92_Aavg['FLNTC']), label='M92', color='tab:green')
plt.legend()
plt.xticks(rotation=45, ha='right')
plt.ylabel(r'W/m$^2$')
plt.savefig(fig_path+'net_clr.png')

# Net clouds
plt.figure()
plt.suptitle("Net flux clouds (all-sky - clear-sky)")
plt.plot(ceres_Aavg.time, ceres_Aavg['toa_net_all_mon'] - ceres_Aavg['toa_net_clr_c_mon'], label='CERES', color='tab:blue')
plt.plot(ceres_Aavg.time, ceres_Aavg['toa_net_all_mon'] - ceres_Aavg['toa_net_clr_c_mon']+se_toa_net_cre, color='tab:blue', linestyle='--')
plt.plot(ceres_Aavg.time, ceres_Aavg['toa_net_all_mon'] - ceres_Aavg['toa_net_clr_c_mon']-se_toa_net_cre, color='tab:blue', linestyle='--')
plt.plot(A21.time, ((A21_Aavg['FSNT']-A21_Aavg['FLNT']) - (A21_Aavg['FSNTC']-A21_Aavg['FLNTC'])), label='A21', color='tab:orange')
plt.plot(M92.time, ((M92_Aavg['FSNT']-M92_Aavg['FLNT'])- (M92_Aavg['FSNTC']-M92_Aavg['FLNTC'])), label='M92', color='tab:green')
plt.legend()
plt.xticks(rotation=45, ha='right')
plt.ylabel(r'W/m$^2$')
plt.savefig(fig_path+'net_clouds.png')

### SURFACE FLUXES
plt.figure()
plt.suptitle("Surface Net LW Flux")
plt.plot(ceres_Aavg.time, ceres_Aavg['sfc_net_lw_all_mon'], label='CERES', color='tab:blue')
plt.plot(ceres_Aavg.time, ceres_Aavg['sfc_net_lw_all_mon']+se_sfc_lwn_all, color='tab:blue', linestyle='--')
plt.plot(ceres_Aavg.time, ceres_Aavg['sfc_net_lw_all_mon']-se_sfc_lwn_all, color='tab:blue', linestyle='--')
plt.plot(A21.time, -A21_Aavg['FLNS'], label='A21', color='tab:orange')
plt.plot(M92.time, -M92_Aavg['FLNS'], label='M92', color='tab:green')
plt.legend()
plt.xticks(rotation=45, ha='right')
plt.ylabel(r'W/m$^2$')
plt.savefig(fig_path+'sfc_lw.png')

plt.figure()
plt.suptitle("Surface Net LW Flux Clear-sky")
plt.plot(ceres_Aavg.time, ceres_Aavg['sfc_net_lw_clr_t_mon'], label='CERES', color='tab:blue')
plt.plot(ceres_Aavg.time, ceres_Aavg['sfc_net_lw_clr_t_mon']+se_sfc_lwn_all, color='tab:blue', linestyle='--')
plt.plot(ceres_Aavg.time, ceres_Aavg['sfc_net_lw_clr_t_mon']-se_sfc_lwn_all, color='tab:blue', linestyle='--')
plt.plot(A21.time, -A21_Aavg['FLNSC'], label='A21', color='tab:orange')
plt.plot(M92.time, -M92_Aavg['FLNSC'], label='M92', color='tab:green')
plt.legend()
plt.xticks(rotation=45, ha='right')
plt.ylabel(r'W/m$^2$')
plt.savefig(fig_path+'sfc_lw_clr.png')

plt.figure()
plt.suptitle("Surface Net SW Flux")
plt.plot(ceres_Aavg.time, ceres_Aavg['sfc_net_sw_all_mon'], label='CERES', color='tab:blue')
plt.plot(ceres_Aavg.time, ceres_Aavg['sfc_net_sw_all_mon']+se_sfc_swn_all, color='tab:blue', linestyle='--')
plt.plot(ceres_Aavg.time, ceres_Aavg['sfc_net_sw_all_mon']-se_sfc_swn_all, color='tab:blue', linestyle='--')
plt.plot(A21.time, A21_Aavg['FSNS'], label='A21', color='tab:orange')
plt.plot(M92.time, M92_Aavg['FSNS'], label='M92', color='tab:green')
plt.legend()
plt.xticks(rotation=45, ha='right')
plt.ylabel(r'W/m$^2$')
plt.savefig(fig_path+'sfc_sw_cre.png')


plt.figure()
plt.suptitle("Surface Net SW CRE")
plt.plot(ceres_Aavg.time, ceres_Aavg['sfc_cre_net_sw_mon'], label='CERES', color='tab:blue')
plt.plot(ceres_Aavg.time, ceres_Aavg['sfc_cre_net_sw_mon']+se_sfc_swn_cre, color='tab:blue', linestyle='--')
plt.plot(ceres_Aavg.time, ceres_Aavg['sfc_cre_net_sw_mon']-se_sfc_swn_cre, color='tab:blue', linestyle='--')
plt.plot(A21.time, A21_Aavg['SWCFS'], label='A21', color='tab:orange')
plt.plot(M92.time, M92_Aavg['SWCFS'], label='M92', color='tab:green')
plt.legend()
plt.xticks(rotation=45, ha='right')
plt.ylabel(r'W/m$^2$')
plt.savefig(fig_path+'sfc_sw_cre.png')

plt.figure()
plt.suptitle("Surface Net LW CRE")
plt.plot(ceres_Aavg.time, ceres_Aavg['sfc_cre_net_lw_mon'], label='CERES', color='tab:blue')
plt.plot(ceres_Aavg.time, ceres_Aavg['sfc_cre_net_lw_mon']+se_sfc_lwn_cre, color='tab:blue', linestyle='--')
plt.plot(ceres_Aavg.time, ceres_Aavg['sfc_cre_net_lw_mon']-se_sfc_lwn_cre, color='tab:blue', linestyle='--')
plt.plot(A21.time, A21_Aavg['LWCFS'], label='A21', color='tab:orange')
plt.plot(M92.time, M92_Aavg['LWCFS'], label='M92', color='tab:green')
plt.legend()
plt.xticks(rotation=45, ha='right')
plt.ylabel(r'W/m$^2$')
plt.savefig(fig_path+'sfc_lw_cre.png')

plt.figure()
plt.suptitle("Surface Net total CRE")
plt.plot(ceres_Aavg.time, ceres_Aavg['sfc_cre_net_tot_mon'], label='CERES', color='tab:blue')
plt.plot(ceres_Aavg.time, ceres_Aavg['sfc_cre_net_tot_mon']+se_sfc_net_cre, color='tab:blue', linestyle='--')
plt.plot(ceres_Aavg.time, ceres_Aavg['sfc_cre_net_tot_mon']-se_sfc_net_cre, color='tab:blue', linestyle='--')
plt.plot(A21.time, (A21_Aavg['SWCFS']+A21_Aavg['LWCFS']), label='A21', color='tab:orange')
plt.plot(M92.time, (M92_Aavg['SWCFS']+M92_Aavg['LWCFS']), label='M92', color='tab:green')
plt.legend()
plt.xticks(rotation=45, ha='right')
plt.ylabel(r'W/m$^2$')
plt.savefig(fig_path+'sfc_tot_cre.png')



### CLOUD CHARACTERISTICS
plt.figure()
plt.suptitle("Total cloud fraction")
plt.plot(ceres_Aavg.time, ceres_Aavg['cldarea_total_daynight_mon'], label='CERES', color='tab:blue')
plt.plot(A21.time, A21_Aavg['CLDTOT']*100, label='A21', color='tab:orange')
plt.plot(M92.time, M92_Aavg['CLDTOT']*100, label='M92', color='tab:green')
plt.legend()
plt.xticks(rotation=45, ha='right')
plt.ylabel('%')
plt.savefig(fig_path+'cldfrc_tot.png')

plt.show()

