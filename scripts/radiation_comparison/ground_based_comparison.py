import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

# Data path
ALE_path = "/home/astridbg/Documents/data/stations/ALE_basic_rad_2004-2014/datasets/"
ALE_fname = "ALE_basic_rad_2004-2014.csv"
ALE_lat = 82.5; ALE_lon = 62.3

NYA_path = "/home/astridbg/Documents/data/stations/NYA_radiation_2006-05_etseq/datasets/"
NYA_fname = "NYA_radiation_2006-2023.csv"
NYA_lat = 78.9227; NYA_lon = 11.927300

BAR_path = "/home/astridbg/Documents/data/stations/BAR_radiation_1992-01_etseq/datasets/"
BAR_fname = "BAR_radiation_1992-2022.csv"
BAR_lat = 71.17; BAR_lon = 156.47

EUR_path = "/home/astridbg/Documents/data/stations/EUR_basic_rad_2007-2011//datasets/"
EUR_fname = "EUR_radiation_2007-2011.csv"
EUR_lat = 78.59; EUR_lon = 85.49

# Model path
noresm_path = "/home/astridbg/Documents/data/noresm_postprocessed/"

# Paths to figures
fig_path = "/home/astridbg/Documents/master/figures/radiation_comparison/groundbased_comparison/"

# Model variables to consider
model_vars = ['FSUTOA', 'FLUT', 'FLUTC', 'FLNT', 'FLNTC', 'FSNT', 'FSNTC', 'SWCFS', 'LWCFS', 'FSNS', 'FLNS', 'FSNSC', 'FLNSC', 'CLDTOT']

# Read data
ALE_df = pd.read_csv(ALE_path + ALE_fname, index_col=0)
NYA_df = pd.read_csv(NYA_path + NYA_fname, index_col=0)
BAR_df = pd.read_csv(BAR_path + BAR_fname, index_col=0)
EUR_df = pd.read_csv(EUR_path + EUR_fname, index_col=0)
A21 = xr.open_mfdataset([noresm_path+var+'_andenes21_20220222_2007-04-15_2010-03-15.nc' for var in model_vars])
M92 = xr.open_mfdataset([noresm_path+var+'_meyers92_20220210_2007-04-15_2010-03-15.nc' for var in model_vars])
print(A21.FLUT)
quit()


### ALERT
"""
# Solar flux

plt.figure()
plt.suptitle("Alert: Net incoming solar flux at surface")
#plt.plot(A21.time, ALE_df['SWD [W/m**2]'].loc["2007-04":"2010-03"], label='Alert station', color='black')
(ALE_df['SWD [W/m**2]'] - ALE_df['SWU [W/m**2]']).loc["2007-04":"2010-03"].plot(label='Alert station', color='black')
#plt.plot(ceres_Aavg.time, ceres_Aavg['toa_lw_all_mon']+se_toa_lw_all, color='tab:blue', linestyle='--')
#plt.plot(ceres_Aavg.time, ceres_Aavg['toa_lw_all_mon']-se_toa_lw_all, color='tab:blue', linestyle='--')
plt.plot(ALE_df.index[32:68],A21['FSNS'].sel(lat=ALE_lat, lon=ALE_lon, method='nearest'), label='A21', color='tab:orange' )
plt.plot(ALE_df.index[32:68],M92['FSNS'].sel(lat=ALE_lat, lon=ALE_lon, method='nearest'), label='M92', color='tab:green' )
plt.legend()
plt.xticks(rotation=45, ha='right')
plt.ylabel(r'W/m$^2$')
plt.savefig(fig_path+'alert_sw_net.png')

# Longwave flux

plt.figure()
plt.suptitle("Alert: Net outgoing longwave flux at surface")
(ALE_df['LWU [W/m**2]'] - ALE_df['LWD [W/m**2]']).loc["2007-04":"2010-03"].plot(label='Alert station', color='black')
#plt.plot(ceres_Aavg.time, ceres_Aavg['toa_lw_all_mon']+se_toa_lw_all, color='tab:blue', linestyle='--')
#plt.plot(ceres_Aavg.time, ceres_Aavg['toa_lw_all_mon']-se_toa_lw_all, color='tab:blue', linestyle='--')
plt.plot(ALE_df.index[32:68],A21['FLNS'].sel(lat=ALE_lat, lon=ALE_lon, method='nearest'), label='A21', color='tab:orange' )
plt.plot(ALE_df.index[32:68],M92['FLNS'].sel(lat=ALE_lat, lon=ALE_lon, method='nearest'), label='M92', color='tab:green' )
plt.legend()
plt.xticks(rotation=45, ha='right')
plt.ylabel(r'W/m$^2$')
plt.savefig(fig_path+'alert_lw_net.png')

# Net flux at suface

plt.figure()
plt.suptitle("Alert: Net flux at surface")
(ALE_df['SWD [W/m**2]'] + ALE_df['LWD [W/m**2]'] - ALE_df['SWU [W/m**2]'] - ALE_df['LWU [W/m**2]']).loc["2007-04":"2010-03"].plot(label='Alert station', color='black')
#plt.plot(ceres_Aavg.time, ceres_Aavg['toa_lw_all_mon']+se_toa_lw_all, color='tab:blue', linestyle='--')
#plt.plot(ceres_Aavg.time, ceres_Aavg['toa_lw_all_mon']-se_toa_lw_all, color='tab:blue', linestyle='--')
plt.plot(ALE_df.index[32:68],(A21['FSNS'] - A21['FLNS']).sel(lat=ALE_lat, lon=ALE_lon, method='nearest'), label='A21', color='tab:orange' )
plt.plot(ALE_df.index[32:68],(M92['FSNS'] - M92['FLNS']).sel(lat=ALE_lat, lon=ALE_lon, method='nearest'), label='M92', color='tab:green' )
plt.legend()
plt.xticks(rotation=45, ha='right')
plt.ylabel(r'W/m$^2$')
plt.savefig(fig_path+'alert_net_flux.png')

### NY-ÅLESUND

# Solar flux

plt.figure()
plt.suptitle("Ny-Ålesund: Net incoming solar flux at surface")
#plt.plot(A21.time, ALE_df['SWD [W/m**2]'].loc["2007-04":"2010-03"], label='Alert station', color='black')
(NYA_df['SWD [W/m**2]'] - NYA_df['SWU [W/m**2]']).loc["2007-04":"2010-03"].plot(label='NyA station', color='black')
#plt.plot(ceres_Aavg.time, ceres_Aavg['toa_lw_all_mon']+se_toa_lw_all, color='tab:blue', linestyle='--')
#plt.plot(ceres_Aavg.time, ceres_Aavg['toa_lw_all_mon']-se_toa_lw_all, color='tab:blue', linestyle='--')
plt.plot(NYA_df.index[32:68],A21['FSNS'].sel(lat=NYA_lat, lon=NYA_lon, method='nearest'), label='A21', color='tab:orange' )
plt.plot(NYA_df.index[32:68],M92['FSNS'].sel(lat=NYA_lat, lon=NYA_lon, method='nearest'), label='M92', color='tab:green' )
plt.legend()
plt.xticks(rotation=45, ha='right')
plt.ylabel(r'W/m$^2$')
plt.savefig(fig_path+'nya_sw_net.png')

# Longwave flux

plt.figure()
plt.suptitle("Ny-Ålesund: Net outgoing longwave flux at surface")
(NYA_df['LWU [W/m**2]'] - NYA_df['LWD [W/m**2]']).loc["2007-04":"2010-03"].plot(label='NyA station', color='black')
#plt.plot(ceres_Aavg.time, ceres_Aavg['toa_lw_all_mon']+se_toa_lw_all, color='tab:blue', linestyle='--')
#plt.plot(ceres_Aavg.time, ceres_Aavg['toa_lw_all_mon']-se_toa_lw_all, color='tab:blue', linestyle='--')
plt.plot(NYA_df.index[32:68],A21['FLNS'].sel(lat=NYA_lat, lon=NYA_lon, method='nearest'), label='A21', color='tab:orange' )
plt.plot(NYA_df.index[32:68],M92['FLNS'].sel(lat=NYA_lat, lon=NYA_lon, method='nearest'), label='M92', color='tab:green' )
plt.legend()
plt.xticks(rotation=45, ha='right')
plt.ylabel(r'W/m$^2$')
plt.savefig(fig_path+'nya_lw_net.png')

# Net flux at suface

plt.figure()
plt.suptitle("Ny-Ålesund: Net flux at surface")
(NYA_df['SWD [W/m**2]'] + NYA_df['LWD [W/m**2]'] - NYA_df['SWU [W/m**2]'] - NYA_df['LWU [W/m**2]']).loc["2007-04":"2010-03"].plot(label='NyA station', color='black')
#plt.plot(ceres_Aavg.time, ceres_Aavg['toa_lw_all_mon']+se_toa_lw_all, color='tab:blue', linestyle='--')
#plt.plot(ceres_Aavg.time, ceres_Aavg['toa_lw_all_mon']-se_toa_lw_all, color='tab:blue', linestyle='--')
plt.plot(NYA_df.index[32:68],(A21['FSNS'] - A21['FLNS']).sel(lat=NYA_lat, lon=NYA_lon, method='nearest'), label='A21', color='tab:orange' )
plt.plot(NYA_df.index[32:68],(M92['FSNS'] - M92['FLNS']).sel(lat=NYA_lat, lon=NYA_lon, method='nearest'), label='M92', color='tab:green' )
plt.legend()
plt.xticks(rotation=45, ha='right')
plt.ylabel(r'W/m$^2$')
plt.savefig(fig_path+'nya_net_flux.png')


### BARROW


# Solar flux

plt.figure()
plt.suptitle("Utqiaġvik (Barrow): Net incoming solar flux at surface")
#plt.plot(A21.time, ALE_df['SWD [W/m**2]'].loc["2007-04":"2010-03"], label='Alert station', color='black')
(BAR_df['SWD [W/m**2]'] - BAR_df['SWU [W/m**2]']).loc["2007-04":"2010-03"].plot(label='Observed', color='black')
#plt.plot(ceres_Aavg.time, ceres_Aavg['toa_lw_all_mon']+se_toa_lw_all, color='tab:blue', linestyle='--')
#plt.plot(ceres_Aavg.time, ceres_Aavg['toa_lw_all_mon']-se_toa_lw_all, color='tab:blue', linestyle='--')
plt.plot(np.arange(36),A21['FSNS'].sel(lat=BAR_lat, lon=BAR_lon, method='nearest'), label='A21', color='tab:orange' )
plt.plot(np.arange(36),M92['FSNS'].sel(lat=BAR_lat, lon=BAR_lon, method='nearest'), label='M92', color='tab:green' )
plt.legend()
plt.xticks(rotation=45, ha='right')
plt.ylabel(r'W/m$^2$')
plt.savefig(fig_path+'bar_sw_net.png')

# Longwave flux

plt.figure()
plt.suptitle("Utqiaġvik (Barrow): Net outgoing longwave flux at surface")
(BAR_df['LWU [W/m**2]'] - BAR_df['LWD [W/m**2]']).loc["2007-04":"2010-03"].plot(label='Observed', color='black')
#plt.plot(ceres_Aavg.time, ceres_Aavg['toa_lw_all_mon']+se_toa_lw_all, color='tab:blue', linestyle='--')
#plt.plot(ceres_Aavg.time, ceres_Aavg['toa_lw_all_mon']-se_toa_lw_all, color='tab:blue', linestyle='--')
plt.plot(np.arange(36),A21['FLNS'].sel(lat=BAR_lat, lon=BAR_lon, method='nearest'), label='A21', color='tab:orange' )
plt.plot(np.arange(36),M92['FLNS'].sel(lat=BAR_lat, lon=BAR_lon, method='nearest'), label='M92', color='tab:green' )
plt.legend()
plt.xticks(rotation=45, ha='right')
plt.ylabel(r'W/m$^2$')
plt.savefig(fig_path+'bar_lw_net.png')

# Net flux at suface

plt.figure()
plt.suptitle("Utqiaġvik (Barrow): Net flux at surface")
(BAR_df['SWD [W/m**2]'] + BAR_df['LWD [W/m**2]'] - BAR_df['SWU [W/m**2]'] - BAR_df['LWU [W/m**2]']).loc["2007-04":"2010-03"].plot(label='Observed', color='black')
#plt.plot(ceres_Aavg.time, ceres_Aavg['toa_lw_all_mon']+se_toa_lw_all, color='tab:blue', linestyle='--')
#plt.plot(ceres_Aavg.time, ceres_Aavg['toa_lw_all_mon']-se_toa_lw_all, color='tab:blue', linestyle='--')
plt.plot(np.arange(36),(A21['FSNS'] - A21['FLNS']).sel(lat=BAR_lat, lon=BAR_lon, method='nearest'), label='A21', color='tab:orange' )
plt.plot(np.arange(36),(M92['FSNS'] - M92['FLNS']).sel(lat=BAR_lat, lon=BAR_lon, method='nearest'), label='M92', color='tab:green' )
plt.legend()
plt.xticks(rotation=45, ha='right')
plt.ylabel(r'W/m$^2$')
plt.savefig(fig_path+'bar_net_flux.png')
"""
### EUREKA
#Index(['Height [m]', 'SWD [W/m**2]', 'SWD std dev [W/m**2]',
#       'SWD min [W/m**2]', 'SWD max [W/m**2]', 'DIR [W/m**2]',
#       'DIR std dev [W/m**2]', 'DIR min [W/m**2]', 'DIR max [W/m**2]',
#       'DIF [W/m**2]', 'DIF std dev [W/m**2]', 'DIF min [W/m**2]',
#       'DIF max [W/m**2]', 'LWD [W/m**2]', 'LWD std dev [W/m**2]',
#       'LWD min [W/m**2]', 'LWD max [W/m**2]'],
#      dtype='object')

# Solar flux


plt.figure()
plt.suptitle("Eureka: Only incoming solar flux at surface")
#plt.plot(A21.time, ALE_df['SWD [W/m**2]'].loc["2007-04":"2010-03"], label='Alert station', color='black')
#(EUR_df['SWD [W/m**2]'] - EUR_df['SWU [W/m**2]']).loc["2007-04":"2010-03"].plot(label='Observed', color='black')
EUR_df['SWD [W/m**2]'].loc["2007-09":"2010-03"].plot(label='Observed', color='black')
#plt.plot(ceres_Aavg.time, ceres_Aavg['toa_lw_all_mon']+se_toa_lw_all, color='tab:blue', linestyle='--')
#plt.plot(ceres_Aavg.time, ceres_Aavg['toa_lw_all_mon']-se_toa_lw_all, color='tab:blue', linestyle='--')
plt.plot(np.arange(36),A21['FSNS'].sel(lat=EUR_lat, lon=EUR_lon, method='nearest'), label='A21', color='tab:orange' )
plt.plot(np.arange(36),M92['FSNS'].sel(lat=EUR_lat, lon=EUR_lon, method='nearest'), label='M92', color='tab:green' )
plt.legend()
plt.xticks(rotation=45, ha='right')
plt.ylabel(r'W/m$^2$')
plt.savefig(fig_path+'eur_sw_net.png')
"""
# Longwave flux

plt.figure()
plt.suptitle("Eureka: Net outgoing longwave flux at surface")
(NYA_df['LWU [W/m**2]'] - NYA_df['LWD [W/m**2]']).loc["2007-04":"2010-03"].plot(label='Observed', color='black')
#plt.plot(ceres_Aavg.time, ceres_Aavg['toa_lw_all_mon']+se_toa_lw_all, color='tab:blue', linestyle='--')
#plt.plot(ceres_Aavg.time, ceres_Aavg['toa_lw_all_mon']-se_toa_lw_all, color='tab:blue', linestyle='--')
plt.plot(np.arange(36),A21['FLNS'].sel(lat=NYA_lat, lon=NYA_lon, method='nearest'), label='A21', color='tab:orange' )
plt.plot(np.arange(36),M92['FLNS'].sel(lat=NYA_lat, lon=NYA_lon, method='nearest'), label='M92', color='tab:green' )
plt.legend()
plt.xticks(rotation=45, ha='right')
plt.ylabel(r'W/m$^2$')
plt.savefig(fig_path+'eur_lw_net.png')

# Net flux at suface

plt.figure()
plt.suptitle("Eureka: Net flux at surface")
(NYA_df['SWD [W/m**2]'] + NYA_df['LWD [W/m**2]'] - NYA_df['SWU [W/m**2]'] - NYA_df['LWU [W/m**2]']).loc["2007-04":"2010-03"].plot(label='Observed', color='black')
#plt.plot(ceres_Aavg.time, ceres_Aavg['toa_lw_all_mon']+se_toa_lw_all, color='tab:blue', linestyle='--')
#plt.plot(ceres_Aavg.time, ceres_Aavg['toa_lw_all_mon']-se_toa_lw_all, color='tab:blue', linestyle='--')
plt.plot(np.arange(36),(A21['FSNS'] - A21['FLNS']).sel(lat=NYA_lat, lon=NYA_lon, method='nearest'), label='A21', color='tab:orange' )
plt.plot(np.arange(36),(M92['FSNS'] - M92['FLNS']).sel(lat=NYA_lat, lon=NYA_lon, method='nearest'), label='M92', color='tab:green' )
plt.legend()
plt.xticks(rotation=45, ha='right')
plt.ylabel(r'W/m$^2$')
plt.savefig(fig_path+'eur_net_flux.png')
"""