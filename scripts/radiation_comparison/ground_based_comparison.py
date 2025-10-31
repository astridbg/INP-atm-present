import xarray as xr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter
# Set font style to match latex document----------
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size':20})
# ------------------------------------------------
from functions import *

model_rpath="/projects/NS9600K/astridbg/data/model/noresm_postprocessed/"
obs_rpath="/projects/NS9600K/astridbg/data/observations/radiation/stations/"
wpath="/projects/NS9600K/astridbg/INP-atm-present/figures/radiation_comparison/groundbased_comparison/"
#wpath = "/projects/NS9600K/astridbg/INP-atm-present/figures/radiation_comparison/groundbased_comparison/20240522/"

station_paths = ["ALE_basic_rad_2004-2014/datasets/ALE_basic_rad_2004-2014.csv",
                 "BAR_radiation_1992-01_etseq/datasets/BAR_radiation_1992-2022.csv",
                 "NYA_radiation_2006-05_etseq/datasets/NYA_radiation_2006-2023.csv"]

#-------------------------------------
# Station coordinates and names
#-------------------------------------

ALE_lat = 82.5; ALE_lon = 360-62.3
NYA_lat = 78.9227; NYA_lon = 11.927300
BAR_lat = 71.17; BAR_lon = 360-156.47

station_coords = [[ALE_lat, ALE_lon], 
                  [BAR_lat, BAR_lon], 
                  [NYA_lat, NYA_lon]]
station_names = ["Alert", "Utqiaġvik", "Ny-Ålesund"]
#-------------------------------------
# Read data
#-------------------------------------

obs_dfs = []
for i in range(len(station_paths)):
    df = pd.read_csv(obs_rpath+station_paths[i], index_col=0)
    obs_dfs.append(df)

# Model variables to consider
#model_vars = ['SWCFS', 'LWCFS', 'FSNS', 'FLNS', 'FSNSC', 'FLNSC', 'CLDTOT', 'FLDS', 'FSDS']
model_vars = ['FSNS', 'FLNS']
#A21 = xr.open_mfdataset([model_rpath+var+'_andenes21_20220222_2007-04-15_2010-03-15.nc' for var in model_vars])
#M92 = xr.open_mfdataset([model_rpath+var+'_meyers92_20220210_2007-04-15_2010-03-15.nc' for var in model_vars])
A21 = xr.open_mfdataset([model_rpath+var+'_A21_20240612_2007-04-15_2010-03-15.nc' for var in model_vars])
#M92 = xr.open_mfdataset([model_rpath+var+'_M92_20240612_2007-04-15_2010-03-15.nc' for var in model_vars])
S23 = xr.open_mfdataset([model_rpath+var+'_NorESM2_3_slf_output_2007-04-15_2010-03-15.nc' for var in model_vars])

#-------------------------------------
# Three subplots in one
#-------------------------------------

# Net Surface Flux

fig, axs = plt.subplots(3, 1, sharex=True, figsize=[13,12])
for i, df, coord, name in zip(range(3), obs_dfs, station_coords, station_names): 
    #Plot observations
    obs = (df['SWD [W/m**2]'] + df['LWD [W/m**2]'] - df['SWU [W/m**2]'] - df['LWU [W/m**2]']).loc["2007-04":"2010-03"]
    obs.plot(label='Observed', color='black', ax=axs[i], linewidth=2)
    #Plot uncertainty
    axs[i].fill_between(np.arange(36), obs - obs*(0.095), obs + obs*(0.095), color='black', alpha=0.1)
    #Plot model
    axs[i].plot(np.arange(36),(A21['FSNS'] - A21['FLNS']).sel(lat=coord[0], lon=coord[1], method='nearest'), label='A21', color='tab:orange', linewidth=2 )
#    axs[i].plot(np.arange(36),(M92['FSNS'] - M92['FLNS']).sel(lat=coord[0], lon=coord[1], method='nearest'), label='M92', color='tab:blue',linewidth=2, linestyle='--' )
    axs[i].plot(np.arange(36),(S23['FSNS'] - S23['FLNS']).sel(lat=coord[0], lon=coord[1], method='nearest'), label='NorESM2.3', color='tab:blue',linewidth=2, linestyle='--' )

    axs[i].set_title(name)
    if i==0:
        axs[i].legend(loc='upper center', bbox_to_anchor=(0.5, 1.45),ncol=3)
    axs[i].set_ylabel(r'W/m$^2$')
    axs[i].set_xlabel('')
    axs[i].set_xticks(np.arange(1, 38, 4))
    axs[i].grid(alpha=0.5)
    axs[i].tick_params('x',labelrotation=45)

#fig.savefig(wpath+'pdf/stations_all_net.pdf', bbox_inches="tight")
#fig.savefig(wpath+'png/stations_all_net.png', bbox_inches="tight")
fig.savefig(wpath+'png/stations_all_net_A21_S23.png', bbox_inches="tight")

plt.clf()

#-------------------------------------
# Monthly averages
#-------------------------------------

A21_month = A21.groupby("time.month").mean("time")
#M92_month = M92.groupby("time.month").mean("time")
S23_month = S23.groupby("time.month").mean("time")


obs = obs_dfs[0].loc["2007-04":"2010-03"]
obs.index = pd.to_datetime(obs.index)

monthsn = np.arange(1,13,1)
months = ['Jan', 'Feb', 'Mar', 'Apr', 'May','Jun','Jul','Aug','Sep','Oct','Nov', 'Dec']

# Net Surface Flux
labels = ['(a)', '(b)', '(c)']
fig, axs = plt.subplots(3, 1, sharex=True, figsize=[13,12])
for i, df, coord, name, label in zip(range(3), obs_dfs, station_coords, station_names, labels): 
    #Plot observations
    df = df.loc["2007-04":"2010-03"]
    df.index = pd.to_datetime(df.index)
    df = df.groupby(df.index.month).mean()
    obs = (df['SWD [W/m**2]'] + df['LWD [W/m**2]'] - df['SWU [W/m**2]'] - df['LWU [W/m**2]'])
    obs.plot(label='Observed', color='black', ax=axs[i], linewidth=2)
    #Plot uncertainty
    axs[i].fill_between(np.arange(1,13,1), obs - obs*(0.095), obs + obs*(0.095), color='black', alpha=0.1)
    #Plot model
    axs[i].plot(np.arange(1,13,1),(A21_month['FSNS'] - A21_month['FLNS']).sel(lat=coord[0], lon=coord[1], method='nearest'), label='A21', color='tab:orange', linewidth=2 )
    #axs[i].plot(np.arange(1,13,1),(M92_month['FSNS'] - M92_month['FLNS']).sel(lat=coord[0], lon=coord[1], method='nearest'), label='M92', color='tab:blue',linewidth=2, linestyle='--' )
    axs[i].plot(np.arange(1,13,1),(S23_month['FSNS'] - S23_month['FLNS']).sel(lat=coord[0], lon=coord[1], method='nearest'), label='NorESM2.3', color='tab:blue',linewidth=2, linestyle='--' )
    axs[i].set_title(name)
    if i==0:
        axs[i].legend(loc='upper center', bbox_to_anchor=(0.5, 1.45),ncol=3)
    axs[i].annotate(label,fontsize=20,
        xy=(0, 1), xycoords='axes fraction',
        xytext=(-30, 30), textcoords='offset points',
        ha='left', va='top')
    axs[i].set_ylabel(r'W/m$^2$')
    axs[i].set_ylim([-55,160])
    axs[i].set_xlabel('')
    axs[i].set_xticks(np.arange(1, 13, 1),months)
    axs[i].grid(alpha=0.5)
    axs[i].tick_params('x',labelrotation=45)

#fig.savefig(wpath+'pdf/stations_all_net_monthly.pdf', bbox_inches="tight")
fig.savefig(wpath+'png/stations_all_net_monthly_A21_NorESM2_3.png', bbox_inches="tight")

plt.clf()

"""
# Downwelling Longwave Flux at Surface

fig, axs = plt.subplots(3, 1, sharex=True, figsize=[13,12])
for i, df, coord, name in zip(range(3), obs_dfs, station_coords, station_names): 
    #Plot observations
    obs = df['LWD [W/m**2]'].loc["2007-04":"2010-03"]
    obs.plot(label='Observed', color='black', ax=axs[i], linewidth=2)
    #Plot uncertainty
    axs[i].fill_between(np.arange(36), obs - obs*(0.02), obs + obs*(0.02), color='black', alpha=0.1)
    #Plot model
    axs[i].plot(np.arange(36),M92['FLDS'].sel(lat=coord[0], lon=coord[1], method='nearest'), label='M92', color='tab:blue',linewidth=2 )
    axs[i].plot(np.arange(36),A21['FLDS'].sel(lat=coord[0], lon=coord[1], method='nearest'), label='A21', color='tab:orange',linewidth=2 )

    axs[i].set_title(name)
    if i==0:
        axs[i].legend()
    axs[i].set_ylabel(r'W/m$^2$')
    axs[i].set_xlabel('')
    axs[i].set_xticks(np.arange(1, 38, 4))
    axs[i].grid(alpha=0.5)
axs[2].tick_params(axis='x', labelrotation=45)
fig.savefig(wpath+'pdf/stations_lwd.pdf', bbox_inches="tight")
fig.savefig(wpath+'png/stations_lwd.png', bbox_inches="tight")

plt.clf()

# Upwelling Longwave Flux at Surface

fig, axs = plt.subplots(3, 1, sharex=True, figsize=[13,12])
for i, df, coord, name in zip(range(3), obs_dfs, station_coords, station_names): 
    #Plot observations
    obs = df['LWU [W/m**2]'].loc["2007-04":"2010-03"]
    obs.plot(label='Observed', color='black', ax=axs[i], linewidth=2)
    #Plot uncertainty
    axs[i].fill_between(np.arange(36), obs - obs*(0.02), obs + obs*(0.02), color='black', alpha=0.1)
    #Plot model
    axs[i].plot(np.arange(36),(M92['FLNS']+M92['FLDS']).sel(lat=coord[0], lon=coord[1], method='nearest'), label='M92', color='tab:blue',linewidth=2 )
    axs[i].plot(np.arange(36),(A21['FLNS']+A21['FLDS']).sel(lat=coord[0], lon=coord[1], method='nearest'), label='A21', color='tab:orange',linewidth=2 )

    axs[i].set_title(name)
    if i==0:
        axs[i].legend()
    axs[i].set_ylabel(r'W/m$^2$')
    axs[i].set_xlabel('')
    axs[i].set_xticks(np.arange(1, 38, 4))
    axs[i].grid(alpha=0.5)
axs[2].tick_params(axis='x', labelrotation=45)
fig.savefig(wpath+'pdf/stations_lwu.pdf', bbox_inches="tight")
fig.savefig(wpath+'png/stations_lwu.png', bbox_inches="tight")

plt.clf()

# Net Longwave Flux

fig, axs = plt.subplots(3, 1, sharex=True, figsize=[13,12])
for i, df, coord, name in zip(range(3), obs_dfs, station_coords, station_names): 
    #Plot observations
    obs = (df['LWD [W/m**2]'] - df['LWU [W/m**2]']).loc["2007-04":"2010-03"]
    obs.plot(label='Observed', color='black', ax=axs[i], linewidth=2)
    #Plot uncertainty
    axs[i].fill_between(np.arange(36), obs - obs*(0.095), obs + obs*(0.095), color='black', alpha=0.1)
    #Plot model
    axs[i].plot(np.arange(36),(-M92['FLNS']).sel(lat=coord[0], lon=coord[1], method='nearest'), label='M92', color='tab:blue',linewidth=2 )
    axs[i].plot(np.arange(36),(-A21['FLNS']).sel(lat=coord[0], lon=coord[1], method='nearest'), label='A21', color='tab:orange', linewidth=2 )

    axs[i].set_title(name)
    if i==0:
        axs[i].legend()
    axs[i].set_ylabel(r'W/m$^2$')
    axs[i].set_xlabel('')
    axs[i].set_xticks(np.arange(1, 38, 4))
    axs[i].grid(alpha=0.5)
    axs[i].tick_params('x',labelrotation=45)

fig.savefig(wpath+'pdf/stations_lw_net.pdf', bbox_inches="tight")
fig.savefig(wpath+'png/stations_lw_net.png', bbox_inches="tight")

plt.clf()


# Downwelling Shortwave Flux at Surface

fig, axs = plt.subplots(3, 1, sharex=True, figsize=[13,12])
for i, df, coord, name in zip(range(3), obs_dfs, station_coords, station_names): 
    #Plot observations
    obs = df['SWD [W/m**2]'].loc["2007-04":"2010-03"]
    obs.plot(label='Observed', color='black', ax=axs[i], linewidth=2)
    #Plot uncertainty
    axs[i].fill_between(np.arange(36), obs - obs*(0.02), obs + obs*(0.02), color='black', alpha=0.1)
    #Plot model
    axs[i].plot(np.arange(36),M92['FSDS'].sel(lat=coord[0], lon=coord[1], method='nearest'), label='M92', color='tab:blue',linewidth=2 )
    axs[i].plot(np.arange(36),A21['FSDS'].sel(lat=coord[0], lon=coord[1], method='nearest'), label='A21', color='tab:orange',linewidth=2 )

    axs[i].set_title(name)
    if i==0:
        axs[i].legend()
    axs[i].set_ylabel(r'W/m$^2$')
    axs[i].set_xticks(np.arange(1, 38, 4))
    axs[i].set_xlabel('')
    axs[i].grid(alpha=0.5)
axs[2].tick_params(axis='x', labelrotation=45)
fig.savefig(wpath+'pdf/stations_swd.pdf', bbox_inches="tight")
fig.savefig(wpath+'png/stations_swd.png', bbox_inches="tight")

plt.clf()

# Upwelling Shortwave Flux at Surface

fig, axs = plt.subplots(3, 1, sharex=True, figsize=[13,12])
for i, df, coord, name in zip(range(3), obs_dfs, station_coords, station_names): 
    #Plot observations
    obs = df['SWU [W/m**2]'].loc["2007-04":"2010-03"]
    obs.plot(label='Observed', color='black', ax=axs[i], linewidth=2)
    #Plot uncertainty
    axs[i].fill_between(np.arange(36), obs - obs*(0.02), obs + obs*(0.02), color='black', alpha=0.1)
    #Plot model
    axs[i].plot(np.arange(36),-(M92['FSNS']-M92['FSDS']).sel(lat=coord[0], lon=coord[1], method='nearest'), label='M92', color='tab:blue',linewidth=2 )
    axs[i].plot(np.arange(36),-(A21['FSNS']-A21['FSDS']).sel(lat=coord[0], lon=coord[1], method='nearest'), label='A21', color='tab:orange',linewidth=2 )

    axs[i].set_title(name)
    if i==0:
        axs[i].legend()
    axs[i].set_ylabel(r'W/m$^2$')
    axs[i].set_xlabel('')
    axs[i].set_xticks(np.arange(1, 38, 4))
    axs[i].grid(alpha=0.5)
axs[2].tick_params(axis='x', labelrotation=45)
fig.savefig(wpath+'pdf/stations_swu.pdf', bbox_inches="tight")
fig.savefig(wpath+'png/stations_swu.png', bbox_inches="tight")

plt.clf()

# Surface Temperature

fig, axs = plt.subplots(3, 1, sharex=True, figsize=[13,12])
for i, df, coord, name in zip(range(3), obs_dfs, station_coords, station_names): 
    #Plot observations
    obs = df['T2 [°C]'].loc["2007-04":"2010-03"]
    obs.plot(label='Observed', color='black', ax=axs[i], linewidth=2)
    #Plot uncertainty
    axs[i].fill_between(np.arange(36), obs - obs*(0.02), obs + obs*(0.02), color='black', alpha=0.1)
    #Plot model
    axs[i].plot(np.arange(36),M92['TREFHT'].sel(lat=coord[0], lon=coord[1], method='nearest'), label='M92', color='tab:blue',linewidth=2 )
    axs[i].plot(np.arange(36),A21['TREFHT'].sel(lat=coord[0], lon=coord[1], method='nearest'), label='A21', color='tab:orange',linewidth=2 )

    axs[i].set_title(name)
    if i==0:
        axs[i].legend()
    axs[i].set_ylabel(r'$^{\circ}$C')
    axs[i].set_xlabel('')
    axs[i].set_xticks(np.arange(1, 38, 4))
    axs[i].grid(alpha=0.5)
#axs[2].set_xticks(np.arange(1, 38, 4))#, ['2007-05', '2007-09', '2008-01', '2008-05', '2008-09', '2009-01', '2009-05', '2009-09', '2010-01', '2010-05'])
axs[2].tick_params(axis='x', labelrotation=45)
fig.savefig(wpath+'pdf/stations_T2.pdf', bbox_inches="tight")
fig.savefig(wpath+'png/stations_T2.png', bbox_inches="tight")

plt.clf()


#-------------------------------------
# Plots for individual stations
#-------------------------------------
for df, coord, name in zip(obs_dfs, station_coords, station_names):
    
    # Net Shortwave Flux

    plt.figure(figsize=[13,4])
    #plt.suptitle(name + ": Net Incoming Shortwave Flux at Surface")
    #Plot observations
    obs = (df['SWD [W/m**2]'] - df['SWU [W/m**2]']).loc["2007-04":"2010-03"]
    obs.plot(label='Observed', color='black')
    #Plot uncertainty
    plt.plot(np.arange(36), obs + obs*(0.055), color='black', linestyle='--')
    plt.plot(np.arange(36), obs - obs*(0.055), color='black', linestyle='--')
    #Plot model
    plt.plot(np.arange(36),M92['FSNS'].sel(lat=coord[0], lon=coord[1], method='nearest'), label='M92', color='tab:blue' )
    plt.plot(np.arange(36),A21['FSNS'].sel(lat=coord[0], lon=coord[1], method='nearest'), label='A21', color='tab:orange' )
    
    plt.legend()
    plt.xticks(rotation=45, ha='right')
    plt.ylabel(r'W/m$^2$')
    plt.xlabel('')
    plt.grid(alpha=0.5)

    plt.savefig(wpath+'pdf/'+name+'_sw_net.pdf', bbox_inches="tight")
    plt.savefig(wpath+'png/'+name+'_sw_net.png', bbox_inches="tight")

    plt.clf()

    # Net Longwave Flux

    plt.figure(figsize=[13,4])
    plt.suptitle(name + ": Net Outgoing Longwave Flux at Surface")
    #Plot observations
    obs = (df['LWU [W/m**2]'] - df['LWD [W/m**2]']).loc["2007-04":"2010-03"]
    obs.plot(label='Observed', color='black')
    #Plot uncertainty
    plt.plot(np.arange(36), obs + obs*(0.04), color='black', linestyle='--')
    plt.plot(np.arange(36), obs - obs*(0.04), color='black', linestyle='--')
    #Plot model
    plt.plot(np.arange(36),M92['FLNS'].sel(lat=coord[0], lon=coord[1], method='nearest'), label='M92', color='tab:blue' )
    plt.plot(np.arange(36),A21['FLNS'].sel(lat=coord[0], lon=coord[1], method='nearest'), label='A21', color='tab:orange' )
    
    plt.legend()
    plt.xticks(rotation=45, ha='right')
    plt.ylabel(r'W/m$^2$')
    plt.xlabel('')
    plt.grid(alpha=0.5)

    plt.savefig(wpath+'pdf/'+name+'_lw_net.pdf', bbox_inches="tight")
    plt.savefig(wpath+'png/'+name+'_lw_net.png', bbox_inches="tight")

    plt.clf()

    # Downwelling Longwave Flux

    plt.figure(figsize=[13,4])
    plt.suptitle(name + ": Downwelling Longwave Flux at Surface")
    #Plot observations
    obs = df['LWD [W/m**2]'].loc["2007-04":"2010-03"]
    obs.plot(label='Observed', color='black')
    #Plot uncertainty
    plt.plot(np.arange(36), obs + obs*(0.04), color='black', linestyle='--')
    plt.plot(np.arange(36), obs - obs*(0.04), color='black', linestyle='--')
    #Plot model
    plt.plot(np.arange(36),M92['FLDS'].sel(lat=coord[0], lon=coord[1], method='nearest'), label='M92', color='tab:blue' )
    plt.plot(np.arange(36),A21['FLDS'].sel(lat=coord[0], lon=coord[1], method='nearest'), label='A21', color='tab:orange' )
    
    plt.legend()
    plt.xticks(rotation=45, ha='right')
    plt.ylabel(r'W/m$^2$')
    plt.xlabel('')
    plt.grid(alpha=0.5)

    plt.savefig(wpath+'pdf/'+name+'_lwd.pdf', bbox_inches="tight")
    plt.savefig(wpath+'png/'+name+'_lwd.png', bbox_inches="tight")

    plt.clf()
    
    # Net Surface Flux

    plt.figure(figsize=[13,4])
    #plt.suptitle(name + ": Net Incoming Flux at Surface")
    #Plot observations
    obs = (df['SWD [W/m**2]'] + df['LWD [W/m**2]'] - df['SWU [W/m**2]'] - df['LWU [W/m**2]']).loc["2007-04":"2010-03"]
    obs.plot(label='Observed', color='black')
    #Plot uncertainty
    plt.fill_between(np.arange(36), obs - obs*(0.095), obs + obs*(0.095), color='black', alpha=0.1)
    #Plot model
    plt.plot(np.arange(36),(M92['FSNS'] - M92['FLNS']).sel(lat=coord[0], lon=coord[1], method='nearest'), label='M92', color='tab:blue' )
    #plt.plot(np.arange(36),(A21['FSNS'] - A21['FLNS']).sel(lat=coord[0], lon=coord[1], method='nearest'), label='A21', color='tab:orange' )
    
    #plt.legend()
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.25),ncol=3)
    plt.xticks(rotation=45, ha='right')
    plt.ylabel(r'W/m$^2$')
    plt.xlabel('')
    plt.ylim([-70, 130])
    plt.grid(alpha=0.5)

    #plt.savefig(wpath+'pdf/'+name+'_all_net.pdf', bbox_inches="tight")
    plt.savefig(wpath+'png/'+name+'_all_net_M92.png', bbox_inches="tight")

    plt.clf()
    """


