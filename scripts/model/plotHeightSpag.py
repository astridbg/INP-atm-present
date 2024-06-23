import xarray as xr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter
from cmcrameri import cm
# Set font style to match latex document----------
#plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size':18})
# ------------------------------------------------
import datetime
import cartopy.crs as ccrs
from functions import *

rpath="/projects/NS9600K/astridbg/data/model/noresm_postprocessed/"
wpath="/projects/NS9600K/astridbg/INP-atm-present/figures/model/heightspag/"

# Default cases----------------
#case1 = "def_20210126"; case1nm = "CAM6"
case1 = "meyers92_20220210"; case1nm = "M92"
# Modified cases---------------
#case2 = "meyers92_20220210"; case2nm = "CAM5"
case2 = "andenes21_20220222"; case2nm = "A21"
#------------------------------	
date1 = "2007-04-15_2010-03-15"
date2 = "2007-04-15_2010-03-15"

#------------------------------
# Areas to analyse
#------------------------------

NYA = [78.9227, 11.9273] # Ny-Ålesund
ALE = [82.5, 297.6] # Alert
BAR = [71.17, 203.55] # Barrow
npole = [[0,360],[85,90]] # North Pole

#------------------------------
# Add open ocean mask
#------------------------------
ocean_mask = False
if ocean_mask:
    ds_ocn = xr.open_dataset(rpath+"OCNFRAC"+"_"+case2+"_"+date2+".nc")
    open_sea = ds_ocn > 0.85
    open_sea = open_sea.mean("time")
    open_sea = open_sea >= 0.5


#------------------------------
# Pairs of two-dimensional fields
#------------------------------
# variables = [[3D-var, 2D-var]]
var_level = 850
#variables = [["CLOUD","CLDTOT"], ["CLOUD", "CLDLOW"]]
variables = [["CLOUD","CLDTOT"]]

#------------------------------
# Plotting monthly mean averages
#------------------------------
for pair in variables:
    
    ds1 = xr.open_mfdataset([rpath+pair[0]+"_"+case1+"_"+date1+".nc", rpath+pair[1]+"_"+case1+"_"+date1+".nc"])
    ds2 = xr.open_mfdataset([rpath+pair[0]+"_"+case2+"_"+date2+".nc", rpath+pair[1]+"_"+case2+"_"+date2+".nc"])
    
    # Get start and end date of period
    date_start = str(ds1.time[0].values).split(" ")[0]
    date_end = str(ds1.time[-1].values).split(" ")[0]

    ### PROFILE

    ds1_ann = ds1.mean("time")
    ds2_ann = ds2.mean("time")

    # Select level
    ds1_level = ds1_ann.sel(lev=var_level, method="nearest")
    ds2_level = ds2_ann.sel(lev=var_level, method="nearest")
    lev_name = str(int(ds1_level.lev.values))

    # Get relative difference between cases time averaged over the whole period
    var = pair[0]
    diff = ds2_level[var]-ds1_level[var]
    reldiff = diff/ds1_level[var].where(ds1_level[var]!=0)*100*np.sign(ds1_level[var])
    diff_used = reldiff
    lev_extent = round(max(abs(np.nanmin(diff_used.sel(lat=slice(66.5,90)).values)), 
                            abs(np.nanmax(diff_used.sel(lat=slice(66.5,90)).values))))
    if lev_extent < 0.004:
            lev_extent = 0.004
    
    levels = np.linspace(-lev_extent,lev_extent,25)

    # Make horizontal average for the Arctic
    ds1_arct_height = computeWeightedMean(ds1_ann[var].sel(lat=slice(66.5,90),lev=slice(350, 1000)))
    ds2_arct_height = computeWeightedMean(ds2_ann[var].sel(lat=slice(66.5,90),lev=slice(350, 1000)))

    height_levels = ds1.lev.sel(lev=slice(350, 1000)).values

    # SPAGETTI PLOT

    # Get monthly mean
    ds1m = ds1.groupby("time.month").mean("time")
    ds2m = ds2.groupby("time.month").mean("time")
    
    # Make time array
    months = []
    for i in range(len(ds1m.month.values)):
        datetime_object = datetime.datetime.strptime(str(ds1m.month.values[i]), "%m")
        months.append(datetime_object.strftime("%b"))
    
    # Get spatial average over Arctic
    ds1_arct = computeWeightedMean(ds1m.sel(lat=slice(66.5,90)))
    ds2_arct = computeWeightedMean(ds2m.sel(lat=slice(66.5,90)))
    print("Arctic")
    print(pair[0])
    print((ds2_arct[pair[0]]-ds1_arct[pair[0]]).mean("month").values)
    print(pair[1])
    print((ds2_arct[pair[1]]-ds1_arct[pair[1]]).mean("month").values)

    # Get values over Ny-Ålesund
    ds1_NYA = ds1m.sel(lat=NYA[0], lon=NYA[1], method="nearest")
    ds2_NYA = ds2m.sel(lat=NYA[0], lon=NYA[1], method="nearest")

    # Get values over Alert
    ds1_ALE = ds1m.sel(lat=ALE[0], lon=ALE[1], method="nearest")
    ds2_ALE = ds2m.sel(lat=ALE[0], lon=ALE[1], method="nearest")

    # Get values over Barrow
    ds1_BAR = ds1m.sel(lat=BAR[0], lon=BAR[1], method="nearest")
    ds2_BAR = ds2m.sel(lat=BAR[0], lon=BAR[1], method="nearest")
        
    # Get spatial average over North Pole
    ds1_npol = computeWeightedMean(ds1m.sel(lon=slice(npole[0][0],npole[0][1]),
                                lat=slice(npole[1][0],npole[1][1])))
    ds2_npol = computeWeightedMean(ds2m.sel(lon=slice(npole[0][0],npole[0][1]),
                                lat=slice(npole[1][0],npole[1][1])))
    

   # fig,axs = plt.subplots(ncols=2,nrows=2, gridspec_kw={'width_ratios': [3, 1], 'height_ratios': [2, 1]}, figsize=[13,9],dpi=300)
    fig,axs = plt.subplot_mosaic([['ax1', 'ax1', 'ax2', 'ax2'], ['ax3', 'ax3', 'ax3', 'ax4']], gridspec_kw={'height_ratios': [2, 1]},figsize=[13,9], dpi=300)
    
    left_label = ["(a)", "(c)"]
    right_label = ["(b)", "(d)"]


    axs['ax1'].plot(ds1_arct_height, height_levels, label=case1nm, color="tab:blue",linestyle="--",linewidth=2)
    axs['ax1'].plot(ds2_arct_height, height_levels, label=case2nm, color="tab:orange",linewidth=2)
    axs['ax1'].hlines(ds1_level.lev.values, axs['ax1'].get_xlim()[0],axs['ax1'].get_xlim()[1], color="black",linestyle=":")
    if var == "NIMEY" or var == "AWNICC" or var == "AWNI" or var == "NIMEYCC" or var == "NIMEYCLD" or var == "AWNICLD" or var == "AWNINONIMEY":
        axs['ax1'].xscale("log")
    axs['ax1'].set_ylabel("hPa")
    axs['ax1'].set_xlabel(ds1[var].units)
    axs['ax1'].legend(loc="upper left",frameon=True,title="Arctic averages",fancybox=False)
    axs['ax1'].grid(alpha=0.5)
    axs['ax1'].invert_yaxis()
    axs['ax1'].annotate(right_label[0],fontsize=20,
        xy=(0, 1), xycoords='axes fraction',
        xytext=(-30, 30), textcoords='offset points',
        ha='left', va='top')

    axs['ax2'].remove()
    ax2 = fig.add_subplot(2, 2, 2, projection=ccrs.Orthographic(0, 90))
    polarCentral_set_latlim([66.5,90], ax2)
    map = diff_used.plot.pcolormesh(ax=ax2, transform=ccrs.PlateCarree(), 
                                            cmap="coolwarm",#cmap=plt.cm.get_cmap('Blues').reversed(), #cmap="coolwarm"
                                            levels=levels,
                                            add_colorbar=False)
    if ocean_mask:
        ax2.contourf(open_sea.lon, open_sea.lat, open_sea["OCNFRAC"], transform=ccrs.PlateCarree(), colors='none',hatches=['..'],levels=[.5, 1.5])
    ax2.coastlines()
    ax2.set_title("Level = "+lev_name+" hPa")  
    cb_ax = fig.add_axes([0.5, 0.44, 0.4, 0.04])

    cbar = plt.colorbar(map, cax=cb_ax, spacing = 'uniform', extend='both', orientation='horizontal', fraction=0.046, pad=0.06)
    #cbar.ax.tick_params(labelsize=18)
    cbar.ax.set_xlabel("%")#, fontsize=23)
    if lev_extent >= 4:
        cbar.ax.xaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}')) # No decimal places        
    elif 0.4 <= lev_extent < 4:
        cbar.ax.xaxis.set_major_formatter(StrMethodFormatter('{x:,.1f}')) # One decimal place
    elif 0.04 <= lev_extent < 0.4:
        cbar.ax.xaxis.set_major_formatter(StrMethodFormatter('{x:,.2f}')) # Two decimal places     
    elif 0.004 <= lev_extent < 0.04:
        cbar.ax.xaxis.set_major_formatter(StrMethodFormatter('{x:,.3f}')) # Three decimal places
    ax2.annotate(left_label[0],fontsize=20,
        xy=(0, 1), xycoords='axes fraction',
        xytext=(-30, 30), textcoords='offset points',
        ha='left', va='top')
    


    i=1
    axs['ax3'].plot(months, ds2_arct[pair[1]]-ds1_arct[pair[1]], color=cm.romaO.colors[0], label="Arctic")
    axs['ax3'].plot(months, ds2_NYA[pair[1]]-ds1_NYA[pair[1]], color=cm.romaO.colors[49], label="Ny-Ålesund")
    axs['ax3'].plot(months, ds2_ALE[pair[1]]-ds1_ALE[pair[1]], color=cm.romaO.colors[99], label="Alert")
    axs['ax3'].plot(months, ds2_BAR[pair[1]]-ds1_BAR[pair[1]], color=cm.romaO.colors[149], label="Utqiagvik")
    axs['ax3'].plot(months, ds2_npol[pair[1]]-ds1_npol[pair[1]], color=cm.romaO.colors[199], label="North Pole")
    axs['ax3'].set_ylabel(r"$\Delta$"+ds1[pair[1]].units)
    axs['ax3'].tick_params(axis="x",labelsize=20)
    axs['ax3'].grid(alpha=0.5)
    axs['ax3'].annotate(left_label[i],fontsize=20,
            xy=(0, 1), xycoords='axes fraction',
            xytext=(-30, 30), textcoords='offset points',
            ha='left', va='top')

    if pair[i] == "TREFHT":
        # Get average absolute change
        
        abs_arct = ds2_arct[pair[i]]-ds1_arct[pair[i]]
        abs_nya = ds2_NYA[pair[i]]-ds1_NYA[pair[i]]
        abs_ale = ds2_ALE[pair[i]]-ds1_ALE[pair[i]]
        abs_bar = ds2_BAR[pair[i]]-ds1_BAR[pair[i]]
        abs_npol = ds2_npol[pair[i]]-ds1_npol[pair[i]]
    
    
        change_all = pd.DataFrame({"Arctic":abs_arct,"Ny-Ålesund":abs_nya,"Alert":abs_ale,"Utqiagvik":abs_bar,"North Pole":abs_npol})
        bplot= axs['ax4'].boxplot(change_all,patch_artist=True,medianprops={"color":"black"})
        axs['ax4'].set_ylabel(r"$\Delta$"+ds1[pair[i]].units)
    
    elif pair[i] == "TGCLDLWP" or pair[i] == "CLDLWEM":

        # Get average relative change
        # Shrink y axis due to extreme values

        rel_arct = ((ds2_arct[pair[i]]-ds1_arct[pair[i]])/ds1_arct[pair[i]].where(ds1_arct[pair[i]]!=0)).fillna(0)*100*np.sign(ds1_arct[pair[i]])
        rel_nya = ((ds2_NYA[pair[i]]-ds1_NYA[pair[i]])/ds1_NYA[pair[i]].where(ds1_NYA[pair[i]]!=0)).fillna(0)*100*np.sign(ds1_NYA[pair[i]])
        rel_ale = ((ds2_ALE[pair[i]]-ds1_ALE[pair[i]])/ds1_ALE[pair[i]].where(ds1_ALE[pair[i]]!=0)).fillna(0)*100*np.sign(ds1_ALE[pair[i]])
        rel_bar = ((ds2_BAR[pair[i]]-ds1_BAR[pair[i]])/ds1_BAR[pair[i]].where(ds1_BAR[pair[i]]!=0)).fillna(0)*100*np.sign(ds1_BAR[pair[i]])
        rel_npol = ((ds2_npol[pair[i]]-ds1_npol[pair[i]])/ds1_npol[pair[i]].where(ds1_npol[pair[i]]!=0)).fillna(0)*100*np.sign(ds1_npol[pair[i]])
    
    
        rel_all = pd.DataFrame({"Arctic":rel_arct,"Ny-Ålesund":rel_nya,"Alert":rel_ale,"Utqiagvik":rel_bar,"North Pole":rel_npol})
        bplot=axs['ax4'].boxplot(rel_all,patch_artist=True,medianprops={"color":"black"},showfliers=False)
        axs['ax4'].set_ylabel("% change")
        axs['ax4'].set_yscale("log")


    else:
        # Get average relative change
    
        rel_arct = ((ds2_arct[pair[i]]-ds1_arct[pair[i]])/ds1_arct[pair[i]].where(ds1_arct[pair[i]]!=0)).fillna(0)*100*np.sign(ds1_arct[pair[i]])
        rel_nya = ((ds2_NYA[pair[i]]-ds1_NYA[pair[i]])/ds1_NYA[pair[i]].where(ds1_NYA[pair[i]]!=0)).fillna(0)*100*np.sign(ds1_NYA[pair[i]])
        rel_ale = ((ds2_ALE[pair[i]]-ds1_ALE[pair[i]])/ds1_ALE[pair[i]].where(ds1_ALE[pair[i]]!=0)).fillna(0)*100*np.sign(ds1_ALE[pair[i]])
        rel_bar = ((ds2_BAR[pair[i]]-ds1_BAR[pair[i]])/ds1_BAR[pair[i]].where(ds1_BAR[pair[i]]!=0)).fillna(0)*100*np.sign(ds1_BAR[pair[i]])
        rel_npol = ((ds2_npol[pair[i]]-ds1_npol[pair[i]])/ds1_npol[pair[i]].where(ds1_npol[pair[i]]!=0)).fillna(0)*100*np.sign(ds1_npol[pair[i]])

    
        rel_all = pd.DataFrame({"Arctic":rel_arct,"Ny-Ålesund":rel_nya,"Alert":rel_ale,"Utqiagvik":rel_bar,"North Pole":rel_npol})
        bplot=axs['ax4'].boxplot(rel_all,patch_artist=True,medianprops={"color":"black"},showfliers=False)
        axs['ax4'].set_ylabel("% change")
    
    colors=cm.romaO.colors[::50][:5]
    for patch, color in zip(bplot['boxes'], colors):
        patch.set_facecolor(color)
    axs['ax4'].grid(alpha=0.5)
    axs['ax4'].set_xticklabels([])
    axs['ax4'].annotate(right_label[i],fontsize=20,
            xy=(0, 1), xycoords='axes fraction',
            xytext=(-30, 30), textcoords='offset points',
            ha='left', va='top')

    # Shrink current axis's height by 15% on the bottom
    box = axs['ax3'].get_position()
    axs['ax3'].set_position([box.x0, box.y0 + box.height * 0.15,
                box.width, box.height * 0.85])
    box = axs['ax4'].get_position()
    axs['ax4'].set_position([box.x0, box.y0 + box.height * 0.15,
                box.width, box.height * 0.85])

    # Put a legend below current axis
    axs['ax3'].legend(loc='upper center', bbox_to_anchor=(0.75, -0.12), ncol=5, columnspacing=0.5, handlelength=1,handletextpad=0.4)
	
    plt.grid(alpha=0.5)
    plt.savefig(wpath+"pdf/"+pair[0]+"_"+pair[1]+"_avg_"+case1+"_"+case2+".pdf",bbox_inches="tight")
    plt.savefig(wpath+"png/"+pair[0]+"_"+pair[1]+"_avg_"+case1+"_"+case2+".png",bbox_inches="tight")
    plt.clf()
