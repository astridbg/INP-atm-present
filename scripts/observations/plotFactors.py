# Author: Tim Carlsen
# Modified by: Astrid Bragstad Gjelsvik

import numpy as np
import pandas as pd
import xarray as xr
import glob
import datetime
import matplotlib.pyplot as plt
import matplotlib as mpl
import functions
mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size':20})


path_islas = "/projects/NS9600K/data/islas/"
path_cor = "/projects/NS9600K/astridbg/data/observations/Coriolis_postprocessed/"
path_station = "/projects/NS9600K/astridbg/data/observations/SN87110/"
path_temp = "/projects/NS9600K/astridbg/data/observations/Inlet_Temp/"
path_aero = "/projects/NS9600K/data/islas/MetOne/GT-526S/MetOne-20210314131343/2021/03/"
path_sea = "/projects/NS9600K/astridbg/data/observations/Sea/"
wpath = "/projects/NS9600K/astridbg/master/figures/observations/"

# -----------------------------
# Coriolis data
# -----------------------------

# Read in Coriolis INP concentrations

nucleiT = pd.read_csv(path_cor+"Coriolis_nucleiT_cal.csv",index_col=0)
nCor = len(nucleiT.iloc[0,:])

# Read in Coriolis log file

df_cor = pd.read_csv(path_islas+"coriolis_log_all.csv", skiprows = 1)

t_start = pd.DataFrame(df_cor['Date'] + ' ' + df_cor['Start (UTC)'], columns = ['t_start'])
t_start = t_start.set_index('t_start')
t_start.index = pd.to_datetime(t_start.index, format = '%d-%m-%Y %H:%M')

t_end = pd.DataFrame(df_cor['Date'] + ' ' + df_cor['End (UTC)'], columns = ['t_end'])
t_end = t_end.set_index('t_end')
t_end.index = pd.to_datetime(t_end.index, format = '%d-%m-%Y %H:%M')


# Set time index to middle to sampling period

t_diff = datetime.timedelta(minutes=20)
t_middle = t_start.index + t_diff

df_frzT = nucleiT.iloc[1:-1].transpose()
df_frzT["mean_time"] = t_middle
df_frzT = df_frzT.set_index("mean_time")

# Get temperature at which fifty percent for all wells had frozen

index_50 = int(94/2)
t50 = df_frzT.iloc[:,index_50]


# -----------------------------
# Aerosol data
# -----------------------------

# Get pressure data

ds_pres = xr.open_dataset(path_station+"air_pressure_at_sea_level_202103.nc")
df_pres = ds_pres.to_dataframe()

# Get inlet temperature data

files = sorted(glob.glob(path_temp+"*.txt"))

count = 0
for f in files:

    df1 = pd.read_csv(f, delimiter = ',', parse_dates = ['Time'], index_col = ['Time'])
    if count == 0:
        df_temp = df1
    else:
        df_temp = pd.concat([df_temp, df1], axis = 0)

    count += 1

# Get OPC data

file_list = sorted(glob.glob(path_aero+"*/*.CSV"))
count = 0
for f in file_list:
    df1 = pd.read_csv(f, delimiter = ';', decimal = ',', parse_dates = ['Time'], index_col = ['Time'], keep_date_col = True)
    if count == 0:
        df_opc = df1
    else:
        df_opc = pd.concat([df_opc, df1])
        
    count += 1

# Convert OPC data to standard litre

p_std = 1013.25
T_std = 273.15

df_opc["Count2 (/std_L)"] = df_opc["Count2 (/L)"]

for i in range(len(df_opc["Count2 (/L)"])):
   p = df_pres["air_pressure_at_sea_level"].iloc[df_pres.index.get_loc(df_opc.index[i], method="nearest")]
   T = df_temp["Temperature(C)"].iloc[df_temp.index.get_loc(df_opc.index[i], method="nearest")]
   df_opc["Count2 (/std_L)"][i] = df_opc["Count2 (/L)"][i] * p_std/p * (273.15 + T)/T_std


# -----------------------------
# Wind data
# -----------------------------

# Get wind data

ds_wisp = xr.open_dataset(path_station+"wind_speed_202103.nc")
df_wisp = ds_wisp.to_dataframe()

ds_widir = xr.open_dataset(path_station+"wind_from_direction_202103.nc")
df_widir = ds_widir.to_dataframe()

# -----------------------------
# Sea data
# -----------------------------

# Read in sea sample log files

df_shore = pd.read_csv(path_sea+"sea_log_all.csv", skiprows = 1)
#df_open = pd.read_csv(path_sea+"bsea_log_all.csv", skiprows = 1)

t_shore = pd.DataFrame(df_shore['Date'] + ' ' + df_shore['Collection Time'], columns = ['t_collection'])
t_shore = t_shore.set_index('t_collection')
t_shore.index = pd.to_datetime(t_shore.index, format = '%d-%m-%Y %H:%M')

#t_open = pd.DataFrame(df_open['Date'] + ' ' + df_open['Collection Time'], columns = ['t_collection'])
# Fill in assumed time of collection for BSEA003
#new_time = str(t_open.iloc[2]).split(" ")[4] + " " + "14:00"
#t_open.iloc[2] = new_time
#t_open = t_open.set_index('t_collection')
#t_open.index = pd.to_datetime(t_open.index, format = '%d-%m-%Y %H:%M')

# Get T50 and TOC data

T_shore = pd.read_csv(path_sea+"Sea_nucleiT_cal.csv",index_col=0)
T50_shore = T_shore.iloc[48,:]
t_shore["T 50 % frozen (degC)"] = np.array(T50_shore)

TOC_shore = pd.read_csv(path_sea+"Sea_TOC.csv", header=None)
t_shore["TOC (mg/L)"] = np.array(TOC_shore)

#T_open = pd.read_csv(path_sea+"BSea_nucleiT_cal.csv",index_col=0)
#T50_open = T_open.iloc[48,:]
#t_open["T 50 % frozen (degC)"] = np.array(T50_open)

#TOC_open = pd.read_csv(path_sea+"BSea_TOC.csv", header=None)
#t_open["TOC (mg/L)"] = np.array(TOC_open)

# Combine sea sample data
df_sea = t_shore

#df_sea = df_sea.append(t_open)
#print(df_sea)
#df_sea = df_sea.sort_index()
#print(df_sea)

# Interpolate TOC data

df_sea_interpol = df_sea.copy()
for i in t_start.index.values[7:50]:
    df_sea_interpol.loc[i,:] = np.nan
df_sea_interpol = df_sea_interpol.sort_index().interpolate(method='time')


# -----------------------------
# Averages over Coriolis period
# -----------------------------

opc_all = []
wisp_all = []
widir_all = []
toc_all = []

# For full-length datasets
i = 1
for cor in t_start.index:

    print(cor.date(),", Coriolis sample: ",i)

    time_opc = df_opc[str(cor.date())].index.hour
    index_cor_opc = np.where(np.logical_or(time_opc == cor.hour, time_opc == int(cor.hour + cor.minute/60. + 40./60.)))
    opc = np.nanmean(np.array(df_opc['Count2 (/std_L)'][str(cor.date())])[index_cor_opc])

    time_wisp = df_wisp[str(cor.date())].index.hour + df_wisp[str(cor.date())].index.minute/60. + df_wisp[str(cor.date())].index.second/3600.
    index_cor_wisp = np.where((time_wisp >= cor.hour+cor.minute/60.) & (time_wisp <= cor.hour + cor.minute/60. + 40./60.))
    wisp = np.nanmean(np.array(df_wisp['wind_speed'][str(cor.date())])[index_cor_wisp])

    time_widir = df_widir[str(cor.date())].index.hour + df_widir[str(cor.date())].index.minute/60. + df_widir[str(cor.date())].index.second/3600.
    index_cor_widir = np.where((time_widir >= cor.hour+cor.minute/60.) & (time_widir <= cor.hour + cor.minute/60. + 40./60.))
    widir = np.nanmean(np.array(df_widir['wind_from_direction'][str(cor.date())])[index_cor_widir])

    opc_all = np.append(opc_all, opc)
    wisp_all = np.append(wisp_all, wisp)
    widir_all = np.append(widir_all, widir)
    
    i += 1

# For interpolated dataset
i = 7
for cor in t_start.index[7:50]:

    print(cor.date(),", Coriolis sample: ",i)

    time_toc = df_sea_interpol[str(cor.date())].index.hour + df_sea_interpol[str(cor.date())].index.minute/60. + df_sea_interpol[str(cor.date())].index.second/3600.
    index_cor_toc = np.where((time_toc >= cor.hour+cor.minute/60.) & (time_toc <= cor.hour + cor.minute/60. + 40./60.))
    toc = np.nanmean(np.array(df_sea_interpol['TOC (mg/L)'][str(cor.date())])[index_cor_toc])

    toc_all = np.append(toc_all, toc)

    i += 1

# Bin average wind

widir_all_binned = []

for i in range(len(wisp_all)):
    widir = widir_all[i]
    if 300 <= widir <= 360 or widir <= 45:
        widir_all_binned = np.append(widir_all_binned, 4)
    if 210 <= widir < 300:
        widir_all_binned = np.append(widir_all_binned, 3)
    if 175 <= widir < 210:
        widir_all_binned = np.append(widir_all_binned, 2)
    if 45 < widir < 175:
        widir_all_binned = np.append(widir_all_binned, 1)

df_widir_binned = pd.DataFrame({"Dir":widir_all_binned,"FrzT50":t50})
df_widir_binned["Easterly"] = df_widir_binned["FrzT50"].where(df_widir_binned["Dir"] == 1)
df_widir_binned["Southerly"] = df_widir_binned["FrzT50"].where(df_widir_binned["Dir"] == 2)
df_widir_binned["Southwesterly"] = df_widir_binned["FrzT50"].where(df_widir_binned["Dir"] == 3)
df_widir_binned["Arctic"] = df_widir_binned["FrzT50"].where(df_widir_binned["Dir"] == 4)

df_widir_binned = df_widir_binned.drop(["Dir","FrzT50"],axis=1)

# -----------------------------
# Plotting routine
# -----------------------------

fig, axs = plt.subplots(5, 2, gridspec_kw={'width_ratios': [4, 1]}, figsize=(12,15), dpi=300, constrained_layout=True)

df_frzT.T.boxplot(
        positions=mpl.dates.date2num(df_frzT.index),
        widths=0.1, flierprops = dict(marker='.',markeredgecolor="tab:blue"), ax=axs[0,0])
locator = mpl.dates.AutoDateLocator(minticks=10, maxticks=15)
axs[0,0].xaxis.set_major_locator(locator)
axs[0,0].set_xticklabels([])
xlims = mpl.dates.num2date(axs[0,0].get_xlim())
xticks = mpl.dates.num2date(axs[0,0].get_xticks())
axs[0,0].set_ylabel("$^{\circ}$C")
axs[0,0].set_title("Freezing temperatures of INPs")
axs[0,0].annotate("(a)",fontsize=25,
                xy=(0, 1), xycoords='axes fraction',
                xytext=(-30, 20), textcoords='offset points',
                ha='left', va='top')

df_opc["Count2 (/std_L)"].plot(ax=axs[1,0],label="Particles $\geq 0.5 \mu$m",zorder=1)
axs[1,0].scatter(df_frzT.index,opc_all,color="orange",zorder=2)
axs[1,0].set_xbound(xlims[0],xlims[1])
axs[1,0].set_xticks(xticks)
axs[1,0].set_xticklabels([])
axs[1,0].grid()
axs[1,0].set_yscale("log")
#axs[1,0].set_ylim([5,1e+6])
axs[1,0].xaxis.set_minor_locator(mpl.ticker.NullLocator())
axs[1,0].set_ylabel("#/L$_{std}$")
axs[1,0].set_xlabel(None)
axs[1,0].set_title("Concentration of particles $\geq 0.5 \mu$m")
#axs[1,0].legend(loc="upper left",frameon=False)
axs[1,0].annotate("(b)",fontsize=25,
                xy=(0, 1), xycoords='axes fraction',
                xytext=(-30, 20), textcoords='offset points',
                ha='left', va='top')

df_wisp["wind_speed"].plot(ax=axs[2,0], label="Wind speed",zorder=1)
axs[2,0].scatter(df_frzT.index,wisp_all,color="orange",zorder=2)
axs[2,0].set_xbound(xlims[0],xlims[1])
axs[2,0].set_xticks(xticks)
axs[2,0].set_xticklabels([])
axs[2,0].set_title("Wind speed")
#axs[2,0].legend(loc="upper left",frameon=False)
axs[2,0].grid()
#axs[2,0].set_ylim([0,25])
axs[2,0].set_ylabel("m/s")
axs[2,0].set_xlabel(None)
axs[2,0].annotate("(c)",fontsize=25,
                xy=(0, 1), xycoords='axes fraction',
                xytext=(-30, 20), textcoords='offset points',
                ha='left', va='top')

df_widir["wind_from_direction"].plot(ax=axs[3,0],label="Wind direction",zorder=1)
axs[3,0].scatter(df_frzT.index,widir_all,color="orange",zorder=2)
axs[3,0].set_xbound(xlims[0],xlims[1])
axs[3,0].set_title("Origin direction of wind")
#axs[3,0].legend(loc="upper left",frameon=False)
axs[3,0].set_xticks(xticks)
axs[3,0].set_xticklabels([])
axs[3,0].grid()
axs[3,0].set_ylabel("degrees")
axs[3,0].set_xlabel(None)
axs[3,0].annotate("(d)",fontsize=25,
                xy=(0, 1), xycoords='axes fraction',
                xytext=(-30, 20), textcoords='offset points',
                ha='left', va='top')

df_sea["TOC (mg/L)"].plot(ax=axs[4,0],linestyle="--",color="seagreen",zorder=2)
axs[4,0].scatter(df_sea.index, df_sea["TOC (mg/L)"],color="seagreen",zorder=3)
axs[4,0].scatter(df_frzT.index[7:50],toc_all,color="orange",zorder=1)
axs[4,0].set_xbound(xlims[0],xlims[1])
axs[4,0].set_xticks(xticks)
#ax_twin = axs[4,0].twinx()
#df_sea["T 50 % frozen (degC)"].plot(ax=ax_twin, linestyle="--")
#ax_twin.scatter(df_sea.index, df_sea["T 50 % frozen (degC)"])
#ax_twin.set_xbound(xlims[0],xlims[1])
#ax_twin.set_xticks(xticks)
axs[4,0].grid()
axs[4,0].set_ylabel("mg/L")
#axs[4,0].yaxis.label.set_color("seagreen")
#axs[4,0].tick_params(axis='y', colors='seagreen')
#ax_twin.set_ylabel("$^{\circ}$C")
#ax_twin.yaxis.label.set_color("tab:blue")
#ax_twin.tick_params(axis='y', colors='tab:blue')
#axs[4,0].set_title("TOC and temperature at 50 % activated INPs in SML")
axs[4,0].set_title("TOC in surface sea water")
axs[4,0].set_xlabel(None)
axs[4,0].annotate("(e)",fontsize=25,
                xy=(0, 1), xycoords='axes fraction',
                xytext=(-30, 20), textcoords='offset points',
                ha='left', va='top')

axs[0,1].scatter(opc_all,wisp_all,color="orange")
axs[0,1].grid(alpha=0.5)
axs[0,1].set_ylabel("Wind speed (m/s)")
axs[0,1].set_xlabel("Particles $\geq 0.5 \mu$m (#/L$_{std}$)")
axs[0,1].set_xscale("log")
axs[0,1].annotate("R: %.2f, R$^2$: %.2f" %(functions.r(opc_all,wisp_all), functions.rsquared(opc_all,wisp_all)),
                xy=(0, 1), xycoords='axes fraction',
                xytext=(5, 16), textcoords='offset points',
                ha='left', va='top')
axs[0,1].annotate("(f)",fontsize=25,
                xy=(0, 1), xycoords='axes fraction',
                xytext=(-30, 20), textcoords='offset points',
                ha='left', va='top')


axs[1,1].scatter(t50,opc_all,color="orange")
axs[1,1].grid(alpha=0.5)
axs[1,1].set_ylabel("Particles $\geq 0.5 \mu$m (#/L$_{std}$)")
axs[1,1].set_yscale("log")
axs[1,1].annotate("R: %.2f, R$^2$: %.2f" %(functions.r(t50,opc_all),functions.rsquared(t50,opc_all)),
                xy=(0, 1), xycoords='axes fraction',
                xytext=(5, 16), textcoords='offset points',
                ha='left', va='top')
axs[1,1].annotate("(g)",fontsize=25,
                xy=(0, 1), xycoords='axes fraction',
                xytext=(-30, 20), textcoords='offset points',
                ha='left', va='top')

axs[2,1].scatter(t50,wisp_all,color="orange")
axs[2,1].grid(alpha=0.5)
axs[2,1].set_ylabel("Wind speed (m/s)")
axs[2,1].annotate("R: %.2f, R$^2$: %.2f" %(functions.r(t50,wisp_all),functions.rsquared(t50,wisp_all)), 
                xy=(0, 1), xycoords='axes fraction',
                xytext=(5, 16), textcoords='offset points',
                ha='left', va='top')
axs[2,1].annotate("(h)",fontsize=25,
                xy=(0, 1), xycoords='axes fraction',
                xytext=(-30, 20), textcoords='offset points',
                ha='left', va='top')

#axs[3,1].scatter(t50,widir_all,color="orange")
#axs[3,1].grid(alpha=0.5)
#axs[3,1].set_ylabel("Wind origin direction")
#axs[3,1].annotate("R: %.2f, R$^2$: %.2f" %(functions.r(t50,widir_all),functions.rsquared(t50,widir_all)), 
#                xy=(0, 1), xycoords='axes fraction',
#                xytext=(5, 16), textcoords='offset points',
#                ha='left', va='top')

df_widir_binned.boxplot(ax=axs[3,1],vert=False,boxprops={"color":"orange"},whiskerprops={"color":"orange"},medianprops={"color":"black"},)
axs[3,1].grid(alpha=0.5)
#axs[3,1].set_ylabel("Wind origin direction")
axs[3,1].annotate("(i)",fontsize=25,
                xy=(0, 1), xycoords='axes fraction',
                xytext=(-30, 20), textcoords='offset points',
                ha='left', va='top')
axtwin = axs[3,1].twinx()
axtwin.set_ybound(axs[3,1].get_ylim())
axtwin.set_yticks(axs[3,1].get_yticks())
axtwin.set_yticklabels(["N=%d"%len(np.where(widir_all_binned==1)[0]),
                        "N=%d"%len(np.where(widir_all_binned==2)[0]),
                        "N=%d"%len(np.where(widir_all_binned==3)[0]),
                        "N=%d"%len(np.where(widir_all_binned==4)[0])])

axs[4,1].scatter(t50[7:50],toc_all,color="orange")
axs[4,1].scatter(t50[7:50].values[np.where(widir_all_binned[7:50]==3)],toc_all[np.where(widir_all_binned[7:50]==3)],color="springgreen")
axs[4,1].scatter(t50[7:50].values[np.where(widir_all_binned[7:50]==4)],toc_all[np.where(widir_all_binned[7:50]==4)],color="deepskyblue")
axs[4,1].grid(alpha=0.5)
axs[4,1].set_ylabel("TOC (mg/L)")
axs[4,1].set_xlabel("Temperature at 50 % \n frozen fraction ($^{\circ}$C)")
axs[4,1].annotate("R: %.2f, R$^2$: %.2f" %(functions.r(t50[7:50],toc_all),functions.rsquared(t50[7:50],toc_all)), 
                xy=(0, 1), xycoords='axes fraction',
                xytext=(5, 16), textcoords='offset points',
                ha='left', va='top')
axs[4,1].annotate("(j)",fontsize=25,
                xy=(0, 1), xycoords='axes fraction',
                xytext=(-30, 20), textcoords='offset points',
                ha='left', va='top')

plt.savefig(wpath+"factors_toc.pdf", bbox_inches="tight")
plt.savefig(wpath+"factors_toc.png", bbox_inches="tight")

print("Southwesterly")
print(functions.r(t50[7:50].values[np.where(widir_all_binned[7:50]==3)],toc_all[np.where(widir_all_binned[7:50]==3)]),
        functions.rsquared(t50[7:50].values[np.where(widir_all_binned[7:50]==3)],toc_all[np.where(widir_all_binned[7:50]==3)]))
print("Arctic")
print(functions.r(t50[7:50].values[np.where(widir_all_binned[7:50]==4)],toc_all[np.where(widir_all_binned[7:50]==4)]),
        functions.rsquared(t50[7:50].values[np.where(widir_all_binned[7:50]==4)],toc_all[np.where(widir_all_binned[7:50]==4)]))

print("All")
print(functions.r(t50[7:50].values[np.where(widir_all_binned[7:50]>=3)],toc_all[np.where(widir_all_binned[7:50]>=3)]),
        functions.rsquared(t50[7:50].values[np.where(widir_all_binned[7:50]>=3)],toc_all[np.where(widir_all_binned[7:50]>=3)]))


