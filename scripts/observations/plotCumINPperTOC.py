import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size':20})
import pandas as pd
from sklearn.linear_model import LinearRegression
from Meyers import meyers

path = "/projects/NS9600K/astridbg/data/observations/Sea/"
fnameFrzT_shore = "Sea_nucleiT_cal.csv"
fnameConc_shore = "Sea_nucleiOut.csv"
fnameTOC_shore = "Sea_TOC.csv"
fnameFrzT_open = "BSea_nucleiT_cal.csv"
fnameConc_open = "BSea_nucleiOut.csv"
fnameTOC_open = "BSea_TOC.csv"
wpath = "/projects/NS9600K/astridbg/master/figures/observations/"

nucleiT_shore = pd.read_csv(path+fnameFrzT_shore, index_col=0)
nucleiOut_shore = pd.read_csv(path+fnameConc_shore)
TOC_shore = pd.read_csv(path+fnameTOC_shore,header=None)
nucleiT_open = pd.read_csv(path+fnameFrzT_open, index_col=0)
nucleiOut_open = pd.read_csv(path+fnameConc_open)
TOC_open = pd.read_csv(path+fnameTOC_open,header=None)

# Convert TOC from mg/L to g/L

TOC_shore = TOC_shore * 0.001
TOC_open = TOC_open * 0.001

# Exclude the first and last well, as the concentration values here are not representative of reality

nucleiT_shore = nucleiT_shore.iloc[1:-1,:] 
nucleiOut_shore = nucleiOut_shore.iloc[1:-1,:]
nucleiT_open = nucleiT_open.iloc[1:-1,:]    
nucleiOut_open = nucleiOut_open.iloc[1:-1,:]

n_shore = len(nucleiOut_shore.iloc[0,:])
n_open = len(nucleiOut_open.iloc[0,:])

# Calculate cumulative INPs per gram of TOC

INPperTOC_shore = nucleiOut_shore
INPperTOC_open = nucleiOut_open

print("shore")
for i in range(n_shore):
    print(i)
    INPperTOC_shore.iloc[:,i] = INPperTOC_shore.iloc[:,i] / np.float64(TOC_shore.iloc[i])

print("open")
for i in range(n_open):
    print(i)
    INPperTOC_open.iloc[:,i] = INPperTOC_open.iloc[:,i] / np.float64(TOC_open.iloc[i])

# Shape data for parameterisation

X_shore = nucleiT_shore.to_numpy() 
Y_shore = INPperTOC_shore.to_numpy()
X_open = nucleiT_open.to_numpy()              
Y_open = INPperTOC_open.to_numpy()

X = np.concatenate((X_shore, X_open), axis=1)
Y = np.concatenate((Y_shore, Y_open), axis=1)

X = X.flatten()
Y = Y.flatten()

linreg = np.polyfit(X,np.log(Y), 1)
slope = linreg[0]
intercept = linreg[1]

slope_Wilson = -0.4459
intercept_Wilson = 11.2186

plt.figure(figsize=(8,6),dpi=300)
plt.title("Surface sea water at Andenes 18.03$-$30.03 2021",fontsize=22)
plt.grid(alpha=0.5)
alpha=1
"""
for i in range(n_shore):
    if i == 0:
        plt.scatter(nucleiT_shore.iloc[:,i], INPperTOC_shore.iloc[:,i],alpha = alpha, color="none", edgecolor="seagreen", label="Shore water")
    else:
        plt.scatter(nucleiT_shore.iloc[:,i], INPperTOC_shore.iloc[:,i],alpha = alpha, color="none", edgecolor="seagreen")
    alpha -= 0.0

for i in range(n_open):
    if i == 0:
        plt.scatter(nucleiT_open.iloc[:,i], INPperTOC_open.iloc[:,i],alpha = alpha, color="none", edgecolor="mediumturquoise", label="Open sea")
    else:
        plt.scatter(nucleiT_open.iloc[:,i], INPperTOC_open.iloc[:,i],alpha = alpha, color="none", edgecolor="mediumturquoise")
    alpha -= 0.0
"""
plt.yscale("log")
plt.ylim(10**(3.5),10**(11))
plt.xlim(-25,-2)
x = np.linspace(-25,-2,100)
#plt.plot(x, np.exp(intercept + slope*x), linewidth=2, color="black",
#        label="exp("+str(round(intercept,3))+" - "+str(round(np.sign(slope)*slope,3))+r"$\times T$)")
plt.plot(x, 1e-1*np.exp(intercept_Wilson + slope_Wilson*x),linewidth=2,color="red",linestyle="--")
plt.plot(x, 1e+1*np.exp(intercept_Wilson + slope_Wilson*x),linewidth=2,color="red",linestyle="--")
plt.plot(x, np.exp(intercept_Wilson + slope_Wilson*x),linewidth=2,color="red",
        label="Wilson et.al.")#": exp("+str(round(intercept_Wilson,3))+" - "+str(round(np.sign(slope_Wilson)*slope_Wilson,3))+r"$\times T$)")

# Shrink current axis's height by 10% on the bottom
ax = plt.gca()
box = ax.get_position()
ax.set_position([box.x0, box.y0 + box.height * 0.1,
                 box.width, box.height * 0.9])

# Put a legend below current axis
ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=3,borderpad=0.3, columnspacing=0.3, handletextpad=0.2)

plt.xlabel(r"Temperature $T$ [$^{\circ}$C]")
plt.ylabel(r"Cumulative INPs per gram TOC")
#plt.savefig(wpath+"INPpergTOC.pdf",bbox_inches="tight")
plt.savefig(wpath+"INPpergTOC_wilson.png",bbox_inches="tight")
