import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size':14})
import pandas as pd

path1 = "/projects/NS9600K/astridbg/data/observations/Coriolis_preprocessed/"
path2 = "/projects/NS9600K/astridbg/data/observations/Coriolis_postprocessed/"
fname1 = "Coriolis_nucleiT.csv"
fname2 = "Coriolis_nucleiOut.csv"
fname3 = "Coriolis_nucleiT_cal.csv"
fname4 = "Coriolis_nucleiOut_std.csv"
wpath = "/projects/NS9600K/astridbg/master/figures/observations/PostProcessing/"


nucleiT = pd.read_csv(path1+fname1)
nucleiOut = pd.read_csv(path1+fname2)
nucleiT_cal = pd.read_csv(path2+fname3,index_col=0)
nucleiOut_std = pd.read_csv(path2+fname4,index_col=0)

nCor = len(nucleiT.iloc[0,:])

plt.figure(figsize=(8,6),dpi=300)
plt.title("INP concentrations at Andenes 15.03 - 30.03 2021")
plt.grid()
alpha=1
for i in range(nCor):
    if i==0:
        plt.scatter(nucleiT.iloc[:,i],nucleiOut.iloc[:,i], 
                alpha = alpha, color="none", edgecolor="red",label='Uncorrected')
        plt.scatter(nucleiT_cal.iloc[:,i],nucleiOut_std.iloc[:,i], 
                alpha = alpha, color="none", edgecolor="black",label='Corrected')
    else:
        plt.scatter(nucleiT.iloc[:,i],nucleiOut.iloc[:,i], 
                alpha = alpha, color="none", edgecolor="red")
        plt.scatter(nucleiT_cal.iloc[:,i],nucleiOut_std.iloc[:,i], 
                alpha = alpha, color="none", edgecolor="black")
    alpha -= 0.01
plt.yscale("log")
plt.ylim(10**(-4),10**(-0.5))
plt.xlabel(r"Temperature $T$ [$^{\circ}$C]")
plt.ylabel(r"$n_S(T)$ [l$^{-1}$]")
plt.legend()
plt.savefig(wpath+"INPconc_compar.pdf",bbox_inches="tight")
plt.show()


