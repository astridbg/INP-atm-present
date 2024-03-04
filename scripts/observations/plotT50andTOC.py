import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size':20})
import pandas as pd
import functions

path = "/projects/NS9600K/astridbg/data/observations/Sea/"
wpath = "/projects/NS9600K/astridbg/master/figures/observations/"
fname1 = "Sea_nucleiT_cal.csv"
fname2 = "Sea_TOC.csv"
fname3 = "BSea_nucleiT_cal.csv"
fname4 = "BSea_TOC.csv"

s_temp = pd.read_csv(path+fname1, index_col=0)
s_toc = pd.read_csv(path+fname2, header=None)

b_temp = pd.read_csv(path+fname3, index_col=0)
b_toc = pd.read_csv(path+fname4, header=None)

s_t50 = s_temp.iloc[48,:]
b_t50 = b_temp.iloc[48,:]

plt.figure(figsize=(8,6),dpi=300)
plt.title("Surface sea water at Andenes 18.03$-$30.03 2021",fontsize=22)
for i in range(len(s_t50)):
    if i == 0:
        plt.scatter(s_t50[i],s_toc.iloc[i],s=100,color="seagreen",alpha=0.8, 
                label="Shore water, R: %.2f" %functions.r(s_t50.values,np.transpose(s_toc.values)[0]))
    else:
        plt.scatter(s_t50[i],s_toc.iloc[i],s=100,color="seagreen",alpha=0.8)
for i in range(len(b_t50)):
    if i == 0:
        plt.scatter(b_t50[i],b_toc.iloc[i],s=100,color="mediumturquoise",alpha=0.8, 
                label="Open sea, R: %.2f" %functions.r(b_t50.values,np.transpose(b_toc.values)[0]))
    else:
        plt.scatter(b_t50[i],b_toc.iloc[i],s=100,color="mediumturquoise",alpha=0.8)
plt.grid(alpha=0.5)
plt.xlabel(r"Temperature at 50 % frozen fraction [$^{\circ}$C]")
plt.ylabel(r"Total Organic Carbon [mg/L]")

# Shrink current axis's height by 10% on the bottom
ax = plt.gca()
box = ax.get_position()
ax.set_position([box.x0, box.y0 + box.height * 0.1,
                 box.width, box.height * 0.9])

# Put a legend below current axis
ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=2,borderpad=0.3, columnspacing=0.3, handletextpad=0.2)

all_t50 = np.append(s_t50.values,b_t50.values)
all_toc = np.append(np.transpose(s_toc.values)[0],np.transpose(b_toc.values)[0])
ax.annotate("R: %.2f, R$^2$: %.2f" %(functions.r(all_t50,all_toc),functions.rsquared(all_t50,all_toc)),
                xy=(0, 1), xycoords='axes fraction',
                xytext=(5, -5), textcoords='offset points',
                ha='left', va='top')

plt.savefig(wpath+"sea_t50toc.pdf",bbox_inches="tight")
plt.savefig(wpath+"sea_t50toc.png",bbox_inches="tight")
