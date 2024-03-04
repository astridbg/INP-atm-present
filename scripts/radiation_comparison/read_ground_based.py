import numpy as np
import pandas as pd
import glob
import os
import re

### PANGAEA IDs ###
# Alert  932867
# Barrow  959215
# Eureka 933177
# Ny-Ålesund 914927
# Cape Baranova (no collection) e.g. 892289 - EXCLDUED!! Only 2016

#import pangaeapy.pandataset as pngd
#ds = pngd.PanDataSet(892289)
#print(ds.title)
#print(ds.data.head())


data_path = "/home/astridbg/Documents/data/stations/"
#station_path = "ALE_basic_rad_2004-2014/datasets/"
#station_path = "EUR_basic_rad_2007-2011/datasets/"
station_path = "BAR_radiation_1992-01_etseq/datasets/"
#station_path = "NYA_radiation_2006-05_etseq/datasets/"
path = data_path + station_path
#savename = "ALE_basic_rad_2004-2014.csv"
#savename = "NYA_radiation_2006-2023.csv"
#savename = "EUR_radiation_2007-2011.csv"
savename = "BAR_radiation_1992-2022.csv"
filenames = glob.glob(path + "*.tab")

i=1
for filename in sorted(filenames):
    start_line=0
    res = False
    with open(filename,'r') as f:
        while not res:
            b=f.readline()
            res = re.match(r"Date*",b)
            if res:
                break 
            start_line = start_line+1    
    month = filename.split('_')[-1].split('.')[0]
    print(month,start_line)
    df = pd.read_csv(filename, header=start_line, sep='\t', index_col=0)
    df = df.mean(axis=0)
    df = pd.DataFrame(list(zip(df.index, df.values))).transpose()
    df.columns = df.iloc[0]
    df = df.drop(0)
    df["Time"] = [month]
    if i>1:
        dfs = pd.concat([dfs, df], ignore_index=True, sort=False)
    else:
        dfs=df
    print(i)
    i=i+1

dfs = dfs.set_index("Time")
print(dfs)
dfs.to_csv(path + savename)