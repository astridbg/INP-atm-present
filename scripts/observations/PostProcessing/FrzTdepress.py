"""
Author: Astrid Bragstad Gjelsvik
"""

import numpy as np
import pandas as pd
import glob

path_preproc = "/projects/NS9600K/astridbg/data/observations/Sea/"
wpath = "/projects/NS9600K/astridbg/data/observations/Sea/"
fname = "BSea_nucleiT.csv"

# Read in Coriolis freezing temperatures

nucleiT = pd.read_csv(path_preproc+fname)


# Adjust freezing temperatures with calibration values

slope = 0.89
intercept = 0.30

nucleiT_cal = nucleiT.iloc[:,:] * slope + intercept


print(nucleiT)
# Adjust freezing temperatures with freezing point depression

depression = 1.8 # For 35 ppt salinity https://www.metoffice.gov.uk/research/climate/cryosphere-oceans/sea-ice/overview

nucleiT_corr = nucleiT_cal + depression

nucleiT_corr.to_csv(wpath+"BSea_nucleiT_cal.csv")


print(nucleiT_cal)
