import numpy as np
import xarray as xr 


path = "/home/astridbg/Documents/data/"
filename = "predictedRAD_UNml_1901-01-15to2019-12-15_r0.nc"

dataset = xr.open_dataset(path+filename)
