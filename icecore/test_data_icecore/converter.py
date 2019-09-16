import pandas as pd
import numpy as np

a = np.arange(1000,2005,1)
b = np.load("accumulation.npy")
c = np.load("total_icecore_depth.npy")
d = np.load("depth_horizons.npy")
e = np.load("d18O_P.npy")
f = np.load("temperature.npy")

df = pd.DataFrame({"TIME" : a, 
	"TEMPERATURE": f, 
	"ACCUMULATION" : b, 
	"DEPTH": c, 
	"DEPTH HORIZONS": d, 
	"D18O": e})
df.to_csv("submission2.csv", index=False)