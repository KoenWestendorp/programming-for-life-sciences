import pandas as pd
import matplotlib.pyplot as plt
#from pandas.compat import StringIO
import numpy as np

file_path = "./chol-z.dat"
df = pd.read_csv(file_path, skiprows=2, header=None, dtype=float, sep=r'\s+')

time =   df[0]  # time (ps)
center = df[1]  # z-position of bilayer center (nm)
#chols =  df[[2]] # z-position of cholesterol's head oxygen (nm). One column per cholesterol  

# plt.plot(time, center)
# plt.plot(time, df[2])
# plt.plot(time, df[3])
# plt.show()

dist_from_center = df[2] - center

M = np.array(dist_from_center)
M[M > 0]

plt.plot(dist_from_center[dist_from_center >= 0])
plt.plot(dist_from_center[dist_from_center < 0])

plt.show()
