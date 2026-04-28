import glob
import numpy as np
import matplotlib.pyplot as plt
import os




data_dir = "data"                                   # same directory as main.c output
pattern  = os.path.join(data_dir, "correlation_vs_time_Temp:*.txt")
files    = sorted(glob.glob(pattern))

if not files:
    raise FileNotFoundError(
        f"No files matched '{pattern}'.")
for i, fpath in enumerate(files):
    plt.figure()
    data = np.loadtxt(fpath, comments="#")
    time  = data[:, 0]
    corr  = data[:, 1]
    Temp  = data[0, 2]
    plt.plot(time, corr, label = rf"$\beta$= {Temp:.3f}")
    plt.legend()
    plt.show()
