import glob
import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib.cm as cm



data_dir = "data"                                   # same directory as main.c output
pattern  = os.path.join(data_dir, "energy_vs_time_dt_temp_*.txt")
files    = sorted(glob.glob(pattern), reverse=True)

if not files:
    raise FileNotFoundError(
        f"No files matched '{pattern}'.")

colors = ["#049711", "#e9880ae7", "#2c2ea0f1", "#ff0e0e"]

plt.figure()
for i, fpath in enumerate(files):
    data = np.loadtxt(fpath, comments="#")
    time   = data[:, 0]
    energy = data[:, 1]
    dt     = data[0, 2]    

    plt.plot(time, energy, label = rf'$\Delta t = {dt:.0e}$',color=colors[i], alpha = 0.6+ 0.1*i)
plt.ylim((0,150))
plt.ylabel(r'$E/\epsilon  $  (reduced energy)',fontsize = 10)
plt.xlabel(r"$t \sqrt{\frac{\epsilon}{\sigma m}}  $ (reduced time)", fontsize = 10)
plt.legend(loc='upper right')
plt.savefig('energy_vs_dt', dpi=200)
plt.show()



