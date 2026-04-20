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



data_dir = "data"                                   # same directory as main.c output
pattern  = os.path.join(data_dir, "energy_vs_time_Temp_*.txt")
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
    Temp    = data[0, 2]    

    plt.plot(time, energy, label = rf'$T = {Temp:.1f}$',color=colors[i], alpha = 0.6+ 0.1*i)
plt.ylim((0,350))
plt.ylabel(r'$E/\epsilon  $  (reduced energy)',fontsize = 10)
plt.xlabel(r"$t \sqrt{\frac{\epsilon}{\sigma m}}  $ (reduced time)", fontsize = 10)
plt.legend(loc='upper right')
plt.savefig('energy_vs_T', dpi=200)
plt.show()

from scipy.signal import savgol_filter

plt.figure()
for i, fpath in enumerate(files):
    data = np.loadtxt(fpath, comments="#")
    time   = data[:, 0]
    energy = data[:, 1]
    energy_smooth = savgol_filter(energy, window_length=500, polyorder=3) #local smoothing algoritm
    Temp    = data[0, 2]    

    plt.plot(time, energy_smooth, label = rf'$T = {Temp:.1f}$',color=colors[i], alpha = 0.6+ 0.1*i)
plt.axvline(0.1, 0, 250, color = 'k')
plt.ylim((0,250))
plt.ylabel(r'$E/\epsilon  $  (reduced energy)',fontsize = 10)
plt.xlabel(r"$t \sqrt{\frac{\epsilon}{\sigma m}}  $ (reduced time)", fontsize = 10)
plt.legend(loc='upper right')
plt.savefig('smooth_energy_vs_T', dpi=200)
plt.show()

from scipy.optimize import curve_fit

data_dir = "data"
pattern  = os.path.join(data_dir, "MSD.txt")
files    = sorted(glob.glob(pattern))




for fpath in files:
    data = np.loadtxt(fpath, comments="#")
    time = data[:, 0]
    msd  = data[:, 1]

fit_time = time[1000:] # cut off the equilibrium time
fit_msd = msd[1000:]
def func(x, a, b,):
    return a + 6*b*x 

popt, pcov = curve_fit(func, fit_time, fit_msd)
print(popt)
plt.figure()
plt.plot(time, msd, label = "MSD")
#plt.plot(fit_time, func(fit_time, popt[0], popt[1]), label = f"6Dt, D = {popt[1]:.2f}")
plt.legend()
plt.ylabel(r'$MSD/\sigma^2  $  (reduced MSD)',fontsize = 10)
plt.xlabel(r"$t \sqrt{\frac{\epsilon}{\sigma m}}  $ (reduced time)", fontsize = 10)
plt.xlim((0,2))
plt.savefig('MSD_115', dpi=200)
plt.show()

