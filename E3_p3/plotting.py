import glob
import numpy as np
import matplotlib.pyplot as plt
import os




data_dir = "data"                                   # same directory as main.c output
pattern  = os.path.join(data_dir, "energy_magnet_vs_time_beta:*.txt")
files    = sorted(glob.glob(pattern))

if not files:
    raise FileNotFoundError(
        f"No files matched '{pattern}'.")

plt.figure()
for i, fpath in enumerate(files):
    data = np.loadtxt(fpath, comments="#")
    time   = data[:, 0]
    energy = data[:, 1]
    magnet = data[:, 2]
    beta   = data[0, 3]
    plt.plot(time, energy, label = rf"$\beta$= {beta:.3f}")
plt.legend()
plt.xlabel("Steps")
plt.ylabel(rf"$E^*: (E/k_BT)$")
plt.savefig('plots/Energy_vs_time', dpi=200)
plt.show()



plt.figure() 
for i, fpath in enumerate(files):
    data = np.loadtxt(fpath, comments="#")
    time   = data[:, 0]
    energy = data[:, 1]
    magnet = data[:, 2]
    beta   = data[0, 3]
    plt.plot(time, magnet, label = rf"$\beta$= {beta:.3f}")

plt.xlabel("Steps")
plt.ylabel("magnetization")
plt.legend()
plt.savefig('plots/Magnetization_vs_time', dpi=200)
plt.show()

fpath_avg = "data/averages.txt"
data = np.loadtxt(fpath_avg, comments="#")
temp    = data[:, 0]
magnet  = data[:, 1]
energy  = data[:, 2]
energy2 = data[:, 3]
Cv      = data[:, 4]



Cv_minimum_index= np.argmax(Cv) 
T_crit = temp[Cv_minimum_index]
print(f"the critical temprature = {T_crit}")


plt.figure()
plt.plot(temp, magnet) 
plt.axvline(T_crit, color='k', linestyle='--')
plt.xlabel(r"$T^*: k_bT/J$")
plt.ylabel(r"$\langle M \rangle / N$")
plt.savefig('plots/Magnetization_vs_temp', dpi=200)
plt.show()

plt.figure()
plt.plot(temp, energy)
plt.axvline(T_crit, color='k', linestyle='--')
plt.ylabel(r"$\langle E \rangle /(N \cdot J)$")
plt.xlabel(r"$T^*: k_bT/J$")
plt.savefig('plots/energy_vs_temp', dpi=200)
plt.show()
    
plt.figure()
plt.plot(temp, energy2)
plt.axvline(T_crit, color='k', linestyle='--')
plt.xlabel(r"$T^*: k_bT/J$")
plt.ylabel("energy2")
plt.savefig('plots/energy2_vs_temp', dpi=200)
plt.show()



plt.figure()
plt.plot(temp, Cv)
plt.axvline(T_crit, color='k', linestyle='--')
plt.ylabel(r"$ \beta^2 \left( \langle E^2 \rangle - \langle E \rangle^2 \right)/N  $")
plt.xlabel(r"$T^*: k_bT/J$")
plt.savefig('plots/Cv_vs_temp', dpi=200)
plt.show()
