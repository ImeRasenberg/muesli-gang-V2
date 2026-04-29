#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 21:49:51 2026

@author: ime_rasenberg
"""

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



# fpath_avg = "data/critical_averages.txt"
# data = np.loadtxt(fpath_avg, comments="#")
# temp    = data[:, 0]
# magnet  = data[:, 1]
# energy  = data[:, 2]
# energy2 = data[:, 3]
# Cv      = data[:, 4]



# Cv_minimum_index= np.argmax(Cv) 
# T_crit = temp[Cv_minimum_index]
# print(f"the critical temprature = {T_crit}")


# plt.figure()
# plt.plot(temp, magnet) 
# plt.axvline(T_crit, color='k', linestyle='--')
# plt.xlabel(r"$T^*: k_bT/J$")
# plt.ylabel(r"$\langle M \rangle / N$")
# plt.savefig('plots/crit_Magnetization_vs_temp', dpi=200)
# plt.show()

# plt.figure()
# plt.plot(temp, energy)
# plt.axvline(T_crit, color='k', linestyle='--')
# plt.ylabel(r"$\langle E \rangle /(N \cdot J)$")
# plt.xlabel(r"$T^*: k_bT/J$")
# plt.savefig('plots/crit_energy_vs_temp', dpi=200)
# plt.show()
    
# plt.figure()
# plt.plot(temp, energy2)
# plt.axvline(T_crit, color='k', linestyle='--')
# plt.xlabel(r"$T^*: k_bT/J$")
# plt.ylabel("energy2")
# plt.savefig('plots/crit_energy2_vs_temp', dpi=200)
# plt.show()


# from scipy.signal import savgol_filter
# Cv_smooth = savgol_filter(Cv, window_length=5, polyorder=1) #local smoothing algoritm


# plt.figure()
# plt.plot(temp, Cv)
# #plt.plot(temp, Cv_smooth)
# plt.axvline(T_crit, color='k', linestyle='--')
# plt.ylabel(r"$ \beta^2 \left( \langle E^2 \rangle - \langle E \rangle^2 \right)/N  $")
# plt.xlabel(r"$T^*: k_bT/J$")
# plt.savefig('plots/crit_Cv_vs_temp', dpi=200)
# plt.show()



#%%
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def heat_capacity_model(T, Tc, alpha,a):
    val = np.abs(T - Tc)**(-alpha)-a
    return np.clip(val, a_min=None, a_max=3)



lower_bounds = [min(temp), 0.05,0]
upper_bounds = [max(temp), 2,10]

popt, pcov = curve_fit(heat_capacity_model, temp, Cv, 
    p0=[5.1, 0.8,0.5], bounds=(lower_bounds, upper_bounds),
    maxfev=10000
)
perr = np.sqrt(np.diag(pcov))

    
temp_fine = np.linspace(min(temp), max(temp), 100)

fit_line = heat_capacity_model(temp_fine, *popt)
fit_line[np.abs(temp_fine - popt[0]) < 0.01] = np.nan 

plt.figure(figsize=(8, 6))

plt.scatter(temp, Cv, color='royalblue', s=15, alpha=0.7, label='Simulation Data')
plt.plot(temp_fine, fit_line, color='crimson', lw=2, 
         label=f'Fit: $\\alpha$={popt[1]:.3f}±{perr[1]:.3f}')

plt.axvline(popt[0], color='black', linestyle='--', alpha=0.6, 
            label=f'Est. $T_c$={popt[0]:.3f}±{perr[0]:.3f}')

plt.title("Heat Capacity Fit: 2D Ising Model", fontsize=14)
plt.ylabel(r"$\beta^2 (\langle E^2 \rangle - \langle E \rangle^2)/N$", fontsize=12)
plt.xlabel(r"$T^*: k_bT/J$", fontsize=12)
plt.legend(frameon=True, loc='upper right')
plt.grid(True, linestyle=':', alpha=0.6)

plt.tight_layout()
plt.savefig('plots/Cv_fit_analysis.png', dpi=200)
plt.show()



#%%
files = glob.glob('data2/correlation_vs_time_Temp*.txt')

plt.figure(figsize=(10, 6))
files.sort()
N = 2

data_stored = {}
for file_path in files:
    try:
        data = np.loadtxt(file_path)
        
        if data.size == 0:
            continue

        # Extract columns
        steps = data[N:, 0]
        correlation = data[N:, 1]
        
        # Get temperature from the first row of the 3rd column
        # (Assuming T is constant throughout the file)
        temp_val = data[0, 2]

        if temp_val not in data_stored.keys():
            data_stored[temp_val] = {"N":1,
                                     1:{"steps":steps,
                                        "correlation":correlation
                                        }
                                     }
        else:
            data_stored[temp_val][ data_stored[temp_val]["N"]+1 ] = {"steps":steps,
                                                                   "correlation":correlation
                                                                   }
            data_stored[temp_val]["N"] +=1
        
    except Exception as e:
        print(f"Error processing {file_path}: {e}")






for key in data_stored.keys():
    steps *= 0
    correlation *= 0
    temp_data = data_stored[key]
    N=temp_data["N"]
    for key2 in temp_data:
        if key2 == "N":
            continue
        else:
            steps += temp_data[key2]["steps"]/N
            correlation += temp_data[key2]["correlation"]/N


    plt.plot(steps, correlation, label=f'T = {key}')
    plt.xlabel('Simulation Step')
    plt.ylabel('Correlation Function')
    
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    
    # plt.savefig('plots/correlation_plot.png')
    plt.show()

#%%


import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Define the model function
def exponential_decay(x, a, tau, c):
    return a * np.exp(-x / tau) + c

# Lists to store data for the summary plot
T_results = []
tau_results = []
tau_errors = []

# Sorting keys ensures the final T vs Tau plot is ordered correctly
sorted_keys = sorted(data_stored.keys(), key=lambda x: float(x))

for key in sorted_keys:
    steps = 0
    correlation = 0
    temp_data = data_stored[key]
    N = temp_data["N"]
    
    for key2 in temp_data:
        if key2 == "N":
            continue
        steps += temp_data[key2]["steps"] / N
        correlation += temp_data[key2]["correlation"] / N

    try:
        # Initial guesses
        initial_guess = [correlation[0] - correlation[-1], 10, correlation[-1]]
        
        # Adjusting sigma decay: 
        # For T=5.5, we use a larger denominator (150.0) so sigma grows slower,
        # giving more weight to later points compared to the default (50.0).
        sigma_scale = 200.0 if float(key) in [5.25, 5.5] else 50.0
        sigma = np.exp(steps / sigma_scale) 

        popt, pcov = curve_fit(
            exponential_decay, 
            steps, 
            correlation, 
            p0=initial_guess,
            sigma=sigma
        )
        
        perr = np.sqrt(np.diag(pcov))
        tau_val, tau_err = popt[1], perr[1]

        # Store results for the summary plot
        T_results.append(float(key))
        tau_results.append(tau_val)
        tau_errors.append(tau_err)
        
        # Plotting individual fits
        fit_curve = exponential_decay(steps, *popt)
        line, = plt.plot(steps, correlation, 'o', markersize=2, alpha=0.2)
        plt.plot(steps, fit_curve, color=line.get_color(), 
                 label=f'$T = {key}$ ($\\tau = {tau_val:.2f} \\pm {tau_err:.2f}$)')
        
        
        plt.xlabel('Simulation Step')
        plt.ylabel('Correlation Function')
        plt.legend(fontsize='small', ncol=2)
        plt.grid(True, linestyle='--', alpha=0.6)
        plt.tight_layout()
        plt.show()
                 
    except Exception as e:
        print(f"Fit failed for T={key}: {e}")

N = 2
# --- New Summary Plot: Tau vs Temperature ---
plt.figure(figsize=(8, 5))
plt.errorbar(T_results[2:], tau_results[2:], yerr=tau_errors[2:], 
             fmt='o-', capsize=5, markersize=6)#, ecolor='red', color='black')
plt.axvline(T_crit, color='k', linestyle='--')
plt.xlabel('Temperature ($T$)')
plt.ylabel('Decay Time ($\\tau$)')
plt.title('Correlation Decay Time vs Temperature')
plt.grid(True, linestyle=':', alpha=0.7)
plt.tight_layout()
plt.show()