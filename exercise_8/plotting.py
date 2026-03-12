#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 17:05:32 2026

@author: ime_rasenberg
"""

import pandas as pd
import matplotlib.pyplot as plt

# 1. Load the data
# Assuming your C code output: fprintf(fp, "%d\t%f\t%f\n", step, ms.average_pressure, ms.mu_excess);
try:
    data = pd.read_csv('data/beta_0.900000___T_0.500000/measurements.dat', sep='\t', names=['Step', 'Pressure', 'Mu_Excess', "delta"])
except FileNotFoundError:
    print("Error: measurements.dat not found. Run your C simulation first!")
    exit()

# 2. Create the figure
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

# --- Plot Pressure ---
ax1.plot(data['Step'], data['Pressure'], alpha=0.3, color='royalblue', label='Raw Pressure')
ax1.plot(data['Step'], data['Pressure'].rolling(window=50).mean(), color='navy', label='Moving Avg (50)')
ax1.set_ylabel('Pressure ($P$)')
ax1.set_title('NVT Monte Carlo Simulation Results')
ax1.legend()
ax1.grid(True, linestyle='--', alpha=0.6)

# --- Plot Chemical Potential ---
ax2.plot(data['Step'], data['Mu_Excess'], alpha=0.3, color='tomato', label='Raw $\mu_{ex}$')
ax2.plot(data['Step'], data['Mu_Excess'].rolling(window=50).mean(), color='darkred', label='Moving Avg (50)')
ax2.set_ylabel('Excess Chemical Potential ($\mu_{ex}$)')
ax2.set_xlabel('Monte Carlo Step')
ax2.legend()
ax2.grid(True, linestyle='--', alpha=0.6)

plt.tight_layout()
plt.savefig('simulation_results.png', dpi=300)
plt.show()

#%%


# import os
# import glob
# import re
# import pathlib as p

# files =glob.glob("data/beta_*___T_*/measurement.dat")

# pattern = re.compile(r"beta_(.+)___T_(.+)")
# results = []

# for folder in p(".").iterdir():
#     match= pattern.match(folder.name)
#     if match and folder.is_dir():
#         dat_file = folder/"measurements.dat"
#         if dat_file.exists():
#             beta,T = match.group(1)




# Import required libraries
import os
import numpy as np
import matplotlib.pyplot as plt
import re


base_path = "./data"
pattern = re.compile(r"beta_(.+)___T_(.+)")


data_dict = {}

# Iterate over subfolders
for subfolder in os.listdir(base_path):
    sub_path = os.path.join(base_path, subfolder)
    info_path = os.path.join(sub_path, "measurements.dat")
    
    match= pattern.match(subfolder)
    if not match:
        print("didnt find a match")
        continue
    rho = float(match.group(1))
    T = float(match.group(2))
    
    if os.path.isdir(sub_path) and os.path.isfile(info_path):
        with open(info_path, "r") as f:
            lines = f.readlines()
        
        
        # Remaining lines contain the inf matrix
        data = []
        for line in lines:
            values = list(map(float, line.split()))
            if any(values):  # skip the last zero-line
                data.append(values)
        
        data = np.array(data)
        
        data_dict[subfolder] = {
            "rho":rho,
            "T":T,
            "Step": data[:, 0],
            "Pressure": data[:, 1],
            "Mu_Excess": data[:, 2],
        }

new_dict = {}
count=0

start_a = 2000
for key1 in data_dict:
    
    data = data_dict[key1]
    new_dict[count]={}
    
    
    mask = data["Step"]>start_a & np.isinf(data["Step"])
    new_dict[count]["Pressure"] = np.average(data["Pressure"][mask])
    new_dict[count]["Mu_Excess"] = np.average(data["Mu_Excess"][mask])
    new_dict[count]["rho"] = data["rho"]
    new_dict[count]["T"] = data["T"]
    

    count+=1
    
    # # 2. Create the figure
    # fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

    # # --- Plot Pressure ---
    # ax1.plot(data['Step'], data['Pressure'], alpha=0.3, color='royalblue', label='Raw Pressure')
    # # ax1.plot(data['Step'], data['Pressure'].rolling(window=50).mean(), color='navy', label='Moving Avg (50)')
    # ax1.set_ylabel('Pressure ($P$)')
    # ax1.set_title('NVT Monte Carlo Simulation Results')
    # ax1.legend()
    # ax1.grid(True, linestyle='--', alpha=0.6)

    # # --- Plot Chemical Potential ---
    # ax2.plot(data['Step'], data['Mu_Excess'], alpha=0.3, color='tomato', label='Raw $\mu_{ex}$')
    # # ax2.plot(data['Step'], data['Mu_Excess'].rolling(window=50).mean(), color='darkred', label='Moving Avg (50)')
    # ax2.set_ylabel('Excess Chemical Potential ($\mu_{ex}$)')
    # ax2.set_xlabel('Monte Carlo Step')
    # ax2.legend()
    # ax2.grid(True, linestyle='--', alpha=0.6)

    # plt.tight_layout()
    # plt.savefig('simulation_results.png', dpi=300)
    # plt.show()
    

#%%

import os
import numpy as np
import matplotlib.pyplot as plt
import re
from collections import defaultdict

base_path = "./data"
# Regex to match folder names like 'beta_0.1___T_1.0'
pattern = re.compile(r"beta_(.+)___T_(.+)")

# Dictionary to hold raw data
data_dict = {}

# 1. Load Data
if os.path.exists(base_path):
    for subfolder in os.listdir(base_path):
        sub_path = os.path.join(base_path, subfolder)
        info_path = os.path.join(sub_path, "measurements.dat")
        
        match = pattern.match(subfolder)
        if not match:
            continue
            
        rho = float(match.group(1))
        T = float(match.group(2))
        
        if os.path.isdir(sub_path) and os.path.isfile(info_path):
            try:
                # Load data, skipping the last zero-line if necessary
                data = np.loadtxt(info_path)
                # Filter out lines that are all zeros if the file ends with a dummy line
                data = data[~np.all(data == 0, axis=1)]
                
                data_dict[subfolder] = {
                    "rho": rho,
                    "T": T,
                    "Step": data[:, 0],
                    "Pressure": data[:, 1],
                    "Mu_Excess": data[:, 2],
                }
            except Exception as e:
                print(f"Error reading {info_path}: {e}")

# 2. Process Data (Averages and Errors)
# We group results by Temperature for plotting
results = defaultdict(list)
start_a = 3000

for key in data_dict:
    d = data_dict[key]
    
    # Corrected Mask: Use parentheses and check for finite values
    # We exclude the equilibration phase (Step > start_a)
    mask = (d["Step"] > start_a) & np.isfinite(d["Pressure"]) & np.isfinite(d["Mu_Excess"])
    
    if np.any(mask):
        p_slice = d["Pressure"][mask]
        mu_slice = d["Mu_Excess"][mask]
        
        # Calculate Average and Standard Error (Standard Deviation / sqrt(N))
        avg_p = np.mean(p_slice)
        err_p = np.std(p_slice)#/ np.sqrt(len(p_slice))
        
        avg_mu = np.mean(mu_slice)
        err_mu = np.std(mu_slice)#/ np.sqrt(len(mu_slice))
        
        results[d["T"]].append({
            "rho": d["rho"],
            "P": avg_p/d["T"],
            "P_err": err_p,
            "Mu": avg_mu,
            "Mu_err": err_mu
        })
        

# 3. Plotting
fig, (ax1) = plt.subplots(1,1, figsize=(7, 5))

# Sort temperatures to have a consistent legend
sorted_temps = sorted(results.keys())

for T in sorted_temps:
    # Sort by rho so lines connect properly
    data_list = sorted(results[T], key=lambda x: x["rho"])
    
    rhos = [d["rho"] for d in data_list]
    ps = [d["P"] for d in data_list]
    p_errs = [d["P_err"] for d in data_list]
    mus = [d["Mu"] for d in data_list]
    mu_errs = [d["Mu_err"] for d in data_list]
    
    # Plot Pressure
    ax1.errorbar(rhos, ps, yerr=p_errs, fmt='-o', capsize=3, label=f"$T^* = {T}$")
    
    # Plot Chemical Potential
    # ax2.errorbar(rhos, mus, yerr=mu_errs, fmt='-s', capsize=3, label=f"$T^* = {T}$")

# Formatting Pressure Plot
ax1.set_xlabel(r"$\rho \sigma^3$")
ax1.set_ylabel(r"$\beta P \sigma^3$")
# ax1.set_title("Pressure vs Density")
ax1.grid(True, linestyle='--', alpha=0.7)
ax1.legend()



plt.tight_layout()
plt.savefig("pressure.png", dpi=300)
plt.show()

fig, (ax1) = plt.subplots(1,1, figsize=(7, 5))
for T in sorted_temps:
    # Sort by rho so lines connect properly
    data_list = sorted(results[T], key=lambda x: x["rho"])
    
    rhos = [d["rho"] for d in data_list]
    ps = [d["P"] for d in data_list]
    p_errs = [d["P_err"] for d in data_list]
    mus = [d["Mu"] for d in data_list]
    mu_errs = [d["Mu_err"] for d in data_list]
    
    # Plot Pressure
    # ax1.errorbar(rhos, ps, yerr=p_errs, fmt='-o', capsize=3, label=f"$T^* = {T}$")
    
    # Plot Chemical Potential
    ax1.errorbar(rhos, mus, yerr=mu_errs, fmt='-s', capsize=3, label=f"$T^* = {T}$")

# Formatting Pressure Plot
# Formatting Mu Plot
ax1.set_xlabel(r"$\rho \sigma^3$")
ax1.set_ylabel(r"${\mu_{ex}}/{ \epsilon} $")
# ax1.set_title("Excess $\mu$ vs Density")
ax1.grid(True, linestyle='--', alpha=0.7)
ax1.legend()



plt.tight_layout()
plt.savefig("Excess Chemical Potential.png", dpi=300)
plt.show()

fig, (ax1) = plt.subplots(1,1, figsize=(7, 5))
for T in sorted_temps:
    # Sort by rho so lines connect properly
    data_list = sorted(results[T], key=lambda x: x["rho"])
    
    # Plot Pressure
    # ax1.errorbar(rhos, ps, yerr=p_errs, fmt='-o', capsize=3, label=f"$T^* = {T}$")
    rhos = np.array([d["rho"] for d in data_list])
    mus = np.array([d["Mu"] for d in data_list])
    mu_errs = np.array([d["Mu_err"] for d in data_list])
    
    # Now this comparison works!
    mask = rhos < 0.85
    
    # Apply the mask to all arrays
    ax1.errorbar(rhos[mask], mus[mask], yerr=mu_errs[mask], 
                 fmt='-s', capsize=3, label=f"$T^* = {T}$")

# Formatting Mu Plot
ax1.set_xlabel(r"$\rho \sigma^3$")
ax1.set_ylabel(r"${\mu_{ex}}/{ \epsilon} $")
# ax1.set_title("Excess $\mu$ vs Density")
ax1.grid(True, linestyle='--', alpha=0.7)
ax1.legend()



plt.tight_layout()
plt.savefig("Excess Chemical Potential.png", dpi=300)
plt.show()

#%%
import numpy as np
import matplotlib.pyplot as plt

# Plot 1: Full Range
fig, ax1 = plt.subplots(1, 1, figsize=(7, 5))

for T in sorted_temps:
    data_list = sorted(results[T], key=lambda x: x["rho"])
    
    rhos = np.array([d["rho"] for d in data_list])
    mu_ex = np.array([d["Mu"] for d in data_list]) # Assuming your 'Mu' key is mu_ex
    mu_errs = np.array([d["Mu_err"] for d in data_list])
    
    # Calculate Ideal Gas part: mu_ig = T * ln(rho)
    # Note: We skip rho=0 to avoid log(0) issues
    mu_ig = T * np.log(rhos)
    
    # Total Chemical Potential
    mu_tot = mu_ig + mu_ex
    
    # Plot Total Chemical Potential (Markers + Solid Line)
    line = ax1.errorbar(rhos, mu_tot, yerr=mu_errs, fmt='-s', capsize=3, label=f"$T^* = {T}$ (Total)")
    
    # Plot Ideal Gas Reference (Dashed Line, same color as the total)
    ax1.plot(rhos, mu_ig, '--', color=line[0].get_color(), alpha=0.6, label=f"$T^* = {T}$ (Ideal)")

# Formatting
ax1.set_xlabel(r"$\rho \sigma^3$")
ax1.set_ylabel(r"$\mu / \epsilon$")
ax1.grid(True, linestyle='--', alpha=0.7)
# Placing legend outside or shrinking it might be needed if it gets too crowded
ax1.legend(fontsize='small', ncol=2)

plt.tight_layout()
plt.savefig("Total_Chemical_Potential_Full.png", dpi=300)
plt.show()

# Plot 2: Zoomed Range (rho < 0.85)
fig, ax1 = plt.subplots(1, 1, figsize=(7, 5))

for T in sorted_temps:
    data_list = sorted(results[T], key=lambda x: x["rho"])
    
    rhos = np.array([d["rho"] for d in data_list])
    mu_ex = np.array([d["Mu"] for d in data_list])
    mu_errs = np.array([d["Mu_err"] for d in data_list])
    
    # Masking
    mask = rhos < 0.85
    r_m = rhos[mask]
    
    mu_ig = T * np.log(r_m)
    mu_tot = mu_ig + mu_ex[mask]
    
    # Plotting
    line = ax1.errorbar(r_m, mu_tot+mu_ig, yerr=mu_errs[mask], fmt='-s', capsize=3, label=f"$T^* = {T}$")


ax1.set_xlabel(r"$\rho \sigma^3$")
ax1.set_ylabel(r"$\mu / \epsilon$")
ax1.grid(True, linestyle='--', alpha=0.7)
ax1.legend()

plt.tight_layout()
plt.savefig("Total_Chemical_Potential_Zoom.png", dpi=300)
plt.show()