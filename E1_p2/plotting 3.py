#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 07:43:27 2026

@author: ime_rasenberg
"""

import matplotlib.pyplot as plt
import glob
import numpy as np

files = glob.glob("data_3/velocity_correlation_*.txt")

if not files:
    print("No energy files found! Run your C program first.")
    

# Data containers
plot1_data = [] # List of (dt, normalized_steps, energies)
dt_values = []
abs_sum_values = []

for file in files:
    steps = []
    energies = []
    
    with open(file, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip(): 
                continue
            parts = line.split()
            # parts[0] is step, parts[1] is energy
            steps.append(float(parts[0]))
            energies.append(float(parts[1]))
    
    steps_arr = np.array(steps)
    energies_arr = np.array(energies)
    
    
    if len(steps_arr) > 1:
        # Calculate dt and normalized steps
        dt = np.average(np.diff(steps_arr))
        normalized_steps = steps_arr / dt
        plot1_data.append((dt, normalized_steps, energies_arr))
        

# --- FIGURE 1: Energy Evolution ---
plt.figure(figsize=(10, 6))

# Sort by dt descending (Big to Small)
plot1_data.sort(key=lambda x: x[0], reverse=True)

for dt, x_vals, y_vals in plot1_data:
    plt.plot(x_vals, y_vals, label=f"dt = {dt:.1e}")
    

plt.xlabel(r'Unitles time $\tau$')
plt.ylabel(r'$\langle v(0)v(t) \rangle$')
# plt.title('Energy Evolution per Iteration Step')
# plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()

y = y_vals.copy() * 0
count =0
for dt, x_vals, y_vals in plot1_data:
    count+=1
    y += y_vals

y /= count

plt.figure(figsize=(10, 6))

# Sort by dt descending (Big to Small)
plot1_data.sort(key=lambda x: x[0], reverse=True)

plt.plot(x_vals, y, label=f"dt = {dt:.1e}")
    

plt.xlabel(r'Unitles time $\tau$')
plt.ylabel(r'$\langle v(0)v(t) \rangle$')
# plt.title('Energy Evolution per Iteration Step')
# plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()

