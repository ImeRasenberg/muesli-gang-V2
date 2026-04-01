#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 12:03:27 2026

@author: ime_rasenberg
"""

import matplotlib.pyplot as plt
import glob
import numpy as np

files = glob.glob("data_2/energy_vs_time_dt_temp_*.txt")

if not files:
    print("No energy files found! Run your C program first.")


# Data containers
plot1_data = [] # List of (dt, normalized_steps, energies)
dt_values = []
average = []

for file in files:
    steps = []
    energies = []
    
    with open(file, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip(): 
                continue
            parts = line.split()
            # parts[0] is step, parts[1] is energy
            steps.append(float(parts[3]))
            energies.append(float(parts[1])+float(parts[2]))
    
    steps_arr = np.array(steps)
    energies_arr = np.array(energies)
    
    if len(steps_arr) > 1:
        # Calculate dt and normalized steps
        dt = np.average(steps_arr)
        normalized_steps = steps_arr / dt
        plot1_data.append((dt, normalized_steps, energies_arr))
        
        max_var = np.average(energies_arr)
        
        dt_values.append(dt)
        average.append(max_var)

# Convert to numpy arrays (good practice for plotting)
betas = np.array(dt_values)
average = np.array(average)

# Sort by dt (important for clean lines)
sort_idx = np.argsort(dt_values)
betas = betas[sort_idx]
average = average[sort_idx]

# Log-log plot
plt.figure()
plt.scatter(1/betas * 1.5, average*1.5, marker='o', linestyle='-')

plt.xlabel(r' $\beta^{-1}/\epsilon$')
plt.ylabel(r'Average Energy/$\epsilon$')
# plt.title('Log-Log Plot of Average Energy vs dt')
plt.grid(True, which="both", ls="--")

plt.show()


# Linear fit: y = a*x + b
x = 1/betas
y = average

coeffs = np.polyfit(x, y, 1)  # degree 1 polynomial
a, b = coeffs

print(f"Slope (a) = {a}")
print(f"Intercept (b) = {b}")

# Plot fitted line
x_fit = np.linspace(min(x), max(x), 200)
y_fit = a * x_fit + b

plt.plot(x_fit, y_fit, label=f'fit: y={a:.3e}x + {b:.3e}')
plt.legend()

