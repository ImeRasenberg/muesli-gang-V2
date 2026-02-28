#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 26 12:46:56 2026

@author: ime_rasenberg
"""

# Import required libraries
import os
import numpy as np
import matplotlib.pyplot as plt

base_path = "./data"

data_dict = {}

# Iterate over subfolders
for subfolder in os.listdir(base_path):
    sub_path = os.path.join(base_path, subfolder)
    info_path = os.path.join(sub_path, "info.dat")
    
    if os.path.isdir(sub_path) and os.path.isfile(info_path):
        with open(info_path, "r") as f:
            lines = f.readlines()
        
        # First line contains metadata
        header = list(map(float, lines[0].split()))
        delta, dV_m, betaP, n_particles, size_average = header
        
        # Remaining lines contain the inf matrix
        data = []
        for line in lines[1:]:
            values = list(map(float, line.split()))
            if any(values):  # skip the last zero-line
                data.append(values)
        
        data = np.array(data)
        
        data_dict[subfolder] = {
            "n_particles": int(n_particles),
            "delta": delta,
            "dV_m": dV_m,
            "betaP": betaP,
            "size_average": size_average,
            "step": data[:, 0],
            "volume": data[:, 1],
            "acceptance_vol": data[:, 2],
            "converged_vol": data[:, 3],
            "acceptance_move": data[:, 4],
            "converged_move": data[:, 5],
        }

# # Plot volume vs step for the first subfolder (sorted for consistency)
# first_key = sorted(data_dict.keys())[0]
# first_data = data_dict[first_key]

# plt.figure()
# plt.plot(first_data["step"], first_data["volume"])
# plt.xlabel("Step")
# plt.ylabel("Volume")
# plt.title(f"Volume vs Step ({first_key})")
# plt.show()

# data_dict

plt.figure()

for key, dataset in sorted(data_dict.items()):
    
    betaP = dataset["betaP"]
    sigma = dataset["size_average"]
    
    label_value = betaP * sigma**3
    
    plt.plot(
        dataset["step"] / sigma**3,
        dataset["volume"],
        label=rf'$\beta P \sigma^3 = {label_value:.3f}$'
    )

plt.xlabel(r'Step Number')
plt.ylabel(r'$V/\sigma^3$')

plt.legend()
plt.tight_layout()
plt.show()

# This code assumes `data_dict` is already defined in memory
# with the structure previously constructed.

import numpy as np
import matplotlib.pyplot as plt

x_values = []
y_values = []

for key, dataset in sorted(data_dict.items()):
    
    betaP = dataset["betaP"]
    sigma = dataset["size_average"]
    
    # Compute βP σ^3 (legend quantity)
    x_val = betaP * sigma**3
    
    # Select steps >= 15000
    steps = dataset["step"]
    volumes = dataset["volume"]
    
    Volume_particle = (sigma/2)**3 * 4/3* np.pi*dataset["n_particles"]
    
    mask = steps >= 15000
    avg_volume = np.mean(volumes[mask])
    
    x_values.append(x_val)
    y_values.append(Volume_particle / avg_volume)

# Convert to arrays and sort by x
x_values = np.array(x_values)
y_values = np.array(y_values)

order = np.argsort(x_values)
x_values = x_values[order]
y_values = y_values[order]

# Single plot (no custom colors, no subplots)
plt.figure()
plt.plot(x_values, y_values)

plt.xlabel(r'$\beta P \sigma^3$')
plt.ylabel(r'$\langle \rho \rangle_{eq} $')
# the average of the volume after 15000 steps as a fu

plt.tight_layout()
plt.show()

# Return numerical values in case you want them
list(zip(x_values, y_values))


import numpy as np
import matplotlib.pyplot as plt

# ... [Keep your existing data extraction code here] ...

# 1. Generate the Carnahan-Starling Theoretical Curve
eta_theory = np.linspace(0.1, 0.7, 100) # Theory is accurate up to ~0.5
# CS Equation rearranged for x-axis (beta * P * sigma^3)
# Note: rho*sigma^3 = 6*eta/pi
x_theory = (1 + eta_theory + eta_theory**2 - eta_theory**3) / (1 - eta_theory)**3

# 2. Plotting
plt.figure(figsize=(8, 5))

# Your Simulation Data
plt.plot(x_values, y_values, 'bo-', label='MC Simulation')

# Carnahan-Starling Theory
plt.plot(x_theory, eta_theory, 'r--', label='Carnahan-Starling Theory')

plt.xlabel(r'$\beta P \sigma^3$')
plt.ylabel(r'Packing Fraction $\eta$')
plt.title('Hard Sphere Equation of State')
plt.legend()
plt.grid(True, alpha=0.3)

plt.tight_layout()
plt.show()

