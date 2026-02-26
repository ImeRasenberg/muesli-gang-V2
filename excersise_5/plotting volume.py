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