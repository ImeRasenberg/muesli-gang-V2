#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 26 12:46:56 2026

@author: ime_rasenberg
"""

import os
import re
import matplotlib.pyplot as plt

# Path to the subfolder
data_dir = "data"
results = []

# Loop through files in the subfolder
for filename in os.listdir(data_dir):
    if filename.endswith(".dat"):
        
        # Extract number from filename (e.g., coords_step0001000.dat -> 1000)
        match = re.search(r'step(\d+)', filename)
        if match:
            step = int(match.group(1))
            file_path = os.path.join(data_dir, filename)
            
            with open(file_path, 'r') as f:
                # 1. Skip n_particles line (fscanf(read_cords, "%d\n", &n_particles))
                f.readline()
                
                # 2. Read the first dimension line (box[0])
                # Format: dummy box_value
                line = f.readline()
                if line:
                    parts = line.split()
                    if len(parts) >= 2:
                        box_0 = float(parts[1])
                        volume = box_0 ** 3
                        results.append((step, volume))

# Sort by step number so the plot is chronological
results.sort()

# Extract for plotting
steps, volumes = zip(*results)

# Plotting
plt.figure(figsize=(8, 5))
plt.plot(steps, volumes, marker='s', markersize=4, linestyle='-')
plt.xlabel('Step Number')
plt.ylabel('Volume')
plt.title('Volume vs. Step')
plt.grid(True)
plt.show()