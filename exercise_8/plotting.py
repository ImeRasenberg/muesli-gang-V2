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
    data = pd.read_csv('measurements.dat', sep='\t', names=['Step', 'Pressure', 'Mu_Excess'])
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