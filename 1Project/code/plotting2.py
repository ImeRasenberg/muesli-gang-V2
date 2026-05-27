#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 26 14:30:58 2026

@author: ime_rasenberg
"""

import numpy as np
import matplotlib.pyplot as plt
import os

filepath="Data/Spin_orientiation.txt"
subsample=1

if not os.path.exists(filepath):
    print(f"Error: Could not find the file at '{filepath}'.")
    print("Make sure you have run your C program and it successfully saved the data.")
    
 # 1. Load data from file
with open(filepath, 'r') as f:
     N = int(f.readline().strip()) # First line is grid size N
     
     # Read the 3 rows corresponding to X, Y, and Z components
     x_lines = f.readline().strip().split('\t')
     y_lines = f.readline().strip().split('\t')
     z_lines = f.readline().strip().split('\t')

 # Convert to numpy arrays and reshape back to (N, N)
 # Using float64 to match C's double precision
S_x = np.array(x_lines, dtype=np.float64).reshape((N, N))
S_y = np.array(y_lines, dtype=np.float64).reshape((N, N))
S_z = np.array(z_lines, dtype=np.float64).reshape((N, N))

# 2. Create coordinates for the grid
X, Y = np.meshgrid(np.arange(N), np.arange(N))

 # 3. Subsample the data so the plot isn't overcrowded
 # (Showing 10,000 arrows at once makes it impossible to see individual directions)
X_sub = X[::subsample, ::subsample]
Y_sub = Y[::subsample, ::subsample]
U_sub = S_x[::subsample, ::subsample]
V_sub = S_y[::subsample, ::subsample]
W_sub = S_z[::subsample, ::subsample] # Z-component used for coloring

 # 4. Plotting
plt.figure(figsize=(10, 8))
 
 # quiver plots the 2D arrows (U, V). 
 # We map the W (Z-component) to a colormap to visualize the 3rd dimension!
quiver = plt.quiver(X_sub, Y_sub, U_sub, V_sub, W_sub, 
                     cmap='viridis', angles='xy', scale_units='xy', scale=0.8, pivot='mid')
 
 # Add a colorbar to show what the arrow colors mean
cbar = plt.colorbar(quiver)
cbar.set_label('Spin Z-component ($S_z$)', rotation=270, labelpad=15)

plt.title(f"2D Spin Lattice (Subsampled 1:{subsample})")
plt.xlabel("X Lattice Site")
plt.ylabel("Y Lattice Site")
plt.xlim(-1, N)
plt.ylim(-1, N)
plt.gca().set_aspect('equal') # Keep the grid square
plt.grid(True, which='both', linestyle=':', alpha=0.5)

plt.show()

#%%

filepath2="Data/Energy_steps.txt"

if not os.path.exists(filepath2):
    print(f"Error: Could not find the file at '{filepath}'.")
    print("Make sure you have run your C program and it successfully saved the data.")
    # return

data = np.genfromtxt(filepath2, skip_header=1)

# Convert to numpy arrays and reshape back to (N, N)
# Using float64 to match C's double precision
steps = np.array(data[:,0])
Energy = np.array(data[:,1])
acceptance = np.array(data[:,2])
beta = np.array(data[:,3])
Q = np.array(data[:, 4])
# S_z = np.array(z_lines, dtype=np.float64).reshape((N, N))

  
# 4. Plotting
plt.figure(figsize=(10, 8))
plt.plot(steps, Energy)
plt.show()

plt.figure(figsize=(10, 8))
plt.plot(steps, acceptance)
plt.show()

plt.figure(figsize=(10, 8))
plt.plot(steps, Q)
plt.show()

#%%

def get_triangle_charge(S1, S2, S3):
    """
    Calculates the signed solid angle (topological charge contribution) 
    for a triplet of spin vectors mapping onto a unit sphere surface.
    """
    # Triple product: S1 . (S2 x S3)
    num = np.sum(S1 * np.cross(S2, S3), axis=-1)
    
    # Denominator: 1 + S1.S2 + S2.S3 + S3.S1
    den = 1.0 + np.sum(S1*S2, axis=-1) + np.sum(S2*S3, axis=-1) + np.sum(S3*S1, axis=-1)
    
    # Omega = 2 * atan2(num, den)
    omega = 2.0 * np.arctan2(num, den)
    
    # Return charge density (Omega / 4*pi)
    return omega / (4.0 * np.pi)


# Package into an (N, N, 3) array to make handling spin vectors easier
spins = np.stack([S_x, S_y, S_z], axis=-1)

# 2. Prepare an array to hold the local charge for each grid site (i, j)
charge_density = np.zeros((N, N))

# 3. Loop over the lattice with periodic boundary conditions
for i in range(N):
    for j in range(N):
        ip = (i + 1) % N
        jp = (j + 1) % N

        # Get the 4 corners of the current square cell
        S_ij   = spins[i, j]
        S_ipj  = spins[ip, j]
        S_ijp  = spins[i, jp]
        S_ipjp = spins[ip, jp]

        # Triangle 1: (i,j) -> (i+1, j) -> (i, j+1)
        q1 = get_triangle_charge(S_ij, S_ipj, S_ijp)

        # Triangle 2: (i+1, j) -> (i+1, j+1) -> (i, j+1)
        q2 = get_triangle_charge(S_ipj, S_ipjp, S_ijp)

        # Assign the sum of both triangles to this local grid site coordinate
        charge_density[i, j] = q1 + q2

# 4. Total global wrapping number (Sum of all tiles)
total_Q = np.sum(charge_density)

print("=" * 50)
print(f"GLOBAL WRAPPING NUMBER (Q): {total_Q:+.5f}")
print(f"Estimated Skyrmion Count:   {round(abs(total_Q))}")
print("=" * 50)

# 5. Plotting the Local Density Map
plt.figure(figsize=(9, 7))

# imshow creates a clean pixel-by-pixel heat map of our grid calculations
# We use 'seismic' (Blue-White-Red) colormap because it perfectly flags 
# positive charge (red), negative charge (blue), and empty space (white).
vmax = np.max(np.abs(charge_density)) # Scale dynamically, baseline 0.1
img = plt.imshow(charge_density, cmap='seismic', origin='lower', 
                 vmin=-vmax, vmax=vmax, extent=[-0.5, N-0.5, -0.5, N-0.5])

cbar = plt.colorbar(img)
cbar.set_label('Local Topological Charge Density ($q_{ij}$)', rotation=270, labelpad=15)

plt.title(f"Local Topological Charge Density Map\nGlobal Wrapping Number $Q$ = {total_Q:+.4f}")
plt.xlabel("X Lattice Site")
plt.ylabel("Y Lattice Site")
plt.grid(False) # Turn off standard lines to clearly observe the pixels
plt.show()


from scipy.ndimage import gaussian_filter

smooth = gaussian_filter(charge_density, sigma=2, mode = "wrap")

plt.figure(figsize=(9, 7))

vmax2 = np.max(np.abs(smooth)) # Scale dynamically, baseline 0.1
img = plt.imshow(smooth, cmap='seismic', origin='lower', 
                 vmin=-vmax2, vmax=vmax2, extent=[-0.5, N-0.5, -0.5, N-0.5])

cbar = plt.colorbar(img)
cbar.set_label('Local Topological Charge Density ($q_{ij}$)', rotation=270, labelpad=15)

plt.title(f"Local Topological Charge Density Map\nGlobal Wrapping Number $Q$ = {total_Q:+.4f}")
plt.xlabel("X Lattice Site")
plt.ylabel("Y Lattice Site")
plt.grid(False) # Turn off standard lines to clearly observe the pixels
plt.show()




