"""
Created on Thu May 21 12:53:39 2026

@author: ime_rasenberg
"""

import numpy as np
import matplotlib.pyplot as plt
import os

def plot_spin_lattice(filepath="Data/Spin_orientiation.txt", subsample=1):
    """
    Reads the spin state from the C output and plots the 2D plane with 3D arrows.
    The arrows are colored by their Z-component to help visualize out-of-plane tilting.
    """
    if not os.path.exists(filepath):
        print(f"Error: Could not find the file at '{filepath}'.")
        print("Make sure you have run your C program and it successfully saved the data.")
        return

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


def Energy(filepath="Data/Energy_steps.txt"):
    """
    Reads the spin state from the C output and plots the 2D plane with 3D arrows.
    The arrows are colored by their Z-component to help visualize out-of-plane tilting.
    """
    filepath="Data/Energy_steps.txt"
    
    if not os.path.exists(filepath):
        print(f"Error: Could not find the file at '{filepath}'.")
        print("Make sure you have run your C program and it successfully saved the data.")
        # return
    
    data = np.genfromtxt(filepath, skip_header=1)
    
    # Convert to numpy arrays and reshape back to (N, N)
    # Using float64 to match C's double precision
    steps = np.array(data[:,0])
    Energy = np.array(data[:,1])
    acceptance = np.array(data[:,2])
    # S_z = np.array(z_lines, dtype=np.float64).reshape((N, N))
    
      
    # 4. Plotting
    plt.figure(figsize=(10, 8))
    
    plt.plot(steps, Energy)
    
    plt.show()
    plt.figure(figsize=(10, 8))
    
    plt.plot(steps, acceptance)
    
    plt.show()



if __name__ == "__main__":
    # Adjust 'subsample' to see more or fewer arrows. 
    # subsample=1 shows every single spin.
    plot_spin_lattice(filepath="Data/Spin_orientiation.txt", subsample=1)
    Energy()