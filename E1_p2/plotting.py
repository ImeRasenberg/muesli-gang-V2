import matplotlib.pyplot as plt
import glob
import numpy as np

files = glob.glob("data/energy_vs_time_dt_*.txt")

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
            energies.append(float(parts[2])+float(parts[1]))
    
    steps_arr = np.array(steps)
    energies_arr = np.array(energies)
    
    if len(steps_arr) > 1:
        # Calculate dt and normalized steps
        dt = np.average(np.diff(steps_arr))
        normalized_steps = steps_arr / dt
        plot1_data.append((dt, normalized_steps, energies_arr))
        
        # max_var = abs(max(energies_arr[1000:])- min(energies_arr[1000:]))
        if (dt>0.9e-3):
            max_var = np.std(energies_arr)
        else:
            max_var = np.std(energies_arr[500:])
        
        
        dt_values.append(dt)
        abs_sum_values.append(max_var)

# --- FIGURE 1: Energy Evolution ---
plt.figure(figsize=(10, 6))

# Sort by dt descending (Big to Small)
plot1_data.sort(key=lambda x: x[0], reverse=True)

for dt, x_vals, y_vals in plot1_data:
    if (dt>0.9e-3) and (dt<1.1e-3):
        plt.plot(x_vals, y_vals, label=f"dt = {dt:.1e}")

plt.xlabel(r'Unitles time $\tau$')
plt.ylabel(r'$E_{tot}/\beta$')
# plt.title('Energy Evolution per Iteration Step')
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()

# --- FIGURE 2: Derivative Analysis (Log-Log) ---
plt.figure(figsize=(10, 6))

dt_values = np.array(dt_values)
abs_sum_values = np.array(abs_sum_values)

# Sort by dt ascending for a clean line connection
sort_idx = np.argsort(dt_values)

plt.loglog(dt_values[sort_idx], abs_sum_values[sort_idx], 'o-', color='tab:red', markersize=8)

plt.xlabel('Time Step ($dt$)')
plt.ylabel(r'$\Delta E_{tot}/\beta$')
# plt.title('Total Absolute Variation vs. Time Step (Log-Log Scale)')
plt.grid(True, which="both", ls="-", alpha=0.2)
plt.tight_layout()

# Show both windows
plt.show()
