import matplotlib.pyplot as plt
import glob
import numpy as np

def plot_analysis():
    files = glob.glob("data/energy_vs_time_dt_*.txt")
    
    if not files:
        print("No energy files found! Run your C program first.")
        return

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
                energies.append(float(parts[2]))
        
        steps_arr = np.array(steps)
        energies_arr = np.array(energies)
        
        if len(steps_arr) > 1:
            # Calculate dt and normalized steps
            dt = np.average(np.diff(steps_arr))
            normalized_steps = steps_arr / dt
            plot1_data.append((dt, normalized_steps, energies_arr))
            
            # Calculate absolute sum of derivatives: sum| (E_n+1 - E_n) / dt |
            delta_e = np.diff(energies_arr)
            # Dividing by dt gives the derivative dE/dt
            derivative = delta_e
            total_variation = np.sum(np.abs(derivative))
            
            dt_values.append(dt)
            abs_sum_values.append(total_variation)

    # --- FIGURE 1: Energy Evolution ---
    plt.figure(figsize=(10, 6))
    
    # Sort by dt descending (Big to Small)
    plot1_data.sort(key=lambda x: x[0], reverse=True)

    for dt, x_vals, y_vals in plot1_data:
        plt.plot(x_vals, y_vals, label=f"dt = {dt:.1e}")

    plt.xlabel(r'Unitles time $\tau$')
    plt.ylabel(r'$E_{pot}/\beta$')
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
    plt.ylabel(r'Sum of Absolute Derivatives potential energy $\sum |\frac{\Delta E}{\Delta t}|/\beta$')
    # plt.title('Total Absolute Variation vs. Time Step (Log-Log Scale)')
    plt.grid(True, which="both", ls="-", alpha=0.2)
    plt.tight_layout()

    # Show both windows
    plt.show()

if __name__ == "__main__":
    plot_analysis()