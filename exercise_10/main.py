import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter

# load data
data = np.loadtxt("sudoku.dat")

# columns
step = data[:,0]
energy = data[:,1]
temperature = data[:,2]
energy_smooth = savgol_filter(energy, window_length=150, polyorder=3)

# plot energy vs step
plt.figure()
plt.plot(step, energy, 'k-', label = 'energy')
plt.plot(step, energy_smooth, 'r-', label="smoothed energy")
plt.xlabel("Monte Carlo step")
plt.ylabel(r"$ \frac{1}{ \beta_0} $ Energy (number of conflicts)")
plt.xlim((0, step[-1]*1.02))
plt.ylim((0, energy[0]*1.02))
plt.grid()
plt.legend()
plt.savefig("sudoku.png", dpi=500)
plt.show()

plt.figure()
plt.plot(step, temperature, 'b-')
plt.xlabel("Monte Carlo step")
plt.ylabel(r"T (temprature)")
plt.xlim((0, step[-1]*1.02))
plt.ylim((temperature[-1]*0.98, temperature[0]*1.02))
plt.grid()
plt.savefig("sudoku_temp.png", dpi=500)
plt.show()
