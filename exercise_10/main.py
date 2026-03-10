import numpy as np
import matplotlib.pyplot as plt

# load data
data = np.loadtxt("sudoku.dat")

# columns
step = data[:,0]
energy = data[:,1]
temperature = data[:,2]

alpha = energy[0]/step[-1]

# plot energy vs step
plt.figure()
plt.plot(step, energy, 'k-')
#plt.plot(step, energy[0]- alpha*step, 'r-')
plt.xlabel("Monte Carlo step")
plt.ylabel("Energy (number of conflicts)")
plt.xlim((0, step[-1]*1.02))
plt.ylim((0, energy[0]*1.02))
plt.grid()
plt.savefig("sudoku.png", dpi=500)
plt.show()

