"""importing and plotting the g(r) data """

import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt('gr.dat', skip_header=1)


r = data[:, 0]
gr = data[:, 1]


plt.plot(r, gr, color = 'black', label='Radial Distribution Function')
plt.axhline(1.0, color='red', linestyle='--', alpha=0.5, label='Ideal Gas Limit (1.0)')
plt.xlabel('Distance (r)')
plt.ylabel('g(r)')
plt.title('Radial Distribution Function g(r)')
plt.ylim((0,np.max(gr)*1.05))
plt.xlim((0,np.max(r)))
plt.legend()
plt.grid(True, linestyle=':', alpha=0.7)


plt.savefig('gr_plot.png')
