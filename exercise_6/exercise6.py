


import numpy as np
import matplotlib.pyplot as plt
import os


all_files = os.listdir('.')
file_list = sorted([f for f in all_files if f.startswith('gr') and f.endswith('.dat')])

if not file_list:
    raise FileNotFoundError("No files matching 'gr*.dat' found in the current directory.")

print(f"Found {len(file_list)} files: {file_list}")


first_data = np.genfromtxt(file_list[0], skip_header=1)
r = first_data[:, 0]
n_points = len(r)


g_all = np.zeros((len(file_list), n_points))


for i, filename in enumerate(file_list):
    data = np.genfromtxt(filename, skip_header=1)
    g_all[i, :] = data[:, 1]

g_mean = np.mean(g_all, axis=0)
g_std = np.std(g_all, axis=0, ddof=1)  



# Plot the averaged radial distribution function with shaded error region
plt.figure(figsize=(8, 5))
plt.plot(r, g_mean, color='black', label='g(r)')
plt.errorbar(r, g_mean, yerr=g_std, fmt='k.')
plt.axhline(1.0, color='red', linestyle='--', alpha=0.5, label='Ideal Gas Limit (1.0)')
plt.xlabel(r'$r/ \sigma$')
plt.ylabel('g(r)')
plt.xlim(0, np.max(r))
plt.ylim(0, np.max(g_mean + g_std) * 1.05) 
plt.legend()
plt.grid(True, linestyle=':', alpha=0.7)
plt.savefig('gr_fcc.png', dpi=300)
plt.show()

import numpy as np
import matplotlib.pyplot as plt



data = np.genfromtxt('gr_fcc.dat', skip_header=1)
r = data[:, 0]
g = data[:, 1]


plt.figure(figsize=(8, 5))
plt.plot(r, g, color='black', label='g(r)')
plt.axhline(1.0, color='red', linestyle='--', alpha=0.5, label='Ideal Gas Limit (1.0)')
plt.xlabel(r'$r/ \sigma$')
plt.ylabel('g(r)')
plt.xlim(0, np.max(r))
plt.ylim(0, np.max(g) * 1.05) 
plt.legend()
plt.grid(True, linestyle=':', alpha=0.7)
plt.savefig('gr_fcc.png', dpi=300)
plt.show()
