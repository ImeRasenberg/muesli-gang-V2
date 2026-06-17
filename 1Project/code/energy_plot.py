import numpy as np
import matplotlib.pyplot as plt

# Read the energy log file
data = []
with open("data2/energy_log.txt", "r") as f:
    for line in f:
        if line.startswith("#"):
            continue          # skip comment/header lines
        parts = line.strip().split()
        if len(parts) == 3:
            step, beta, energy = parts
            data.append((int(step), float(beta), float(energy)))

# Convert to numpy arrays for easier plotting
data = np.array(data, dtype=[("step", int), ("beta", float), ("energy", float)])
steps = data["step"]
betas = data["beta"]
energies = data["energy"]

plt.figure(figsize=(12, 6))
plt.plot(steps, energies)
plt.xlabel("Monte Carlo step")
plt.ylabel(r"Energy")
plt.title("Energy evolution during Monte Carlo simulation")
plt.grid(True, alpha=0.3)
plt.tight_layout()
#plt.savefig("energy_evolution.png", dpi=150)
plt.show()


energies = energies *betas

# Plot energy vs step, colored by beta
plt.figure(figsize=(12, 6))
plt.plot(steps, energies)
plt.xlabel("Monte Carlo step")
plt.ylabel(r"Energy $\beta$")
plt.title("Energy evolution during Monte Carlo simulation")
plt.grid(True, alpha=0.3)
plt.tight_layout()
#plt.savefig("energy_evolution.png", dpi=150)
plt.show()
