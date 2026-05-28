import json
import matplotlib.pyplot as plt

json_path = "Data/skyrmions.json"
output_plot_path = "skyrmions_plot.png"

# Read the raw file content
with open(json_path, "r") as f:
    content = f.read().strip()

# 1. CONVERT ALL RAW OBJECTS INTO ONE BIG DICTIONARY
# We will use the 'step' as the key so everything is combined into one object
full_dataset_dict = {}

decoder = json.JSONDecoder()
pos = 0

while pos < len(content):
    while pos < len(content) and content[pos].isspace():
        pos += 1
    if pos >= len(content):
        break
        
    try:
        obj, idx = decoder.raw_decode(content[pos:])
        pos += idx  
        
        # Store the entire original object inside our big dictionary
        step_key = obj["step"]
        full_dataset_dict[step_key] = obj
        
    except json.JSONDecodeError as e:
        print(f"Error parsing JSON near position {pos}: {e}")
        break

# 2. READ OUT AND PROCESS DATA FROM THE TOTAL DICTIONARY
steps = []
n_plus_counts = []
n_minus_counts = []

# Sort by step keys to ensure the plot is chronological
for step in sorted(full_dataset_dict.keys()):
    data_entry = full_dataset_dict[step]
    
    steps.append(step)
    n_plus_counts.append(len(data_entry["N+"]))
    n_minus_counts.append(len(data_entry["N-"]))


# --- Generate the plot ---
fig, ax = plt.subplots(figsize=(10, 5))

ax.plot(steps, n_plus_counts, label=r'$N_+$ (Positive Peaks)', color='#d62728', linewidth=2, marker='o', markersize=4)
ax.plot(steps, n_minus_counts, label=r'$N_-$ (Negative Peaks)', color='#1f77b4', linewidth=2, marker='x', markersize=4)

# Format labels, title, and layout
ax.set_xlabel("Simulation Steps", fontsize=12)
ax.set_ylabel("Peak Count", fontsize=12)
ax.set_title("Number of Skyrmion Peaks ($N_+$ and $N_-$) vs. Simulation Steps", fontsize=14, fontweight='bold')
ax.grid(True, linestyle='--', alpha=0.6)
ax.legend(fontsize=11, loc='upper right')

plt.tight_layout()

# Save the resulting visualization
plt.savefig(output_plot_path, dpi=300)
print(f"Successfully processed {len(steps)} data points from the total dictionary. Plot saved to '{output_plot_path}'.")

#%%

import json
import numpy as np
import matplotlib.pyplot as plt
import imageio.v2 as imageio

json_path = "Data/skyrmions.json"
gif_path = "skyrmions_evolution.gif"

# ----------------------------
# PARAMETERS
# ----------------------------
N = 40
sigma = 3
grid_x, grid_y = np.meshgrid(np.arange(N), np.arange(N), indexing="ij")

def gaussian(x, y, x0, y0, sigma):
    return np.exp(-((x - x0)**2 + (y - y0)**2) / (2 * sigma**2))

# ----------------------------
# PARSE STREAMED JSON OBJECTS
# ----------------------------
with open(json_path, "r") as f:
    content = f.read().strip()

decoder = json.JSONDecoder()
pos = 0
data = {}

while pos < len(content):
    while pos < len(content) and content[pos].isspace():
        pos += 1
    if pos >= len(content):
        break

    obj, idx = decoder.raw_decode(content[pos:])
    pos += idx
    data[obj["step"]] = obj

steps = sorted(data.keys())

# ----------------------------
# BUILD FRAMES
# ----------------------------
frames = []

for step in steps:
    if step<600000:
        continue
    
    entry = data[step]
    

    field = np.zeros((N, N), dtype=float)

    # Positive skyrmions
    for _, (y,x) in entry["N+"].items():
        field += gaussian(grid_x, grid_y, x, y, sigma)

    # Negative skyrmions
    for _, (y,x) in entry["N-"].items():
        field -= gaussian(grid_x, grid_y, x, y, sigma)

    vmax = np.max(np.abs(field)) if np.max(np.abs(field)) > 0 else 1.0

    fig, ax = plt.subplots(figsize=(6, 5), dpi=70)

    im = ax.imshow(
        field.T,
        cmap="seismic",
        origin="lower",
        vmin=-vmax,
        vmax=vmax
    )

    ax.set_title(f"Step {step} | Q = {entry['Q']}")
    ax.set_xlabel("x")
    ax.set_ylabel("y")

    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label("Topological density (arb. units)")

    plt.tight_layout()

    # Convert figure to image array
    fig.canvas.draw()

    w, h = fig.canvas.get_width_height()
    img = np.frombuffer(fig.canvas.buffer_rgba(), dtype=np.uint8)
    img = img.reshape((h, w, 4))  # RGBA
    
    img = img[:, :, :3]  # drop alpha channel

    frames.append(img)
    plt.close(fig)

# ----------------------------
# SAVE GIF
# ----------------------------
imageio.mimsave(gif_path, frames, fps=5)

print(f"GIF saved to {gif_path} with {len(frames)} frames.")
