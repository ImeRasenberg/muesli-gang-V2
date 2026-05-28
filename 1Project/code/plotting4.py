import os
import json

data_folder = "Data"

master_dict = {}
decoder = json.JSONDecoder()

# ==========================================================
# loop through all json files
# ==========================================================
filepaths = [
    os.path.join(data_folder, f)
    for f in os.listdir(data_folder)
    if os.path.isfile(os.path.join(data_folder, f))
    and f.endswith(".json")
]

print(f"Found {len(filepaths)} JSON files")

for filepath in filepaths:

    filename = os.path.basename(filepath)

    # ------------------------------------------------------
    # scrape metadata dynamically from filename
    # ------------------------------------------------------
    name_no_ext = os.path.splitext(filename)[0]
    parts = name_no_ext.split("__")

    metadata = {}

    for part in parts:
        if "_" in part:
            key, value = part.split("_", 1)
            metadata[key] = value

    # skip malformed filenames
    if not all(k in metadata for k in ["D", "Hz", "I"]):
        print(f"Skipping {filename}")
        continue

    try:
        J = float(metadata["D"])
        Hz = float(metadata["Hz"])
        I = int(metadata["I"])
    except ValueError:
        print(f"Malformed metadata in {filename}")
        continue

    print(f"Loading {filename}")

    # ------------------------------------------------------
    # initialize nested dictionary
    # master_dict[J][Hz][I]
    # ------------------------------------------------------
    if J not in master_dict:
        master_dict[J] = {}

    if Hz not in master_dict[J]:
        master_dict[J][Hz] = {}

    master_dict[J][Hz][I] = {}

    # ------------------------------------------------------
    # load concatenated json objects
    # ------------------------------------------------------
    with open(filepath, "r") as f:
        content = f.read().strip()

    pos = 0

    while pos < len(content):

        while pos < len(content) and content[pos].isspace():
            pos += 1

        if pos >= len(content):
            break

        try:
            obj, idx = decoder.raw_decode(content[pos:])
            pos += idx

            step = obj["step"]

            # store original object
            master_dict[J][Hz][I][step] = obj

        except json.JSONDecodeError as e:
            print(f"Error in {filename} near position {pos}: {e}")
            break

print("Finished loading.")

#%%

import numpy as np
import matplotlib.pyplot as plt


N_LAST = 30

# -------------------------------------------------
# collect unique parameter values
# -------------------------------------------------
D_vals = sorted(master_dict.keys())

Hz_vals = sorted({
    hz
    for D in master_dict
    for hz in master_dict[D]
})

I_vals = sorted({
    I
    for D in master_dict
    for hz in master_dict[D]
    for I in master_dict[D][hz]
})


for I in I_vals:

    # matrices for plotting
    plus_grid = np.full((len(Hz_vals), len(D_vals)), np.nan)
    minus_grid = np.full((len(Hz_vals), len(D_vals)), np.nan)

    for ix, D in enumerate(D_vals):
        for iy, Hz in enumerate(Hz_vals):

            if Hz not in master_dict[D]:
                continue

            if I not in master_dict[D][Hz]:
                continue

            step_dict = master_dict[D][Hz][I]

            if len(step_dict) == 0:
                continue

            # --------------------------------------
            # largest N_LAST steps
            # --------------------------------------
            largest_steps = sorted(step_dict.keys())[-N_LAST:]

            nplus_lengths = []
            nminus_lengths = []

            for step in largest_steps:

                obj = step_dict[step]

                if "N+" in obj:
                    nplus_lengths.append(len(obj["N+"]))

                if "N-" in obj:
                    nminus_lengths.append(len(obj["N-"]))

            # average over last N steps
            if len(nplus_lengths) > 0:
                plus_grid[iy, ix] = np.mean(nplus_lengths)

            if len(nminus_lengths) > 0:
                minus_grid[iy, ix] = np.mean(nminus_lengths)

    # ==================================================
    # plotting
    # ==================================================
    fig, ax = plt.subplots(figsize=(8, 6))

    im = ax.imshow(
        plus_grid,
        origin="lower",
        aspect="auto",
        extent=[
            min(D_vals),
            max(D_vals),
            min(Hz_vals),
            max(Hz_vals),
        ]
    )

    plt.colorbar(im, ax=ax, label="mean len(N+)")

    ax.set_xlabel("D")
    ax.set_ylabel("Hz")
    ax.set_title(f"N+ size (I={I}) averaged over last {N_LAST} steps")

    plt.show()

    fig, ax = plt.subplots(figsize=(8, 6))

    im = ax.imshow(
        minus_grid,
        origin="lower",
        aspect="auto",
        extent=[
            min(D_vals),
            max(D_vals),
            min(Hz_vals),
            max(Hz_vals),
        ]
    )

    plt.colorbar(im, ax=ax, label="mean len(N-)")

    ax.set_xlabel("D")
    ax.set_ylabel("Hz")
    ax.set_title(f"N- size (I={I}) averaged over last {N_LAST} steps")

    plt.show()
    
    fig, ax = plt.subplots(figsize=(8, 6))

    im = ax.imshow(
        minus_grid-plus_grid,
        origin="lower",
        aspect="auto",
        extent=[
            min(D_vals),
            max(D_vals),
            min(Hz_vals),
            max(Hz_vals),
        ]
    )

    plt.colorbar(im, ax=ax, label="mean len(N-)")

    ax.set_xlabel("D")
    ax.set_ylabel("Hz")
    ax.set_title(f"N- - N+ size (I={I}) averaged over last {N_LAST} steps")

    plt.show()

