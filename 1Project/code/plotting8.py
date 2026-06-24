#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 18 10:43:04 2026

@author: ime_rasenberg
"""

import os
import json

data_folder = "Data"


here = "../documentation/sections/Graphs/"


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
from collections import defaultdict

# ==========================================================
# CONFIG — set the J and Hz you want to plot
# (script picks the closest match present in master_dict)
# ==========================================================
TARGET_J  = 1.5
TARGET_Hz = 1.1

# ----------------------------------------------------------
# Find closest (J, Hz) combination in master_dict
# ----------------------------------------------------------
all_keys = [
    (J, Hz)
    for J in master_dict
    for Hz in master_dict[J]
]

closest = min(
    all_keys,
    key=lambda k: (k[0] - TARGET_J) ** 2 + (k[1] - TARGET_Hz) ** 2,
)
J_sel, Hz_sel = closest
print(f"Plotting J={J_sel}, Hz={Hz_sel}")

# ----------------------------------------------------------
# Collect (beta -> list of values) over all I and steps
# ----------------------------------------------------------
beta_data = defaultdict(list)

for I, steps in master_dict[J_sel][Hz_sel].items():
    for step, obj in steps.items():
        beta = obj["beta"]
        beta_data[beta].append(obj["Q"])

betas  = np.array(sorted(beta_data.keys()))
means  = np.array([np.mean(beta_data[b]) for b in betas])
stds   = np.array([np.std(beta_data[b], ddof=1) if len(beta_data[b]) > 1
                   else 0.0
                   for b in betas])
counts = np.array([len(beta_data[b]) for b in betas])

# ----------------------------------------------------------
# Plot
# ----------------------------------------------------------
fig, ax = plt.subplots(figsize=(8, 5))

ax.errorbar(
    1 / betas, means, yerr=stds,
    marker="o", linewidth=1.8, markersize=5,
    capsize=4, capthick=1.2, elinewidth=1.2,
    color="#3266ad", label="mean ± 1 std",
)

ax.set_xlabel(r"$J\beta^{-1}$", fontsize=18)
ax.set_ylabel(r"$Q$", fontsize=18)
# ax.set_title(
#     rf"$Q$ vs $\beta^{{-1}}$  (J={J_sel}, Hz={Hz_sel})",
#     fontsize=13,
# )
ax.legend(fontsize=11)
ax.grid(True, linestyle="--", linewidth=0.5, alpha=0.6)

# # annotate sample counts (keyed on 1/beta x-position)
# for b, n in zip(betas, counts):
#     ax.annotate(
#         f"n={n}",
#         xy=(1 / b, ax.get_ylim()[1]),
#         xytext=(0, 4),
#         textcoords="offset points",
#         ha="center", va="bottom",
#         fontsize=7, color="gray",
#     )

plt.tight_layout()
plt.savefig(here + "QT.png", dpi=150)
plt.show()
print("Saved Q_vs_beta.png")



#%%

import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict

critical_points = []

for D in master_dict:
    for Hz in master_dict[D]:

        # ------------------------------------------
        # collect averages vs beta
        # ------------------------------------------
        beta_data = defaultdict(list)

        for I, steps in master_dict[D][Hz].items():
            for step, obj in steps.items():
                beta_data[obj["beta"]].append(obj["Q"])

        if len(beta_data) < 2:
            continue

        betas = np.array(sorted(beta_data.keys()))

        means = np.array([
            np.mean(beta_data[b])
            for b in betas
        ])

        # ------------------------------------------
        # reference value:
        # largest beta = lowest temperature
        # ------------------------------------------
        Q0 = abs(means[-1])

        if Q0 == 0:
            continue

        threshold = 0.3 * Q0

        # ------------------------------------------
        # find first crossing
        # ------------------------------------------
        crossing = None

        for beta, q in zip(betas[::-1], means[::-1]):
            if abs(q) < threshold:
                crossing = D / beta      # J β^{-1}
                break

        if crossing is not None:
            critical_points.append((D, Hz, crossing))


D_vals = sorted({x[0] for x in critical_points})
Hz_vals = sorted({x[1] for x in critical_points})

mtx = np.full((len(Hz_vals), len(D_vals)), np.nan)

for D, Hz, crit in critical_points:
    i = Hz_vals.index(Hz)
    j = D_vals.index(D)

    mtx[i, j] = crit
#%%
mask = np.loadtxt("mask.txt", dtype=int).astype(bool)

fig, ax = plt.subplots(figsize=(8, 6))

im = ax.imshow(
    mtx,
    origin="lower",
    aspect="auto",
    extent=[
        min(D_vals),
        max(D_vals),
        min(Hz_vals),
        max(Hz_vals),
    ]
)

ax.set_xlabel(r"$D/J$")
ax.set_ylabel(r"$H/J$")
# ax.set_title(r"$T_c^*$")

cbar = plt.colorbar(im, ax=ax)
cbar.set_label(r"$T_c^*$")

plt.tight_layout()
plt.show()

#%%

mask = np.loadtxt("mask.txt", dtype=int).astype(bool)

# Hide values where mask is False
mtx_masked = np.ma.masked_where(~mask, mtx)

fig, ax = plt.subplots(figsize=(8, 6))

im = ax.imshow(
    mtx_masked,
    origin="lower",
    aspect="auto",
    extent=[
        min(D_vals),
        max(D_vals),
        min(Hz_vals),
        max(Hz_vals),
    ]
)

ax.set_xlabel(r"$D/J$",size=18)
ax.set_ylabel(r"$H/J$",size =18)

cbar = plt.colorbar(im, ax=ax)
cbar.set_label(r"$T_c^*$",size =18)

plt.tight_layout()
plt.savefig(here + "Tc2.png")
plt.show()