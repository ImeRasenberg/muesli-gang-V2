#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 19:57:33 2026

@author: ime_rasenberg
"""

import os
import json

master_dict = {}
decoder = json.JSONDecoder()

# ==========================================================
# find all beta folders
# ==========================================================
beta_folders = [
    f for f in os.listdir(".")
    if os.path.isdir(f) and f.startswith("Data ")
]

print(f"Found {len(beta_folders)} beta folders")

for data_folder in beta_folders:

    # --------------------------------------
    # extract beta from folder name
    # --------------------------------------
    try:
        beta = float(data_folder.split(" ", 1)[1])
    except ValueError:
        print(f"Skipping folder {data_folder}")
        continue

    print(f"\nLoading beta = {beta}")

    if beta not in master_dict:
        master_dict[beta] = {}

    # --------------------------------------
    # collect json files
    # --------------------------------------
    filepaths = [
        os.path.join(data_folder, f)
        for f in os.listdir(data_folder)
        if os.path.isfile(os.path.join(data_folder, f))
        and f.endswith(".json")
    ]

    print(f"Found {len(filepaths)} JSON files")

    for filepath in filepaths:

        filename = os.path.basename(filepath)

        # --------------------------------------
        # scrape metadata
        # --------------------------------------
        name_no_ext = os.path.splitext(filename)[0]
        parts = name_no_ext.split("__")

        metadata = {}

        for part in parts:
            if "_" in part:
                key, value = part.split("_", 1)
                metadata[key] = value

        if not all(k in metadata for k in ["D", "Hz", "I"]):
            print(f"Skipping {filename}")
            continue

        try:
            D = float(metadata["D"])
            Hz = float(metadata["Hz"])
            I = int(metadata["I"])
        except ValueError:
            print(f"Malformed metadata in {filename}")
            continue

        # --------------------------------------
        # initialise nested structure
        # master_dict[beta][D][Hz][I]
        # --------------------------------------
        if D not in master_dict[beta]:
            master_dict[beta][D] = {}

        if Hz not in master_dict[beta][D]:
            master_dict[beta][D][Hz] = {}

        master_dict[beta][D][Hz][I] = {}

        # --------------------------------------
        # load concatenated json objects
        # --------------------------------------
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

                master_dict[beta][D][Hz][I][step] = obj

            except json.JSONDecodeError as e:
                print(f"Error in {filename} near position {pos}: {e}")
                break

print("Finished loading.")

#%%
import numpy as np
import matplotlib.pyplot as plt

N_LAST = 30

beta_vals = sorted(master_dict.keys())

for beta in beta_vals:
    beta_dict = master_dict[beta]

    # -------------------------------------------------
    # collect unique parameter values for this beta
    # -------------------------------------------------
    D_vals = sorted(beta_dict.keys())

    Hz_vals = sorted({
        hz
        for D in beta_dict
        for hz in beta_dict[D]
    })

    I_vals = sorted({
        I
        for D in beta_dict
        for hz in beta_dict[D]
        for I in beta_dict[D][hz]
    })

    tp = np.array([])
    pp = np.array([])
    tn = np.array([])
    pn = np.array([])
    Q  = np.array([])

    for I in I_vals:

        plus_grid  = np.full((len(Hz_vals), len(D_vals)), np.nan)
        pp_grid    = np.full((len(Hz_vals), len(D_vals)), np.nan)
        minus_grid = np.full((len(Hz_vals), len(D_vals)), np.nan)
        pn_grid    = np.full((len(Hz_vals), len(D_vals)), np.nan)
        Q_grid     = np.full((len(Hz_vals), len(D_vals)), np.nan)

        for ix, D in enumerate(D_vals):
            for iy, Hz in enumerate(Hz_vals):

                if Hz not in beta_dict[D]:
                    continue
                if I not in beta_dict[D][Hz]:
                    continue

                step_dict = beta_dict[D][Hz][I]
                if len(step_dict) == 0:
                    continue

                largest_steps = sorted(step_dict.keys())[-N_LAST:]

                nplus_lengths, pp_val = [], []
                nminus_lengths, pn_val = [], []
                Q_val = []

                for step in largest_steps:
                    obj = step_dict[step]
                    Q_val.append(obj["Q"])
                    if "N+" in obj:
                        nplus_lengths.append(len(obj["N+"]))
                        pp_val.append(obj["max_sum"])
                    if "N-" in obj:
                        nminus_lengths.append(len(obj["N-"]))
                        pn_val.append(obj["min_sum"])

                if nplus_lengths:
                    plus_grid[iy, ix]  = np.mean(nplus_lengths)
                    pp_grid[iy, ix]    = np.mean(pp_val)
                if nminus_lengths:
                    minus_grid[iy, ix] = np.mean(nminus_lengths)
                    pn_grid[iy, ix]    = np.mean(pn_val)
                if Q_val:
                    Q_grid[iy, ix]     = np.mean(Q_val)

        if len(tp) == 0:
            tp = plus_grid.copy()
            tn = minus_grid.copy()
            pp = pp_grid.copy()
            pn = pn_grid.copy()
            Q  = Q_grid.copy()
        else:
            tp += plus_grid
            tn += minus_grid
            pp += pp_grid
            pn += pn_grid
            Q  += Q_grid

    n_I = len(I_vals)
    extent = [min(D_vals), max(D_vals), min(Hz_vals), max(Hz_vals)]

    def _imshow(data, title, cbar_label=""):
        fig, ax = plt.subplots(figsize=(8, 6))
        im = ax.imshow(data, origin="lower", aspect="auto", extent=extent)
        plt.colorbar(im, ax=ax, label=cbar_label)
        ax.set_xlabel("D/J")
        ax.set_ylabel("H/J")
        ax.set_title(f"[β={beta}]  {title}")
        plt.show()

    # _imshow((tn - tp) / n_I,           "N- − N+")
    # _imshow(-(Q) / n_I,                "Q measured")
    # _imshow((tn - tp) / n_I + Q / n_I, "Q − (N- − N+)")
    # _imshow((pp) / n_I / tp,           "positive peak height")
    # _imshow(-(pn) / n_I / tn,          "negative peak height")

    # ratio plot with contour
    cut_off = 1.5
    mtx = -(pn) / n_I / tn / ((pp) / n_I / tp)

    fig, ax = plt.subplots(figsize=(8, 6))
    im = ax.imshow(mtx, origin="lower", aspect="auto", extent=extent)
    plt.colorbar(im, ax=ax, label="neg / pos peak height")
    ax.contour(mtx, levels=[cut_off], colors="red", linewidths=1.5, extent=extent)
    ax.set_xlabel("D/J")
    ax.set_ylabel("H/J")
    ax.set_title(f"[β={beta}]  negative peak height / positive peak height")
    plt.show()

    # --------------------------------------------------
    # lattice / spin phase plots
    # --------------------------------------------------
    N_LAST_SPIN = 1

    x = np.array([])
    y = np.array([])
    z = np.array([])

    for I in I_vals:

        sx = np.full((len(Hz_vals), len(D_vals)), np.nan)
        sy = sx.copy()
        sz = sx.copy()

        for ix, D in enumerate(D_vals):
            for iy, Hz in enumerate(Hz_vals):

                if Hz not in beta_dict[D]:
                    continue
                if I not in beta_dict[D][Hz]:
                    continue

                step_dict = beta_dict[D][Hz][I]
                if len(step_dict) == 0:
                    continue

                largest_steps = sorted(step_dict.keys())[-N_LAST_SPIN:]

                sxl, syl, szl = [], [], []
                for step in largest_steps:
                    obj = step_dict[step]
                    sxl.append(obj["spins"][0])
                    syl.append(obj["spins"][1])
                    szl.append(obj["spins"][2])

                if sxl:
                    sx[iy, ix] = np.mean(sxl)
                    sy[iy, ix] = np.mean(syl)
                    sz[iy, ix] = np.mean(szl)

        if len(x) == 0:
            x = sx.copy()
            y = sy.copy()
            z = sz.copy()
        else:
            x += sx
            y += sy
            z += sz

    fig, ax = plt.subplots(figsize=(8, 6))
    im = ax.imshow(z / 40**2 / n_I, origin="lower", aspect="auto", extent=extent)
    plt.colorbar(im, ax=ax, label="⟨Sz⟩ / N")
    ax.contour(mtx, levels=[cut_off], colors="red", linewidths=1.5, extent=extent)
    ax.set_xlabel("D/J")
    ax.set_ylabel("H/J")
    ax.set_title(fr"[β={beta}]  Average Magnetisation $\hat{{z}}$")
    plt.show()
    
#%%
# ==========================================================
# Q vs beta for the grid point closest to a target [H, D]
# ==========================================================

target_H = 1.75   # <-- set your target H/J
target_D = 1.7   # <-- set your target D/J

# find closest available H and D values (global across all betas)
all_D_vals  = sorted({ D  for b in master_dict for D  in master_dict[b] })
all_Hz_vals = sorted({ Hz for b in master_dict for D  in master_dict[b] for Hz in master_dict[b][D] })

closest_D  = min(all_D_vals,  key=lambda d: abs(d  - target_D))
closest_Hz = min(all_Hz_vals, key=lambda h: abs(h  - target_H))

print(f"Target  : H={target_H}, D={target_D}")
print(f"Snapped : H={closest_Hz}, D={closest_D}")

# for each beta, average Q over all I and the last N_LAST steps
N_LAST_Q = 30

beta_vals = sorted(master_dict.keys())
Q_vs_beta = []

for beta in beta_vals:
    beta_dict = master_dict[beta]

    if closest_D  not in beta_dict:                  continue
    if closest_Hz not in beta_dict[closest_D]:       continue

    I_dict = beta_dict[closest_D][closest_Hz]
    q_all  = []

    for I, step_dict in I_dict.items():
        if len(step_dict) == 0:
            continue
        largest_steps = sorted(step_dict.keys())[-N_LAST_Q:]
        for step in largest_steps:
            obj = step_dict[step]
            if "Q" in obj:
                q_all.append(obj["Q"])

    if q_all:
        # Q_vs_beta.append((beta, np.mean(q_all)))
        q_arr = np.array(q_all)
        Q_vs_beta.append((beta, np.mean(q_arr), np.std(q_arr) / np.sqrt(len(q_arr))))

if not Q_vs_beta:
    print("No data found for this point.")
else:
    betas, Qs, errs = zip(*Q_vs_beta)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.errorbar(betas, Qs, yerr=errs, marker="o", linewidth=1.5, capsize=4, capthick=1.5)
    ax.set_xlabel(r"$\beta$")
    ax.set_ylabel("Q")
    ax.set_title(rf"Q vs $\beta$  —  H/J={closest_Hz}, D/J={closest_D}")
    ax.grid(True, linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.show()