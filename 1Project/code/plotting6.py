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



here = "../documentation/sections/Graphs/"
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
        # plt.savefig(here+f"")
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
    cbar = plt.colorbar(im, ax=ax,)
    cbar.set_label("neg / pos peak height", size = 18)
    ax.contour(mtx, levels=[cut_off], colors="red", linewidths=1.5, extent=extent)
    ax.set_xlabel("D/J", size = 18)
    ax.set_ylabel("H/J", size = 18)
    plt.savefig(here+f"HD_B{beta}_PH_NH.png")
    # ax.set_title(f"[β={beta}]  negative peak height / positive peak height")
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
    cbar = plt.colorbar(im, ax=ax)
    ax.contour(mtx, levels=[cut_off], colors="red", linewidths=1.5, extent=extent)
    cbar.set_label(r"$\langle S_z \rangle / N$",size = 18)
    ax.set_xlabel("D/J", size = 18)
    ax.set_ylabel("H/J", size = 18)
    plt.savefig(here+f"HD_B{beta}_AM.png")
    # ax.set_title(fr"[β={beta}]  Average Magnetisation $\hat{{z}}$")
    plt.show()
    
#%%
# ==========================================================
# Q vs beta for the grid point closest to a target [H, D]
# ==========================================================

target_H = 0.25   # <-- set your target H/J
target_D = 0.7  # <-- set your target D/J

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
    ax.errorbar(1/np.array(betas), Qs, yerr=errs, marker="o", linewidth=1.5, capsize=4, capthick=1.5)
    ax.set_xlabel(r"$J\beta^{-1}$",size = 18)
    ax.set_ylabel("Q", size = 18)
    # ax.set_title(rf"Q vs $\beta$  —  H/J={closest_Hz}, D/J={closest_D}")
    ax.grid(True, linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.savefig(here + "QB3.png")
    plt.show()
#%%

# ==========================================================
# T_c map: temperature at which Q decays to 30% of its
# low-T (high-beta) value, for every (H, D) grid point
# ==========================================================

N_LAST_Q = 30

beta_vals   = sorted(master_dict.keys())
all_D_vals  = sorted({ D  for b in master_dict for D  in master_dict[b] })
all_Hz_vals = sorted({ Hz for b in master_dict for D  in master_dict[b]
                        for Hz in master_dict[b][D] })

Tc_grid = np.full((len(all_Hz_vals), len(all_D_vals)), np.nan)

for ix, D in enumerate(all_D_vals):
    for iy, Hz in enumerate(all_Hz_vals):

        # --- collect (beta, mean_Q) for this grid point ---
        Q_vs_beta = []
        for beta in beta_vals:
            beta_dict = master_dict[beta]
            if D   not in beta_dict:             continue
            if Hz  not in beta_dict[D]:          continue

            q_all = []
            for I, step_dict in beta_dict[D][Hz].items():
                if not step_dict:
                    continue
                for step in sorted(step_dict.keys())[-N_LAST_Q:]:
                    obj = step_dict[step]
                    if "Q" in obj:
                        q_all.append(obj["Q"])

            if q_all:
                Q_vs_beta.append((beta, np.mean(q_all)))

        if len(Q_vs_beta) < 2:
            continue

        # sort by descending beta (= ascending T)
        Q_vs_beta.sort(key=lambda x: -x[0])
        betas_arr = np.array([b for b, _ in Q_vs_beta])
        Qs_arr    = np.array([q for _, q in Q_vs_beta])

        # reference value: average of the 3 highest-beta points
        n_ref   = min(3, len(Qs_arr))
        Q_ref   = np.mean(Qs_arr[:n_ref])

        if Q_ref == 0:
            continue

        # threshold: 30 % of reference, respecting sign
        # "decayed to 30%" means |Q| has shrunk, i.e. Q moved toward zero
        Q_thresh = Q_ref * 0.30

        # find first index (going toward high T) where Q crosses threshold
        # and *stays* there (all subsequent points also past threshold)
        Tc = np.nan
        for i in range(1, len(Qs_arr)):
            # check whether Q has passed the 30% threshold
            if Q_ref > 0:
                crossed = Qs_arr[i] < Q_thresh
            else:
                crossed = Qs_arr[i] > Q_thresh   # Q_ref negative: toward 0 means >

            if crossed:
                # verify it stays past threshold for all remaining points
                remaining = Qs_arr[i:]
                if Q_ref > 0:
                    stays = np.all(remaining < Q_thresh)
                else:
                    stays = np.all(remaining > Q_thresh)

                if stays:
                    # interpolate between i-1 and i for a smoother estimate
                    b0, q0 = betas_arr[i-1], Qs_arr[i-1]
                    b1, q1 = betas_arr[i],   Qs_arr[i]
                    if q1 != q0:
                        beta_c = b0 + (Q_thresh - q0) * (b1 - b0) / (q1 - q0)
                    else:
                        beta_c = (b0 + b1) / 2
                    Tc = 1.0 / beta_c
                    break

        Tc_grid[iy, ix] = Tc

# ==========================================================
# plot T_c as a function of H and D
# ==========================================================
extent = [min(all_D_vals), max(all_D_vals), min(all_Hz_vals), max(all_Hz_vals)]

fig, ax = plt.subplots(figsize=(8, 6))
im = ax.imshow(Tc_grid, origin="lower", aspect="auto", extent=extent)
cbr = plt.colorbar(im, ax=ax )
cbr.set_label(r"$T_c \ (J/k_B)$",size = 18)

ax.set_xlabel("D/J",size = 18)
ax.set_ylabel("H/J",size = 18)
# ax.set_title(r"$T_c$ — temperature where $Q$ decays to 30% of low-$T$ value")
plt.tight_layout()
plt.show()


# ==========================================================
# mask T_c to the domain where mtx > cut_off at beta=4
# ==========================================================

cut_off = 1.5
target_beta = 4.0
closest_beta = min(master_dict.keys(), key=lambda b: abs(b - target_beta))
print(f"Using beta={closest_beta} for domain mask")

beta_dict = master_dict[closest_beta]
D_vals_b  = sorted(beta_dict.keys())
Hz_vals_b = sorted({ Hz for D in beta_dict for Hz in beta_dict[D] })

I_vals_b  = sorted({ I for D in beta_dict for Hz in beta_dict[D]
                       for I in beta_dict[D][Hz] })

tp_b = np.zeros((len(Hz_vals_b), len(D_vals_b)))
pp_b = np.zeros((len(Hz_vals_b), len(D_vals_b)))
tn_b = np.zeros((len(Hz_vals_b), len(D_vals_b)))
pn_b = np.zeros((len(Hz_vals_b), len(D_vals_b)))

for I in I_vals_b:
    for ix, D in enumerate(D_vals_b):
        for iy, Hz in enumerate(Hz_vals_b):
            if Hz not in beta_dict[D]:          continue
            if I  not in beta_dict[D][Hz]:      continue
            step_dict = beta_dict[D][Hz][I]
            if not step_dict:                   continue

            largest_steps = sorted(step_dict.keys())[-N_LAST_Q:]
            npl, ppl, nml, pnl = [], [], [], []

            for step in largest_steps:
                obj = step_dict[step]
                if "N+" in obj:
                    npl.append(len(obj["N+"]))
                    ppl.append(obj["max_sum"])
                if "N-" in obj:
                    nml.append(len(obj["N-"]))
                    pnl.append(obj["min_sum"])

            if npl: tp_b[iy, ix] += np.mean(npl);  pp_b[iy, ix] += np.mean(ppl)
            if nml: tn_b[iy, ix] += np.mean(nml);  pn_b[iy, ix] += np.mean(pnl)

n_I_b = len(I_vals_b)
mtx_b = -(pn_b / n_I_b / tn_b) / (pp_b / n_I_b / tp_b)

# build mask on the same (all_Hz_vals, all_D_vals) grid as Tc_grid
domain_mask = np.full((len(all_Hz_vals), len(all_D_vals)), False)

for ix, D in enumerate(all_D_vals):
    for iy, Hz in enumerate(all_Hz_vals):
        if D  not in D_vals_b:  continue
        if Hz not in Hz_vals_b: continue
        bix = D_vals_b.index(D)
        biy = Hz_vals_b.index(Hz)
        if mtx_b[biy, bix] > cut_off:
            domain_mask[iy, ix] = True

Tc_masked = np.where(domain_mask, Tc_grid, np.nan)

# ==========================================================
# plot masked T_c
# ==========================================================
extent = [min(all_D_vals), max(all_D_vals), min(all_Hz_vals), max(all_Hz_vals)]

fig, ax = plt.subplots(figsize=(8, 6))
im = ax.imshow(Tc_masked, origin="lower", aspect="auto", extent=extent)
cbar = plt.colorbar(im, ax=ax)
cbar.set_label(r"$T_f^*$", size = 18)
ax.set_xlabel("D/J", size = 18)
ax.set_ylabel("H/J", size = 18)
# ax.set_title(r"$T_c^*$")
plt.tight_layout()
plt.show()

#%% exponential fit instead


from scipy.optimize import curve_fit

def exp_decay(T, a, b, tau):
    return a + b * np.exp(-T / tau)

N_LAST_Q = 30

beta_vals   = sorted(master_dict.keys())
all_D_vals  = sorted({D for b in master_dict for D in master_dict[b]})
all_Hz_vals = sorted({Hz for b in master_dict
                         for D in master_dict[b]
                         for Hz in master_dict[b][D]})

tau_grid = np.full((len(all_Hz_vals), len(all_D_vals)), np.nan)

for ix, D in enumerate(all_D_vals):
    for iy, Hz in enumerate(all_Hz_vals):

        # --------------------------------------------------
        # collect (T, <Q>)
        # --------------------------------------------------
        Q_vs_T = []

        for beta in beta_vals:

            if D not in master_dict[beta]:
                continue
            if Hz not in master_dict[beta][D]:
                continue

            q_all = []

            for I, step_dict in master_dict[beta][D][Hz].items():

                if not step_dict:
                    continue

                for step in sorted(step_dict.keys())[-N_LAST_Q:]:

                    obj = step_dict[step]

                    if "Q" in obj:
                        q_all.append(obj["Q"])

            if q_all:
                Q_vs_T.append((1.0 / beta, np.mean(q_all)))

        if len(Q_vs_T) < 4:
            continue

        # --------------------------------------------------
        # arrays sorted by temperature
        # --------------------------------------------------
        Q_vs_T.sort(key=lambda x: x[0])

        T_arr = np.array([x[0] for x in Q_vs_T])
        Q_arr = np.array([x[1] for x in Q_vs_T])

        T_span = T_arr.max() - T_arr.min()

        if T_span <= 0:
            continue

        # --------------------------------------------------
        # fit
        # --------------------------------------------------
        a0 = Q_arr[-1]
        b0 = Q_arr[0] - a0
        tau0 = T_span / 2

        try:

            popt, _ = curve_fit(
                exp_decay,
                T_arr,
                Q_arr,
                p0=[a0, b0, tau0],
                bounds=(
                    [-np.inf, -np.inf, 1e-12],
                    [ np.inf,  np.inf, np.inf]
                ),
                maxfev=10000
            )

            a_fit, b_fit, tau_fit = popt

            # reject fits whose decay length exceeds
            # the measured temperature window
            if tau_fit > T_span*1.5:
                continue

            tau_grid[iy, ix] = tau_fit

        except (RuntimeError, ValueError):
            continue
        
fig, ax = plt.subplots(figsize=(8, 6))

im = ax.imshow(
    tau_grid,
    origin="lower",
    aspect="auto",
    extent=extent
)

cbar = plt.colorbar(im, ax=ax)
cbar.set_label(r"$T_c \ (J/k_B)$", size=18)

ax.set_xlabel("D/J", size=18)
ax.set_ylabel("H/J", size=18)

plt.tight_layout()
# plt.savefig("Tc_tau.png", dpi=300, bbox_inches="tight")
plt.show()

tau_masked = np.where(domain_mask, tau_grid, np.nan)


fig, ax = plt.subplots(figsize=(8, 6))

im = ax.imshow(
    tau_masked,
    origin="lower",
    aspect="auto",
    extent=extent
)

cbar = plt.colorbar(im, ax=ax)
cbar.set_label(r"$T_f^*$", size=18)

ax.set_xlabel("D/J", size=18)
ax.set_ylabel("H/J", size=18)

plt.tight_layout()
plt.savefig(here + "Tc_masked_tau.png", dpi=300, bbox_inches="tight")
plt.show()
#%%
#%%
# ==========================================================
# Q vs 1/beta for the selected point, with exponential fit
# ==========================================================

from scipy.optimize import curve_fit

def exp_decay(T, a, b, tau):
    return a + b * np.exp(-T / tau)

N_LAST_Q = 30

# reuse closest_D, closest_Hz from the Q vs beta cell above
Q_vs_T_sel = []

for beta in sorted(master_dict.keys()):
    beta_dict = master_dict[beta]
    if closest_D  not in beta_dict:             continue
    if closest_Hz not in beta_dict[closest_D]:  continue

    q_all = []
    for I, step_dict in beta_dict[closest_D][closest_Hz].items():
        if not step_dict: continue
        for step in sorted(step_dict.keys())[-N_LAST_Q:]:
            obj = step_dict[step]
            if "Q" in obj:
                q_all.append(obj["Q"])

    if q_all:
        q_arr = np.array(q_all)
        Q_vs_T_sel.append((1.0 / beta, np.mean(q_arr), np.std(q_arr) / np.sqrt(len(q_arr))))

Q_vs_T_sel.sort(key=lambda x: x[0])
T_sel  = np.array([x[0] for x in Q_vs_T_sel])
Q_sel  = np.array([x[1] for x in Q_vs_T_sel])
err_sel = np.array([x[2] for x in Q_vs_T_sel])

fig, ax = plt.subplots(figsize=(8, 5))
ax.errorbar(T_sel, Q_sel, yerr=err_sel, marker="o", linewidth=1.5,
            capsize=4, capthick=1.5, label="mean ± s.e.m.")

try:
    T_span = T_sel.max() - T_sel.min()
    popt, _ = curve_fit(
        exp_decay, T_sel, Q_sel,
        p0=[Q_sel[-1], Q_sel[0] - Q_sel[-1], T_span / 2],
        bounds=([-np.inf, -np.inf, 1e-12], [np.inf, np.inf, np.inf]),
        maxfev=10000
    )
    T_fit = np.linspace(T_sel.min(), T_sel.max(), 300)
    ax.plot(T_fit, exp_decay(T_fit, *popt), color="tomato", linewidth=2,
            label=rf"fit: $T_f^*={popt[2]:.3f}$")
    print(f"Fit (H={closest_Hz}, D={closest_D}): a={popt[0]:.4f}, b={popt[1]:.4f}, Tf*={popt[2]:.4f}")
except RuntimeError:
    print("Fit failed for selected point.")

ax.set_xlabel(r"$J\beta^{-1}$", size=18)
ax.set_ylabel("Q", size=18)
ax.legend(fontsize=11)
ax.grid(True, linestyle="--", alpha=0.5)
plt.tight_layout()
plt.savefig(here + "QB3_fit.png", dpi=150)
plt.show()

# ==========================================================
# T_f* map from exponential fit, unmasked
# ==========================================================

Tf_grid = np.full((len(all_Hz_vals), len(all_D_vals)), np.nan)

for ix, D in enumerate(all_D_vals):
    for iy, Hz in enumerate(all_Hz_vals):

        Q_vs_T = []
        for beta in beta_vals:
            if D  not in master_dict[beta]:          continue
            if Hz not in master_dict[beta][D]:       continue

            q_all = []
            for I, step_dict in master_dict[beta][D][Hz].items():
                if not step_dict: continue
                for step in sorted(step_dict.keys())[-N_LAST_Q:]:
                    obj = step_dict[step]
                    if "Q" in obj:
                        q_all.append(obj["Q"])

            if q_all:
                Q_vs_T.append((1.0 / beta, np.mean(q_all)))

        if len(Q_vs_T) < 4:
            continue

        Q_vs_T.sort(key=lambda x: x[0])
        T_arr = np.array([x[0] for x in Q_vs_T])
        Q_arr = np.array([x[1] for x in Q_vs_T])
        T_span = T_arr.max() - T_arr.min()

        if T_span <= 0:
            continue

        try:
            popt, _ = curve_fit(
                exp_decay, T_arr, Q_arr,
                p0=[Q_arr[-1], Q_arr[0] - Q_arr[-1], T_span / 2],
                bounds=([-np.inf, -np.inf, 1e-12], [np.inf, np.inf, np.inf]),
                maxfev=10000
            )
            if popt[2] <= 1.5 * T_arr.max():
                Tf_grid[iy, ix] = popt[2]
        except (RuntimeError, ValueError):
            continue

extent = [min(all_D_vals), max(all_D_vals), min(all_Hz_vals), max(all_Hz_vals)]

fig, ax = plt.subplots(figsize=(8, 6))
im = ax.imshow(Tf_grid, origin="lower", aspect="auto", extent=extent)
cbar = plt.colorbar(im, ax=ax)
cbar.set_label(r"$T_f^*\ (J/k_B)$", size=18)
ax.set_xlabel("D/J", size=18)
ax.set_ylabel("H/J", size=18)
plt.tight_layout()
plt.show()

# ==========================================================
# T_f* map masked to skyrmion domain
# ==========================================================

Tf_masked = np.where(domain_mask, Tf_grid, np.nan)

fig, ax = plt.subplots(figsize=(8, 6))
im = ax.imshow(Tf_masked, origin="lower", aspect="auto", extent=extent)
cbar = plt.colorbar(im, ax=ax)
cbar.set_label(r"$T_f^*$", size=18)
ax.set_xlabel("D/J", size=18)
ax.set_ylabel("H/J", size=18)
plt.tight_layout()
plt.savefig(here + "Tf_masked.png", dpi=300, bbox_inches="tight")
plt.show()
