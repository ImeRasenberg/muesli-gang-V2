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


tp = np.array([])
pp = np.array([])
tn = np.array([])
pn = np.array([])
Q = np.array([])


for I in I_vals:

    # matrices for plotting
    plus_grid = np.full((len(Hz_vals), len(D_vals)), np.nan)
    pp_grid = np.full((len(Hz_vals), len(D_vals)), np.nan)
    minus_grid = np.full((len(Hz_vals), len(D_vals)), np.nan)
    pn_grid = np.full((len(Hz_vals), len(D_vals)), np.nan)
    Q_grid = np.full((len(Hz_vals), len(D_vals)), np.nan)
    
    
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
            pp_val = []
            nminus_lengths = []
            pn_val = []
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

            # average over last N steps
            if len(nplus_lengths) > 0:
                plus_grid[iy, ix] = np.mean(nplus_lengths)
                pp_grid[iy, ix] = np.mean(pp_val)

            if len(nminus_lengths) > 0:
                minus_grid[iy, ix] = np.mean(nminus_lengths)
                pn_grid[iy, ix] = np.mean(pn_val)
            
            if len(Q_val)>0:
                Q_grid[iy, ix] = np.mean(Q_val)

    # ==================================================
    # plotting
    # ==================================================
    
    if len(tp) == 0:
        tp = plus_grid.copy()
        tn = minus_grid.copy() 
        pp = pp_grid.copy()
        pn = pn_grid.copy()
        Q = Q_grid.copy()
    else:
        tp += plus_grid
        tn += minus_grid
        pp += pp_grid
        pn += pn_grid
        Q += Q_grid


#%%
fig, ax = plt.subplots(figsize=(8, 6))

im = ax.imshow(
    (tn-tp)/len(I_vals),
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

# ======================================= #

fig, ax = plt.subplots(figsize=(8, 6))

im = ax.imshow(
    -(Q)/len(I_vals),
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
ax.set_title(f"Q measured")

plt.show()

# ======================================= #

fig, ax = plt.subplots(figsize=(8, 6))

im = ax.imshow(
    (pp)/len(I_vals)/tp,
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
ax.set_title(f"positie peaks hight")

plt.show()

# ======================================= #


fig, ax = plt.subplots(figsize=(8, 6))

im = ax.imshow(
    -(pn)/len(I_vals)/tn,
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
ax.set_title(f"negative peak hight")

plt.show()


#%%
cut_off = 1.5
mtx =  -(pn)/len(I_vals)/tn / ((pp)/len(I_vals)/tp)
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

plt.colorbar(im, ax=ax, label="mean len(N-)")
ax.contour(mtx, levels = [cut_off] , colors="red", linewidths=1.5,
           extent=[
               min(D_vals),
               max(D_vals),
               min(Hz_vals),
               max(Hz_vals),
           ]
           )
ax.set_xlabel("D")
ax.set_ylabel("Hz")
ax.set_title(f"negative peak hight / positive peak hight")

plt.show()
#%% looking at the lattice phase

import numpy as np
import matplotlib.pyplot as plt


N_LAST =1

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


x = np.array([])
y = np.array([])
z = np.array([])


for I in I_vals:

    sx = np.full((len(Hz_vals), len(D_vals)), np.nan)
    sy = sx.copy()
    sz = sx.copy()
    
    for ix, D in enumerate(D_vals):
        for iy, Hz in enumerate(Hz_vals):

            if Hz not in master_dict[D]:
                continue

            if I not in master_dict[D][Hz]:
                continue

            if mtx[iy][ix] > cut_off:
                # print(f"skipping H{Hz}, D{D}")
                continue
        
            # if ix + iy  < 1:
            #     print(f"skipping H{Hz}, D{\D}")
            #     continue
            
            
            step_dict = master_dict[D][Hz][I]

            if len(step_dict) == 0:
                continue


            largest_steps = sorted(step_dict.keys())[-N_LAST:]

            sxl = []
            syl = []
            szl = []
            

            for step in largest_steps:

                obj = step_dict[step]
                
                sxl.append(obj["spins"][0])
                syl.append(obj["spins"][1])
                szl.append(obj["spins"][2])

            
            if len(Q_val)>0:
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

im = ax.imshow(
    z/40**2/len(I_vals),
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
ax.set_title(f"negative peak hight")

plt.show()

#%%
import numpy as np
import matplotlib.pyplot as plt

# ==========================================================
# parameters
# ==========================================================
L = 40
sigma = 1.0
N_LAST = 1

# ==========================================================
# real-space grid
# ==========================================================
X, Y = np.meshgrid(np.arange(L), np.arange(L), indexing="ij")


def add_gaussian(field, x, y, amp):
    dx = X - x
    dy = Y - y
    field += amp * np.exp(-(dx**2 + dy**2) / (2 * sigma**2))


# ==========================================================
# particle extraction (robust)
# ==========================================================
def extract_particles(obj):
    out = []

    for key in ["N+", "N-"]:
        if key not in obj:
            continue

        entries = obj[key]
        iterable = entries.values() if isinstance(entries, dict) else entries

        for e in iterable:
            if e is None:
                continue

            if len(e) >= 3:
                x, y, amp = e[0], e[1], e[2]
            elif len(e) == 2:
                x, y = e
                amp = 1.0
            else:
                continue

            out.append((x, y, amp))

    return out


# ==========================================================
# density
# ==========================================================
def build_density(step_dict):
    field = np.zeros((L, L))

    for step in sorted(step_dict.keys())[-N_LAST:]:
        obj = step_dict[step]

        for x, y, amp in extract_particles(obj):
            add_gaussian(field, x, y, amp)

    return field


# ==========================================================
# FFT + TRUE k-grid
# ==========================================================
def compute_fft(field):
    fft = np.fft.fftshift(np.fft.fft2(field))

    # physical k-grid (angular wave numbers)
    k_vals = np.fft.fftshift(np.fft.fftfreq(L, d=1.0)) * 2 * np.pi

    KX, KY = np.meshgrid(k_vals, k_vals, indexing="ij")

    return np.abs(fft), KX, KY


def low_k_mask(fft_field, KX, KY, kmax=2.0):
    K = np.sqrt(KX**2 + KY**2)
    mask = K < kmax

    out = np.copy(fft_field)
    out[~mask] = 0.0
    return out


# ==========================================================
# run first valid system
# ==========================================================
for D in master_dict:
    for Hz in master_dict[D]:
        for I in master_dict[D][Hz]:

            step_dict = master_dict[D][Hz][I]
            if not step_dict:
                continue

            print(f"Using D={D}, Hz={Hz}, I={I}")

            # ---------------------------
            # real space
            # ---------------------------
            field = build_density(step_dict)

            # ---------------------------
            # k space
            # ---------------------------
            fft_field, KX, KY = compute_fft(field)
            fft_lowk = low_k_mask(fft_field, KX, KY, kmax=2.0)

            K = np.sqrt(KX**2 + KY**2)

            # ==================================================
            # plotting
            # ==================================================
            fig, ax = plt.subplots(1, 2, figsize=(11, 4))

            # real space
            im0 = ax[0].imshow(field, origin="lower", cmap="inferno")
            ax[0].set_title("Real-space density")
            plt.colorbar(im0, ax=ax[0])

            # k-space (NOW CORRECT AXES)
            im1 = ax[1].pcolormesh(
                KX,
                KY,
                np.log1p(fft_lowk),
                shading="auto",
                cmap="viridis"
            )

            ax[1].set_title("Low-k structure factor")
            ax[1].set_xlabel(r"$k_x$ (rad / site)")
            ax[1].set_ylabel(r"$k_y$ (rad / site)")

            plt.colorbar(im1, ax=ax[1])

            plt.tight_layout()
            plt.show()

            raise SystemExit