import numpy as np
import matplotlib.pyplot as plt

# ==========================================================
# parameters
# ==========================================================
L = 40
sigma = 2.0
N_LAST = 1

k_min = 0.1
k_max = 2.0

# ==========================================================
# real-space grid
# ==========================================================
X, Y = np.meshgrid(np.arange(L), np.arange(L), indexing="ij")


# ==========================================================
# periodic Gaussian (torus geometry)
# ==========================================================
def add_gaussian(field, x, y, amp):
    dx = np.abs(X - x)
    dy = np.abs(Y - y)

    dx = np.minimum(dx, L - dx)
    dy = np.minimum(dy, L - dy)

    field += amp * np.exp(-(dx**2 + dy**2) / (2 * sigma**2))


# ==========================================================
# particle extraction
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
# FFT + k-grid
# ==========================================================
def compute_fft(field):
    fft = np.fft.fftshift(np.fft.fft2(field))

    k_vals = np.fft.fftshift(np.fft.fftfreq(L, d=1.0)) * 2 * np.pi
    KX, KY = np.meshgrid(k_vals, k_vals, indexing="ij")

    return np.abs(fft), KX, KY


# ==========================================================
# annular filter: k_min < |k| < k_max
# ==========================================================
def k_window(fft_field, KX, KY, k_min=0.1, k_max=2.0):
    K = np.sqrt(KX**2 + KY**2)

    mask = (K >= k_min) & (K <= k_max)

    out = np.zeros_like(fft_field)
    out[mask] = fft_field[mask]

    return out


# ==========================================================
# run first valid configuration
# ==========================================================
for D in master_dict:
    for Hz in master_dict[D]:
        for I in master_dict[D][Hz]:

            step_dict = master_dict[D][Hz][I]
            if not step_dict:
                continue

            print(f"Using D={D}, Hz={Hz}, I={I}")

            # --------------------------
            # real space
            # --------------------------
            field = build_density(step_dict)

            # --------------------------
            # FFT
            # --------------------------
            fft_field, KX, KY = compute_fft(field)

            fft_filtered = k_window(
                fft_field,
                KX,
                KY,
                k_min=k_min,
                k_max=k_max
            )

            # ==================================================
            # plot
            # ==================================================
            fig, ax = plt.subplots(1, 2, figsize=(11, 4))

            im0 = ax[0].imshow(field, origin="lower", cmap="viridis")
            ax[0].set_title("Periodic real-space density")
            plt.colorbar(im0, ax=ax[0])

            # restrict displayed k-range to [-2, 2]
            im1 = ax[1].pcolormesh(
                KX,
                KY,
                np.log1p(fft_filtered),
                shading="auto",
                cmap="inferno",
                vmin=0,
                vmax=np.log1p(np.max(fft_filtered) + 1e-12)
            )

            ax[1].set_xlim(-2, 2)
            ax[1].set_ylim(-2, 2)

            ax[1].set_title("Filtered structure factor (0.1 < |k| < 2)")
            ax[1].set_xlabel(r"$k_x$ (rad / site)")
            ax[1].set_ylabel(r"$k_y$ (rad / site)")

            plt.colorbar(im1, ax=ax[1])

            plt.tight_layout()
            plt.show()

            raise SystemExit
            
#%%

import numpy as np

# ==========================================================
# order parameters from FFT
# ==========================================================
def compute_order_parameters(fft_field, KX, KY):
    K = np.sqrt(KX**2 + KY**2)

    mask = (K > 0.1) & (K < 2.0)
    S = np.copy(fft_field) * mask

    if np.sum(S) == 0:
        return 0, 0, 0

    # --------------------------------------------------
    # find dominant ring peak correctly
    # --------------------------------------------------
    idx = np.unravel_index(np.argmax(S), S.shape)

    k0 = K[idx]

    # --------------------------------------------------
    # shell selection
    # --------------------------------------------------
    shell = np.abs(K - k0) < 0.2

    ky = KY[shell]
    kx = KX[shell]
    s  = S[shell]

    if len(s) == 0:
        return 0, 0, 0

    theta = np.arctan2(ky, kx)

    w = s / (np.sum(s) + 1e-12)

    psi_stripe = np.abs(np.sum(w * np.exp(1j * theta)))
    psi_square = np.abs(np.sum(w * np.exp(1j * 4 * theta)))
    psi_hex    = np.abs(np.sum(w * np.exp(1j * 6 * theta)))

    return psi_stripe, psi_square, psi_hex


# ==========================================================
# phase classification
# ==========================================================
def classify(psi_s, psi_sq, psi_h):
    if psi_h > 0.4:
        return "Skyrmion lattice (hex)"
    elif psi_sq > 0.4:
        return "Square phase"
    elif psi_s > 0.5:
        return "Stripe phase"
    else:
        return "Disordered"


# ==========================================================
# MAIN LOOP (replace plotting section only)
# ==========================================================
for D in master_dict:
    for Hz in master_dict[D]:

        # pick only ONE I (first valid)
        chosen = None

        for I in master_dict[D][Hz]:
            if master_dict[D][Hz][I]:
                chosen = I
                break

        if chosen is None:
            continue

        step_dict = master_dict[D][Hz][chosen]

        field = build_density(step_dict)
        fft_field, KX, KY = compute_fft(field)

        psi_s, psi_sq, psi_h = compute_order_parameters(fft_field, KX, KY)

        # --------------------------
        # phase classification
        # --------------------------
        if psi_h > 0.4:
            phase = "Skyrmion lattice (hex)"
        elif psi_sq > 0.4:
            phase = "Square phase"
        elif psi_s > 0.5:
            phase = "Stripe phase"
        else:
            phase = "Disordered"

        # ==================================================
        # plot
        # ==================================================
        fig, ax = plt.subplots(1, 2, figsize=(10, 4))

        im0 = ax[0].imshow(field, origin="lower", cmap="viridis")
        ax[0].set_title("Real space")

        im1 = ax[1].pcolormesh(
            KX, KY,
            np.log1p(fft_field),
            shading="auto",
            cmap="inferno"
        )

        ax[1].set_xlim(-2, 2)
        ax[1].set_ylim(-2, 2)
        ax[1].set_title("FFT (structure factor)")

        title = f"D={D}, Hz={Hz}  →  {phase}"
        fig.suptitle(title, fontweight="bold")

        plt.tight_layout()
        plt.show()
#%%
import numpy as np

# ==========================================================
# classification + OPs (same physics as before)
# ==========================================================
def compute_order_parameters(fft_field, KX, KY):
    K = np.sqrt(KX**2 + KY**2)

    mask = (K > 0.1) & (K < 2.0)
    S = np.copy(fft_field) * mask

    if np.sum(S) == 0:
        return 0, 0, 0

    idx = np.unravel_index(np.argmax(S), S.shape)
    k0 = K[idx]

    shell = np.abs(K - k0) < 0.2

    ky = KY[shell]
    kx = KX[shell]
    s = S[shell]

    if len(s) == 0:
        return 0, 0, 0

    theta = np.arctan2(ky, kx)
    w = s / (np.sum(s) + 1e-12)

    psi_stripe = np.abs(np.sum(w * np.exp(1j * theta)))
    psi_square = np.abs(np.sum(w * np.exp(1j * 4 * theta)))
    psi_hex = np.abs(np.sum(w * np.exp(1j * 6 * theta)))

    return psi_stripe, psi_square, psi_hex


def classify(psi_s, psi_sq, psi_h):
    if psi_h > 0.4:
        return "hex"
    elif psi_sq > 0.4:
        return "square"
    elif psi_s > 0.5:
        return "stripe"
    else:
        return "disorder"


# ==========================================================
# containers for counts
# ==========================================================
D_list = sorted(master_dict.keys())
Hz_list = sorted({hz for D in master_dict for hz in master_dict[D]})

counts_disorder = np.zeros((len(D_list), len(Hz_list)))
counts_stripe   = np.zeros((len(D_list), len(Hz_list)))
counts_square   = np.zeros((len(D_list), len(Hz_list)))
counts_hex      = np.zeros((len(D_list), len(Hz_list)))


# ==========================================================
# main scan
# ==========================================================
for iD, D in enumerate(D_list):
    for iH, Hz in enumerate(Hz_list):

        if Hz not in master_dict[D]:
            continue

        for I in master_dict[D][Hz]:

            step_dict = master_dict[D][Hz][I]
            if not step_dict:
                continue

            field = build_density(step_dict)
            fft_field, KX, KY = compute_fft(field)

            psi_s, psi_sq, psi_h = compute_order_parameters(fft_field, KX, KY)
            phase = classify(psi_s, psi_sq, psi_h)

            if phase == "hex":
                counts_hex[iD, iH] += 1
            elif phase == "square":
                counts_square[iD, iH] += 1
            elif phase == "stripe":
                counts_stripe[iD, iH] += 1
            else:
                counts_disorder[iD, iH] += 1


# ==========================================================
# plotting helper
# ==========================================================
def plot_map(data, title):
    plt.figure(figsize=(6, 5))
    plt.imshow(
        data,
        origin="lower",
        aspect="auto",
        cmap="viridis"
    )
    plt.colorbar(label="count")
    plt.xticks(range(len(Hz_list)), np.round(Hz_list, 2), rotation=90)
    plt.yticks(range(len(D_list)), np.round(D_list, 2))
    plt.title(title)
    plt.xlabel("Hz")
    plt.ylabel("D")
    plt.tight_layout()
    plt.show()


# ==========================================================
# 4 final phase maps
# ==========================================================
plot_map(counts_disorder, "Disorder phase count")
plot_map(counts_stripe,   "Stripe phase count")
plot_map(counts_square,   "Square phase count")
plot_map(counts_hex,      "Hexagonal (skyrmion lattice) count")