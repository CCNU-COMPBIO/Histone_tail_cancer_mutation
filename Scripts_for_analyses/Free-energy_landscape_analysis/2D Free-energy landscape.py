# ===============================================================
# Free-energy landscape (Python) with K-means cluster mapping
# - Reaction coordinates: NOC (x-axis) and EED (y-axis)
# - Free energy reported as dimensionless: ΔG/(kBT) = -ln(P/Pmax)
# - Number of clusters selected automatically using the elbow method
# ===============================================================

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.colors import LinearSegmentedColormap
from scipy.ndimage import gaussian_filter, label
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler

# ---------------- Global plotting style ----------------
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = ["Arial"]
plt.rcParams["mathtext.fontset"] = "stix"

# ---------------- User-configurable parameters ----------------
excel_file = "NOC_EED.xlsx"
sheet_names = ["H2B WT"]

nbins = 65
k_min = 2
k_max = 30

tick_label_size = 25
border_width = 2.0
tick_line_length = 0
contour_levels_num = 6
n_ticks_x = 5
n_ticks_y = 6
cbar_label_size = 20
cluster_point_size = 80

# Jitter to reduce aliasing from integer-valued NOC (for visualization)
jitter_amp = 1.0      # ±1.0
jitter_step = 0.01    # discretization

# ---------------- Elbow-method helper ----------------
def find_elbow_k(ks, inertias):
    """Pick elbow point using maximum distance to the line connecting endpoints."""
    ks = np.array(ks, dtype=float)
    inertias = np.array(inertias, dtype=float)

    x1, y1 = ks[0], inertias[0]
    x2, y2 = ks[-1], inertias[-1]

    denom = np.sqrt((y2 - y1) ** 2 + (x2 - x1) ** 2)
    if denom == 0:
        return int(ks[0])

    dists = []
    for x0, y0 in zip(ks, inertias):
        d = abs((y2 - y1) * x0 - (x2 - x1) * y0 + x2 * y1 - y2 * x1) / denom
        dists.append(d)

    return int(ks[np.argmax(dists)])

# ---------------- Main ----------------
if not os.path.exists(excel_file):
    raise FileNotFoundError(f"File not found: {excel_file}")

for sheet_name in sheet_names:
    print(f"\nProcessing {excel_file} - {sheet_name}...")

    try:
        data = pd.read_excel(excel_file, sheet_name=sheet_name)
    except Exception as e:
        print(f"Error reading sheet '{sheet_name}': {e}")
        continue

    if ("EED" not in data.columns) or ("NOC" not in data.columns):
        print(f"Columns 'EED' or 'NOC' missing in {sheet_name}. Skipping...")
        continue

    # Reaction coordinates
    NOC = data["NOC"].values.astype(float)  # x-axis
    EED = data["EED"].values.astype(float)  # y-axis

    # Add discrete jitter to NOC for visualization only (reduces vertical striping)
    max_steps = int(round(jitter_amp / jitter_step))
    jitter = np.random.randint(-max_steps, max_steps + 1, size=NOC.shape) * jitter_step
    NOC_jittered = NOC + jitter

    # ------------------------------------------------------------
    # Manual ranges used for different systems (keep as reference)
    # ------------------------------------------------------------
    # H2B:
    xlim_manual = [35, 135]
    ylim_manual = [0, 60]

    # H2A_C:
    # xlim_manual = [15, 65]
    # ylim_manual = [0, 35]

    # H2A_N:
    # xlim_manual = [14, 76]
    # ylim_manual = [0, 40]

    # H3:
    # xlim_manual = [75, 210]
    # ylim_manual = [0, 60]

    # H4:
    # xlim_manual = [30, 125]
    # ylim_manual = [0, 50]

    # Histogram bin edges
    xedges = np.linspace(xlim_manual[0], xlim_manual[1], nbins + 1)
    yedges = np.linspace(ylim_manual[0], ylim_manual[1], nbins + 1)

    # 2D histogram: x=NOC, y=EED
    H, _, _ = np.histogram2d(NOC_jittered, EED, bins=[xedges, yedges])

    total_counts = np.sum(H)
    if total_counts == 0:
        print(f"No points within the specified ranges for {sheet_name}. Skipping...")
        continue

    # Occupancy probability per bin: Pij = Hij / sum_{i,j} Hij
    P = H / total_counts

    # Dimensionless free energy: ΔG/(kBT) = -ln(P/Pmax)
    free_energy = np.zeros_like(P, dtype=float)
    mask = H > 0
    Pmax = np.max(P[mask])
    free_energy[mask] = -np.log(P[mask] / Pmax)

    fe_max_val = np.max(free_energy[mask]) if np.any(mask) else 1.0
    if fe_max_val == 0:
        fe_max_val = 1.0

    # Assign max value to zero-occupancy bins for visualization only
    free_energy[~mask] = fe_max_val

    # Colormap
    color_stops = [
        [0.3, 0.1, 0.1],
        [0.3, 0.2, 0.7],
        [0.2, 0.6, 0.8],
        [0.2, 0.8, 0.8],
        [0.2, 0.8, 0.2],
        [0.8, 0.8, 0.1],
    ]
    cm = LinearSegmentedColormap.from_list("custom_cmap", color_stops, N=1024)

    # Plot grid (bin centers)
    xcenters = 0.5 * (xedges[:-1] + xedges[1:])
    ycenters = 0.5 * (yedges[:-1] + yedges[1:])
    X, Y = np.meshgrid(xcenters, ycenters)

    fig, ax = plt.subplots(figsize=(10, 8), dpi=100)
    mesh = ax.pcolormesh(
        X, Y, free_energy.T, cmap=cm, shading="gouraud", vmin=0, vmax=fe_max_val
    )

    # Contours of free energy
    levels_all = np.linspace(0, fe_max_val, contour_levels_num)
    if len(np.unique(levels_all)) > 1:
        ax.contour(
            X, Y, free_energy.T, levels=levels_all,
            colors="k", linewidths=0.6, alpha=1, linestyles="solid"
        )

    # Smooth outer boundary contour of occupied region (optional visual aid)
    mask_grid = (H.T > 0)
    labeled_array, num_features = label(mask_grid)
    if num_features > 0:
        component_sizes = np.bincount(labeled_array.ravel())
        min_component_size = 20
        filtered_mask = np.zeros_like(mask_grid, dtype=bool)
        for i in range(1, num_features + 1):
            if component_sizes[i] >= min_component_size:
                filtered_mask[labeled_array == i] = True
    else:
        filtered_mask = mask_grid

    smoothed_mask = gaussian_filter(filtered_mask.astype(float), sigma=0.5)
    if np.max(smoothed_mask) > 0.1:
        ax.contour(X, Y, smoothed_mask, levels=[0.1], colors="k", linewidths=0.6)

    # ---------------- K-means clustering (elbow-selected k) ----------------
    valid_mask = (
        (NOC_jittered >= xlim_manual[0]) & (NOC_jittered <= xlim_manual[1]) &
        (EED >= ylim_manual[0]) & (EED <= ylim_manual[1])
    )

    Xf = NOC_jittered[valid_mask]
    Yf = EED[valid_mask]

    if len(Xf) >= k_max:
        features = np.column_stack((Xf, Yf))
        scaler = StandardScaler()
        features_scaled = scaler.fit_transform(features)

        ks_range = range(k_min, k_max + 1)
        inertias = []
        for k in ks_range:
            km = KMeans(n_clusters=k, random_state=42, n_init=10)
            km.fit(features_scaled)
            inertias.append(km.inertia_)

        best_k = find_elbow_k(list(ks_range), inertias)
        print(f"  Elbow-selected k = {best_k}")

        km_final = KMeans(n_clusters=best_k, random_state=42, n_init=10)
        km_final.fit(features_scaled)
        centers = scaler.inverse_transform(km_final.cluster_centers_)

        labels = km_final.labels_
        cluster_counts = np.bincount(labels, minlength=best_k)

        s_min, s_max = 30, 350
        if np.max(cluster_counts) > np.min(cluster_counts):
            norm = (cluster_counts - np.min(cluster_counts)) / (np.max(cluster_counts) - np.min(cluster_counts))
            point_sizes = s_min + norm * (s_max - s_min)
        else:
            point_sizes = np.full(best_k, cluster_point_size)

        ax.scatter(
            centers[:, 0], centers[:, 1],
            c="black", marker="o", s=point_sizes, zorder=100
        )
    else:
        print(
            f"Skipping K-means for {sheet_name}: not enough points ({len(Xf)}) to search up to k={k_max}."
        )

    # Axes and labels
    ax.set_xlabel("Number of Contacts", fontweight="bold", fontsize=25)
    ax.set_ylabel("End to End Distance (Å)", fontweight="bold", fontsize=25)
    ax.set_title(f"{sheet_name}", fontweight="bold", fontsize=35)

    ax.set_xlim(xlim_manual)
    ax.set_ylim(ylim_manual)

    ax.xaxis.set_major_locator(ticker.MaxNLocator(n_ticks_x, prune="both"))
    ax.yaxis.set_major_locator(ticker.MaxNLocator(n_ticks_y, prune="both"))
    ax.tick_params(
        axis="both", which="major",
        labelsize=tick_label_size, width=2.0, length=tick_line_length
    )

    for spine in ax.spines.values():
        spine.set_linewidth(border_width)

    cbar = plt.colorbar(mesh, ax=ax)
    cbar.ax.yaxis.set_major_locator(ticker.MultipleLocator(1))
    cbar.set_label(r"$\Delta G/(k_B T)$", fontsize=cbar_label_size, fontweight="bold", labelpad=10)
    cbar.ax.tick_params(labelsize=tick_label_size, width=2.0, length=tick_line_length)
    cbar.outline.set_linewidth(border_width)

    for lab in ax.get_xticklabels() + ax.get_yticklabels() + cbar.ax.get_yticklabels():
        lab.set_fontname("Arial")
        lab.set_fontweight("bold")

    fig_name = f"FreeEnergy_KMeans_Elbow_{sheet_name}.png"
    fig.savefig(fig_name, dpi=600, bbox_inches="tight")
    print(f"Saved: {fig_name}")

    plt.show()
