import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from scipy.interpolate import interp1d

# ==============================================================
# 1D free-energy profile projected onto the end-to-end distance (EED)
# - Free energy is reported as dimensionless: ΔG/(kBT) = -ln p(x) + const
# - The curve is shifted so that the minimum within the data-supported
#   region is zero.
# - A low-energy accessible window is defined by ΔG/(kBT) ≤ 1.
#   D_EE: window width (right intersection - left intersection)
#   D_MAX: upper boundary of the window (right intersection; per manuscript)
# ==============================================================

# ---------------- User parameters ----------------
excel_file = "NOC_EED.xlsx"
sheet_name = "H3 WT"
col_name = "EED"

level_kbT = 1.0          # threshold: min + 1 (in units of kBT)
n_grid = 1000            # grid resolution for the profile

# Manual axis ranges (set to None for automatic)
#H2B
#xlim_manual = [0, 60]    # EED range (Å)
#ylim_manual = [-1, 6]    # ΔG/(kBT)
#H3
xlim_manual = [0, 60]    
ylim_manual = [-1, 7]  
#H2A_N
#xlim_manual = [0, 40]    
#ylim_manual = [-1, 5] 
#H2A_C
#xlim_manual = [0, 40]    
#ylim_manual = [-1, 6] 
#H4
#xlim_manual = [0, 50]    
#ylim_manual = [-1, 7] 

# Line widths
lineWidthFreeE = 2.5
lineWidthRefLine = 0.6
arrowHeadLineWidth = 2.0
axisLineWidth = 3.0

# Tick spacing (set <=0 for automatic)
xtickSpacing = 10
ytickSpacing = 1

# Fonts
fontName = "Arial"
fontSizeAxes = 25
fontSizeTitle = 30
fontSizeArrowText = 20

arrowColor = [0, 0, 0]

# ---------------- Read EED data ----------------
try:
    df = pd.read_excel(excel_file, sheet_name=sheet_name)
    if col_name not in df.columns:
        raise ValueError(f"Column '{col_name}' not found.")
    EED = df[col_name].dropna().values.astype(float)
    if len(EED) == 0:
        raise ValueError("EED data are empty.")
except Exception as e:
    print(f"Data loading error: {e}")
    try:
        xl = pd.ExcelFile(excel_file)
        print(f"Sheets in {excel_file}:\n{xl.sheet_names}")
    except Exception:
        pass
    sys.exit(1)

# ---------------- x-grid ----------------
data_xmin = np.min(EED)
data_xmax = np.max(EED)

xlim_use = [data_xmin, data_xmax] if xlim_manual is None else xlim_manual
xvals = np.linspace(xlim_use[0], xlim_use[1], n_grid)

valid_mask = (xvals >= data_xmin) & (xvals <= data_xmax)
if not np.any(valid_mask):
    raise ValueError("The x-grid does not overlap with the data range. Check xlim_manual.")

# ---------------- Estimate probability density p(x) ----------------
eps_p = 1e-12
try:
    kde = gaussian_kde(EED)
    p = kde(xvals)
    p = np.maximum(p, eps_p)
    p = p / np.trapz(p, xvals)
except Exception:
    try:
        nbins = max(10, min(100, len(EED) // 2))
        counts, bin_edges = np.histogram(EED, bins=nbins, density=True)
        bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
        f = interp1d(bin_centers, counts, kind="linear", bounds_error=False, fill_value=0.0)
        p = f(xvals)
        p = np.maximum(p, eps_p)
        p = p / np.trapz(p, xvals)
    except Exception:
        mu = np.mean(EED)
        sigma = max(0.1, 0.01 * abs(mu))
        p = np.exp(-(xvals - mu) ** 2 / (2 * sigma**2))
        p = np.maximum(p, eps_p)
        p = p / np.trapz(p, xvals)

p = np.maximum(p, eps_p)

# ---------------- Free energy: ΔG/(kBT) ----------------
freeE = -np.log(p)
freeE = freeE - np.min(freeE[valid_mask])  # set the minimum within the data-supported region to 0

# Do not draw outside the data-supported region
freeE_plot = freeE.copy()
freeE_plot[~valid_mask] = np.nan

# ---------------- Define the low-energy window at min + level_kbT ----------------
fe_min = np.min(freeE[valid_mask])
fe_level = fe_min + level_kbT

# Continuous window around the global minimum (avoids merging multiple basins)
imin = np.where(valid_mask)[0][np.argmin(freeE[valid_mask])]

if freeE[imin] > fe_level:
    left_cross = np.nan
    right_cross = np.nan
    D_EE = np.nan
    D_MAX = np.nan
else:
    ileft = imin
    while ileft > 0 and valid_mask[ileft] and freeE[ileft] <= fe_level:
        ileft -= 1

    iright = imin
    while iright < len(xvals) - 1 and valid_mask[iright] and freeE[iright] <= fe_level:
        iright += 1

    # Left intersection (linear interpolation)
    if ileft < imin:
        x1, x2 = xvals[ileft], xvals[ileft + 1]
        y1, y2 = freeE[ileft], freeE[ileft + 1]
        left_cross = x1 + (fe_level - y1) / (y2 - y1) * (x2 - x1)
    else:
        left_cross = xvals[imin]

    # Right intersection (linear interpolation)
    if iright > imin:
        x1, x2 = xvals[iright - 1], xvals[iright]
        y1, y2 = freeE[iright - 1], freeE[iright]
        right_cross = x1 + (fe_level - y1) / (y2 - y1) * (x2 - x1)
    else:
        right_cross = xvals[imin]

    D_EE = right_cross - left_cross
    D_MAX = right_cross  # per manuscript: upper boundary (right intersection) of the window

# ---------------- Plot ----------------
plt.rcParams["font.family"] = fontName
plt.rcParams["font.weight"] = "bold"

fig, ax = plt.subplots(figsize=(6, 6))
fig.subplots_adjust(left=0.01, bottom=0.1, right=0.9, top=0.9)

ax.plot(xvals, freeE_plot, "b-", linewidth=lineWidthFreeE)
ax.axhline(y=fe_level, color=[0.8, 0, 0], linestyle="--", linewidth=lineWidthRefLine)

# Axis ranges
if ylim_manual is None:
    ypad = 0.05 * np.ptp(freeE[valid_mask])
    ylim_use = [np.min(freeE[valid_mask]) - ypad, np.max(freeE[valid_mask]) + ypad]
else:
    ylim_use = ylim_manual

ax.set_xlim(xlim_use)
ax.set_ylim(ylim_use)

# Ticks
if xtickSpacing > 0:
    ax.set_xticks(np.arange(xlim_use[0], xlim_use[1] + xtickSpacing, xtickSpacing))
if ytickSpacing > 0:
    ax.set_yticks(np.arange(ylim_use[0], ylim_use[1] + ytickSpacing, ytickSpacing))

# Double-headed arrow and annotations
if not (np.isnan(left_cross) or np.isnan(right_cross)):
    ax.annotate(
        "",
        xy=(right_cross, fe_level),
        xytext=(left_cross, fe_level),
        arrowprops=dict(
            arrowstyle="<->,head_length=0.6,head_width=0.4",
            color=arrowColor,
            lw=arrowHeadLineWidth,
            shrinkA=0,
            shrinkB=0,
        ),
    )

    ax.plot(right_cross, fe_level, "ro", markersize=6.5, zorder=3)

    arrow_text_y = fe_level + 0.06 * (ylim_use[1] - ylim_use[0])
    ax.text(
        (left_cross + right_cross) / 2,
        arrow_text_y,
        f"D$_{{EE}}$={D_EE:.2f}",
        horizontalalignment="center",
        fontname=fontName,
        fontsize=fontSizeArrowText,
        fontweight="bold",
    )

    dmax_text_x = right_cross + 0.01 * (xlim_use[1] - xlim_use[0])
    dmax_text_y = fe_level - 0.06 * (ylim_use[1] - ylim_use[0])
    ax.text(
        dmax_text_x,
        dmax_text_y,
        f"D$_{{MAX}}$={D_MAX:.2f}",
        fontname=fontName,
        fontsize=20,
        fontweight="bold",
    )
# Styling
ax.tick_params(axis="both", which="major", labelsize=fontSizeAxes, length=0, width=axisLineWidth)
for spine in ax.spines.values():
    spine.set_linewidth(axisLineWidth)

ax.set_xlabel("End to End Distance (Å)", fontname=fontName, fontsize=fontSizeAxes, fontweight="bold")
ax.set_ylabel(r"$\Delta G(k_B T)$", fontname=fontName, fontsize=fontSizeTitle, fontweight="bold")
ax.set_title(sheet_name, fontname=fontName, fontsize=fontSizeTitle, fontweight="bold")
outfile = f"FreeEnergy_1D_EED_{sheet_name}.png"

plt.savefig(outfile, dpi=600, bbox_inches="tight")
print(f"Saved: {outfile}")
plt.show()
