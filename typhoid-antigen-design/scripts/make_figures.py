"""Regenerate the report figures (matplotlib only).

usage (from the repository root): python3 -I scripts/make_figures.py
needs: numpy, matplotlib; run scripts/gene_prevalence.py first (Figure 1).
Writes figures/fig1_prevalence.png, fig4_epitope_map.png, fig6_specificity_heatmap.png,
fig7_ranking.png, fig8_rank_robustness.png.  Structure panels: scripts/render_structures.py.
"""
import csv, os
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch

plt.rcParams.update({"font.family": "serif", "font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
INK = "#0b0b0b"
COL = {"hlyE": "#2a78d6", "cdtB": "#eb6834", "yncE": "#1baf7a"}
FLAG = "#c2185b"                       # flagged literature segment HlyE 175-188
LEN = {"hlyE": 303, "cdtB": 269, "yncE": 353}
os.makedirs("figures", exist_ok=True)

rows = list(csv.DictReader(open("results/epitopes/epitope_ranking_v3.csv")))
for r in rows:
    r["start"], r["end"] = int(r["start"]), int(r["end"])
    r["np"], r["rk"], r["rkp"] = float(r["score_no_prior"]), int(r["rank_no_prior"]), int(r["rank"])
rows.sort(key=lambda r: (r["rk"], -r["np"]))
is_flag = lambda r: r["protein"] == "hlyE" and r["start"] == 175
seg_color = lambda r: FLAG if is_flag(r) else COL[r["protein"]]

# Figure 1: gene prevalence (unique genomes)
prev = {}
for r in csv.DictReader(open("results/genomes/gene_prevalence_unique.csv")):
    prev.setdefault(r["gene"], {})[r["group"]] = float(r["pct"])
fig, ax = plt.subplots(figsize=(6.2, 2.9))
groups = [("typhoidal", "Typhoidal (Typhi + Paratyphi A)", "#2a78d6"), ("NTS", "Non-typhoidal Salmonella", "#eb6834"),
          ("confounders", "Confounders", "#8a8a85")]
for j, (g, lab, c) in enumerate(groups):
    xs = [i + (j - 1) * 0.28 for i in range(3)]
    vals = [prev[k][g] for k in ("hlyE", "cdtB", "yncE")]
    ax.bar(xs, vals, 0.26, color=c, label=lab)
    for x, v in zip(xs, vals):
        ax.text(x, v + 1.5, f"{v:g}", ha="center", fontsize=8)
ax.set_xticks(range(3)); ax.set_xticklabels(["hlyE", "cdtB", "yncE"]); ax.set_ylim(0, 115)
ax.set_ylabel("Genomes carrying the gene (%)")
ax.legend(frameon=False, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=3)
ax.yaxis.grid(True, color="#e4e3df", lw=.6); ax.set_axisbelow(True)
fig.tight_layout(); fig.savefig("figures/fig1_prevalence.png", dpi=220); plt.close(fig)

# Figure 4: ranked segments along each protein
fig, axs = plt.subplots(3, 1, figsize=(6.2, 2.9))
for ax, p in zip(axs, LEN):
    ax.barh(0, LEN[p], color="#e4e3df", height=.5)
    for r in rows:
        if r["protein"] != p: continue
        ax.barh(0, r["end"] - r["start"] + 1, left=r["start"], color=seg_color(r), height=.5, edgecolor="white", lw=.5)
        if r["klass"] == "Typhi-selective":
            ax.plot((r["start"] + r["end"]) / 2, 0.5, marker="D", color=INK, ms=4)
    ax.set_xlim(0, 360); ax.set_ylim(-.6, .9); ax.set_yticks([]); ax.spines["left"].set_visible(False)
    ax.set_ylabel(p, rotation=0, ha="right", va="center")
    if p != "yncE": ax.set_xticks([]); ax.spines["bottom"].set_visible(False)
axs[2].set_xlabel("Residue position (full-length numbering)")
fig.tight_layout(); fig.savefig("figures/fig4_epitope_map.png", dpi=220); plt.close(fig)

# Figure 6: specificity heatmap from the sliding 9-mer window files
cm = LinearSegmentedColormap.from_list("s", ["#f4f3ef", "#bcd3f2", "#6ea1e0", "#2a5fa8", "#12305c"])
fig, axs = plt.subplots(3, 1, figsize=(6.4, 5.0), gridspec_kw={"hspace": 0.55})
for ax, p in zip(axs, LEN):
    W = list(csv.DictReader(open(f"results/genomes/{p}_windows_by_group.csv")))
    M = np.full((3, LEN[p]), np.nan)
    for w in W:
        c = int(w["Start"]) + 4
        for i, k in enumerate(("nts", "paratyphiA", "confounders")):
            M[i, c - 1] = int(w[k])
    for i in range(3):                  # carry each window value across the residues it is centred on
        last = np.nan
        for j in range(LEN[p]):
            if np.isnan(M[i, j]): M[i, j] = last
            else: last = M[i, j]
    im = ax.imshow(np.ma.masked_invalid(M), aspect="auto", cmap=cm, vmin=0, vmax=6, extent=(0.5, LEN[p] + 0.5, 3, 0),
                   interpolation="nearest")
    ax.set_yticks([.5, 1.5, 2.5]); ax.set_yticklabels(["NTS", "Paratyphi A", "Confounders"], fontsize=8)
    ax.set_xlim(0, LEN[p] + 1); ax.set_title(p, loc="left", fontsize=9.5, pad=9); ax.tick_params(length=2)
    for r in rows:
        if r["protein"] == p:
            ax.plot([r["start"], r["end"]], [-.18, -.18], color=seg_color(r), lw=3.5, solid_capstyle="butt", clip_on=False)
    if p != "yncE": ax.set_xticklabels([])
axs[2].set_xlabel("Residue position (centre of 9-residue window)")
cb = fig.colorbar(im, ax=axs, orientation="vertical", fraction=.025, pad=.02, ticks=[0, 2, 4, 6])
cb.set_label("Minimum mismatches vs Typhi 9-mer", fontsize=8); cb.ax.set_yticklabels(["0", "2", "4", ">=6"])
fig.savefig("figures/fig6_specificity_heatmap.png", dpi=220, bbox_inches="tight"); plt.close(fig)

# Figure 7: primary (no prior-evidence) ranking
fig, ax = plt.subplots(figsize=(6.2, 5.9))
y = list(range(len(rows)))[::-1]
ax.barh(y, [r["np"] for r in rows], color=[seg_color(r) for r in rows], height=.68)
for yy, r in zip(y, rows):
    tag = ("  (Typhi-selective)" if r["klass"] == "Typhi-selective" else "") + ("  (flagged, literature)" if is_flag(r) else "")
    ax.text(r["np"] + .008, yy, f"{r['np']:.3f}{tag}", va="center", fontsize=7.5)
ax.set_yticks(y); ax.set_yticklabels([f"{r['protein']} {r['start']}–{r['end']}" for r in rows], fontsize=8)
ax.set_xlim(0, 1.15); ax.set_xlabel("Composite score without prior-evidence term (0–1)")
ax.legend(handles=[Patch(color=COL[k], label=k) for k in COL], frameon=False, fontsize=8, loc="lower right")
ax.xaxis.grid(True, color="#e4e3df", lw=.6); ax.set_axisbelow(True)
fig.tight_layout(); fig.savefig("figures/fig7_ranking.png", dpi=220); plt.close(fig)

# Figure 8: rank robustness of the prior-weighted ranking
rp = sorted(rows, key=lambda r: (r["rkp"], -float(r["score"])))
fig, ax = plt.subplots(figsize=(6.2, 5.4)); y = list(range(len(rp)))[::-1]
ax.barh([v + .19 for v in y], [float(r["keeps_rank_pm10_pct"]) for r in rp], height=.36, color="#2a78d6", label="Weights ±10%")
ax.barh([v - .19 for v in y], [float(r["keeps_rank_pm20_pct"]) for r in rp], height=.36, color="#eb6834", label="Weights ±20%")
ax.set_yticks(y); ax.set_yticklabels([f"{r['rkp']}. {r['protein']} {r['start']}–{r['end']}" for r in rp], fontsize=7.5)
ax.set_xlim(0, 105); ax.set_xlabel("Runs in which the segment keeps its rank (%, 1000 random weight perturbations)")
ax.legend(frameon=False, fontsize=8, loc="lower center", bbox_to_anchor=(.5, 1.0), ncol=2)
ax.xaxis.grid(True, color="#e4e3df", lw=.6); ax.set_axisbelow(True)
fig.tight_layout(); fig.savefig("figures/fig8_rank_robustness.png", dpi=220); plt.close(fig)
print("wrote", sorted(os.listdir("figures")))
