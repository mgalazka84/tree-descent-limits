"""Draw the broom and its three limiting densities. Requires Matplotlib."""
from math import erf, exp, pi, sqrt
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({
    "font.family": "DejaVu Serif", "font.size": 8.3,
    "mathtext.fontset": "stix", "pdf.fonttype": 42,
    "axes.spines.top": False, "axes.spines.right": False,
})
fig, (tree, density) = plt.subplots(
    1, 2, figsize=(5.55, 1.92), gridspec_kw={"width_ratios": [0.85, 1.4]}
)
fig.subplots_adjust(left=0.015, right=0.983, bottom=0.20, top=0.81, wspace=0.33)
ink, accent, gray = "#23272a", "#215d76", "#565856"
angles = np.radians([107, 137, 168, 199, 230, 259])
for a in angles:
    x, y = np.cos(a), np.sin(a)
    tree.plot([0, x], [0, y], c=ink, lw=0.75)
    tree.scatter([x], [y], c="white", edgecolors=ink, linewidths=0.75, s=14, zorder=3)
tree.plot([0, 1.45], [0, 0], c=ink, lw=0.85)
tree.plot([1.45, 2.14], [0, 0], c=ink, lw=0.9, ls=(0, (1.2, 2.2)))
tree.plot([2.14, 2.67], [0, 0], c=ink, lw=0.85)
tree.scatter([0.72, 1.45, 2.14, 2.67], [0]*4, c="white",
             edgecolors=ink, s=14, linewidths=0.75, zorder=3)
tree.scatter([0], [0], c=accent, edgecolors=ink, s=27, linewidths=0.6, zorder=4)
tree.annotate("root", (0, 0), (0.30, -0.73), ha="left", va="center",
              arrowprops=dict(arrowstyle="-", color=gray, lw=0.55, shrinkB=5))
tree.text(-0.67, 1.30, r"$k$ leaves", ha="center")
tree.text(1.76, 0.43, "path", ha="center")
tree.set(aspect="equal", xlim=(-1.20, 2.94), ylim=(-1.15, 1.56))
tree.axis("off")

x = np.linspace(-3.6, 3.6, 650)
normal = np.array([exp(-v*v/2)/sqrt(2*pi) for v in x])
phi = lambda u: (1 + erf(u/sqrt(2))) / 2
s = sqrt(3)
critical = np.array([sqrt(2)*(phi(sqrt(2)*v+s)-phi(sqrt(2)*v-s))/(2*s) for v in x])
density.plot(x, normal, color=ink, lw=1.1, ls="--", label=r"$c=0$ (normal)")
density.plot(x, critical, color=accent, lw=1.65, label=r"$c=1$")
density.plot([-3.6,-s,-s,s,s,3.6], [0,0,1/(2*s),1/(2*s),0,0],
             color=gray, lw=1.05, ls=":", label=r"$c=\infty$ (uniform)")
density.set(xlim=(-3.55, 3.55), ylim=(-0.009,0.44), xticks=[-3,-2,-1,0,1,2,3],
            yticks=[0,0.2,0.4], yticklabels=["0","0.2","0.4"], xlabel="standardised descents")
density.set_ylabel("density", fontsize=8, labelpad=1.6)
density.legend(loc="lower center", bbox_to_anchor=(0.5,1.03), ncol=3,
               frameon=False, fontsize=7.1, handlelength=1.35,
               handletextpad=0.35, columnspacing=0.9, borderaxespad=0.0)
density.tick_params(length=2.4, width=0.6, labelsize=7.3)
for side in ["left", "bottom"]:
    density.spines[side].set_linewidth(0.6)
fig.savefig(Path(__file__).with_name("broom_limits.pdf"),
            metadata={"Title": "A broom and its limiting descent densities",
                      "CreationDate": None, "ModDate": None})
