"""Draw the finite path/star laws and the deletion example as vector PDFs."""
import csv
from math import cos, pi, sin
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
INK = "#263339"
ACCENT = "#15576b"
PALE = "#dde9eb"
FAINT = "#718085"
plt.rcParams.update({
    "font.family": "DejaVu Serif", "font.size": 8, "mathtext.fontset": "stix",
    "pdf.fonttype": 42, "axes.spines.top": False, "axes.spines.right": False,
})


def output(figure, stem, title):
    figure.savefig(HERE / (stem + ".pdf"),
                   metadata={"Title": title, "CreationDate": None, "ModDate": None})
    plt.close(figure)


def segment(ax, start, end, color=INK, width=0.95, style="-"):
    ax.plot([start[0], end[0]], [start[1], end[1]], color=color, lw=width,
            linestyle=style, solid_capstyle="round", zorder=1)


def node(ax, pos, label, filled=False, size=73, font=7.2):
    ax.scatter(*pos, s=size, facecolor=ACCENT if filled else "white",
               edgecolor=INK, lw=0.8, zorder=3)
    ax.text(*pos, str(label), ha="center", va="center", fontsize=font,
            color="white" if filled else INK, zorder=4)


def path_star():
    with (HERE / "path_star_n7.csv").open(newline="") as data:
        rows = list(csv.DictReader(data))
    ranks = [4, 1, 6, 2, 7, 3, 5]
    parent_sets = [[(i, i+1) for i in range(6)], [(0, i) for i in range(1, 7)]]
    assert all(sum(ranks[j] > ranks[i] for i, j in pairs) == 3 for pairs in parent_sets)
    fig = plt.figure(figsize=(6, 2.76))
    grid = fig.add_gridspec(2, 2, height_ratios=(1.26, 1), hspace=0.24, wspace=0.18,
                            left=0.105, right=0.975, top=0.90, bottom=0.16)
    for col, (title, count_key, pairs) in enumerate(zip(
        ["Path", "Star"], ["path_labellings", "star_labellings"], parent_sets
    )):
        tree, pmf = fig.add_subplot(grid[0, col]), fig.add_subplot(grid[1, col])
        if col == 0:
            pos = [(i, 0) for i in range(7)]
            tree.set_xlim(-0.55, 6.55)
            tree.set_ylim(-1.29, 1.29)
        else:
            angles = [pi + 2*pi*i/6 for i in range(6)]
            pos = [(0, 0)] + [(1.4*cos(a), 1.0*sin(a)) for a in angles]
            tree.set_xlim(-2.75, 2.75)
            tree.set_ylim(-1.28, 1.28)
        tree.set_title(title, fontsize=9.1, color=INK, pad=0)
        tree.axis("off")
        for i, j in pairs:
            scored = ranks[j] > ranks[i]
            segment(tree, pos[i], pos[j], ACCENT if scored else FAINT,
                    2.3 if scored else 0.85)
        for i, (xy, rank) in enumerate(zip(pos, ranks)):
            node(tree, xy, rank, filled=(i == 0), size=104, font=7.5)
        k = [int(row["descents"]) for row in rows]
        probs = [int(row[count_key])/int(row["total_labellings"]) for row in rows]
        # Stems with a dot keep even the path endpoint mass 1/5040 visible.
        pmf.vlines(k, 0, probs, lw=1.7, color=ACCENT)
        pmf.scatter(k, probs, s=11.5, c=ACCENT, zorder=3, clip_on=False)
        pmf.set(xticks=k, yticks=[0, .25, .5], yticklabels=["0", ".25", ".50"],
                xlim=(-.36, 6.36), ylim=(-.012, .525), xlabel="Score")
        pmf.xaxis.labelpad = 1
        if col == 0:
            pmf.set_ylabel("Probability", labelpad=2.0)
        pmf.tick_params(length=2.4, width=.6, labelsize=7.1, pad=2)
        for side in ("left", "bottom"):
            pmf.spines[side].set_linewidth(.6)
    output(fig, "path_star_exact", "Seven positions: exact path and star descent probabilities")


def cut_example():
    children = {"r": ("h", "a", "b", "c"), "h": ("d", "e", "f"),
                "a": ("x", "y"), "b": ("z",), "d": ("w",)}
    hubs = {v for v, kids in children.items() if len(kids) > 2}
    edges = [(p, v) for p, kids in children.items() for v in kids]
    retained = [(p, v) for p, v in edges if p not in hubs and v not in hubs]
    assert hubs == {"r", "h"} and len(edges) == 11 and len(retained) == 4
    original = {"r": (2.9, 3), "h": (1.2, 2), "a": (3.15, 2), "b": (4.42, 2),
                "c": (5.4, 2), "d": (0.25, 1), "e": (1.18, 1), "f": (2.08, 1),
                "x": (2.78, 1), "y": (3.60, 1), "z": (4.51, 1), "w": (0.25, 0)}
    forest = {"a": (.96, 2.05), "b": (2.64, 2.05), "d": (3.92, 2.05),
              "x": (.45, 1.10), "y": (1.48, 1.10), "z": (2.64, 1.10),
              "w": (3.92, 1.10), "c": (1.30, .13), "e": (2.45, .13), "f": (3.60, .13)}
    fig, axes = plt.subplots(1, 2, figsize=(6, 1.85), gridspec_kw={"width_ratios":[1.22, 1]})
    fig.subplots_adjust(left=.01, right=.99, bottom=.18, top=.82, wspace=.13)
    for ax, title, pos, shown in zip(axes, ["Original tree", "After deleting r and h"],
                                    [original, forest], [edges, retained]):
        ax.set_title(title, fontsize=9.1, pad=2, color=INK)
        ax.axis("off")
        ax.set(xlim=(-.12, 5.78) if pos is original else (0, 4.56), ylim=(-.4, 3.38))
        for p, v in shown:
            cut = p in hubs or v in hubs
            segment(ax, pos[p], pos[v], FAINT if cut else ACCENT,
                    .95 if cut else 1.4, (0, (3, 2)) if cut else "-")
        for name, xy in pos.items():
            node(ax, xy, name, filled=(name in hubs), size=75, font=7)
    fig.legend(handles=[Line2D([0],[0],color=FAINT,lw=.95,ls=(0,(3,2)),label="Deleted edge"),
                        Line2D([0],[0],color=ACCENT,lw=1.4,label="Retained edge")],
               loc="lower center", bbox_to_anchor=(.5,.005), ncol=2,
               frameon=False, fontsize=7.3, handlelength=2)
    output(fig, "cut_example", "Removing the two vertices with more than two children")


if __name__ == "__main__":
    path_star()
    cut_example()
