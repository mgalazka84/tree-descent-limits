# Descents in randomly labelled trees

Computations and figures for *Gaussian and non-Gaussian limits for descents in randomly labelled trees* by Marek Gałązka and Hanna Wdowicka.

Shuffle distinct numbers onto the vertices of a rooted tree. A descent occurs on an edge when the child receives the larger number. Figure 1 compares a path and a star with seven vertices. Exactly one of the 5,040 labellings gives six descents on the path. On the star, all 720 labellings with the smallest number at the centre give six descents.

## Reproduce the exact checks

From the repository root, run the following commands with Python 3.10 or later. Both scripts use the standard library.

```sh
python verification/check_exact.py
python verification/export_path_star.py > figures/path_star_n7.csv
```

The verification compares direct enumeration with a subset recurrence on all 200 rooted unlabelled tree shapes with at most eight vertices (4,894,311 shape–labelling pairs). It checks the mean, variance and symmetry of the count, the mean-square deletion error in 85 cases, and the exact star Wasserstein formula against an independent integral for 131 star sizes. Running the first command reproduces [the recorded output](results/exact_checks.txt). Each row of [the figure data](figures/path_star_n7.csv) gives the number of labellings with a specified count, along with the common denominator.

## Rebuild the figures

The vector PDFs can be opened directly. To regenerate them, use Python 3.11 or later and install the plotting packages listed in `requirements.txt`. Run from the repository root.

```sh
python -m pip install -r requirements.txt
python figures/draw_tree_examples.py
python figures/draw_broom.py
```

| Figure | Files | What it shows |
| --- | --- | --- |
| 1 | [Path and star](figures/path_star_exact.pdf) | One labelling on each tree and the exact probability of every score over all labellings. |
| 2 | [Removing two vertices](figures/cut_example.pdf) | The edges that disappear when adjacent vertices with four and three children are removed, and the resulting forest. |
| 3 | [Broom limits](figures/broom_limits.pdf) | A tree with leaves at the end of a path and three limiting densities calculated from the theorem. |
