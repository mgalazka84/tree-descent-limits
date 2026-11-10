"""Print the exact counts of descents for the seven-vertex path and star."""
import csv
import sys
from math import factorial

from check_exact import polynomial_by_subsets


def main():
    path = (-1, 0, 1, 2, 3, 4, 5)
    star = (-1, 0, 0, 0, 0, 0, 0)
    counts_path = polynomial_by_subsets(path)
    counts_star = polynomial_by_subsets(star)
    total = factorial(7)
    assert sum(counts_path) == sum(counts_star) == total
    assert counts_path[-1] == 1 and counts_star == (factorial(6),) * 7
    writer = csv.writer(sys.stdout, lineterminator="\n")
    writer.writerow(["descents", "path_labellings", "star_labellings", "total_labellings"])
    for score, (a, b) in enumerate(zip(counts_path, counts_star)):
        writer.writerow([score, a, b, total])


if __name__ == "__main__":
    main()
