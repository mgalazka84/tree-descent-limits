"""Exact checks for descents in uniform random labellings of small rooted trees.

Run with Python 3.10 or newer. Only the Python standard library is used.
A tree is a sorted tuple of the rooted trees hanging from its root.
"""
from collections import Counter, defaultdict
from decimal import Decimal, localcontext
from functools import cache
from itertools import permutations
from math import factorial


@cache
def size(tree):
    return 1 + sum(map(size, tree))


def rooted_shapes(max_size):
    shapes = {1: [()]}
    for n in range(2, max_size + 1):
        available = [t for k in range(1, n) for t in shapes[k]]

        def fill(left, start, selected):
            if left == 0:
                yield tuple(sorted(selected))
                return
            for i in range(start, len(available)):
                child = available[i]
                k = size(child)
                if k > left:
                    break
                yield from fill(left - k, i, selected + [child])

        shapes[n] = sorted(fill(n - 1, 0, []))
    return shapes


def parents(tree):
    pred = []

    def visit(branch, parent):
        vertex = len(pred)
        pred.append(parent)
        for child in branch:
            visit(child, vertex)

    visit(tree, -1)
    return tuple(pred)


def polynomial_by_subsets(pred):
    """Build permutations in order of increasing mark, using a subset DP."""
    n = len(pred)
    dp = [{} for _ in range(1 << n)]
    dp[0] = {0: 1}
    for mask, counts in enumerate(dp):
        for v, p in enumerate(pred):
            bit = 1 << v
            if mask & bit:
                continue
            increment = int(p >= 0 and bool(mask & (1 << p)))
            output = dp[mask | bit]
            for k, count in counts.items():
                output[k + increment] = output.get(k + increment, 0) + count
    return tuple(dp[-1].get(k, 0) for k in range(n))


def verify_cut_moment(pred, degrees, threshold):
    """Integrate the linear mark term conditional on each permutation of ranks."""
    n = len(pred)
    hubs = [v for v, d in enumerate(degrees) if d > threshold]
    hub_set = set(hubs)
    cut = [(p, v) for v, p in enumerate(pred)
           if p >= 0 and (v in hub_set or p in hub_set)]
    signed = [0] * n
    for p, v in cut:
        signed[p] += 1
        signed[v] -= 1
    for v in hubs:
        signed[v] -= degrees[v]
    formula = len(cut) + sum(b * b for b in signed)
    # If C = 2*(cut descents)-number of cut edges, the centred cut score
    # is sqrt(3)*C. Given rank r, a uniform mark has mean r/(n+1).
    total = 0
    for ranks in permutations(range(1, n + 1)):
        c = 2 * sum(ranks[v] > ranks[p] for p, v in cut) - len(cut)
        linear = sum(degrees[v] * (n + 1 - 2 * ranks[v]) for v in hubs)
        total += 3 * c * c * (n + 1) - 6 * c * linear
    denom = factorial(n) * (n + 1)
    assert total + denom * sum(degrees[v] ** 2 for v in hubs) == denom * formula
    assert formula * threshold <= 5 * sum(d * d for d in degrees)


def star_wasserstein_checks():
    """Compare the closed form with the independent integral of |F_star - F_uniform|."""
    leaf_counts = list(range(1, 129)) + [255, 511, 997]
    with localcontext() as context:
        context.prec = 65
        half = Decimal(3).sqrt()
        for k in leaf_counts:
            n = Decimal(k + 1)
            root = (n * n - 1).sqrt()
            step = 2 * half / root
            atoms = [(Decimal(i) - Decimal(k) / 2) * step for i in range(k + 1)]
            knots = [-half, *atoms, half]
            area = Decimal(0)
            for i, (left, right) in enumerate(zip(knots, knots[1:])):
                left_value = (left + half) / (2 * half) - Decimal(i) / n
                right_value = (right + half) / (2 * half) - Decimal(i) / n
                if left_value * right_value < 0:
                    average = (left_value**2 + right_value**2) / (
                        2 * (abs(left_value) + abs(right_value))
                    )
                else:
                    average = (abs(left_value) + abs(right_value)) / 2
                area += (right - left) * average
            closed_form = half / (2 * n) * (1 + (n - root)**2 / 3)
            assert abs(area - closed_form) < Decimal("1e-55"), k
    return len(leaf_counts)


def run():
    shape_count = label_count = cut_count = 0
    first_pair = None
    counts_by_size = []
    for n, shapes in rooted_shapes(8).items():
        counts_by_size.append(len(shapes))
        examples = defaultdict(dict)
        for tree in shapes:
            pred = parents(tree)
            edges = [(p, v) for v, p in enumerate(pred) if p >= 0]
            degree_counter = Counter(pred[1:])
            degrees = [degree_counter[v] for v in range(n)]
            a = sum(d * d for d in degrees)
            frequency = [0] * n
            for ranks in permutations(range(n)):
                frequency[sum(ranks[v] > ranks[p] for p, v in edges)] += 1
            polynomial = tuple(frequency)
            assert polynomial == polynomial_by_subsets(pred)
            assert sum(polynomial) == factorial(n)
            assert polynomial == polynomial[::-1]
            assert 2 * sum(k * v for k, v in enumerate(polynomial)) == (n - 1) * factorial(n)
            scaled_variance = 3 * sum((2*k - (n-1))**2 * v for k, v in enumerate(polynomial))
            assert scaled_variance == (a + 2 * degrees[0]) * factorial(n)
            if n <= 6 and n >= 2:
                for threshold in range(1, max(degrees) + 1):
                    verify_cut_moment(pred, degrees, threshold)
                    cut_count += 1
            if first_pair is None:
                group = examples[degrees[0], tuple(sorted(degrees))]
                if group and polynomial not in group:
                    old_poly, old_pred = next(iter(group.items()))
                    first_pair = (n, old_pred, old_poly, pred, polynomial)
                group[polynomial] = pred
            shape_count += 1
            label_count += factorial(n)
    assert counts_by_size == [1, 1, 2, 4, 9, 20, 48, 115]
    print("Rooted unlabeled shapes by size:", counts_by_size)
    print("Rooted shapes:", shape_count)
    print("Shape-labelling pairs:", label_count)
    print("Exact deletion-moment checks through size six:", cut_count)
    print("First equal-profile pair (size, parents, coefficients, parents, coefficients):")
    print(first_pair)
    print("Star Wasserstein formulas checked by CDF integration:", star_wasserstein_checks())


if __name__ == "__main__":
    run()
