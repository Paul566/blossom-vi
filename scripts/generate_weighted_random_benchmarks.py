#!/usr/bin/env python3

import math
from pathlib import Path

import networkx as nx
import numpy as np
from scipy.spatial import Delaunay


NUM_RANDOM_INSTANCES = 10
MAX_RANDOM_EDGES = 10_000_00
NUM_DELAUNAY_INSTANCES = 30
MAX_DELAUNAY_NODES = 3_000_000
NUM_GEOMETRIC_INSTANCES = 20
MAX_GEOMETRIC_NODES = 100_000
FIXED_RECTANGLE_SIDE = 1_000_000.0
WEIGHT_LOW = -1_000_000
WEIGHT_HIGH_EXCLUSIVE = 1_000_001


class WeightedGraph:
    def __init__(self, num_nodes: int, us: np.ndarray, vs: np.ndarray):
        self.num_nodes = num_nodes
        self.us = us
        self.vs = vs
        self.weights = np.zeros(len(us), dtype=np.int32)


def cumulative_edges_before(n: int, u: np.ndarray) -> np.ndarray:
    return u * (2 * n - u - 1) // 2


def ranks_to_edges(n: int, ranks: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    row_starts = cumulative_edges_before(n, np.arange(n, dtype=np.int64))
    u = np.searchsorted(row_starts, ranks, side="right") - 1
    offsets = ranks - row_starts[u]
    v = u + 1 + offsets
    return u.astype(np.int32), v.astype(np.int32)


class RandomGraphGenerator:
    def __init__(self, seed: int):
        self.rng = np.random.default_rng(seed)

    def matching_plus_random_graph(self, num_nodes: int, num_edges: int) -> WeightedGraph:
        if num_nodes % 2 != 0:
            raise ValueError("num_nodes must be even to contain a perfect matching")
        min_edges_needed = num_nodes // 2
        if num_edges < min_edges_needed:
            raise ValueError("num_edges must be at least num_nodes / 2")

        permutation = self.rng.permutation(num_nodes).astype(np.int32)
        matching_us = np.minimum(permutation[0::2], permutation[1::2])
        matching_vs = np.maximum(permutation[0::2], permutation[1::2])
        matching_ranks = (
            cumulative_edges_before(num_nodes, matching_us.astype(np.int64))
            + (matching_vs.astype(np.int64) - matching_us.astype(np.int64) - 1)
        )

        remaining_edges = num_edges - min_edges_needed
        if remaining_edges > 0:
            extra_ranks = np.empty(0, dtype=np.int64)
            sorted_matching_ranks = np.sort(matching_ranks)
            while len(extra_ranks) < remaining_edges:
                batch_size = max(2 * (remaining_edges - len(extra_ranks)), 1024)
                raw_us = self.rng.integers(num_nodes, size=batch_size, dtype=np.int32)
                raw_vs = self.rng.integers(num_nodes, size=batch_size, dtype=np.int32)
                valid_mask = raw_us != raw_vs
                us = np.minimum(raw_us[valid_mask], raw_vs[valid_mask])
                vs = np.maximum(raw_us[valid_mask], raw_vs[valid_mask])
                candidate_ranks = (
                    cumulative_edges_before(num_nodes, us.astype(np.int64))
                    + (vs.astype(np.int64) - us.astype(np.int64) - 1)
                )
                candidate_ranks = np.unique(candidate_ranks)
                keep_mask = ~np.isin(candidate_ranks, sorted_matching_ranks, assume_unique=True)
                candidate_ranks = candidate_ranks[keep_mask]
                extra_ranks = np.unique(np.concatenate((extra_ranks, candidate_ranks)))

            ranks = np.sort(np.concatenate((sorted_matching_ranks, extra_ranks[:remaining_edges])))
        else:
            ranks = np.sort(matching_ranks)

        us, vs = ranks_to_edges(num_nodes, ranks)
        return WeightedGraph(num_nodes, us, vs)


def make_graph_weighted(graph: WeightedGraph) -> None:
    graph.weights = np.random.randint(
        WEIGHT_LOW,
        WEIGHT_HIGH_EXCLUSIVE,
        size=len(graph.us),
        dtype=np.int32,
    )


def export_graph(graph: WeightedGraph, path: Path) -> None:
    data = np.column_stack((graph.us, graph.vs, graph.weights))
    with path.open("w") as f:
        f.write(f"{graph.num_nodes} {len(graph.us)}\n")
        np.savetxt(f, data, fmt="%d %d %d")


def save_random_graph(output_path: Path, graph_generator: RandomGraphGenerator, num_nodes: int, num_edges: int) -> None:
    graph = graph_generator.matching_plus_random_graph(num_nodes, num_edges)
    make_graph_weighted(graph)
    export_graph(graph, output_path)


def dense_num_nodes(num_edges: int) -> int:
    num_nodes = max(2, round(math.sqrt(20 * num_edges)))
    if num_nodes % 2 != 0:
        num_nodes += 1
    while num_nodes * (num_nodes - 1) // 2 < num_edges:
        num_nodes += 2
    return num_nodes


def sparse_num_nodes(num_edges: int) -> int:
    num_nodes = max(2, math.ceil(num_edges / 5))
    if num_nodes % 2 != 0:
        num_nodes += 1
    while num_nodes * (num_nodes - 1) // 2 < num_edges:
        num_nodes += 2
    return num_nodes


def delaunay_graph(points: np.ndarray) -> WeightedGraph:
    tri = Delaunay(points)
    simplices = tri.simplices.astype(np.int32, copy=False)

    edges = np.vstack(
        (
            simplices[:, [0, 1]],
            simplices[:, [1, 2]],
            simplices[:, [2, 0]],
        )
    )
    edges.sort(axis=1)
    edges = np.unique(edges, axis=0)

    deltas = points[edges[:, 0]] - points[edges[:, 1]]
    weights = np.rint(np.linalg.norm(deltas, axis=1)).astype(np.int32)

    graph = WeightedGraph(len(points), edges[:, 0].copy(), edges[:, 1].copy())
    graph.weights = weights
    return graph


def fixed_rectangle_points(num_nodes: int, rng: np.random.Generator) -> np.ndarray:
    return rng.uniform(0.0, FIXED_RECTANGLE_SIDE, size=(num_nodes, 2))


def sqrt_n_rectangle_points(num_nodes: int, rng: np.random.Generator) -> np.ndarray:
    side = math.sqrt(num_nodes)
    return rng.uniform(0.0, side, size=(num_nodes, 2))


def geometric_radius(num_nodes: int) -> float:
    return math.sqrt(10.0 / (math.pi * num_nodes))


def random_perfect_matching_edges(num_nodes: int, rng: np.random.Generator) -> np.ndarray:
    if num_nodes % 2 != 0:
        raise ValueError("num_nodes must be even to contain a perfect matching")

    permutation = rng.permutation(num_nodes).astype(np.int32)
    edges = np.column_stack((permutation[0::2], permutation[1::2]))
    edges.sort(axis=1)
    return edges


def scaling_factor_from_geometric_edges(
    geometric_distances: np.ndarray,
    weight_scale_fn,
) -> float:
    if len(geometric_distances) == 0:
        return 0.0

    scaled_geometric_distances = weight_scale_fn(geometric_distances)
    nonzero_mask = geometric_distances != 0
    if np.any(nonzero_mask):
        idx = int(np.flatnonzero(nonzero_mask)[0])
        return float(scaled_geometric_distances[idx] / geometric_distances[idx])
    return 0.0


def geometric_graph(
    num_nodes: int,
    radius: float,
    seed: int,
    weight_scale_fn,
) -> WeightedGraph:
    rng = np.random.default_rng(seed)
    graph_nx = nx.random_geometric_graph(num_nodes, radius, dim=2, seed=seed)
    positions = np.array([graph_nx.nodes[i]["pos"] for i in range(num_nodes)]) * np.sqrt(num_nodes)

    geometric_edges = np.array(graph_nx.edges(), dtype=np.int32)
    if len(geometric_edges) == 0:
        geometric_edges = np.empty((0, 2), dtype=np.int32)
        scale_factor = 0.0
    else:
        geometric_edges.sort(axis=1)
        geometric_deltas = positions[geometric_edges[:, 0]] - positions[geometric_edges[:, 1]]
        geometric_distances = np.linalg.norm(geometric_deltas, axis=1)
        scale_factor = 1.
        if weight_scale_fn is not None:
            scale_factor = scaling_factor_from_geometric_edges(geometric_distances, weight_scale_fn)

    matching_edges = random_perfect_matching_edges(num_nodes, rng)
    edges = np.vstack((geometric_edges, matching_edges))
    edges = np.unique(edges, axis=0)

    deltas = positions[edges[:, 0]] - positions[edges[:, 1]]
    distances = np.linalg.norm(deltas, axis=1)
    scaled_distances = distances * scale_factor

    graph = WeightedGraph(num_nodes, edges[:, 0].copy(), edges[:, 1].copy())
    graph.weights = np.rint(scaled_distances).astype(np.int32)
    return graph


def scale_weights_to_max(distances: np.ndarray) -> np.ndarray:
    max_distance = np.max(distances)
    if max_distance == 0:
        return np.zeros_like(distances)
    return distances * (1_000_000.0 / max_distance)


# def scale_weights_to_mean(distances: np.ndarray) -> np.ndarray:
#     mean_distance = np.mean(distances)
#     if mean_distance == 0:
#         return np.zeros_like(distances)
#     return distances * (5.0 / mean_distance)


def generate_random_family(
    tests_weighted_dir: Path,
    family_name: str,
    num_nodes_fn,
    graph_generator: RandomGraphGenerator,
) -> None:
    family_dir = tests_weighted_dir / family_name
    family_dir.mkdir(parents=True, exist_ok=True)

    for idx in range(1, NUM_RANDOM_INSTANCES + 1):
        num_edges = idx * MAX_RANDOM_EDGES // NUM_RANDOM_INSTANCES
        num_nodes = num_nodes_fn(num_edges)
        filename = f"random-{num_nodes}-{num_edges}"
        family_output = family_dir / filename

        save_random_graph(family_output, graph_generator, num_nodes, num_edges)

        mean_degree = 2 * num_edges / num_nodes
        print(
            f"{family_name}: wrote {filename} "
            f"(n={num_nodes}, m={num_edges}, mean_degree={mean_degree:.2f})"
        )


def generate_delaunay_family(
    tests_weighted_dir: Path,
    family_name: str,
    points_fn,
    family_seed: int,
) -> None:
    family_dir = tests_weighted_dir / family_name
    family_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(family_seed)

    for idx in range(1, NUM_DELAUNAY_INSTANCES + 1):
        num_nodes = idx * MAX_DELAUNAY_NODES // NUM_DELAUNAY_INSTANCES
        points = points_fn(num_nodes, rng)
        graph = delaunay_graph(points)
        filename = family_dir / f"delaunay-{num_nodes}-{len(graph.us)}"
        export_graph(graph, filename)
        print(f"{family_name}: wrote {filename.name} (n={num_nodes}, m={len(graph.us)})")


def generate_geometric_family(
    tests_weighted_dir: Path,
    family_name: str,
    weight_scale_fn,
    family_seed: int,
) -> None:
    family_dir = tests_weighted_dir / family_name
    family_dir.mkdir(parents=True, exist_ok=True)

    for idx in range(1, NUM_GEOMETRIC_INSTANCES + 1):
        num_nodes = idx * MAX_GEOMETRIC_NODES // NUM_GEOMETRIC_INSTANCES
        radius = geometric_radius(num_nodes)
        graph = geometric_graph(num_nodes, radius, family_seed + idx, weight_scale_fn)
        filename = family_dir / f"geometric-{num_nodes}-{len(graph.us)}"
        export_graph(graph, filename)
        mean_degree = 2 * len(graph.us) / num_nodes if num_nodes > 0 else 0.0
        print(
            f"{family_name}: wrote {filename.name} "
            f"(n={num_nodes}, m={len(graph.us)}, mean_degree={mean_degree:.2f})"
        )


def main() -> None:
    repo_root = Path(__file__).resolve().parents[1]
    tests_weighted_dir = repo_root / "tests-weighted"

    np.random.seed(0)
    graph_generator = RandomGraphGenerator(seed=0)

    # generate_random_family(tests_weighted_dir, "dense-random", dense_num_nodes, graph_generator)
    # generate_random_family(tests_weighted_dir, "sparse-random", sparse_num_nodes, graph_generator)
    # generate_delaunay_family(
    #     tests_weighted_dir,
    #     "delaunay-big-weights",
    #     fixed_rectangle_points,
    #     family_seed=10_000,
    # )
    # generate_delaunay_family(
    #     tests_weighted_dir,
    #     "delaunay-small-weights",
    #     sqrt_n_rectangle_points,
    #     family_seed=20_000,
    # )
    # generate_geometric_family(
    #     tests_weighted_dir,
    #     "geometric-big-weights",
    #     scale_weights_to_max,
    #     family_seed=30_000,
    # )
    generate_geometric_family(
        tests_weighted_dir,
        "geometric-small-weights",
        None,
        family_seed=40_000,
    )


if __name__ == "__main__":
    main()
