#!/usr/bin/env python3

import json
import math
from pathlib import Path

import networkx as nx
import numpy as np
from scipy.spatial import ConvexHull, Delaunay


NUM_RANDOM_INSTANCES = 10
MAX_RANDOM_EDGES = 10_000_00
NUM_DELAUNAY_INSTANCES = 30
MAX_DELAUNAY_NODES = 3_000_000
NUM_GEOMETRIC_INSTANCES = 20
MAX_GEOMETRIC_NODES = 100_000
NUM_SPHERE_MAXCUT_INSTANCES = 20
MAX_SPHERE_MAXCUT_NODES = 1000_000
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


def read_tsp_point_cloud(path: Path) -> np.ndarray:
    points = []
    in_section = False

    with path.open() as f:
        for line in f:
            line = line.strip()

            if line.startswith("NODE_COORD_SECTION"):
                in_section = True
                continue

            if not in_section:
                continue

            if line.startswith("EOF"):
                break

            parts = line.split()
            if len(parts) >= 3:
                points.append((float(parts[1]), float(parts[2])))

    return np.asarray(points, dtype=np.float64)


def sample_points_on_sphere(num_nodes: int, radius: float, rng: np.random.Generator) -> np.ndarray:
    points = rng.normal(size=(num_nodes, 3))
    norms = np.linalg.norm(points, axis=1, keepdims=True)
    return radius * points / norms


def sphere_radius_fixed(_: int) -> float:
    return 1_000_000.0


def sphere_radius_sqrt_n(num_nodes: int) -> float:
    return math.sqrt(num_nodes)


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
    rectangle_size: float,
) -> WeightedGraph:
    rng = np.random.default_rng(seed)
    graph_nx = nx.random_geometric_graph(num_nodes, radius, dim=2, seed=seed)
    positions = np.array([graph_nx.nodes[i]["pos"] for i in range(num_nodes)]) * rectangle_size

    geometric_edges = np.array(graph_nx.edges(), dtype=np.int32)
    geometric_edges.sort(axis=1)

    matching_edges = random_perfect_matching_edges(num_nodes, rng)
    edges = np.vstack((geometric_edges, matching_edges))
    edges = np.unique(edges, axis=0)

    deltas = positions[edges[:, 0]] - positions[edges[:, 1]]
    distances = np.linalg.norm(deltas, axis=1)

    graph = WeightedGraph(num_nodes, edges[:, 0].copy(), edges[:, 1].copy())
    graph.weights = np.rint(distances).astype(np.int32)
    return graph


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


def generate_tsp_family(tests_weighted_dir: Path) -> None:
    family_dir = tests_weighted_dir / "tsp"
    family_dir.mkdir(parents=True, exist_ok=True)
    tsp_dir = Path(__file__).resolve().parents[1] / "tests-aux" / "vlsi_tsp" / "vlsi"

    for tsp_path in sorted(tsp_dir.glob("*.tsp")):
        points = read_tsp_point_cloud(tsp_path)
        if len(points) < 50_000:
            continue
        if len(points) % 2 != 0:
            points = points[:-1]

        graph = delaunay_graph(points)
        filename = family_dir / f"{tsp_path.stem}-{len(graph.us)}"
        export_graph(graph, filename)
        print(f"tsp: wrote {filename.name} (n={graph.num_nodes}, m={len(graph.us)})")


def generate_geometric_family(
    tests_weighted_dir: Path,
    family_name: str,
    rectangle_size_fn,
    family_seed: int,
) -> None:
    family_dir = tests_weighted_dir / family_name
    family_dir.mkdir(parents=True, exist_ok=True)

    for idx in range(1, NUM_GEOMETRIC_INSTANCES + 1):
        num_nodes = idx * MAX_GEOMETRIC_NODES // NUM_GEOMETRIC_INSTANCES
        radius = geometric_radius(num_nodes)
        graph = geometric_graph(num_nodes, radius, family_seed + idx, rectangle_size_fn(num_nodes))
        filename = family_dir / f"geometric-{num_nodes}-{len(graph.us)}"
        export_graph(graph, filename)
        mean_degree = 2 * len(graph.us) / num_nodes if num_nodes > 0 else 0.0
        print(
            f"{family_name}: wrote {filename.name} "
            f"(n={num_nodes}, m={len(graph.us)}, mean_degree={mean_degree:.2f})"
        )


def sphere_triangulation_faces(points: np.ndarray) -> np.ndarray:
    hull = ConvexHull(points)
    faces = hull.simplices.astype(np.int32, copy=False)

    face_normals = np.cross(
        points[faces[:, 1]] - points[faces[:, 0]],
        points[faces[:, 2]] - points[faces[:, 0]],
    )
    outward_mask = np.einsum("ij,ij->i", face_normals, points[faces[:, 0]]) < 0
    faces[outward_mask] = faces[outward_mask][:, [0, 2, 1]]
    return faces


def sphere_max_cut_instance_data(
    points: np.ndarray,
    rng: np.random.Generator,
    max_abs_weight: int,
) -> dict:
    faces = sphere_triangulation_faces(points)
    num_faces = len(faces)

    city_base = 3 * np.arange(num_faces, dtype=np.int32)
    zero_us = np.concatenate((city_base, city_base + 1, city_base + 2))
    zero_vs = np.concatenate((city_base + 1, city_base + 2, city_base))
    zero_weights = np.zeros(3 * num_faces, dtype=np.int32)

    edge_to_face_side: dict[tuple[int, int], list[tuple[int, int]]] = {}
    for face_idx, face in enumerate(faces):
        for side_idx in range(3):
            u = int(face[side_idx])
            v = int(face[(side_idx + 1) % 3])
            key = (u, v) if u < v else (v, u)
            edge_to_face_side.setdefault(key, []).append((face_idx, side_idx))

    cross_us = []
    cross_vs = []
    cross_weights = []
    cross_face_pairs = []
    triangulation_edges = []
    triangulation_edge_weights = []
    for (u, v), incidences in edge_to_face_side.items():
        if len(incidences) != 2:
            raise ValueError(f"expected exactly two incident faces for edge {(u, v)}, got {len(incidences)}")

        (face_a, side_a), (face_b, side_b) = incidences
        weight = int(rng.integers(0, max_abs_weight + 1))
        cross_us.append(3 * face_a + side_a)
        cross_vs.append(3 * face_b + side_b)
        cross_weights.append(weight)
        cross_face_pairs.append((face_a, face_b))
        triangulation_edges.append((u, v))
        triangulation_edge_weights.append(weight)

    cross_us_arr = np.asarray(cross_us, dtype=np.int32)
    cross_vs_arr = np.asarray(cross_vs, dtype=np.int32)
    cross_weights_arr = np.asarray(cross_weights, dtype=np.int32)

    graph = WeightedGraph(
        3 * num_faces,
        np.concatenate((zero_us, cross_us_arr)),
        np.concatenate((zero_vs, cross_vs_arr)),
    )
    graph.weights = np.concatenate((zero_weights, cross_weights_arr))
    return {
        "graph": graph,
        "faces": faces,
        "triangulation_edges": np.asarray(triangulation_edges, dtype=np.int32),
        "triangulation_edge_weights": np.asarray(triangulation_edge_weights, dtype=np.int32),
        "cross_face_pairs": np.asarray(cross_face_pairs, dtype=np.int32),
    }


def max_cut_instance_from_sphere_triangulation(
    points: np.ndarray,
    rng: np.random.Generator,
    max_abs_weight: int,
) -> WeightedGraph:
    return sphere_max_cut_instance_data(points, rng, max_abs_weight)["graph"]


def export_sphere_max_cut_example(
    output_path: Path,
    num_points: int = 10,
    radius: float = 1_000_000.0,
    seed: int = 123,
    max_abs_weight: int = 1_000_000,
) -> None:
    rng = np.random.default_rng(seed)
    points = sample_points_on_sphere(num_points, radius, rng)
    data = sphere_max_cut_instance_data(points, rng, max_abs_weight)
    graph = data["graph"]
    faces = data["faces"]
    triangulation_edges = data["triangulation_edges"]
    triangulation_edge_weights = data["triangulation_edge_weights"]
    cross_face_pairs = data["cross_face_pairs"]

    payload = {
        "num_points": num_points,
        "radius": radius,
        "seed": seed,
        "max_abs_weight": max_abs_weight,
        "points": points.tolist(),
        "faces": faces.tolist(),
        "triangulation_edges": triangulation_edges.tolist(),
        "triangulation_edge_weights": triangulation_edge_weights.tolist(),
        "instance_num_nodes": int(graph.num_nodes),
        "instance_edges": [
            [int(u), int(v), int(w)]
            for u, v, w in zip(graph.us, graph.vs, graph.weights)
        ],
        "cross_face_pairs": cross_face_pairs.tolist(),
    }
    output_path.write_text(json.dumps(payload, indent=2))


def generate_sphere_max_cut_family(
    tests_weighted_dir: Path,
    family_name: str,
    radius_fn,
    max_abs_weight: int,
    family_seed: int,
) -> None:
    family_dir = tests_weighted_dir / family_name
    family_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(family_seed)

    for idx in range(1, NUM_SPHERE_MAXCUT_INSTANCES + 1):
        num_points = idx * MAX_SPHERE_MAXCUT_NODES // NUM_SPHERE_MAXCUT_INSTANCES
        radius = radius_fn(num_points)
        points = sample_points_on_sphere(num_points, radius, rng)
        graph = max_cut_instance_from_sphere_triangulation(points, rng, max_abs_weight)
        filename = family_dir / f"sphere-maxcut-{num_points}-{len(graph.us)}"
        export_graph(graph, filename)
        print(
            f"{family_name}: wrote {filename.name} "
            f"(points={num_points}, nodes={graph.num_nodes}, m={len(graph.us)})"
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
    # generate_tsp_family(tests_weighted_dir)
    # generate_geometric_family(
    #     tests_weighted_dir,
    #     "geometric-big-weights",
    #     lambda n : 1_000_000,
    #     family_seed=30_000,
    # )
    # generate_geometric_family(
    #     tests_weighted_dir,
    #     "geometric-small-weights",
    #     lambda n : np.sqrt(n),
    #     family_seed=40_000,
    # )
    generate_sphere_max_cut_family(
        tests_weighted_dir,
        "maxcut-big-weights",
        sphere_radius_fixed,
        1_000_000,
        family_seed=50_000,
    )
    generate_sphere_max_cut_family(
        tests_weighted_dir,
        "maxcut-small-weights",
        sphere_radius_sqrt_n,
        1,
        family_seed=60_000,
    )


if __name__ == "__main__":
    main()
