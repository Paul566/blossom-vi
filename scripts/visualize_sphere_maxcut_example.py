#!/usr/bin/env python3

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import generate_weighted_random_benchmarks as gen


TRIPLET_OFFSET = 0.055
ARC_STEPS = 48
FRONT_ALPHA = 0.85
BACK_ALPHA = 0.18


def load_example(path: Path) -> dict:
    return json.loads(path.read_text())


def normalize_rows(vectors: np.ndarray) -> np.ndarray:
    norms = np.linalg.norm(vectors, axis=1, keepdims=True)
    return vectors / norms


def tangent_basis(direction: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    reference = np.array([0.0, 0.0, 1.0])
    if abs(np.dot(direction, reference)) > 0.9:
        reference = np.array([0.0, 1.0, 0.0])
    basis1 = np.cross(direction, reference)
    basis1 /= np.linalg.norm(basis1)
    basis2 = np.cross(direction, basis1)
    basis2 /= np.linalg.norm(basis2)
    return basis1, basis2


def geodesic_arc(p0: np.ndarray, p1: np.ndarray, radius: float, steps: int = ARC_STEPS) -> np.ndarray:
    u0 = p0 / np.linalg.norm(p0)
    u1 = p1 / np.linalg.norm(p1)
    dot = float(np.clip(np.dot(u0, u1), -1.0, 1.0))

    if dot > 0.9999:
        ts = np.linspace(0.0, 1.0, steps)[:, None]
        points = normalize_rows((1.0 - ts) * u0[None, :] + ts * u1[None, :])
        return radius * points

    theta = np.arccos(dot)
    ts = np.linspace(0.0, 1.0, steps)
    sin_theta = np.sin(theta)
    weights0 = np.sin((1.0 - ts) * theta) / sin_theta
    weights1 = np.sin(ts * theta) / sin_theta
    points = weights0[:, None] * u0[None, :] + weights1[:, None] * u1[None, :]
    return radius * normalize_rows(points)


def build_instance_layout(data: dict) -> dict[int, np.ndarray]:
    faces = np.asarray(data["faces"], dtype=np.int32)
    points = np.asarray(data["points"], dtype=np.float64)
    radius = float(data["radius"])

    layout: dict[int, np.ndarray] = {}
    for face_idx, face in enumerate(faces):
        face_points = points[face]
        direction = normalize_rows(face_points.mean(axis=0, keepdims=True))[0]
        for side_idx in range(3):
            edge_midpoint = 0.5 * (
                face_points[side_idx] + face_points[(side_idx + 1) % 3]
            )
            edge_direction = edge_midpoint / np.linalg.norm(edge_midpoint)
            offset = edge_direction - np.dot(edge_direction, direction) * direction
            offset_norm = np.linalg.norm(offset)
            if offset_norm == 0.0:
                basis1, _ = tangent_basis(direction)
                offset = basis1
            else:
                offset /= offset_norm
            node_direction = direction + TRIPLET_OFFSET * offset
            node_direction /= np.linalg.norm(node_direction)
            node_id = 3 * face_idx + side_idx
            layout[node_id] = radius * node_direction
    return layout


def set_axes_equal(ax, radius: float) -> None:
    limit = radius * 1.15
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    ax.set_zlim(-limit, limit)
    ax.set_box_aspect((1, 1, 1))


def draw_sphere(ax, radius: float) -> None:
    u = np.linspace(0.0, 2.0 * np.pi, 80)
    v = np.linspace(0.0, np.pi, 40)
    x = radius * np.outer(np.cos(u), np.sin(v))
    y = radius * np.outer(np.sin(u), np.sin(v))
    z = radius * np.outer(np.ones_like(u), np.cos(v))
    ax.plot_surface(x, y, z, color="#9ec5c4", alpha=0.2, linewidth=0, shade=False, zorder=0)


def draw_geodesic_edges(ax, points: np.ndarray, edges: np.ndarray, radius: float, color: str, linewidth: float, alpha: float) -> None:
    for u, v in edges:
        arc = geodesic_arc(points[u], points[v], radius)
        ax.plot(arc[:, 0], arc[:, 1], arc[:, 2], color=color, linewidth=linewidth, alpha=alpha)


def view_direction(elev_deg: float, azim_deg: float) -> np.ndarray:
    elev = np.deg2rad(elev_deg)
    azim = np.deg2rad(azim_deg)
    return np.array([
        np.cos(elev) * np.cos(azim),
        np.cos(elev) * np.sin(azim),
        np.sin(elev),
    ])


def side_alpha(point: np.ndarray, camera_dir: np.ndarray) -> float:
    return FRONT_ALPHA if float(np.dot(point, camera_dir)) >= 0.0 else BACK_ALPHA


def draw_dimmed_arc(ax, p0: np.ndarray, p1: np.ndarray, radius: float, color: str, linewidth: float, camera_dir: np.ndarray) -> None:
    arc = geodesic_arc(p0, p1, radius)
    midpoint = arc[len(arc) // 2]
    alpha = side_alpha(midpoint, camera_dir)
    ax.plot(arc[:, 0], arc[:, 1], arc[:, 2], color=color, linewidth=linewidth, alpha=alpha)


def draw_overlay(ax, data: dict) -> None:
    points = np.asarray(data["points"], dtype=np.float64)
    triangulation_edges = np.asarray(data["triangulation_edges"], dtype=np.int32)
    radius = float(data["radius"])
    elev = 22.0
    azim = 34.0
    camera_dir = view_direction(elev, azim)

    draw_sphere(ax, radius)
    for u, v in triangulation_edges:
        draw_dimmed_arc(ax, points[u], points[v], radius, color="#2f6b66", linewidth=1.3, camera_dir=camera_dir)

    point_alphas = np.array([side_alpha(point, camera_dir) for point in points])
    point_colors = [(239 / 255, 111 / 255, 108 / 255, alpha) for alpha in point_alphas]
    ax.scatter(points[:, 0], points[:, 1], points[:, 2], s=55, c=point_colors, depthshade=False)

    instance_layout = build_instance_layout(data)
    zero_edges = [(int(u), int(v)) for u, v, w in data["instance_edges"] if int(w) == 0]
    cross_edges = [(int(u), int(v)) for u, v, w in data["instance_edges"] if int(w) != 0]

    instance_points = np.array([instance_layout[node_id] for node_id in sorted(instance_layout)])
    instance_alphas = np.array([side_alpha(point, camera_dir) for point in instance_points])
    instance_colors = [(242 / 255, 207 / 255, 102 / 255, alpha) for alpha in instance_alphas]
    ax.scatter(
        instance_points[:, 0],
        instance_points[:, 1],
        instance_points[:, 2],
        s=36,
        c=instance_colors,
        edgecolors="black",
        linewidths=0.5,
        depthshade=False,
    )

    for u, v in zero_edges:
        draw_dimmed_arc(
            ax,
            instance_layout[u],
            instance_layout[v],
            radius,
            color="#9aa0a6",
            linewidth=1.0,
            camera_dir=camera_dir,
        )

    for u, v in cross_edges:
        draw_dimmed_arc(
            ax,
            instance_layout[u],
            instance_layout[v],
            radius,
            color="#4c78a8",
            linewidth=1.1,
            camera_dir=camera_dir,
        )

    ax.set_title("Sphere Triangulation And Reduced MWPM Instance")
    ax.set_axis_off()
    set_axes_equal(ax, radius)
    ax.view_init(elev=elev, azim=azim)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        type=Path,
        default=Path(__file__).resolve().parent / "sphere_maxcut_example.json",
        help="Path to exported example JSON.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().parent / "sphere_maxcut_example.png",
        help="Path to output image.",
    )
    parser.add_argument("--force-regenerate", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    if args.force_regenerate or not args.input.exists():
        gen.export_sphere_max_cut_example(args.input)

    data = load_example(args.input)

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection="3d")
    draw_overlay(ax, data)
    fig.tight_layout()
    fig.savefig(args.output, dpi=180, bbox_inches="tight")
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
