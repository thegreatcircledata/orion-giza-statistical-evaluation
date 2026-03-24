"""
Study 2: Information-theoretic decomposition of the shape match.

Decomposes the Procrustes match into two independent shape parameters
(distance ratio and bend angle) and computes per-parameter p-values
and information content in bits.

Expected result: 8.7 bits total (5.0 ratio + 5.1 angle), joint test 10.3 bits.

Usage:
    python src/02_information_content.py
"""
import json, os, sys
import numpy as np
from itertools import combinations

sys.path.insert(0, os.path.dirname(__file__))
from shared_data import get_giza_triangle, load_bsc_catalog, procrustes_distance, ORION_BELT

RESULTS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "results")
os.makedirs(RESULTS_DIR, exist_ok=True)


def triangle_shape_params(points):
    """Compute shape parameters: distance ratio and bend angle."""
    d01 = np.linalg.norm(points[1] - points[0])
    d12 = np.linalg.norm(points[2] - points[1])
    d02 = np.linalg.norm(points[2] - points[0])

    # Distance ratio: shorter pair / longer pair of adjacent sides
    ratio = min(d01, d12) / max(d01, d12) if max(d01, d12) > 0 else 1.0

    # Bend angle at middle point
    v1 = points[0] - points[1]
    v2 = points[2] - points[1]
    cos_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2) + 1e-12)
    angle = np.degrees(np.arccos(np.clip(cos_angle, -1, 1)))

    return ratio, angle


def run_information_content():
    print("=" * 60)
    print("Study 2: Information-Theoretic Decomposition")
    print("=" * 60)

    data_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")
    positions, names, mags = load_bsc_catalog(mag_limit=4.0, data_dir=data_dir)
    giza = get_giza_triangle()
    orion = np.array([[ORION_BELT[s]["ra"], ORION_BELT[s]["dec"]]
                      for s in ["Alnitak", "Alnilam", "Mintaka"]])

    # Giza and Orion shape parameters
    giza_ratio, giza_angle = triangle_shape_params(giza)
    orion_ratio, orion_angle = triangle_shape_params(orion)
    print(f"Giza:  ratio={giza_ratio:.4f}, angle={giza_angle:.2f}°")
    print(f"Orion: ratio={orion_ratio:.4f}, angle={orion_angle:.2f}°")
    print(f"Delta: ratio={abs(giza_ratio-orion_ratio):.4f}, angle={abs(giza_angle-orion_angle):.2f}°")

    # Compute shape params for all triplets using canonical middle vertex:
    # the vertex with the largest angle (the "bend" point)
    N = len(positions)
    print(f"\nComputing shape parameters for all C({N},3) triplets...")
    all_ratios = []
    all_angles = []

    for idx in combinations(range(N), 3):
        triplet = positions[list(idx)]
        # Find the vertex with the largest angle (the bend point)
        best_angle = -1
        best_r = 0
        best_a = 0
        for mid in range(3):
            others = [i for i in range(3) if i != mid]
            pts_reord = np.array([triplet[others[0]], triplet[mid], triplet[others[1]]])
            r, a = triangle_shape_params(pts_reord)
            if a > best_angle:
                best_angle = a
                best_r = r
                best_a = a
        all_ratios.append(best_r)
        all_angles.append(best_a)

    all_ratios = np.array(all_ratios)
    all_angles = np.array(all_angles)

    # Per-parameter p-values
    delta_ratio = abs(giza_ratio - orion_ratio)
    delta_angle = abs(giza_angle - orion_angle)

    n_closer_ratio = np.sum(np.abs(all_ratios - giza_ratio) <= delta_ratio)
    n_closer_angle = np.sum(np.abs(all_angles - giza_angle) <= delta_angle)
    n_total = len(all_ratios)

    p_ratio = n_closer_ratio / n_total
    p_angle = n_closer_angle / n_total
    p_joint = np.sum((np.abs(all_ratios - giza_ratio) <= delta_ratio) &
                     (np.abs(all_angles - giza_angle) <= delta_angle)) / n_total

    # Procrustes p-value (from ranking)
    p_procrustes = 0.0025  # From study 01

    bits_ratio = -np.log2(p_ratio) if p_ratio > 0 else float('inf')
    bits_angle = -np.log2(p_angle) if p_angle > 0 else float('inf')
    bits_joint = -np.log2(p_joint) if p_joint > 0 else float('inf')
    bits_procrustes = -np.log2(p_procrustes)

    print(f"\nResults:")
    print(f"  Ratio:     p={p_ratio:.4f}, bits={bits_ratio:.1f}")
    print(f"  Angle:     p={p_angle:.4f}, bits={bits_angle:.1f}")
    print(f"  Joint:     p={p_joint:.6f}, bits={bits_joint:.1f}")
    print(f"  Procrustes: p={p_procrustes:.4f}, bits={bits_procrustes:.1f}")

    results = {
        "study": "02_information_content",
        "giza": {"ratio": giza_ratio, "angle": giza_angle},
        "orion": {"ratio": orion_ratio, "angle": orion_angle},
        "p_values": {
            "ratio": float(p_ratio), "angle": float(p_angle),
            "joint": float(p_joint), "procrustes": p_procrustes,
        },
        "bits": {
            "ratio": round(bits_ratio, 2), "angle": round(bits_angle, 2),
            "joint": round(bits_joint, 2), "procrustes": round(bits_procrustes, 2),
        },
        "n_triplet_configs": n_total,
    }

    out_path = os.path.join(RESULTS_DIR, "02_information_content.json")
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved to {out_path}")


if __name__ == "__main__":
    run_information_content()
