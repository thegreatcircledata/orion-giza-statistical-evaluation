"""
Study 5: Reflection and chirality analysis.

Tests whether the Giza-Orion Procrustes match requires reflection (mirror),
addressing Fairall's (1999) objection. Runs Procrustes with and without
reflection permitted and compares ranks.

Expected result: No reflection needed (det(R)=+1), unreflected rank better (1,380 vs 1,716).

Usage:
    python src/05_reflection_chirality.py
"""
import json, os, sys
import numpy as np
from itertools import combinations

sys.path.insert(0, os.path.dirname(__file__))
from shared_data import get_giza_triangle, load_bsc_catalog, procrustes_distance, ORION_BELT, fast_exhaustive_ranking


RESULTS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "results")
os.makedirs(RESULTS_DIR, exist_ok=True)


def run_reflection_test():
    print("=" * 60)
    print("Study 5: Reflection and Chirality Analysis")
    print("=" * 60)

    data_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")
    positions, names, mags = load_bsc_catalog(mag_limit=4.0, data_dir=data_dir)
    giza = get_giza_triangle()

    # Use CATALOG coordinates for Orion (same as exhaustive ranking)
    orion_indices = []
    for sname in ["Alnitak", "Alnilam", "Mintaka"]:
        for i, n in enumerate(names):
            if sname.lower() in n.lower():
                orion_indices.append(i)
                break
    orion = positions[orion_indices]

    # Apply equirectangular projection (same as in fast_exhaustive_ranking)
    mean_dec = orion[:, 1].mean()
    orion_proj = orion.copy()
    orion_proj[:, 0] *= np.cos(np.radians(mean_dec))

    # Test Giza-Orion with and without reflection (best of 6 permutations)
    D_reflect = float('inf')
    R_reflect = None
    D_no_reflect = float('inf')
    R_no_reflect = None
    best_perm_reflect = None
    best_perm_no_reflect = None
    for perm in [(0,1,2), (0,2,1), (1,0,2), (1,2,0), (2,0,1), (2,1,0)]:
        op = orion_proj[list(perm)]
        Dr, Rr, _ = procrustes_distance(giza, op, allow_reflection=True)
        Dnr, Rnr, _ = procrustes_distance(giza, op, allow_reflection=False)
        if Dr < D_reflect:
            D_reflect = Dr
            R_reflect = Rr
            best_perm_reflect = perm
        if Dnr < D_no_reflect:
            D_no_reflect = Dnr
            R_no_reflect = Rnr
            best_perm_no_reflect = perm

    det_reflect = np.linalg.det(R_reflect)
    det_no_reflect = np.linalg.det(R_no_reflect)

    # Rotation angle
    angle = np.degrees(np.arctan2(R_no_reflect[1, 0], R_no_reflect[0, 0]))

    print(f"With reflection:    D = {D_reflect:.6e}, det(R) = {det_reflect:.1f}")
    print(f"Without reflection: D = {D_no_reflect:.6e}, det(R) = {det_no_reflect:.1f}")
    print(f"Rotation angle: {angle:.1f}°")
    print(f"Reflection needed: {'Yes' if det_reflect < 0 else 'No'}")

    # Rank both modes across all triplets (vectorized)
    print(f"\nRanking all triplets in both modes (vectorized)...")
    N = len(positions)
    n_triplets = N * (N - 1) * (N - 2) // 6

    print("  Mode 1: with reflection...")
    D_all_reflect, _ = fast_exhaustive_ranking(giza, positions, allow_reflection=True)
    print("  Mode 2: without reflection...")
    D_all_no_reflect, _ = fast_exhaustive_ranking(giza, positions, allow_reflection=False)

    rank_reflect = int(np.sum(D_all_reflect < D_reflect)) + 1
    rank_no_reflect = int(np.sum(D_all_no_reflect < D_no_reflect)) + 1

    # Count how many top-100 (reflected mode) need reflection
    n_top100_need_reflect = 0
    sorted_idx = np.argsort(D_all_reflect)[:100]
    for i in sorted_idx:
        if D_all_reflect[i] < D_all_no_reflect[i]:
            n_top100_need_reflect += 1

    print(f"\nRank (with reflect): {rank_reflect:,}")
    print(f"Rank (no reflect):   {rank_no_reflect:,}")
    print(f"Top-100 needing reflection: {n_top100_need_reflect}/100")

    results = {
        "study": "05_reflection_chirality",
        "n_stars": N,
        "n_triplets": n_triplets,
        "giza_orion": {
            "D_with_reflect": float(D_reflect),
            "D_without_reflect": float(D_no_reflect),
            "det_R_reflect": float(det_reflect),
            "det_R_no_reflect": float(det_no_reflect),
            "rotation_angle_deg": round(angle, 1),
            "reflection_needed": bool(det_reflect < 0),
        },
        "rank_with_reflect": rank_reflect,
        "rank_without_reflect": rank_no_reflect,
        "p_with_reflect": round(rank_reflect / n_triplets, 6),
        "p_without_reflect": round(rank_no_reflect / n_triplets, 6),
        "top100_need_reflection": n_top100_need_reflect,
    }

    out_path = os.path.join(RESULTS_DIR, "05_reflection_chirality.json")
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved to {out_path}")


if __name__ == "__main__":
    run_reflection_test()
