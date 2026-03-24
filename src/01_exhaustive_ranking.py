"""
Study 1: Exhaustive Procrustes triplet ranking.

Computes the Procrustes distance from the Giza pyramid triangle to every
possible stellar triplet from the BSC (V < 4.0, N=162 stars, 695,520 triplets).
Ranks all triplets by ascending D and reports the position of the Orion Belt triplet.

Expected result: Orion Belt ranks ~1,717th (p = 0.0025).

Usage:
    python src/01_exhaustive_ranking.py

Output:
    results/01_exhaustive_ranking.json
"""
import json
import time
import os
import sys
import numpy as np
from itertools import combinations

# Add parent to path for imports
sys.path.insert(0, os.path.dirname(__file__))
from shared_data import (
    get_giza_triangle, procrustes_distance,
    load_bsc_catalog, ORION_BELT, fast_exhaustive_ranking
)

RESULTS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "results")
os.makedirs(RESULTS_DIR, exist_ok=True)


def find_orion_indices(names):
    """Find the indices of the three Orion Belt stars in the catalog."""
    orion_names = {"Alnitak", "Alnilam", "Mintaka"}
    # Also match on common alternates
    alternates = {
        "zet Ori": "Alnitak", "eps Ori": "Alnilam", "del Ori": "Mintaka",
        "Zeta Ori": "Alnitak", "Epsilon Ori": "Alnilam", "Delta Ori": "Mintaka",
    }
    indices = {}
    for i, name in enumerate(names):
        clean = name.strip()
        if clean in orion_names:
            indices[clean] = i
        for alt, proper in alternates.items():
            if alt.lower() in clean.lower():
                indices[proper] = i
    return indices


def run_exhaustive_ranking():
    """Run the full exhaustive Procrustes ranking."""
    print("=" * 60)
    print("Study 1: Exhaustive Procrustes Triplet Ranking")
    print("=" * 60)

    # Load data
    print("\nLoading stellar catalog...")
    data_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")
    positions, names, mags = load_bsc_catalog(mag_limit=4.0, data_dir=data_dir)
    N = len(positions)
    n_triplets = N * (N-1) * (N-2) // 6
    print(f"  N = {N} stars, C({N},3) = {n_triplets:,} triplets")

    # Get Giza triangle
    giza = get_giza_triangle()
    print(f"  Giza triangle loaded (3 points in local meters)")

    # Find Orion Belt indices
    orion_idx = find_orion_indices(names)
    print(f"  Orion Belt stars found at indices: {orion_idx}")
    if len(orion_idx) != 3:
        print("  WARNING: Could not find all 3 Orion Belt stars!")
        print(f"  Found: {orion_idx}")

    # Exhaustive enumeration (vectorized)
    print(f"\nEnumerating all {n_triplets:,} triplets (vectorized)...")
    t0 = time.time()

    all_D, triplet_indices = fast_exhaustive_ranking(giza, positions, allow_reflection=True)

    elapsed = time.time() - t0
    print(f"\nCompleted in {elapsed:.1f}s ({n_triplets/elapsed:.0f} triplets/sec)")

    # Find Orion Belt triplet
    orion_set = set(orion_idx.values())
    orion_D = None
    for k, idx in enumerate(triplet_indices):
        if set(idx) == orion_set:
            orion_D = float(all_D[k])
            break

    # Sort and rank
    sorted_D = np.sort(all_D)

    if orion_D is not None:
        orion_rank = int(np.searchsorted(sorted_D, orion_D)) + 1  # 1-indexed
        p_value = orion_rank / n_triplets
        print(f"\n{'='*60}")
        print(f"RESULTS:")
        print(f"  Orion Belt Procrustes D = {orion_D:.6e}")
        print(f"  Rank: {orion_rank:,} / {n_triplets:,}")
        print(f"  p-value = {p_value:.6f}")
        print(f"  Percentile: {(1-p_value)*100:.2f}th")
        print(f"{'='*60}")
    else:
        print("\nWARNING: Orion Belt triplet not found in enumeration!")
        orion_rank = None

    # Top 20 matches
    top20 = []
    sorted_indices_order = np.argsort(all_D)
    for i in range(min(20, len(sorted_indices_order))):
        idx_combo = triplet_indices[sorted_indices_order[i]]
        top20.append({
            "rank": i + 1,
            "D": float(all_D[sorted_indices_order[i]]),
            "stars": [names[j] for j in idx_combo],
        })

    # Save results
    results = {
        "study": "01_exhaustive_ranking",
        "description": "Exhaustive Procrustes ranking of all stellar triplets vs Giza",
        "n_stars": N,
        "n_triplets": n_triplets,
        "mag_limit": 4.0,
        "allow_reflection": True,
        "orion_belt": {
            "D": float(orion_D) if orion_D else None,
            "rank": orion_rank,
            "p_value": float(orion_rank / n_triplets) if orion_rank else None,
        },
        "top_20": top20,
        "elapsed_seconds": round(elapsed, 1),
    }

    out_path = os.path.join(RESULTS_DIR, "01_exhaustive_ranking.json")
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to {out_path}")

    return results


if __name__ == "__main__":
    run_exhaustive_ranking()
