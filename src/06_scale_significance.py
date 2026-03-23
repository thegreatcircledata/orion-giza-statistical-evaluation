"""
Study 6: Scale significance test.

Tests whether the Procrustes scale factor (ground-to-sky ratio) constitutes
independent evidence of intentional design via Monte Carlo simulation.

Expected result: Scale not significant (p = 0.219).

Usage:
    python src/06_scale_significance.py
"""
import json, os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from shared_data import get_giza_triangle, get_orion_triangle, procrustes_distance

RESULTS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "results")
os.makedirs(RESULTS_DIR, exist_ok=True)

N_TRIALS = 10000
PLATEAU_SIZE_M = 2000  # Approximate Giza plateau buildable area


def run_scale_test():
    print("=" * 60)
    print("Study 6: Scale Significance Test")
    print("=" * 60)

    giza = get_giza_triangle()
    orion = get_orion_triangle()

    # Observed scale
    _, _, obs_scale = procrustes_distance(giza, orion)
    print(f"Observed Procrustes scale: {obs_scale:.1f} m/rad")

    # Monte Carlo: random triangles on a 2km plateau
    print(f"\nMC: {N_TRIALS:,} random triangles on {PLATEAU_SIZE_M}m plateau...")
    np.random.seed(42)
    mc_scales = []

    for _ in range(N_TRIALS):
        # 3 random points on a plateau
        pts = np.random.uniform(0, PLATEAU_SIZE_M, (3, 2))
        _, _, scale = procrustes_distance(pts, orion)
        mc_scales.append(scale)

    mc_scales = np.sort(mc_scales)
    percentile = 100 * np.searchsorted(mc_scales, obs_scale) / N_TRIALS
    p_value = percentile / 100

    print(f"\nResults:")
    print(f"  Observed scale: {obs_scale:.1f} m/rad")
    print(f"  MC median: {np.median(mc_scales):.1f} m/rad")
    print(f"  Percentile: {percentile:.1f}th")
    print(f"  p-value: {p_value:.4f}")
    print(f"  Significant: {'Yes' if p_value < 0.05 else 'No'}")

    results = {
        "study": "06_scale_significance",
        "observed_scale_m_per_rad": round(float(obs_scale), 1),
        "mc_trials": N_TRIALS,
        "plateau_size_m": PLATEAU_SIZE_M,
        "mc_median": round(float(np.median(mc_scales)), 1),
        "percentile": round(percentile, 1),
        "p_value": round(p_value, 4),
    }

    out_path = os.path.join(RESULTS_DIR, "06_scale_significance.json")
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved to {out_path}")


if __name__ == "__main__":
    run_scale_test()
