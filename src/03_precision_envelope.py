"""
Study 3: Builder's precision envelope.

Simulates 10,000 intentional-encoding attempts with naked-eye angular error
(sigma = 0.5 deg, Schaefer 1993) to determine the distribution of Procrustes
distances achievable by intentional construction.

Expected result: Giza sits at 0.8th percentile (anomalously tight).

Usage:
    python src/03_precision_envelope.py
"""
import json, os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from shared_data import get_giza_triangle, get_orion_triangle, procrustes_distance, SIGMA_STAR_DEG

RESULTS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "results")
os.makedirs(RESULTS_DIR, exist_ok=True)

N_TRIALS = 10000


def run_precision_envelope():
    print("=" * 60)
    print("Study 3: Builder's Precision Envelope")
    print("=" * 60)

    giza = get_giza_triangle()
    orion = get_orion_triangle()
    sigma = SIGMA_STAR_DEG

    # Observed D
    D_obs, _, _ = procrustes_distance(giza, orion)
    print(f"Observed Giza-Orion D = {D_obs:.6e}")

    # Simulate intentional encoding attempts
    print(f"\nSimulating {N_TRIALS:,} intentional-encoding attempts (sigma={sigma} deg)...")
    np.random.seed(42)
    D_envelope = []

    for i in range(N_TRIALS):
        # Builder observes Orion with Gaussian noise
        noise = np.random.normal(0, sigma, orion.shape)
        observed = orion + noise

        # Builder constructs ground triangle matching observed positions
        # Then we measure how close it is to the TRUE Orion pattern
        D, _, _ = procrustes_distance(orion, observed)
        D_envelope.append(D)

    D_envelope = np.sort(D_envelope)

    # Where does Giza sit?
    percentile = 100 * np.searchsorted(D_envelope, D_obs) / N_TRIALS
    median = np.median(D_envelope)
    p5 = np.percentile(D_envelope, 5)
    p95 = np.percentile(D_envelope, 95)

    print(f"\nEnvelope statistics:")
    print(f"  Median D = {median:.4f}")
    print(f"  5th percentile = {p5:.4f}")
    print(f"  95th percentile = {p95:.4f}")
    print(f"  Giza percentile = {percentile:.1f}th")

    results = {
        "study": "03_precision_envelope",
        "sigma_star_deg": sigma,
        "n_trials": N_TRIALS,
        "D_observed": float(D_obs),
        "percentile": round(percentile, 1),
        "envelope": {
            "median": float(median),
            "p5": float(p5), "p95": float(p95),
            "p1": float(np.percentile(D_envelope, 1)),
            "p99": float(np.percentile(D_envelope, 99)),
        },
    }

    out_path = os.path.join(RESULTS_DIR, "03_precision_envelope.json")
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved to {out_path}")


if __name__ == "__main__":
    run_precision_envelope()
