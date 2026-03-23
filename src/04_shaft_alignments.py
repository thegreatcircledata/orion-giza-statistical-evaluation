"""
Study 4: Great Pyramid shaft alignment analysis.

Tests all four internal shafts against the bright-star catalog at 2560 BCE
using astropy for precession modeling. Reports individual matches and
epoch convergence.

Expected result: All 4 targets found but epochs span 1,325 years (no convergence).

Usage:
    python src/04_shaft_alignments.py

Requires: astropy
"""
import json, os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from shared_data import SHAFT_ANGLES, SHAFT_CLAIMS, GIZA_LATITUDE, load_bsc_catalog

RESULTS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "results")
os.makedirs(RESULTS_DIR, exist_ok=True)

CONSTRUCTION_EPOCH = -2559  # 2560 BCE in astronomical year numbering
TOLERANCE_DEG = 2.0  # Search window for each shaft


def shaft_target_dec(angle_deg, latitude_deg):
    """Compute target declination for a shaft at given angle from horizontal."""
    # For a south-pointing shaft: dec = -(90 - lat - angle)
    # For a north-pointing shaft: dec = 90 - lat + angle
    # Simplified: dec = lat - (90 - angle) for south, dec = lat + (90 - angle) for north
    # More precisely: altitude = angle, azimuth = 180 (south) or 0 (north)
    # dec = arcsin(sin(lat)*sin(alt) + cos(lat)*cos(alt)*cos(az))
    lat_r = np.radians(latitude_deg)
    alt_r = np.radians(angle_deg)

    # South shafts (azimuth = 180)
    # dec = arcsin(sin(lat)*sin(alt) - cos(lat)*cos(alt))
    # North shafts (azimuth = 0)
    # dec = arcsin(sin(lat)*sin(alt) + cos(lat)*cos(alt))

    return angle_deg, lat_r, alt_r


def run_shaft_analysis():
    print("=" * 60)
    print("Study 4: Shaft Alignment Analysis")
    print("=" * 60)

    try:
        from astropy.coordinates import SkyCoord, FK5
        from astropy.time import Time
        import astropy.units as u
    except ImportError:
        print("ERROR: astropy is required for precession modeling.")
        print("Install with: pip install astropy")
        return

    data_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")
    positions, names, mags = load_bsc_catalog(mag_limit=4.0, data_dir=data_dir)

    # Compute shaft target declinations
    lat = GIZA_LATITUDE
    targets = {}
    for shaft, angle in SHAFT_ANGLES.items():
        if "south" in shaft:
            dec = -(90 - lat - angle)  # Below celestial equator
        else:
            dec = 90 - (90 - lat - angle)  # Near pole
        targets[shaft] = {"angle": angle, "target_dec": dec, "claimed": SHAFT_CLAIMS[shaft]}

    print(f"\nShaft target declinations at lat={lat}°:")
    for shaft, info in targets.items():
        print(f"  {shaft}: angle={info['angle']}°, target dec={info['target_dec']:.2f}°, claimed={info['claimed']}")

    # Propagate star positions to 2560 BCE
    print(f"\nPropagating {len(positions)} stars to {abs(CONSTRUCTION_EPOCH)} BCE...")

    # Use astropy for precession
    epoch_str = f"J{2000 + CONSTRUCTION_EPOCH + 2000}"  # Approximate Julian epoch
    # Better: use Julian date directly
    # 2560 BCE = -2559 CE. J2000 = 2000.0. Delta = -4559 years.
    target_epoch = Time(2000.0 + CONSTRUCTION_EPOCH, format='jyear', scale='tt')

    coords_j2000 = SkyCoord(ra=positions[:, 0]*u.deg, dec=positions[:, 1]*u.deg,
                             frame=FK5(equinox='J2000'))

    # Precess to construction epoch
    target_frame = FK5(equinox=target_epoch)
    coords_epoch = coords_j2000.transform_to(target_frame)

    dec_at_epoch = coords_epoch.dec.deg

    # Find matches for each shaft
    shaft_results = {}
    for shaft, info in targets.items():
        target_dec = info["target_dec"]
        matches = []
        for i in range(len(names)):
            offset = dec_at_epoch[i] - target_dec
            if abs(offset) <= TOLERANCE_DEG:
                matches.append({
                    "name": names[i],
                    "dec_at_epoch": round(float(dec_at_epoch[i]), 3),
                    "offset_deg": round(float(offset), 3),
                    "vmag": float(mags[i]),
                })
        matches.sort(key=lambda x: abs(x["offset_deg"]))

        # Check if claimed star is found
        claimed_found = any(info["claimed"].lower() in m["name"].lower() for m in matches)

        shaft_results[shaft] = {
            "angle": info["angle"],
            "target_dec": round(info["target_dec"], 2),
            "claimed_star": info["claimed"],
            "claimed_found": claimed_found,
            "n_matches": len(matches),
            "best_matches": matches[:5],
        }
        print(f"\n  {shaft}: {len(matches)} stars within ±{TOLERANCE_DEG}°, claimed={info['claimed']} found={claimed_found}")

    # Epoch sweep for best-fit epochs
    print("\nRunning epoch sweep (-3500 to -1500 BCE)...")
    best_epochs = {}
    for shaft, info in targets.items():
        target_dec = info["target_dec"]
        claimed = info["claimed"]
        best_offset = float('inf')
        best_epoch = None

        for year in range(-3500, -1500, 10):
            t = Time(2000.0 + year, format='jyear', scale='tt')
            frame = FK5(equinox=t)
            coords = coords_j2000.transform_to(frame)

            # Find claimed star
            for i, name in enumerate(names):
                if claimed.lower() in name.lower():
                    offset = abs(coords.dec.deg[i] - target_dec)
                    if offset < best_offset:
                        best_offset = offset
                        best_epoch = year
                    break

        if best_epoch is not None:
            best_epochs[shaft] = {"year": best_epoch, "offset_deg": round(best_offset, 4)}
            print(f"  {shaft}: best epoch {abs(best_epoch)} BCE (offset {best_offset:.4f}°)")

    epochs = [v["year"] for v in best_epochs.values()]
    spread = max(epochs) - min(epochs) if epochs else 0

    results = {
        "study": "04_shaft_alignments",
        "construction_epoch": CONSTRUCTION_EPOCH,
        "n_stars": len(positions),
        "tolerance_deg": TOLERANCE_DEG,
        "shafts": shaft_results,
        "epoch_sweep": best_epochs,
        "epoch_spread_years": spread,
    }

    out_path = os.path.join(RESULTS_DIR, "04_shaft_alignments.json")
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved to {out_path}")


if __name__ == "__main__":
    run_shaft_analysis()
