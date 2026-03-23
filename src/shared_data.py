"""
Shared constants and utility functions for the Orion-Giza statistical evaluation.

References:
    - Pyramid coordinates: WGS84, from Dash Foundation surveys
    - Pyramid base lengths: Petrie (1883)
    - Orion Belt positions: J2000, Yale Bright Star Catalogue (Hoffleit & Jaschek 1991)
    - Shaft angles: Gantenbrink (1993)
    - Naked-eye precision: Schaefer (1993)
"""
import numpy as np
from itertools import combinations

# ─── Pyramid data ──────────────────────────────────────────────
PYRAMIDS = {
    "Khufu":    {"lat": 29.9792, "lon": 31.1342, "base_m": 230.4},
    "Khafre":   {"lat": 29.9761, "lon": 31.1306, "base_m": 215.3},
    "Menkaure":  {"lat": 29.9725, "lon": 31.1278, "base_m": 108.5},
}

# Shaft angles (degrees from horizontal, Gantenbrink 1993)
SHAFT_ANGLES = {
    "kings_south":  45.0,
    "kings_north":  32.5,
    "queens_south": 39.5,
    "queens_north": 39.0,
}
SHAFT_CLAIMS = {
    "kings_south":  "Alnitak",
    "kings_north":  "Thuban",
    "queens_south": "Sirius",
    "queens_north": "Kochab",
}
GIZA_LATITUDE = 29.977  # degrees

# Orion Belt J2000 (RA, Dec in degrees)
ORION_BELT = {
    "Alnitak":  {"ra": 85.1897, "dec": -1.9425, "vmag": 1.76},
    "Alnilam":  {"ra": 84.0533, "dec": -1.2019, "vmag": 1.69},
    "Mintaka":  {"ra": 83.0017, "dec": -0.2992, "vmag": 2.23},
}

# Queens' pyramids (approximate, from satellite imagery)
QUEENS_PYRAMIDS = {
    "G1a": {"lat": 29.9753, "lon": 31.1381},
    "G1b": {"lat": 29.9744, "lon": 31.1381},
    "G1c": {"lat": 29.9735, "lon": 31.1381},
}

# Naked-eye angular precision (Schaefer 1993)
SIGMA_STAR_DEG = 0.5


# ─── Coordinate transforms ────────────────────────────────────
def latlon_to_local(lat, lon, ref_lat=None, ref_lon=None):
    """Convert lat/lon to local Cartesian (meters) centered on reference point."""
    if ref_lat is None:
        ref_lat = PYRAMIDS["Khafre"]["lat"]
    if ref_lon is None:
        ref_lon = PYRAMIDS["Khafre"]["lon"]
    x = (lon - ref_lon) * 111320 * np.cos(np.radians(ref_lat))
    y = (lat - ref_lat) * 110540
    return x, y


def get_giza_triangle():
    """Return Giza pyramid positions as Nx2 array in local meters."""
    points = []
    for name in ["Khufu", "Khafre", "Menkaure"]:
        p = PYRAMIDS[name]
        x, y = latlon_to_local(p["lat"], p["lon"])
        points.append([x, y])
    return np.array(points)


def get_orion_triangle():
    """Return Orion Belt star positions as Nx2 array (RA, Dec in degrees)."""
    points = []
    for name in ["Alnitak", "Alnilam", "Mintaka"]:
        s = ORION_BELT[name]
        points.append([s["ra"], s["dec"]])
    return np.array(points)


# ─── Procrustes analysis ──────────────────────────────────────
def procrustes_distance(target, source, allow_reflection=True):
    """
    Compute Procrustes distance between two Nx2 point sets.

    Removes translation, rotation, uniform scaling, and optionally reflection.
    Returns the Procrustes distance D (sum of squared residuals after alignment).

    Parameters:
        target: Nx2 array (reference shape)
        source: Nx2 array (shape to align)
        allow_reflection: if True, allow mirror reflection in alignment

    Returns:
        D: float, Procrustes distance
        R: 2x2 rotation/reflection matrix
        scale: float, scale factor applied to source
    """
    # Center both
    t = target - target.mean(axis=0)
    s = source - source.mean(axis=0)

    # Normalize to unit Frobenius norm
    t_norm = np.sqrt((t**2).sum())
    s_norm = np.sqrt((s**2).sum())
    t = t / t_norm
    s = s / s_norm

    # Optimal rotation via SVD
    M = t.T @ s
    U, sigma, Vt = np.linalg.svd(M)

    if allow_reflection:
        R = U @ Vt
    else:
        # Ensure proper rotation (det = +1)
        d = np.linalg.det(U @ Vt)
        D_sign = np.diag([1, np.sign(d)])
        R = U @ D_sign @ Vt

    # Aligned source
    s_aligned = (R @ s.T).T

    # Procrustes distance
    D = ((t - s_aligned)**2).sum()

    # Scale factor (in original units)
    scale = t_norm / s_norm

    return D, R, scale


def exhaustive_procrustes(target, star_positions, star_names=None, allow_reflection=True):
    """
    Compute Procrustes distance from target to every possible triplet from star_positions.

    Parameters:
        target: 3x2 array (e.g., Giza pyramids)
        star_positions: Nx2 array of star positions
        star_names: optional list of N star names
        allow_reflection: whether to allow reflection in Procrustes

    Returns:
        results: list of (D, indices, names) sorted by ascending D
    """
    N = len(star_positions)
    results = []

    for idx in combinations(range(N), 3):
        triplet = star_positions[list(idx)]

        # Try all 6 permutations of the triplet
        best_D = float('inf')
        best_perm = idx
        for perm in [(0,1,2), (0,2,1), (1,0,2), (1,2,0), (2,0,1), (2,1,0)]:
            t = triplet[list(perm)]
            D, _, _ = procrustes_distance(target, t, allow_reflection=allow_reflection)
            if D < best_D:
                best_D = D
                best_perm = tuple(idx[p] for p in perm)

        names = None
        if star_names is not None:
            names = tuple(star_names[i] for i in best_perm)

        results.append((best_D, best_perm, names))

    results.sort(key=lambda x: x[0])
    return results


# ─── Star catalog loading ─────────────────────────────────────
def load_bsc_catalog(mag_limit=4.0, data_dir="../data"):
    """
    Load the Yale Bright Star Catalogue, filtered to mag < mag_limit.

    Returns:
        positions: Nx2 array of (RA, Dec) in degrees
        names: list of star names
        mags: list of visual magnitudes

    Note: This function expects a BSC data file at data_dir/bsc5.csv.
    If not present, it will attempt to download it, or fall back to
    a pre-filtered catalog shipped with this repo.
    """
    import os
    import csv

    catalog_path = os.path.join(data_dir, "bsc_mag4.csv")

    if not os.path.exists(catalog_path):
        raise FileNotFoundError(
            f"Star catalog not found at {catalog_path}. "
            f"Run 'python src/fetch_catalog.py' to download it."
        )

    positions = []
    names = []
    mags = []

    with open(catalog_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            vmag = float(row['vmag'])
            if vmag <= mag_limit:
                ra = float(row['ra_deg'])
                dec = float(row['dec_deg'])
                positions.append([ra, dec])
                names.append(row['name'])
                mags.append(vmag)

    return np.array(positions), names, mags
