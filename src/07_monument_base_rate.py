"""
Study 7: Cross-site monument-asterism base rate calibration.

Tests 12 famous triangular monument sites worldwide against their culture's
most prominent asterism using identical Procrustes methodology.

Expected result: Giza ranks first (p=0.0025), Thornborough second (p=0.0054).

Usage:
    python src/07_monument_base_rate.py
"""
import json, os, sys, time
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from shared_data import procrustes_distance, load_bsc_catalog, fast_exhaustive_ranking

RESULTS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "results")
os.makedirs(RESULTS_DIR, exist_ok=True)

# 12 monument sites with their proposed stellar counterparts
SITES = [
    {"name": "Giza (Khufu-Khafre-Menkaure)",
     "coords": [(29.9792, 31.1342), (29.9761, 31.1306), (29.9725, 31.1278)],
     "asterism": "Orion Belt", "stars": ["Alnitak", "Alnilam", "Mintaka"]},
    {"name": "Thornborough Henges",
     "coords": [(54.2356, -1.5685), (54.2282, -1.5690), (54.2206, -1.5695)],
     "asterism": "Orion Belt", "stars": ["Alnitak", "Alnilam", "Mintaka"]},
    {"name": "Teotihuacan (Sun-Moon-Feathered Serpent)",
     "coords": [(19.6925, -98.8438), (19.6953, -98.8431), (19.6842, -98.8444)],
     "asterism": "Pleiades", "stars": ["Alcyone", "Atlas", "Electra"]},
    {"name": "Angkor (Angkor Wat-Bayon-Ta Prohm)",
     "coords": [(13.4125, 103.8670), (13.4411, 103.8590), (13.4350, 103.8893)],
     "asterism": "Draco", "stars": ["Eltanin", "Rastaban", "Grumium"]},
    {"name": "Xi'an Pyramids",
     "coords": [(34.3381, 108.9342), (34.3363, 108.9367), (34.3345, 108.9389)],
     "asterism": "Big Dipper", "stars": ["Dubhe", "Merak", "Phecda"]},
    {"name": "Carnac Alignments",
     "coords": [(47.5950, -3.0770), (47.5925, -3.0690), (47.5905, -3.0610)],
     "asterism": "Pleiades", "stars": ["Alcyone", "Atlas", "Electra"]},
    {"name": "Newgrange-Knowth-Dowth",
     "coords": [(53.6947, -6.4755), (53.7003, -6.4908), (53.6992, -6.4500)],
     "asterism": "Orion Belt", "stars": ["Alnitak", "Alnilam", "Mintaka"]},
    {"name": "Avebury-Silbury-West Kennet",
     "coords": [(51.4288, -1.8544), (51.4157, -1.8574), (51.4087, -1.8500)],
     "asterism": "Pleiades", "stars": ["Alcyone", "Atlas", "Electra"]},
    {"name": "Easter Island (3 largest ahu)",
     "coords": [(-27.1248, -109.2772), (-27.1167, -109.3466), (-27.1893, -109.4301)],
     "asterism": "Orion Belt", "stars": ["Alnitak", "Alnilam", "Mintaka"]},
    {"name": "Tikal (Temples I-II-IV)",
     "coords": [(17.2220, -89.6237), (17.2216, -89.6243), (17.2237, -89.6287)],
     "asterism": "Orion Belt", "stars": ["Alnitak", "Alnilam", "Mintaka"]},
    {"name": "Cahokia (Monks-Mound72-Fox)",
     "coords": [(38.6605, -90.0622), (38.6552, -90.0610), (38.6614, -90.0555)],
     "asterism": "Orion Belt", "stars": ["Alnitak", "Alnilam", "Mintaka"]},
    {"name": "Stonehenge-Woodhenge-Durrington",
     "coords": [(51.1789, -1.8262), (51.1891, -1.7859), (51.1936, -1.7853)],
     "asterism": "Pleiades", "stars": ["Alcyone", "Atlas", "Electra"]},
]

# Extra stars not in the mag<4.0 catalog but needed for some asterisms
EXTRA_STARS = {
    "Alcyone":  (56.871, 24.105),
    "Atlas":    (57.291, 24.053),
    "Electra":  (56.219, 24.113),
    "Rastaban": (262.608, 52.301),
    "Grumium":  (268.382, 56.873),
}


def site_to_local(coords):
    """Convert lat/lon coordinates to local meters."""
    lats = np.array([c[0] for c in coords])
    lons = np.array([c[1] for c in coords])
    ref_lat = np.mean(lats)
    ref_lon = np.mean(lons)
    x = (lons - ref_lon) * 111320 * np.cos(np.radians(ref_lat))
    y = (lats - ref_lat) * 110540
    return np.column_stack([x, y])


def run_base_rate():
    print("=" * 60)
    print("Study 7: Cross-Site Monument-Asterism Base Rate")
    print("=" * 60)

    data_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")
    positions, names, mags = load_bsc_catalog(mag_limit=4.0, data_dir=data_dir)
    N = len(positions)
    n_triplets = N * (N - 1) * (N - 2) // 6
    print(f"Catalog: {N} stars, {n_triplets:,} triplets")

    # Build extended catalog for stars not in mag<4.0
    ext_positions = list(positions)
    ext_names = list(names)
    for sname, (ra, dec) in EXTRA_STARS.items():
        if sname not in ext_names:
            ext_positions.append([ra, dec])
            ext_names.append(sname)
    ext_positions = np.array(ext_positions)

    results_list = []
    for site_idx, site in enumerate(SITES):
        print(f"\n[{site_idx+1}/{len(SITES)}] {site['name']} vs {site['asterism']}...")

        # Convert monument to local coordinates
        monument = site_to_local(site["coords"])

        # Find asterism star positions (from extended catalog)
        star_pos = []
        for sname in site["stars"]:
            found = False
            for i, n in enumerate(ext_names):
                if sname.lower() in n.lower():
                    star_pos.append(ext_positions[i])
                    found = True
                    break
            if not found:
                print(f"  WARNING: Star {sname} not found for {site['name']}")

        if len(star_pos) != 3:
            print(f"  SKIPPING: Could not find all 3 stars")
            continue

        star_triangle = np.array(star_pos)

        # Procrustes distance for the claimed match (best of 6 permutations)
        best_D = float('inf')
        for perm in [(0,1,2), (0,2,1), (1,0,2), (1,2,0), (2,0,1), (2,1,0)]:
            p_tri = star_triangle[list(perm)]
            D, _, _ = procrustes_distance(monument, p_tri, allow_reflection=True)
            if D < best_D:
                best_D = D

        # Rank against all triplets from the ORIGINAL catalog (not extended)
        t0 = time.time()
        all_D, _ = fast_exhaustive_ranking(monument, positions, allow_reflection=True)
        rank = int(np.sum(all_D <= best_D))
        p_value = rank / n_triplets
        elapsed = time.time() - t0

        print(f"  D = {best_D:.6f}, rank = {rank:,}/{n_triplets:,}, p = {p_value:.6f} ({elapsed:.1f}s)")

        results_list.append({
            "site": site["name"],
            "asterism": site["asterism"],
            "stars": site["stars"],
            "D": float(best_D),
            "rank": rank,
            "p_value": float(p_value),
        })

    # Sort by p-value
    results_list.sort(key=lambda x: x["p_value"])

    # Benjamini-Hochberg correction
    n_sites = len(results_list)
    for i, r in enumerate(results_list):
        bh_p = r["p_value"] * n_sites / (i + 1)
        r["bh_corrected_p"] = min(bh_p, 1.0)

    # Summary
    n_sig_001 = sum(1 for r in results_list if r["p_value"] < 0.01)
    n_sig_005 = sum(1 for r in results_list if r["p_value"] < 0.05)

    print(f"\n{'='*60}")
    print(f"SUMMARY")
    print(f"  Sites tested: {n_sites}")
    print(f"  p < 0.01: {n_sig_001}/{n_sites}")
    print(f"  p < 0.05: {n_sig_005}/{n_sites}")
    for r in results_list:
        sig = "**" if r["p_value"] < 0.01 else ("*" if r["p_value"] < 0.05 else "")
        print(f"  {r['site']:45s} D={r['D']:.6f} p={r['p_value']:.6f} BH={r['bh_corrected_p']:.4f} {sig}")

    results = {
        "study": "07_monument_base_rate",
        "n_stars": N,
        "n_triplets": n_triplets,
        "n_sites": n_sites,
        "n_significant_001": n_sig_001,
        "n_significant_005": n_sig_005,
        "sites": results_list,
    }

    out_path = os.path.join(RESULTS_DIR, "07_monument_base_rate.json")
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved to {out_path}")


if __name__ == "__main__":
    run_base_rate()
