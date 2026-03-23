"""
Fetch and filter the Yale Bright Star Catalogue (5th edition) for use in the analysis.

Downloads the BSC from VizieR, filters to V < 4.0 and visibility from Giza latitude,
and saves a clean CSV to data/bsc_mag4.csv.

Usage:
    python src/fetch_catalog.py

Output:
    data/bsc_mag4.csv  (columns: hr, name, ra_deg, dec_deg, vmag, spectral_type)
"""
import os
import csv
import urllib.request
import re

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")
OUTPUT = os.path.join(DATA_DIR, "bsc_mag4.csv")

# VizieR TAP query for BSC5 (catalog V/50)
VIZIER_URL = (
    "https://vizier.cds.unistra.fr/viz-bin/ase/tsv?"
    "-source=V/50&-out=HR,Name,RAJ2000,DEJ2000,Vmag,SpType"
    "&Vmag=<4.0&-out.max=unlimited"
)

GIZA_LAT = 29.977  # degrees


def fetch_bsc():
    """Download BSC from VizieR and save filtered catalog."""
    os.makedirs(DATA_DIR, exist_ok=True)

    print("Downloading Yale Bright Star Catalogue from VizieR...")
    try:
        req = urllib.request.Request(VIZIER_URL, headers={"User-Agent": "Python/OrionGiza"})
        with urllib.request.urlopen(req, timeout=30) as resp:
            raw = resp.read().decode("utf-8")
    except Exception as e:
        print(f"Download failed: {e}")
        print("Please download the BSC manually from VizieR (catalog V/50)")
        print(f"and save a CSV with columns: hr,name,ra_deg,dec_deg,vmag,spectral_type")
        print(f"to: {OUTPUT}")
        return False

    # Parse TSV (skip header lines starting with #)
    lines = [l for l in raw.strip().split("\n") if not l.startswith("#") and l.strip()]
    if len(lines) < 2:
        print("Unexpected response format. Attempting alternative parse...")
        return False

    stars = []
    for line in lines[2:]:  # Skip header + separator
        parts = line.split("\t")
        if len(parts) < 5:
            continue
        try:
            hr = parts[0].strip()
            name = parts[1].strip() if len(parts) > 1 else f"HR{hr}"
            ra = float(parts[2].strip())
            dec = float(parts[3].strip())
            vmag = float(parts[4].strip())
            sptype = parts[5].strip() if len(parts) > 5 else ""
        except (ValueError, IndexError):
            continue

        # Filter: visible from Giza (dec > -(90 - lat) = -60.023)
        if dec < -(90 - GIZA_LAT):
            continue
        # Filter: never rises above horizon (circumpolar below)
        # A star is visible if dec > -(90 - lat)
        # Already filtered above

        if not name or name == "":
            name = f"HR{hr}"

        stars.append({
            "hr": hr, "name": name, "ra_deg": f"{ra:.6f}",
            "dec_deg": f"{dec:.6f}", "vmag": f"{vmag:.2f}",
            "spectral_type": sptype
        })

    print(f"Downloaded {len(stars)} stars with V < 4.0 visible from Giza latitude")

    with open(OUTPUT, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["hr", "name", "ra_deg", "dec_deg", "vmag", "spectral_type"])
        writer.writeheader()
        writer.writerows(stars)

    print(f"Saved to {OUTPUT}")
    return True


if __name__ == "__main__":
    fetch_bsc()
