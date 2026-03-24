"""
Fetch and filter the Yale Bright Star Catalogue (5th edition) for use in the analysis.

Downloads the BSC from VizieR via astroquery, filters to V < 4.0 and visibility
from Giza latitude, and saves a clean CSV to data/bsc_mag4.csv.

Usage:
    python src/fetch_catalog.py

Output:
    data/bsc_mag4.csv  (columns: hr, name, ra_deg, dec_deg, vmag, spectral_type)
"""
import os
import csv

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")
OUTPUT = os.path.join(DATA_DIR, "bsc_mag4.csv")
GIZA_LAT = 29.977  # degrees
DEC_LIMIT = -(90 - GIZA_LAT)  # -60.023 degrees


def fetch_bsc():
    """Download BSC from VizieR via astroquery and save filtered catalog."""
    os.makedirs(DATA_DIR, exist_ok=True)

    print("Downloading Yale Bright Star Catalogue from VizieR via astroquery...")
    from astroquery.vizier import Vizier
    from astropy.coordinates import SkyCoord
    import astropy.units as u

    # Query BSC5 (catalog V/50) for all stars with Vmag < 4.0
    v = Vizier(columns=["HR", "Name", "RAJ2000", "DEJ2000", "Vmag", "SpType"],
               row_limit=-1)
    v.column_filters = {"Vmag": "<4.0"}
    result = v.get_catalogs("V/50")

    if not result or len(result) == 0:
        print("ERROR: No results from VizieR query.")
        return False

    table = result[0]
    print(f"VizieR returned {len(table)} rows")

    # Convert sexagesimal RA/Dec to degrees
    coords = SkyCoord(table["RAJ2000"], table["DEJ2000"],
                      unit=(u.hourangle, u.deg), frame="fk5")
    ra_deg = coords.ra.deg
    dec_deg = coords.dec.deg

    stars = []
    for i, row in enumerate(table):
        try:
            hr = str(row["HR"])
            name = str(row["Name"]).strip() if row["Name"] else f"HR{hr}"
            ra = float(ra_deg[i])
            dec = float(dec_deg[i])
            vmag = float(row["Vmag"])
            sptype = str(row["SpType"]).strip() if row["SpType"] else ""
        except (ValueError, KeyError):
            continue

        # Filter: visible from Giza (dec > -60.023°)
        if dec < DEC_LIMIT:
            continue

        if not name or name == "" or name == "--":
            name = f"HR{hr}"

        stars.append({
            "hr": hr, "name": name, "ra_deg": f"{ra:.6f}",
            "dec_deg": f"{dec:.6f}", "vmag": f"{vmag:.2f}",
            "spectral_type": sptype
        })

    print(f"Filtered to {len(stars)} stars with V < 4.0 visible from Giza latitude (dec > {DEC_LIMIT:.1f}°)")

    with open(OUTPUT, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["hr", "name", "ra_deg", "dec_deg", "vmag", "spectral_type"])
        writer.writeheader()
        writer.writerows(stars)

    print(f"Saved to {OUTPUT}")
    return True


if __name__ == "__main__":
    fetch_bsc()
