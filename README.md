# Orion–Giza Statistical Evaluation

**Institute:** [Deep Time Research Institute](https://deeptime-research.org)
**Author:** Elliot Allan — elliot@deeptime-research.org — ORCID [0009-0008-8541-0944](https://orcid.org/0009-0008-8541-0944)

> ⚠️ This repository is under active revision. The analysis described here corresponds to the manuscript currently under review. An earlier version targeted Journal of Astronomical History and Heritage; the current submission is at Journal of Astronomy in Culture.

Reproducibility package for:

> Allan, E. (2026). *A comprehensive statistical evaluation of the Orion–Giza correlation.* Deep Time Research Institute. Submitted to *Journal of Astronomy in Culture*.

## What this is

The first exhaustive statistical evaluation of Bauval's (1989) Orion Correlation Theory, testing whether the spatial arrangement of the three Great Pyramids of Giza encodes the shape of Orion's Belt. We rank the Giza layout against all 695,520 bright-star triplets (V < 4.0) from the Yale Bright Star Catalogue using Procrustes shape analysis.

### Key findings

- **Shape match is significant**: Orion's Belt ranks 1,717th of 695,520 triplets (p = 0.0025, 8.7 bits of information)
- **No reflection needed**: Pure rotation (92.7°), neutralizing Fairall's (1999) mirror objection
- **Every extension fails**: Orientation (p = 0.74), shaft convergence (1,325-year spread), Nile Valley extension (p = 0.50), scale (p = 0.219)
- **Giza is the best among 12 sites**: Cross-site comparison with BH-corrected p = 0.030
- **Precision paradox**: The match sits at the 0.8th percentile of the builder's precision envelope

## Repository structure

```
├── src/                    # Analysis scripts (Python)
│   ├── shared_data.py      # Constants, coordinates, Procrustes implementation
│   ├── fetch_catalog.py    # Download BSC from VizieR
│   ├── 01_exhaustive_ranking.py    # Core triplet ranking
│   ├── 02_information_content.py   # Bits decomposition
│   ├── 03_precision_envelope.py    # Builder's envelope MC
│   ├── 04_shaft_alignments.py      # Shaft analysis with precession
│   ├── 05_reflection_chirality.py  # Reflection/chirality test
│   ├── 06_scale_significance.py    # Scale MC test
│   └── 07_monument_base_rate.py    # Cross-site comparison
├── data/                   # Input data (BSC catalog, fetched on first run)
├── results/                # Output JSON from each study
├── figures/                # Publication-quality figures
├── paper/                  # Manuscript (JAHH format, .docx)
├── requirements.txt
└── README.md
```

## Quick start

```bash
# Clone and install dependencies
git clone https://github.com/thegreatcircledata/orion-giza-statistical-evaluation.git
cd orion-giza-statistical-evaluation
pip install -r requirements.txt

# Fetch the star catalog
python src/fetch_catalog.py

# Run the core analysis (warning: ~30-60 min for exhaustive enumeration)
python src/01_exhaustive_ranking.py

# Run individual studies
python src/02_information_content.py
python src/03_precision_envelope.py
python src/04_shaft_alignments.py      # Requires astropy
python src/05_reflection_chirality.py   # Warning: ~60 min
python src/06_scale_significance.py
python src/07_monument_base_rate.py
```

## Data sources

- **Yale Bright Star Catalogue** (5th ed.): Hoffleit & Jaschek (1991), accessed via VizieR
- **Gaia DR3**: Gaia Collaboration (2023), used for independent cross-check
- **Pyramid coordinates**: WGS84, from Dash Foundation surveys
- **Pyramid dimensions**: Petrie (1883)
- **Shaft angles**: Gantenbrink (1993)

## Methodology

All analyses use **ordinary Procrustes analysis** (OPA), which removes translation, rotation, uniform scaling, and optionally reflection, leaving 2 effective degrees of freedom for a 3-point configuration. The empirical p-value is defined as rank/N where N = C(162,3) = 695,520 total triplets.

See the paper for full methodological details.

## AI disclosure

Claude (Anthropic) was used as a computational research assistant for code generation, Monte Carlo simulation design, and manuscript drafting. All analyses were conceived, directed, and verified by the author. See the paper's acknowledgments for full disclosure.

## License

MIT. The corresponding author is affiliated with the Deep Time Research Institute, a registered 501(c)(3) nonprofit research organisation.

## Citation

```bibtex
@article{allan2026orion,
  author       = {Allan, Elliot},
  title        = {A comprehensive statistical evaluation of the {Orion}--{Giza} correlation},
  journal      = {Journal of Astronomy in Culture},
  year         = {2026},
  note         = {Submitted},
  organization = {Deep Time Research Institute},
  orcid        = {0009-0008-8541-0944}
}
```
