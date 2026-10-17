# Bang-time perturbations and redshift drift in LTB cosmology

[![Paper: preprint v2.0](https://img.shields.io/badge/paper-preprint%20v2.0-1455B5?style=flat-square)](paper.pdf)
[![Models: LTB and Λ-LTB](https://img.shields.io/badge/models-LTB%20%2B%20%CE%9B--LTB-242A35?style=flat-square)](paper.tex)
[![Reproduce: Python 3.10+](https://img.shields.io/badge/reproduce-Python%203.10%2B-3776AB?style=flat-square&logo=python&logoColor=white)](analysis.py)
[![Manuscript: CC BY 4.0](https://img.shields.io/badge/manuscript-CC%20BY%204.0-596579?style=flat-square)](LICENSE)
[![Code: MIT](https://img.shields.io/badge/code-MIT-596579?style=flat-square)](LICENSE)

Reproducible research preprint by Dani Postma (DXN1), 17 October 2026.

The paper derives the neighboring-light-cone drift equation and answers the initial model questions with computed results. The original hump-shaped bang-time profiles develop shell crossings, while monotone profiles avoid them in the expanding parabolic dust and flat matter-plus-Λ families. Their sampled equal-distance drift response changes sign relative to the hump family. A local curvature reconstruction matches the complete distance curve at first order but grows sharply near the EdS angular-distance maximum, where its inverse equation is singular. These are model-family results, not an observational detection or a universal theorem.

## Reproduce

```bash
python -m venv .venv
. .venv/bin/activate
pip install -r requirements.txt
python analysis.py
python viable.py
python lambda_model.py
python constrained.py
pdflatex -interaction=nonstopmode paper.tex
pdflatex -interaction=nonstopmode paper.tex
```

Outputs: `paper.pdf`, machine-readable JSON values, five vector figures, and PNG previews. Python 3.10+ and a LaTeX installation with the packages listed in `paper.tex` are required. Coordinates and responses are dimensionless; the flat Λ model uses `H0=1` and `(Ωm, ΩΛ)=(0.3, 0.7)`.

## Scope

The full distance curve is compensated **locally at first order** in `constrained.py`; that solution becomes uncontrolled near the EdS distance maximum. The flat Λ background is calibrated by dimensionless density parameters but no inhomogeneous model has been fitted to cosmological data. Read the paper's assumptions and observational limits before interpreting a response as a physical signal. The old draft's FLRW endpoint shortcut and unsupported zero-crossing claim are withdrawn.

## Files

- `paper.tex`: manuscript and references.
- `paper.pdf`: compiled preprint.
- `analysis.py`: independent numerical reproduction.
- `viable.py`: monotone bang-time profiles without shell crossings in parabolic dust.
- `lambda_model.py`: exact flat Λ-LTB family and homogeneous-limit checks.
- `constrained.py`: local first-order distance-curve compensation.
- `results/numbers.json`: full grid values and diagnostics.
- `results/fixed_distance_response.pdf`: vector figure used in the paper.
- `results/comparison_and_endpoint_error.pdf`: fixed-observable and endpoint-error comparisons.
- `results/monotone_response.pdf`, `results/lambda_response.pdf`, `results/constrained_response.pdf`: regular-family and reconstruction figures.

License: CC BY 4.0 for the manuscript and MIT for code; see `LICENSE`.
