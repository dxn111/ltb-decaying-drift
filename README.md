# Bang-time perturbations and redshift drift in parabolic LTB

[![Paper: preprint v1.0](https://img.shields.io/badge/paper-preprint%20v1.0-1455B5?style=flat-square)](paper.pdf)
[![Model: parabolic LTB](https://img.shields.io/badge/model-parabolic%20LTB-242A35?style=flat-square)](paper.tex)
[![Reproduce: Python 3.10+](https://img.shields.io/badge/reproduce-Python%203.10%2B-3776AB?style=flat-square&logo=python&logoColor=white)](analysis.py)
[![Manuscript: CC BY 4.0](https://img.shields.io/badge/manuscript-CC%20BY%204.0-596579?style=flat-square)](LICENSE)
[![Code: MIT](https://img.shields.io/badge/code-MIT-596579?style=flat-square)](LICENSE)

Reproducible research preprint by Dani Postma (DXN1), 17 October 2026.

The paper derives the neighboring-light-cone equation for redshift drift, the parameter sensitivity equation, and the derivative at equal luminosity distance. It evaluates two centrally smooth bang-time profiles on an exact parabolic Einstein–de Sitter LTB background. The sampled responses have no zero crossing. This is an illustrative calculation, not a realistic void fit, a unique signature, or an observational detection.

## Reproduce

```bash
python -m venv .venv
. .venv/bin/activate
pip install -r requirements.txt
python analysis.py
pdflatex -interaction=nonstopmode paper.tex
pdflatex -interaction=nonstopmode paper.tex
```

Outputs: `paper.pdf`, `results/numbers.json`, `results/fixed_distance_response.pdf` and `.png`. Python 3.10+ and a LaTeX installation with the packages listed in `paper.tex` are required. All quantities in the numerical example are dimensionless.

## Scope

The entire distance–redshift curve has not been compensated by adjusting the growing-mode data. The two profiles vary the decaying-mode shape on the same homogeneous background. The old draft's FLRW endpoint shortcut and unsupported zero-crossing claim are withdrawn. Read the paper's limitations before interpreting these results as cosmology.

## Files

- `paper.tex`: manuscript and references.
- `paper.pdf`: compiled preprint.
- `analysis.py`: independent numerical reproduction.
- `results/numbers.json`: full grid values and diagnostics.
- `results/fixed_distance_response.pdf`: vector figure used in the paper.

License: CC BY 4.0 for the manuscript and MIT for code; see `LICENSE`.
