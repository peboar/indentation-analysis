# Indentation Analysis

Python tools for linear elastic and viscoelastic spherical indentation analysis.

## Overview

This repository contains numerical and analytical tools for spherical indentation analysis, including linear elastic contact mechanics, viscoelastic material characterization, and numerical fitting of viscoelastic properties.

## Repository Contents

| File | Description |
|---|---|
| `linear_elastic_indentation.py` | Analytical solutions for linear elastic spherical indentation. |
| `numerical_fitted.py` | Numerical fitting of viscoelastic properties using indentation force data. |
| `convergence_terms_per_decade.py` | Investigation of convergence with respect to the number of Prony series terms per decade. |
| `abaqus_indentation_analysis.py` | Extraction and analysis of indentation results from Abaqus output databases (ODB files). |
| `force_indenter.csv` | Indentation force data used by `numerical_fitted.py` and `convergence_terms_per_decade.py`. |

## Requirements

The scripts use the following Python libraries:

- NumPy
- Matplotlib
- SciPy

The Abaqus-related analysis requires an appropriate Abaqus Python environment.

## Usage

The Python scripts can be executed individually for their respective analyses.

The files `numerical_fitted.py` and `convergence_terms_per_decade.py` require `force_indenter.csv` to be available in the expected working directory.

`abaqus_indentation_analysis.py` is intended to be executed within the Abaqus/CAE Python environment.

## Status

This repository is under development. Numerical implementations and results are subject to ongoing validation.

## Author

Developed as part of doctoral research at KTH Royal Institute of Technology.