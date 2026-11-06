# Indentation Analysis

Python tools for linear elastic and viscoelastic spherical indentation analysis.

## Overview

This repository contains analytical and numerical tools for spherical indentation analysis, including linear elastic contact mechanics, viscoelastic material characterization, and numerical implementation of methods proposed by Fadil et al.

## Repository Contents

| File | Description |
|---|---|
| `linear_elastic_indentation.py` | Analytical solutions for linear elastic spherical indentation, including contact radius, contact area, normal force, contact pressure, and radial tensile stress. |
| `fadil_equation_13.py` | Numerical implementation of Fadil et al., Equation (13), including force fitting and shear modulus extraction. |
| `abaqus_fadil_equation_comparison.py` | Extracts indentation results from Abaqus and compares shear modulus predictions from Fadil et al., Equations (13) and (16), against the reference Prony series used to define the viscoelastic material model in Abaqus. |
| `force_indenter.csv` | Indentation force data used by `fadil_equation_13.py`. |

## Requirements

The scripts use the following Python libraries:

- NumPy
- Matplotlib
- SciPy

The Abaqus-related analysis requires an appropriate Abaqus Python environment.

## Assumptions

The numerical implementation of Fadil et al., Equation (13), assumes that the instantaneous shear modulus is known and fixed. The Prony coefficients are determined through constrained least-squares optimization, subject to non-negativity and a sum not exceeding unity.

For the current Abaqus validation, the instantaneous shear modulus is calculated from the Young's modulus and Poisson's ratio used to define the material model.

## Usage

The Python scripts can be executed individually for their respective analyses.

`fadil_equation_13.py` requires `force_indenter.csv`, which should be located in the same directory as the script.

### Abaqus analysis

`abaqus_fadil_equation_comparison.py` is intended to be executed within the Abaqus/CAE Python environment and requires access to the relevant Abaqus output database (ODB).

The script assumes that the indentation analysis consists of two load steps. The first step must be either a `*Visco` step or a static step.

- **Static first step:** Only the data from the second step are considered.
- **Two `*Visco` steps:** If the indentation depth increases monotonically across both steps, the full data from both steps are used for Equation (13).

For Equation (16), only the holding phase is considered, regardless of the load-step configuration.

The extracted shear modulus predictions are compared against the reference Prony series used to define the viscoelastic material model in Abaqus.

## Status

This repository is under development. Numerical implementations and results are subject to ongoing validation.