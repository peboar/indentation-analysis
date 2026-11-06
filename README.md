# Indentation Analysis

Python tools for linear elastic and viscoelastic spherical indentation analysis.

## Overview

This repository contains analytical and numerical tools for spherical indentation analysis, including linear elastic contact mechanics, viscoelastic material characterization, and numerical implementation of methods proposed by [Fadil et al. (2017)](https://www.researchgate.net/publication/320319893_On_the_Measurement_of_two_Independent_Viscoelastic_Functions_with_Instrumented_Indentation_Tests).

## Validation and Results

### Linear Elastic Indentation

The analytical solution for linear elastic spherical indentation is compared against finite element results obtained using Abaqus. The comparison includes contact radius, contact area, normal force, contact pressure, and radial tensile stress.

![Hertz validation](figures/hertz_validation.png)

Mesh convergence is investigated to evaluate the influence of the finite element discretisation on the predicted indentation response and stress distribution.

![Mesh convergence](figures/mesh_convergence.png)

### Viscoelastic Indentation

The shear modulus extraction methods proposed by Fadil et al. (2017), Equations (13) and (16), are evaluated against the reference Prony series used to define the viscoelastic material model in Abaqus.

The comparison is presented using both linear and logarithmic time axes.

![Shear modulus comparison - linear time](figures/shear_modulus_linear.png)

![Shear modulus comparison - logarithmic time](figures/shear_modulus_log.png)

### Relaxation-Time Convergence

The influence of the relaxation-time discretisation on the shear modulus predictions obtained using Equation (13) from Fadil et al. (2017) is investigated to establish a suitable number of Prony terms.

![Relaxation-time convergence](figures/relaxation_time_convergence.png)

## Repository Contents

| File | Description |
|---|---|
| `linear_elastic/linear_elastic_indentation.py` | Analytical solutions for linear elastic spherical indentation, including contact radius, contact area, normal force, contact pressure, and radial tensile stress. |
| `viscoelastic/fadil_equation_13.py` | Numerical implementation of Equation (13) from Fadil et al. (2017), including force fitting and shear modulus extraction. |
| `viscoelastic/fadil_equation_13_relaxation_time_convergence.py` | Investigates the influence of the relaxation-time discretisation on the accuracy and stability of shear modulus predictions obtained using Equation (13). |
| `abaqus/abaqus_fadil_equation_comparison.py` | Extracts indentation results from Abaqus and compares shear modulus predictions obtained using Equations (13) and (16) from Fadil et al. (2017) against the reference Prony series used to define the viscoelastic material model in Abaqus. |
| `viscoelastic/force_indenter.csv` | Indentation force data used by the standalone viscoelastic analysis scripts. |

## Theoretical Background

The viscoelastic indentation methods implemented in this repository are based on the work of Fadil et al. (2017), particularly Equations (13) and (16), presented in:

[Fadil et al. (2017) – On the Measurement of two Independent Viscoelastic Functions with Instrumented Indentation Tests](https://www.researchgate.net/publication/320319893_On_the_Measurement_of_two_Independent_Viscoelastic_Functions_with_Instrumented_Indentation_Tests).

The numerical implementations in this repository follow these formulations, with the implementation details and assumptions described below.

## Requirements

The scripts use the following Python libraries:

- NumPy
- Matplotlib
- SciPy

The Abaqus-related analysis requires an appropriate Abaqus Python environment.

## Numerical Implementation and Assumptions

The numerical implementation of Fadil et al., Equation (13), is formulated as a constrained least-squares problem to determine the Prony coefficients.

### Instantaneous Shear Modulus

The instantaneous shear modulus is assumed to be known and fixed. For the current Abaqus validation, it is calculated from the Young's modulus and Poisson's ratio used to define the material model:

$$
G_0 = \frac{E}{2(1+\nu)}
$$

Only the Prony coefficients are determined through optimization.

### Optimization Procedure

The force response is expressed in matrix form:

$$
\mathbf{A}\mathbf{g} \approx \mathbf{F}_{\mathrm{norm}}
$$

where $\mathbf{A}$ is the force matrix, $\mathbf{g}$ contains the Prony coefficients, and $\mathbf{F}_{\mathrm{norm}}$ contains the normalized force data.

The coefficients are obtained by minimizing the mean squared error (MSE), with an optional smoothness regularization term:

$$
J(\mathbf{g}) =
\frac{1}{N}
\left[
\left\|\mathbf{A}\mathbf{g}-\mathbf{F}_{\mathrm{norm}}\right\|^2
+
\lambda^2\left\|\mathbf{D}\mathbf{g}\right\|^2
\right]
$$

where $N$ is the number of data points, $\lambda$ is the regularization parameter, and $\mathbf{D}$ is a second-order difference matrix used to penalize variations between consecutive Prony coefficients. The regularization matrix acts on the Prony coefficients only; the normalized instantaneous contribution is fixed.

The optimization is performed using the Sequential Least Squares Programming (SLSQP) algorithm in SciPy, subject to the following constraints:

- Prony coefficients must be non-negative.
- The sum of the Prony coefficients must not exceed unity.
- The normalized instantaneous contribution is fixed to unity.

### Regularization

The regularization parameter is currently set to $\lambda=0$, meaning that the optimization minimizes the ordinary mean squared force error without applying a smoothness penalty.

This assumption is considered appropriate for the current validation, where the force data are extracted from Abaqus simulations and measurement noise is negligible. However, numerical errors associated with the finite element analysis may still be present.

For future experimental applications, the regularization parameter will need to be reassessed to account for measurement noise and improve the stability and smoothness of the fitted Prony coefficients.

### Relaxation-Time Discretisation

The relaxation-time convergence study investigates the influence of the number of logarithmic relaxation-time decades and the number of Prony terms per decade on the shear modulus predictions obtained using Equation (13) from Fadil et al. (2017).

The study indicates that two relaxation-time decades provide a suitable balance between accuracy and stability. A discontinuity in the predicted shear modulus $G(t)$ is observed when increasing the number of coefficients from four to five within a single decade.

To provide a safety margin and avoid relying on this potentially unstable discretisation, two relaxation-time decades with two Prony terms per decade are adopted in both the standalone and Abaqus implementations. This results in four Prony terms across the selected relaxation-time range.

## Usage

The Python scripts can be executed individually for their respective analyses.

### Standalone Analysis

The scripts `fadil_equation_13.py` and `fadil_equation_13_relaxation_time_convergence.py` require `force_indenter.csv`, which should be located in the same directory as the scripts.

`fadil_equation_13.py` performs the numerical implementation of Fadil et al., Equation (13), including force fitting and shear modulus extraction.

`fadil_equation_13_relaxation_time_convergence.py` investigates the influence of the relaxation-time discretisation on the predicted shear modulus. It is used to determine the number of relaxation-time decades adopted in the standalone and Abaqus implementations.

### Abaqus Analysis

`abaqus_fadil_equation_comparison.py` is intended to be executed within the Abaqus/CAE Python environment and requires access to the relevant Abaqus output database (ODB).

The script assumes that the indentation analysis consists of two load steps. The first step must be either a `*Visco` step or a static step.

- **Static first step:** Only the data from the second step are considered.
- **Two `*Visco` steps:** If the indentation depth increases monotonically across both steps, the full data from both steps are used for Equation (13).

For Equation (16), only the holding phase is considered, regardless of the load-step configuration.

The extracted shear modulus predictions are compared against the reference Prony series used to define the viscoelastic material model in Abaqus.

## References

Fadil, H., et al. (2017). [On the Measurement of two Independent Viscoelastic Functions with Instrumented Indentation Tests](https://www.researchgate.net/publication/320319893_On_the_Measurement_of_two_Independent_Viscoelastic_Functions_with_Instrumented_Indentation_Tests).

## Status

This repository is under development. The linear elastic implementation has been validated against analytical Hertz solutions. The viscoelastic implementations and their numerical stability are subject to ongoing validation.