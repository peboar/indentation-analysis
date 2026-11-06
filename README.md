# Indentation Analysis

Python tools for analyzing spherical indentation of linear-elastic and linear-viscoelastic materials. The repository includes implementations of analytical and numerical approaches for extracting material properties from indentation data, with particular emphasis on the viscoelastic method proposed by Fadil et al. (2017).

The repository also contains an Abaqus implementation for evaluating the methods using finite element simulations.

## Repository Contents

```text
Indentation-Analysis/
│
├── linear_elastic/
│   └── linear_elastic_indentation.py
│
├── viscoelastic/
│   ├── fadil_equation_13.py
│   ├── fadil_equation_13_relaxation_time_convergence.py
│   └── force_indenter.csv
│
├── abaqus/
│   ├── model/
│   │   └── 02_viscoelastic_axisymmetric.cae
│   │
│   ├── input_files/
│   │   ├── indentation_1s.inp
│   │   ├── indentation_0_1s.inp
│   │   └── indentation_0_01s.inp
│   │
│   └── abaqus_fadil_equation_comparison.py
│
├── figures/
│   ├── axisymmetric_fe_mesh.png
│   ├── shear_modulus_loading_time_1s.png
│   ├── shear_modulus_loading_time_0_1s.png
│   ├── shear_modulus_loading_time_0_01s.png
│   └── relaxation_time_convergence.png
│
└── README.md
```

## Theoretical Background

Theoretical assumption: The derivation of Equations (13) and (16) requires a time-independent Poisson's ratio. Consequently, the methods are applicable under the assumption that Poisson's ratio remains constant throughout the viscoelastic response.
Two approaches are considered:

- **Equation (13):** A numerical method for recovering the viscoelastic relaxation response from indentation data.
- **Equation (16):** A closed-form expression for estimating the shear relaxation modulus during the holding phase of an indentation test.

The methods are implemented independently to allow comparison of their performance and numerical behavior.

**Assumption:** Both Equations (13) and (16) are derived under the assumption of a time-independent Poisson's ratio. Consequently, the implementations in this repository assume that Poisson's ratio remains constant throughout the viscoelastic response.

### Applicability of Equation (13)

Equation (13) requires a monotonically increasing indentation history. Consequently, the numerical implementation is applied to indentation histories consisting of a loading phase followed by a holding phase. The method is not intended for unloading or non-monotonic indentation histories.

## Finite Element Model

### Geometry and Element Formulation

The indentation simulations are performed in Abaqus using an axisymmetric finite element model. A rigid spherical indenter is pressed into a viscoelastic specimen under displacement control. Geometric nonlinearity is enabled.

The finite specimen has a height of 40 mm and a radius of 50 mm. An additional infinite-element region is included to approximate an unbounded specimen domain and reduce the influence of the model boundaries on the indentation response. This region also has a height of 40 mm and a radial extent of 50 mm.

The selected model dimensions and element density provide good agreement with the analytical solution for the quantities investigated. This agreement is observed for the reaction force, contact area, contact pressure, and other stress components. However, the peak radial stress shows a larger discrepancy, indicating that this quantity is more sensitive to the finite element discretization and/or the domain approximation.

The dimensions and positioning of the infinite-element region are defined with respect to the pole, which is assumed to coincide with the contact point between the indenter and the specimen. The pole serves as the origin for the far-field mapping used by the infinite elements. Consequently, the finite and infinite-element regions are positioned and dimensioned consistently with this assumed pole location.

The finite specimen is discretized using four-node axisymmetric hybrid elements (CAX4H), while CINAX4 elements are used in the infinite-element region.

The indenter is modeled as an analytical rigid surface with a spherical profile and a radius of 15.875 mm.

The finite element mesh, including the specimen, infinite-element region, and analytical rigid indenter, is shown below.

![Axisymmetric finite element mesh](figures/axisymmetric_fe_mesh.png)

*Figure 1. Axisymmetric finite element mesh showing the specimen, infinite-element region, and analytical rigid indenter.*

**Important:** When generating a new input file from the Abaqus `.cae` model, manually change the element type declaration for the infinite-element region from `CAX4` to `CINAX4`.

The generated input file initially contains:

```text
*Element, type=CAX4
```

to:

```text
*Element, type=CINAX4
```

This modification is required to obtain the intended infinite-element formulation.

### Contact Definition

Frictionless tangential behavior is assumed. The normal contact constraint is enforced using the direct method.

### Material Properties

The specimen is modeled as a linear viscoelastic material with the following elastic properties:

| Property | Value |
|---|---:|
| Young's modulus, $E$ | 476.0142 MPa |
| Poisson's ratio, $\nu$ | 0.48 |
| Initial shear modulus, $G_0$ | 160.8156 MPa |

The initial shear modulus is calculated from the elastic properties:

$$
G_0=\frac{E}{2(1+\nu)}
$$

The time-dependent shear response is described by an 11-term Prony series:

| $i$ | $g_i$ | $\tau_i$ (s) |
|---:|---:|---:|
| 1 | 0.1621 | 0.0010 |
| 2 | 0.2670 | 0.0036 |
| 3 | 0.1898 | 0.0132 |
| 4 | 0.1440 | 0.0479 |
| 5 | 0.0990 | 0.1741 |
| 6 | 0.0640 | 0.6325 |
| 7 | 0.0375 | 2.2974 |
| 8 | 0.0202 | 8.3453 |
| 9 | 0.0097 | 30.3143 |
| 10 | 0.0042 | 110.1169 |
| 11 | 0.0019 | 400.0000 |

The Prony coefficients describe the normalized shear relaxation function:

$$
G(t)=G_0\left[1-\sum_{i=1}^{N}g_i
\left(1-\exp\left(-\frac{t}{\tau_i}\right)\right)\right]
$$

where $G_0$ is the instantaneous shear modulus, $g_i$ are the Prony coefficients, and $\tau_i$ are the corresponding relaxation times.

### Analysis Steps

The indentation simulations consist of two consecutive viscoelastic steps:

1. **Ramp:** The indenter is displaced linearly to a final indentation depth of 0.3 mm.
2. **Relaxation:** The indentation depth is held constant for a further 10 seconds to capture the time-dependent relaxation response.

Three different ramp durations are investigated:

| Case | Ramp duration (s) | Final indentation depth (mm) | Indentation rate (mm/s) |
|---|---:|---:|---:|
| 1 | 1.0 | 0.3 | 0.3 |
| 2 | 0.1 | 0.3 | 3.0 |
| 3 | 0.01 | 0.3 | 30.0 |

The same final indentation depth is used in all cases, allowing the influence of indentation rate on the extracted viscoelastic response to be investigated.

Both steps use Abaqus' viscoelastic procedure with automatic time incrementation. The ramp duration is the prescribed analysis-step duration and should not be confused with the automatic time-increment size.

### Output Requests

The following quantities are extracted from the Abaqus simulations:

- Indenter reaction force (RF2).
- Indenter displacement (U2).
- Contact area (CAREA).
- Contact force (CFN2).

The reaction force and displacement are recorded as history output at the indenter reference point. Contact area and contact force are recorded at regular time intervals of 0.001 s. Field output is requested for selected nodal, element, and contact quantities.

## Numerical Implementation

### Equation (13): Numerical Method

Equation (13) is implemented as a numerical procedure for recovering the shear relaxation modulus from indentation force and displacement data.

The relaxation response is represented using a discretized Prony series. The corresponding coefficients are determined through constrained least-squares optimization:

$$
\min_{\mathbf{g}}
\left\|\mathbf{A}\mathbf{g}-\mathbf{F}_{\mathrm{norm}}\right\|^2
+
\lambda^2\left\|\mathbf{D}\mathbf{g}\right\|^2
$$

where:

- $\mathbf{A}$ is the system matrix constructed from the indentation history.
- $\mathbf{g}$ contains the unknown Prony coefficients.
- $\mathbf{F}_{\mathrm{norm}}$ is the normalized indentation force.
- $\mathbf{D}$ is a second-difference matrix used for optional smoothness regularization.
- $\lambda$ controls the regularization strength.

The optimization is subject to the constraints:

$$
g_i \geq 0
$$

$$
\sum_i g_i \leq 1
$$

The instantaneous shear modulus is assumed to be known and is calculated from the elastic properties.

The regularization parameter is currently set to $\lambda=0$, meaning that no additional smoothness regularization is applied.

### Equation (16): Closed-Form Method

Equation (16) provides a closed-form expression for estimating the shear relaxation modulus during the holding phase of an indentation test.

For a spherical indenter, the contact radius is calculated as:

$$
a=\sqrt{Rh}
$$

where $R$ is the indenter radius and $h$ is the indentation depth.

The shear relaxation modulus is then estimated as:

$$
G(t)=\frac{3}{8}\frac{(1-\nu)RF(t)}{a^3}
$$

where $F(t)$ is the measured indentation force.

Equation (16) is evaluated using the holding phase only.

### Relaxation-Time Discretization

The numerical implementation of Equation (13) requires a discretization of the relaxation-time spectrum.

A convergence study is performed to investigate the influence of the number of relaxation times on the recovered shear relaxation modulus. The relaxation times are distributed logarithmically over the selected time range.

The number of relaxation times per decade is varied to assess the influence of the discretization on the numerical solution. The selected discretization is intended to provide a balance between numerical stability and sufficient resolution of the relaxation response.

## Results

### Influence of Loading Time

The influence of the indentation loading duration on the extracted shear relaxation modulus is investigated using three different ramp durations: 1.0 s, 0.1 s, and 0.01 s.

The extracted responses from Equations (13) and (16) are compared with the reference Prony-series response for each loading duration.

**1.0 s ramp duration**

![Shear relaxation modulus for a 1.0 s ramp duration](figures/shear_modulus_loading_time_1s.png)

**0.1 s ramp duration**

![Shear relaxation modulus for a 0.1 s ramp duration](figures/shear_modulus_loading_time_0_1s.png)

**0.01 s ramp duration**

![Shear relaxation modulus for a 0.01 s ramp duration](figures/shear_modulus_loading_time_0_01s.png)

### Relaxation-Time Discretization Convergence

The influence of the relaxation-time discretization on the numerical solution of Equation (13) is investigated by varying the number of relaxation times per decade.

The extracted shear relaxation modulus is compared with the reference Prony-series response to evaluate the effect of the discretization.

![Relaxation-time discretization convergence](figures/relaxation_time_convergence.png)

### Quantitative Error Comparison

The normalized root mean square error (NRMSE) is used to quantify the difference between the extracted and reference shear relaxation moduli:

$$
\mathrm{NRMSE}=
\frac{
\sqrt{\mathrm{mean}\left[(G_{\mathrm{extracted}}-G_{\mathrm{reference}})^2\right]}
}{
\max\left(|G_{\mathrm{reference}}|\right)
}
\times 100\%
$$

The error is evaluated for both extraction methods and each loading duration.

| Ramp duration | Equation (13) NRMSE (%) | Equation (16) NRMSE (%) |
|---|---:|---:|
| 1.0 s | TBD | TBD |
| 0.1 s | TBD | TBD |
| 0.01 s | TBD | TBD |

## Requirements

The standalone Python implementations require:

- NumPy
- Matplotlib
- SciPy

The Abaqus implementation requires:

- Abaqus/CAE Python environment
- NumPy
- Matplotlib
- SciPy

## Usage

### Standalone Python Scripts

The standalone scripts can be executed using a Python environment with the required dependencies installed.

The indentation force data are read from `force_indenter.csv`. The file should be located in the same directory as the corresponding script, unless its path is modified.

### Abaqus Implementation

The Abaqus implementation requires an existing output database (`.odb`) containing the relevant indentation history output.

The script is executed within the Abaqus/CAE Python environment.

The implementation assumes two analysis steps. The first step can be either a viscoelastic or static step. If the first step is static, the indentation history is extracted from the second step. If both steps are viscoelastic, the numerical implementation of Equation (13) can use the full monotonically increasing indentation history. Equation (16) is evaluated using the holding phase only.

## References

Fadil, H., et al. (2017). *On the Measurement of Two Independent Viscoelastic Functions with Instrumented Indentation Tests.*

## Status

The linear-elastic implementation has been developed as an initial investigation of the numerical behavior of the indentation model.

The current work focuses on the implementation and evaluation of the viscoelastic indentation methods, including their numerical stability, relaxation-time discretization, and performance under different indentation loading durations.

The Abaqus implementation is under continued development and validation.