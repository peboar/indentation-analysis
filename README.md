# Indentation Analysiss

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

### Spherical Indentation of Linear-Viscoelastic Materials

The relationship between indentation force and indentation depth for a spherical indenter can be described using linear viscoelastic theory. Following the formulation presented by Fadil et al. (2017), and assuming a time-independent Poisson's ratio, the indentation force is expressed as

$$
P(t)=
\frac{8\sqrt{R}}{3(1-\nu)}
\int_0^t G(t-\tau)
\frac{d h^{3/2}(\tau)}{d\tau}\,d\tau.
\tag{1}
$$

Here:

- \(P(t)\) is the indentation force at time \(t\).
- \(R\) is the radius of the spherical indenter.
- \(\nu\) is the Poisson's ratio, assumed to be independent of time.
- \(G(t)\) is the shear relaxation modulus.
- \(h(t)\) is the indentation depth.
- \(\tau\) is the integration variable, representing an earlier time.
- \(t-\tau\) is the elapsed time over which the material has relaxed.

For convenience, the geometric factor is defined as

$$
C=\frac{8\sqrt{R}}{3(1-\nu)}.
\tag{2}
$$

Equation (1) can consequently be written as

$$
P(t)=C\int_0^t G(t-\tau)
\frac{d h^{3/2}(\tau)}{d\tau}\,d\tau.
\tag{3}
$$

The convolution integral accounts for the time-dependent contribution of the material response to the indentation history. Each increment in indentation contributes to the force according to the relaxation modulus evaluated at the elapsed time since that increment occurred.

The formulation assumes linear viscoelasticity and a time-independent Poisson's ratio. These assumptions are important when interpreting the applicability of the resulting material properties.

### Prony-Series Representation of the Relaxation Modulus

For the numerical implementation, the shear relaxation modulus is approximated using a finite Prony series:

$$
G(t)=G_0\left[
1-\sum_{i=1}^{p}g_i
\left(1-e^{-t/\tau_i}\right)
\right].
\tag{4}
$$

Here:

- \(G_0\) is the instantaneous shear modulus.
- \(p\) is the number of Prony terms.
- \(g_i\) is the normalized Prony coefficient associated with term \(i\).
- \(\tau_i\) is the corresponding relaxation time.
- \(t\) is time.

The coefficients \(g_i\) describe the relative contributions of the individual relaxation terms, while the relaxation times \(\tau_i\) determine their characteristic timescales.

The long-term shear modulus is obtained by evaluating Equation (4) as \(t\) approaches infinity:

$$
G_\infty=G_0\left(1-\sum_{i=1}^{p}g_i\right).
\tag{5}
$$

Consequently, physically admissible relaxation coefficients are subject to the constraints

$$
g_i\geq 0,
\qquad
\sum_{i=1}^{p}g_i\leq 1.
\tag{6}
$$

The Prony-series representation is an approximation introduced for the numerical identification of the relaxation response. The number of terms and the associated relaxation times determine the resolution of the approximated relaxation spectrum.

### Numerical Evaluation of the Convolution Integral

Substituting the Prony-series representation from Equation (4) into the force expression in Equation (3) gives

$$
\begin{aligned}
P(t)=CG_0\int_0^t
\Bigg[
1-\sum_{i=1}^{p}g_i
\left(1-e^{-(t-\tau)/\tau_i}\right)
\Bigg]
\frac{d h^{3/2}(\tau)}{d\tau}\,d\tau.
\end{aligned}
\tag{7}
$$

To evaluate the convolution integral efficiently, an individual integral is defined for each Prony term:

$$
J_i(t)=
\int_0^t e^{-(t-\tau)/\tau_i}
\frac{d h^{3/2}(\tau)}{d\tau}\,d\tau.
\tag{8}
$$

Rather than evaluating this integral independently over the complete indentation history at every time point, a recursive expression can be derived.

Consider a discrete sequence of time points \(t_0,t_1,\ldots,t_n\). At time \(t_n\), the integral can be divided into two intervals:

$$
\begin{aligned}
J_i(t_n)
={}&
\int_0^{t_{n-1}}
e^{-(t_n-\tau)/\tau_i}
\frac{d h^{3/2}(\tau)}{d\tau}\,d\tau\\
&+
\int_{t_{n-1}}^{t_n}
e^{-(t_n-\tau)/\tau_i}
\frac{d h^{3/2}(\tau)}{d\tau}\,d\tau.
\end{aligned}
\tag{9}
$$

Using the identity

$$
e^{-(t_n-\tau)/\tau_i}
=
e^{-(t_n-t_{n-1})/\tau_i}
e^{-(t_{n-1}-\tau)/\tau_i},
$$

the first integral can be expressed using the previously evaluated convolution integral:

$$
\begin{aligned}
J_i(t_n)
={}&e^{-(t_n-t_{n-1})/\tau_i}J_i(t_{n-1})\\
&+
\int_{t_{n-1}}^{t_n}
e^{-(t_n-\tau)/\tau_i}
\frac{d h^{3/2}(\tau)}{d\tau}\,d\tau.
\end{aligned}
\tag{10}
$$

#### Piecewise-Linear Approximation

To evaluate the remaining integral, the indentation history is approximated as piecewise linear over each time interval. Specifically, \(h^{3/2}(t)\) is assumed to vary linearly between consecutive sampling points. Its derivative is therefore approximated by

$$
\frac{d h^{3/2}(\tau)}{d\tau}
\approx
\frac{h_n^{3/2}-h_{n-1}^{3/2}}
{t_n-t_{n-1}},
\tag{11}
$$

where \(h_n=h(t_n)\) denotes the indentation depth at time \(t_n\).

Substituting this approximation into Equation (10) gives

$$
\begin{aligned}
J_i(t_n)
={}&e^{-(t_n-t_{n-1})/\tau_i}J_i(t_{n-1})\\
&+
\frac{h_n^{3/2}-h_{n-1}^{3/2}}
{t_n-t_{n-1}}
\int_{t_{n-1}}^{t_n}
e^{-(t_n-\tau)/\tau_i}\,d\tau.
\end{aligned}
\tag{12}
$$

Evaluating the remaining integral results in the recurrence relation

$$
\boxed{
\begin{aligned}
J_i(t_n)
={}&e^{-\Delta t_n/\tau_i}J_i(t_{n-1})\\
&+
\frac{\tau_i}{\Delta t_n}
\left(h_n^{3/2}-h_{n-1}^{3/2}\right)
\left(1-e^{-\Delta t_n/\tau_i}\right),
\end{aligned}
}
\tag{13}
$$

where

$$
\Delta t_n=t_n-t_{n-1}.
\tag{14}
$$

This expression allows the convolution integral to be evaluated sequentially in time. Only the previously calculated value \(J_i(t_{n-1})\) is required to obtain the next value.

The initial condition is

$$
J_i(t_0)=0,
\tag{15}
$$

assuming that the integration starts at zero time.

The exponential integration in Equation (13) is exact under the piecewise-linear approximation of \(h^{3/2}(t)\). However, this does not imply that the actual indentation depth \(h(t)\) varies linearly in the transformed quantity \(h^{3/2}(t)\).

### Force Expression

Returning to Equation (7), the force can be separated into individual contributions:

$$
\begin{aligned}
P(t)=CG_0\Bigg[
&\int_0^t\frac{d h^{3/2}(\tau)}{d\tau}\,d\tau\\
&-\sum_{i=1}^{p}g_i
\int_0^t
\left(1-e^{-(t-\tau)/\tau_i}\right)
\frac{d h^{3/2}(\tau)}{d\tau}\,d\tau
\Bigg].
\end{aligned}
\tag{16}
$$

Using the fundamental theorem of calculus,

$$
\int_0^t\frac{d h^{3/2}(\tau)}{d\tau}\,d\tau
=
h^{3/2}(t)-h^{3/2}(0).
$$

Assuming zero initial indentation, \(h(0)=0\), and applying the definition of \(J_i(t)\), the force becomes

$$
\begin{aligned}
P(t)
=CG_0\Bigg[
h^{3/2}(t)
-\sum_{i=1}^{p}g_i
\left(h^{3/2}(t)-J_i(t)\right)
\Bigg].
\end{aligned}
\tag{17}
$$

Rearranging the expression gives

$$
\boxed{
P(t)=CG_0\left[
h^{3/2}(t)+
\sum_{i=1}^{p}g_i
\left(J_i(t)-h^{3/2}(t)\right)
\right].
}
\tag{18}
$$

This expression is particularly useful because the indentation force is now written linearly in terms of the unknown Prony coefficients \(g_i\).

### Matrix Formulation

Evaluating Equation (18) at the discrete time points \(t_1,\ldots,t_n\) allows the force expression to be written in matrix form:

$$
\mathbf{P}=CG_0\mathbf{A}
\begin{bmatrix}
1\\
g_1\\
\vdots\\
g_p
\end{bmatrix},
\tag{19}
$$

where \(\mathbf{P}\) is the vector of indentation forces and \(\mathbf{A}\) is the system matrix:

$$
\mathbf{A}=
\begin{bmatrix}
h_1^{3/2} & J_1(t_1)-h_1^{3/2} & \cdots & J_p(t_1)-h_1^{3/2}\\
h_2^{3/2} & J_1(t_2)-h_2^{3/2} & \cdots & J_p(t_2)-h_2^{3/2}\\
\vdots & \vdots & \ddots & \vdots\\
h_n^{3/2} & J_1(t_n)-h_n^{3/2} & \cdots & J_p(t_n)-h_n^{3/2}
\end{bmatrix}.
\tag{20}
$$

The first column contains the squared-transformation indentation history \(h_j^{3/2}\), while the remaining columns contain the corresponding convolution integrals shifted by the indentation history.

The coefficient vector contains a fixed first entry equal to one, followed by the unknown Prony coefficients.

For the optimization, it is convenient to separate the known first column from the unknown coefficients. Define

$$
\mathbf{h}=
\begin{bmatrix}
h_1^{3/2}\\
h_2^{3/2}\\
\vdots\\
h_n^{3/2}
\end{bmatrix},
$$

and construct the reduced system matrix

$$
\mathbf{B}=
\begin{bmatrix}
J_1(t_1)-h_1^{3/2} & \cdots & J_p(t_1)-h_1^{3/2}\\
J_1(t_2)-h_2^{3/2} & \cdots & J_p(t_2)-h_2^{3/2}\\
\vdots & \ddots & \vdots\\
J_1(t_n)-h_n^{3/2} & \cdots & J_p(t_n)-h_n^{3/2}
\end{bmatrix}.
$$

The reduced system is then

$$
\underbrace{
\frac{\mathbf{P}}{CG_0}-\mathbf{h}
}_{\mathbf{b}}
=
\mathbf{B}\mathbf{g},
\tag{21}
$$

where

$$
\mathbf{g}=
\begin{bmatrix}
g_1\\
g_2\\
\vdots\\
g_p
\end{bmatrix}
$$

contains the unknown Prony coefficients.

The coefficients are determined by solving a constrained least-squares optimization problem:

$$
\boxed{
\min_{\mathbf{g}}
\left\|\mathbf{B}\mathbf{g}-\mathbf{b}\right\|^2
+
\lambda^2\left\|\mathbf{D}\mathbf{g}\right\|^2
}
\tag{22}
$$

subject to

$$
g_i\geq 0,
\qquad
\sum_{i=1}^{p}g_i\leq 1.
\tag{23}
$$

Here:

- \(\mathbf{B}\) is the reduced system matrix constructed from the indentation history.
- \(\mathbf{b}\) is the normalized force vector after subtracting the instantaneous elastic contribution.
- \(\mathbf{g}\) contains the unknown Prony coefficients.
- \(\mathbf{D}\) is a second-difference matrix used for optional smoothness regularization.
- \(\lambda\) controls the regularization strength.

The instantaneous shear modulus \(G_0\) and the relaxation times \(\tau_i\) are assumed to be known. Consequently, only the Prony coefficients are determined through the optimization.

The regularization parameter is currently set to \(\lambda=0\), meaning that no additional smoothness regularization is applied.

## Applicability of Equation (13)

Equation (13) of Fadil et al. requires a monotonically increasing indentation history. Consequently, the numerical implementation is applied to indentation histories consisting of a loading phase followed by a holding phase.

The method is not intended for unloading or non-monotonic indentation histories.

The recursive convolution formulation presented above is an implementation of the Prony-series-based numerical approach. The equation numbering in this README refers to the equations presented in this document and should not be confused with equation numbers in the original publication.

### Equation (16): Closed-Form Method

Equation (16) in Fadil et al. provides a closed-form expression for estimating the shear relaxation modulus during the holding phase of an indentation test.

For a spherical indenter, the contact radius is calculated as

$$
a=\sqrt{Rh}.
$$

where \(R\) is the indenter radius and \(h\) is the indentation depth.

The shear relaxation modulus is then estimated as

$$
G(t)=\frac{3}{8}\frac{(1-\nu)RF(t)}{a^3}.
$$

where \(F(t)\) is the measured indentation force.

Equation (16) is evaluated using the holding phase only.

### Relaxation-Time Discretization

The numerical implementation of Equation (13) requires a discretization of the relaxation-time spectrum.

The relaxation times are distributed logarithmically over the selected time range. The number of relaxation times is varied to investigate the influence of the discretization on the recovered shear relaxation modulus.

A convergence study is performed to assess the influence of the number of relaxation times per decade on the numerical solution. The selected discretization is intended to provide a balance between numerical stability and sufficient resolution of the relaxation response.

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

Which should be changed to:

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
| Young's modulus, \(E\) | 476.0142 MPa |
| Poisson's ratio, \(\nu\) | 0.48 |
| Initial shear modulus, \(G_0\) | 160.8156 MPa |

The initial shear modulus is calculated from the elastic properties:

$$
G_0=\frac{E}{2(1+\nu)}
$$

The time-dependent shear response is described by an 11-term Prony series:

| \(i\) | \(g_i\) | \(\tau_i\) (s) |
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
G(t)=G_0\left[
1-\sum_{i=1}^{N}g_i
\left(1-\exp\left(-\frac{t}{\tau_i}\right)\right)
\right].
$$

where \(G_0\) is the instantaneous shear modulus, \(g_i\) are the Prony coefficients, and \(\tau_i\) are the corresponding relaxation times.

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

The relaxation response is represented using a discretized Prony series. The corresponding coefficients are determined through constrained least-squares optimization, as described in the theoretical background.

The optimization uses the reduced system

$$
\mathbf{b}=\mathbf{B}\mathbf{g},
$$

where the matrix is constructed from the indentation history and the prescribed relaxation times.

The optimization is subject to the constraints

$$
g_i\geq 0,
\qquad
\sum_{i=1}^{p}g_i\leq 1.
$$

The instantaneous shear modulus is assumed to be known and is calculated from the elastic properties.

The regularization parameter is currently set to \(\lambda=0\), meaning that no additional smoothness regularization is applied.

### Equation (16): Closed-Form Method

Equation (16) provides a closed-form expression for estimating the shear relaxation modulus during the holding phase of an indentation test.

For a spherical indenter, the contact radius is calculated as

$$
a=\sqrt{Rh}.
$$

The shear relaxation modulus is then estimated as

$$
G(t)=\frac{3}{8}\frac{(1-\nu)RF(t)}{a^3}.
$$

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
\sqrt{\frac{1}{n}\sum_{i=1}^{n}
\left(G_{\mathrm{extracted}}(t_i)-G_{\mathrm{reference}}(t_i)\right)^2}
}{
\max_{1\leq i\leq n}|G_{\mathrm{reference}}(t_i)|
}
\times 100\%.
$$

The error is evaluated for both extraction methods and each loading duration.

| Ramp duration | Equation (13) NRMSE (%) | Equation (16) NRMSE (%) |
|---|---:|---:|
| 1.0 s | 0.189 | 14.526 |
| 0.1 s | 0.191 | 3.779 |
| 0.01 s | 0.193 | 1.046 |

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