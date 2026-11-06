from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize


def read_indentation_data(csv_path):
    """Return the indentation times, depths, and normal forces."""
    raw_data = np.genfromtxt(
        csv_path,
        delimiter=",",
        skip_header=1,
        usecols=(0, 1, 2, 3),
    )

    times = raw_data[:, 0] - raw_data[0, 0]
    depths = np.abs(raw_data[:, 2])
    forces = np.abs(raw_data[:, 3])

    is_increasing = np.r_[True, np.diff(times) > 0]

    return (
        times[is_increasing],
        depths[is_increasing],
        forces[is_increasing],
    )


def compute_instantaneous_shear_modulus(youngs_modulus, poisson_ratio):
    """Compute the instantaneous shear modulus."""
    return youngs_modulus / (2.0 * (1.0 + poisson_ratio))


def compute_contact_coefficient(indenter_radius, poisson_ratio):
    """Compute the contact coefficient."""
    return (
        8.0 * np.sqrt(indenter_radius)
        / (3.0 * (1.0 - poisson_ratio))
    )


def get_relaxation_time_range(times, min_relaxation_time=1e-4):
    """Return the shortest and longest relaxation times."""
    shortest = min(np.min(np.diff(times)), min_relaxation_time)
    longest = times[-1] - times[0]

    if longest <= shortest:
        raise ValueError(
            "Test duration is not longer than the minimum relaxation time."
        )

    return shortest, longest


def get_number_of_coefficients(
    times,
    terms_per_decade=5,
    min_relaxation_time=1e-4,
):
    """Compute the number of Prony terms for the selected time range."""
    shortest, longest = get_relaxation_time_range(
        times,
        min_relaxation_time,
    )

    return (
        int(np.ceil(terms_per_decade * np.log10(longest / shortest))) + 1
    )


def get_prony_relaxation_times(
    times,
    number_of_relaxation_times,
    min_relaxation_time=1e-4,
):
    """Generate logarithmically spaced Prony relaxation times."""
    shortest, longest = get_relaxation_time_range(
        times,
        min_relaxation_time,
    )

    return np.logspace(
        np.log10(shortest),
        np.log10(longest),
        number_of_relaxation_times,
    )


def compute_convolution_integral(
    times,
    depth_pow_3_2,
    relaxation_time,
):
    """Compute the exponential convolution integral.

    Exact for piecewise-linear depth**1.5, using a recursive update.
    """
    delta_time = np.diff(times)
    depth_rates = np.diff(depth_pow_3_2) / delta_time

    exponential_decay = np.exp(-delta_time / relaxation_time)

    interval_contribution = (
        relaxation_time
        * depth_rates
        * -np.expm1(-delta_time / relaxation_time)
    )

    convolution_integral = np.zeros_like(times, dtype=float)

    for i in range(1, len(times)):
        convolution_integral[i] = (
            exponential_decay[i - 1] * convolution_integral[i - 1]
            + interval_contribution[i - 1]
        )

    return convolution_integral


def compute_force_matrix(
    times,
    depth_pow_3_2,
    relaxation_times,
    contact_coefficient,
    instantaneous_shear_modulus,
    force_scale,
):
    """Compute the normalised force matrix.

    force_matrix @ [1, g_1, ..., g_n] = F(t) / force_scale

    Column 0 multiplies the fixed instantaneous response coefficient.
    The remaining columns multiply the Prony coefficients.
    """
    convolution_integrals = [
        compute_convolution_integral(
            times,
            depth_pow_3_2,
            relaxation_time,
        )
        for relaxation_time in relaxation_times
    ]

    return (
        contact_coefficient
        * instantaneous_shear_modulus
        / force_scale
        * np.column_stack(
            [
                depth_pow_3_2,
                *[
                    convolution_integral - depth_pow_3_2
                    for convolution_integral in convolution_integrals
                ],
            ]
        )
    )


def compute_regularization_matrix(number_of_terms):
    """Compute the second-order regularisation matrix."""
    regularization_matrix = np.zeros(
        (number_of_terms - 2, number_of_terms + 1)
    )

    for i in range(number_of_terms - 2):
        regularization_matrix[i, i + 1:i + 4] = [1, -2, 1]

    return regularization_matrix


def objective_function(
    parameters,
    force_matrix,
    normalized_forces,
    regularization_matrix,
    regularization_parameter,
):
    """Compute the regularised least-squares objective."""
    residual = force_matrix @ parameters - normalized_forces

    regularization = (
        regularization_parameter
        * (regularization_matrix @ parameters)
    )

    objective = (
        np.sum(residual**2)
        + np.sum(regularization**2)
    )

    return objective / len(force_matrix)


def objective_gradient(
    parameters,
    force_matrix,
    normalized_forces,
    regularization_matrix,
    regularization_parameter,
):
    """Compute the gradient of the regularised least-squares objective."""
    residual = force_matrix @ parameters - normalized_forces

    gradient = (
        2.0 * force_matrix.T @ residual
        + 2.0 * regularization_parameter**2
        * regularization_matrix.T
        @ (regularization_matrix @ parameters)
    )

    return gradient / len(force_matrix)


def compute_coefficients(
    force_matrix,
    normalized_forces,
    regularization_matrix,
    regularization_parameter,
):
    """Solve the constrained regularised least-squares problem."""
    number_of_parameters = force_matrix.shape[1]

    initial_guess = np.zeros(number_of_parameters)
    initial_guess[0] = 1.0
    initial_guess[1:] = 0.5 / (number_of_parameters - 1)

    jacobian_constraint_instantaneous_modulus = np.zeros(
        number_of_parameters
    )
    jacobian_constraint_instantaneous_modulus[0] = 1.0

    jacobian_constraint_prony_coefficients = -np.ones(
        number_of_parameters
    )
    jacobian_constraint_prony_coefficients[0] = 0.0

    result = minimize(
        objective_function,
        initial_guess,
        args=(
            force_matrix,
            normalized_forces,
            regularization_matrix,
            regularization_parameter,
        ),
        jac=objective_gradient,
        method="SLSQP",
        bounds=[(0.0, None)] * number_of_parameters,
        constraints=[
            {
                "type": "eq",
                "fun": lambda parameters: parameters[0] - 1.0,
                "jac": lambda parameters: (
                    jacobian_constraint_instantaneous_modulus
                ),
            },
            {
                "type": "ineq",
                "fun": lambda parameters: (
                    1.0 - np.sum(parameters[1:])
                ),
                "jac": lambda parameters: (
                    jacobian_constraint_prony_coefficients
                ),
            },
        ],
        options={"maxiter": 1000, "ftol": 1e-15},
    )

    if not result.success:
        raise RuntimeError(result.message)

    misfit = np.sqrt(
        np.mean(
            (force_matrix @ result.x - normalized_forces) ** 2
        )
    )

    if misfit > 1e-2:
        raise RuntimeError(
            f"Solver stopped with a poor fit (RMS misfit {misfit:.1e})."
        )

    return result.x


def plot_results(
    times,
    fitted_shear_modulus,
    reference_shear_modulus,
    fitted_forces,
    reference_forces,
    xlog=True,
):
    """Plot shear modulus and force comparisons."""
    if xlog:
        times_modulus = times[1:]
        fitted_modulus = fitted_shear_modulus[1:]
        reference_modulus = reference_shear_modulus[1:]
    else:
        times_modulus = times
        fitted_modulus = fitted_shear_modulus
        reference_modulus = reference_shear_modulus

    spacing = np.geomspace if xlog else np.linspace

    marker_targets = spacing(
        times_modulus[0],
        times_modulus[-1],
        5,
    )

    if xlog:
        marker_distances = np.abs(
            np.log(times_modulus[:, None])
            - np.log(marker_targets)
        )
    else:
        marker_distances = np.abs(
            times_modulus[:, None] - marker_targets
        )

    marker_indices = np.unique(
        marker_distances.argmin(axis=0)
    )

    with plt.rc_context({
        "font.family": "serif",
        "mathtext.fontset": "stix",
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "xtick.minor.visible": True,
        "ytick.minor.visible": True,
        "axes.linewidth": 1.0,
        "font.size": 10,
    }):
        fig, axes = plt.subplots(1, 2, figsize=(10, 4))

        golden_ratio = (1 + np.sqrt(5)) / 2

        for ax in axes:
            ax.set_box_aspect(1 / golden_ratio)

        axes[0].plot(
            times_modulus,
            reference_modulus,
            color="forestgreen",
            linewidth=1.0,
            zorder=2,
            label="Reference Prony series",
        )

        axes[0].plot(
            times_modulus,
            fitted_modulus,
            color="red",
            linestyle="--",
            linewidth=1.0,
            marker="o",
            markersize=4.2,
            mfc="none",
            mew=1.0,
            markevery=list(marker_indices),
            zorder=4,
            label="Fadil et al., Eq. (13)",
        )

        if xlog:
            axes[0].set_xscale("log")

        axes[0].set_xlabel("Time (s)")
        axes[0].set_ylabel("Shear modulus (MPa)")
        axes[0].legend(
            frameon=True,
            fancybox=False,
            edgecolor="black",
        )

        axes[1].plot(
            times,
            reference_forces,
            color="green",
            linewidth=1.0,
            zorder=2,
            label="Reference Prony series",
        )

        axes[1].plot(
            times,
            fitted_forces,
            color="red",
            linestyle="--",
            linewidth=1.0,
            zorder=4,
            label="Fadil et al., Eq. (13)",
        )

        axes[1].set_xlabel("Time (s)")
        axes[1].set_ylabel("Force (N)")
        axes[1].legend(
            frameon=True,
            fancybox=False,
            edgecolor="black",
        )

        for ax in axes:
            ax.tick_params(top=True, right=True, width=1.0)
            ax.set_axisbelow(True)

        fig.tight_layout()
        plt.show()


if __name__ == "__main__":
    # Parameters
    csv_path = Path(__file__).resolve().parent / "force_indenter.csv"

    indenter_radius = 15.875
    youngs_modulus = 476.0142
    poisson_ratio = 0.48

    terms_per_decade = 2
    min_relaxation_time = 1e-4
    regularization_parameter = 0.0

    # Known Abaqus Prony series
    prony_coefficients = np.array([
        0.1621,
        0.2670,
        0.1898,
        0.1440,
        0.0990,
        0.0640,
        0.0375,
        0.0202,
        0.0097,
        0.0042,
        0.0019,
    ])

    relaxation_times = np.array([
        0.0010,
        0.0036,
        0.0132,
        0.0479,
        0.1741,
        0.6325,
        2.2974,
        8.3453,
        30.3143,
        110.1169,
        400.0,
    ])

    times, depths, forces = read_indentation_data(csv_path)

    depth_pow_3_2 = depths**(3.0 / 2.0)
    force_scale = np.max(forces)
    normalized_forces = forces / force_scale

    instantaneous_shear_modulus = compute_instantaneous_shear_modulus(
        youngs_modulus,
        poisson_ratio,
    )

    contact_coefficient = compute_contact_coefficient(
        indenter_radius,
        poisson_ratio,
    )

    number_of_coefficients = get_number_of_coefficients(
        times,
        terms_per_decade,
        min_relaxation_time,
    )

    fitted_relaxation_times = get_prony_relaxation_times(
        times,
        number_of_coefficients,
        min_relaxation_time,
    )

    force_matrix = compute_force_matrix(
        times,
        depth_pow_3_2,
        fitted_relaxation_times,
        contact_coefficient,
        instantaneous_shear_modulus,
        force_scale,
    )

    regularization_matrix = compute_regularization_matrix(
        number_of_coefficients
    )

    fitted_coefficients = compute_coefficients(
        force_matrix,
        normalized_forces,
        regularization_matrix,
        regularization_parameter,
    )

    fitted_prony_coefficients = fitted_coefficients[1:]

    numerical_shear_modulus = instantaneous_shear_modulus * (
        1.0 - np.sum(fitted_prony_coefficients)
        + np.sum(
            fitted_prony_coefficients[:, None]
            * np.exp(
                -times / fitted_relaxation_times[:, None]
            ),
            axis=0,
        )
    )

    fitted_forces = force_scale * (
        force_matrix @ fitted_coefficients
    )

    prony_shear_modulus = instantaneous_shear_modulus * (
        1.0 - np.sum(prony_coefficients)
        + np.sum(
            prony_coefficients[:, None]
            * np.exp(
                -times / relaxation_times[:, None]
            ),
            axis=0,
        )
    )

    residual = force_matrix @ fitted_coefficients - normalized_forces
    mean_square_error = np.mean(residual**2)

    plot_results(
        times,
        numerical_shear_modulus,
        prony_shear_modulus,
        fitted_forces,
        forces,
    )
    plot_results(
        times,
        numerical_shear_modulus,
        prony_shear_modulus,
        fitted_forces,
        forces,
        False,
    )