import numpy as np
import matplotlib.pyplot as plt
from scipy.special import dawsn
from scipy.optimize import minimize


def compute_exact_convolution(
    times,
    relaxation_time,
    maximum_depth,
    ramp_duration,
):
    """Compute the exact exponential convolution for a linear ramp and hold."""

    velocity = (maximum_depth / ramp_duration) ** 1.5

    def compute_ramp_convolution(times):
        dimensionless_time = np.sqrt(times / relaxation_time)

        return (
            1.5
            * velocity
            * relaxation_time**1.5
            * (
                dimensionless_time
                - dawsn(dimensionless_time)
            )
        )

    ramp_end_convolution = compute_ramp_convolution(ramp_duration)

    return np.where(
        times <= ramp_duration,
        compute_ramp_convolution(np.minimum(times, ramp_duration)),
        np.exp(
            -np.maximum(times - ramp_duration, 0.0)
            / relaxation_time
        ) * ramp_end_convolution,
    )


def compute_exact_forces(
    times,
    prony_coefficients,
    relaxation_times,
    instantaneous_shear_modulus,
    contact_coefficient,
    maximum_depth,
    ramp_duration,
):
    """Compute the exact indentation force from the reference Prony series."""

    depth_pow_3_2 = (
        maximum_depth
        * np.minimum(times / ramp_duration, 1.0)
    ) ** 1.5

    forces = (
        instantaneous_shear_modulus
        * (1.0 - np.sum(prony_coefficients))
        * depth_pow_3_2
    )

    for prony_coefficient, relaxation_time in zip(
        prony_coefficients,
        relaxation_times,
    ):
        convolution_integral = compute_exact_convolution(
            times,
            relaxation_time,
            maximum_depth,
            ramp_duration,
        )

        forces += (
            instantaneous_shear_modulus
            * prony_coefficient
            * convolution_integral
        )

    return contact_coefficient * forces


def compute_instantaneous_shear_modulus(
    youngs_modulus,
    poisson_ratio,
):
    """Compute the instantaneous shear modulus."""
    return youngs_modulus / (2.0 * (1.0 + poisson_ratio))


def compute_contact_coefficient(
    indenter_radius,
    poisson_ratio,
):
    """Compute the contact coefficient."""
    return (
        8.0 * np.sqrt(indenter_radius)
        / (3.0 * (1.0 - poisson_ratio))
    )


def get_relaxation_time_range(
    times,
    min_relaxation_time=1e-4,
):
    """Return the shortest and longest relaxation times."""
    shortest = max(
        np.min(np.diff(times)),
        min_relaxation_time,
    )
    longest = times[-1] - times[0]

    if longest <= shortest:
        raise ValueError(
            "Test duration is not longer than the minimum relaxation time."
        )

    return shortest, longest


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

    return result.x


def plot_mean_square_errors(
    number_of_coefficients_list,
    force_mean_square_error_list,
    shear_modulus_mean_square_error_list,
    times_per_decade_list=None,
):
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))

    golden_ratio = (1 + np.sqrt(5)) / 2

    for ax in axes:
        ax.set_box_aspect(1 / golden_ratio)

    axes[0].plot(
        number_of_coefficients_list,
        force_mean_square_error_list,
        color="blue",
        marker="o",
        markerfacecolor="none",
    )

    axes[1].plot(
        number_of_coefficients_list,
        shear_modulus_mean_square_error_list,
        color="blue",
        marker="o",
        markerfacecolor="none",
    )

    axes[0].set_ylabel("Normalized force MSE")
    axes[1].set_ylabel("Normalized shear modulus MSE")

    for ax in axes:
        ax.set_xlabel("Number of Prony terms")

        ax.tick_params(
            axis="both",
            which="both",
            direction="in",
            top=True,
            right=True,
            width=1.0,
        )

        ax.minorticks_on()

        if times_per_decade_list is not None:
            times_per_decade_array = np.asarray(
                times_per_decade_list,
                dtype=float,
            )

            number_of_coefficients_array = np.asarray(
                number_of_coefficients_list,
                dtype=float,
            )

            integer_values = np.arange(
                1,
                int(np.floor(np.max(times_per_decade_array))) + 1,
            )

            tick_positions = np.interp(
                integer_values,
                times_per_decade_array,
                number_of_coefficients_array,
            )

            secondary_ax = ax.secondary_xaxis("top")

            secondary_ax.set_xticks(tick_positions)
            secondary_ax.set_xticklabels(
                [str(value) for value in integer_values]
            )

            secondary_ax.set_xlabel("Prony terms per decade")

            secondary_ax.tick_params(
                axis="x",
                which="major",
                direction="in",
                top=True,
                bottom=False,
                length=4,
                width=1.0,
            )

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    ramp_duration = 1.0
    hold_end_time = 11.0
    maximum_depth = 0.3

    number_of_data_points = 10000

    times = np.linspace(
        0.0,
        hold_end_time,
        number_of_data_points,
    )

    ramp_mask = times <= ramp_duration
    number_of_ramp_points = np.sum(ramp_mask)

    depths = maximum_depth * np.ones_like(times)
    depths[ramp_mask] = np.linspace(
        0.0,
        maximum_depth,
        number_of_ramp_points,
    )

    indenter_radius = 15.875
    youngs_modulus = 476.0142
    poisson_ratio = 0.48

    instantaneous_shear_modulus = compute_instantaneous_shear_modulus(
        youngs_modulus,
        poisson_ratio,
    )

    contact_coefficient = compute_contact_coefficient(
        indenter_radius,
        poisson_ratio,
    )

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

    # Ground truth
    exact_forces = compute_exact_forces(
        times=times,
        prony_coefficients=prony_coefficients,
        relaxation_times=relaxation_times,
        instantaneous_shear_modulus=instantaneous_shear_modulus,
        contact_coefficient=contact_coefficient,
        maximum_depth=maximum_depth,
        ramp_duration=ramp_duration,
    )

    exact_shear_modulus = instantaneous_shear_modulus * (
        1.0 - np.sum(prony_coefficients)
        + np.sum(
            prony_coefficients[:, None]
            * np.exp(
                -times / relaxation_times[:, None]
            ),
            axis=0,
        )
    )

    # Fit
    force_scale = np.max(exact_forces)
    normalized_forces = exact_forces / force_scale
    depth_pow_3_2 = depths**(3.0 / 2.0)

    number_of_coefficients_list = range(3, 31)
    times_per_decade_list = []

    force_mean_square_error_list = []
    shear_modulus_mean_square_error_list = []

    for number_of_coefficients in number_of_coefficients_list:
        fitted_relaxation_times = get_prony_relaxation_times(
            times=times,
            number_of_relaxation_times=number_of_coefficients,
        )

        times_per_decade = (
            (number_of_coefficients - 1)
            / np.log10(
                fitted_relaxation_times[-1]
                / fitted_relaxation_times[0]
            )
        )

        times_per_decade_list.append(times_per_decade)

        force_matrix = compute_force_matrix(
            times=times,
            depth_pow_3_2=depth_pow_3_2,
            relaxation_times=fitted_relaxation_times,
            contact_coefficient=contact_coefficient,
            instantaneous_shear_modulus=instantaneous_shear_modulus,
            force_scale=force_scale,
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

        fitted_forces = (
            force_scale
            * force_matrix
            @ fitted_coefficients
        )

        force_residual = fitted_forces - exact_forces
        shear_modulus_residual = (
            numerical_shear_modulus - exact_shear_modulus
        )

        force_mean_square_error = (
            np.sum(force_residual**2)
            / np.sum(exact_forces**2)
        )

        shear_modulus_mean_square_error = (
            np.sum(shear_modulus_residual**2)
            / np.sum(exact_shear_modulus**2)
        )

        force_mean_square_error_list.append(force_mean_square_error)
        shear_modulus_mean_square_error_list.append(
            shear_modulus_mean_square_error
        )

    print(
        "Prony terms per decade: ",
        np.round(times_per_decade_list, 2),
    )

    plot_mean_square_errors(
        number_of_coefficients_list,
        force_mean_square_error_list,
        shear_modulus_mean_square_error_list,
        times_per_decade_list,
    )