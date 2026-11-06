from abaqus import session
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize


class ViscoelasticIndentation:
    def __init__(
        self,
        indenter_radius,
        poisson_ratio,
        youngs_modulus=None,
        prony_coefficients=None,
        relaxation_times=None,
    ):
        self.indenter_radius = indenter_radius
        self.poisson_ratio = poisson_ratio        
        if youngs_modulus:
            self.instantaneous_shear_modulus = youngs_modulus / (
                2.0 * (1.0 + poisson_ratio)
            )
        
        self.contact_coefficient = (
            8.0 * np.sqrt(indenter_radius) / (3.0 * (1.0 - poisson_ratio))
        )
        self.prony_coefficients = prony_coefficients
        self.relaxation_times = relaxation_times
                
        self.current_viewport_name = session.currentViewportName
        self.current_viewport = session.viewports[self.current_viewport_name]
        
        # Pull displayed object and ensure it contains valid step metadata
        displayed_obj = self.current_viewport.displayedObject
        if hasattr(displayed_obj, 'steps'):
            self.odb = displayed_obj
        else:
            raise TypeError(
                "Active viewport is currently showing a Part or Assembly. "
                "Please switch to the Visualization module before pasting."
            )
            
        self.step_items = self.odb.steps.items()
        self._indentation_depth_times = []
        self.depths = []
        self._normal_force_times = []
        self.forces = []
        self.times = None
        self.closed_form_shear_modulus = None
        self._u2_data = {}
        self._rf2_data = {}
        self.step_data = []
        self.start_time_hold = 0.0
        for step_name, step in self.step_items:
            self.step_data.append(
                (
                    step_name,
                    step.procedure,
                    step.totalTime,
                    step.timePeriod,
                )
            )
            self._u2_data[step_name] = {}
            self._rf2_data[step_name] = {}
        self._extract_history_data()
        self._extract_indentation_data()
        
        if (
            self.depths is not None 
            and self.forces is not None
        ):
            self.closed_form_shear_modulus = (
                self.compute_closed_form_shear_modulus()
            )
            self.numerical_shear_modulus = (
                self.compute_numerical_shear_modulus()
            )
        
        self.prony_shear_modulus = self.compute_prony_shear_modulus()
    
    
    def _extract_history_data(self):
        """Internal helper to extract matching U2 and RF2 data."""
        for step_name, step in self.step_items:
            for region_name, region in step.historyRegions.items():
                if "U2" in region.historyOutputs:
                    self._u2_data[step_name][region_name] = (
                        region.historyOutputs["U2"].data
                    )
                if "RF2" in region.historyOutputs:
                    self._rf2_data[step_name][region_name] = (
                        region.historyOutputs["RF2"].data
                    )
    
    def _extract_indentation_data(self):
        step_data_list = []
        for data in self.step_data:
            if "*visco" in data[1].lower():
                step_data_list.append(data)
        if not step_data_list:
            print("Warning: No *Visco step was found")
            return
        times = []
        for step_data in step_data_list:
            step_name = step_data[0]
            regions_u2 = self._u2_data[step_name].keys()
            if len(regions_u2) > 1:
                print(
                    "Warning: More than one node found for the U2 output. "
                    "Ensure that the correct rigid body node has been selected"
                )
            regions_rf2 = self._rf2_data[step_name].keys()
            if len(regions_rf2) > 1:
                print(
                    "Warning: More than one node found for the RF2 output. "
                    "All forces are summed. Ensure that the boundary sets "
                    "have been defined in the history output"
                )
            # 1. Extract and concatenate U2 data inline
            u2_region = next(iter(regions_u2))
            u2_array = np.array(self._u2_data[step_name][u2_region])
            
            self._indentation_depth_times.append(u2_array[:, 0])
            self.depths.append(abs(u2_array[:, 1]))
            
            # 2. Extract first RF2 region for THIS step
            rf2_regions_list = list(regions_rf2)
            first_rf2_region = rf2_regions_list[0]
            rf2_array = np.array(self._rf2_data[step_name][first_rf2_region])
            
            step_force_times = rf2_array[:, 0]
            step_forces = abs(rf2_array[:, 1])
            if not np.array_equal(u2_array[:, 0], rf2_array[:, 0]):
                print(
                    "Time samples do not match for the displacements and "
                    "normal forces"
                )
            times.append(step_force_times)
            # 3. Sum any remaining RF2 regions for this step
            for rf2_region in rf2_regions_list[1:]:
                additional_rf2_array = np.array(
                    self._rf2_data[step_name][rf2_region]
                )
                step_forces += abs(additional_rf2_array[:, 1])
            
            # 4. Concatenate RF2 data inline
            self._normal_force_times.append(step_force_times)
            self.forces.append(step_forces)
        self.times = times
    
    def _get_indentation_history(self, method="prony"):
        if method == "closed" and len(self.times) > 1:
            indices = [-1]
        elif method in ("numerical", "prony") or len(self.times) == 1:
            indices = range(len(self.times))
        else:
            raise ValueError("Unknown method.")
        
        times = []
        depths = []
        forces = []
        
        time_offset = 0.0
        
        for i in indices:
            if times:
                step_times = self.times[i][1:]
                step_depths = self.depths[i][1:]
                step_forces = self.forces[i][1:]
            else:
                step_times = self.times[i]
                step_depths = self.depths[i]
                step_forces = self.forces[i]
            
            times.extend(step_times + time_offset)
            depths.extend(step_depths)
            forces.extend(step_forces)
            
            time_offset += self.step_data[i][3]
            
        return (
            np.asarray(times),
            np.asarray(depths),
            np.asarray(forces),
        )
    
    
    def _get_relaxation_time_range(self, times, min_relaxation_time=1e-4):
        """Return (shortest, longest) relaxation time for the Prony grid."""
        shortest = min(np.min(np.diff(times)), min_relaxation_time)
        longest = times[-1] - times[0]
        if longest <= shortest:
            raise ValueError(
                "Test duration is not longer than the minimum relaxation time"
            )
        return shortest, longest
    
    
    def _get_number_of_coefficients(
        self, 
        times, 
        terms_per_decade=2, 
        min_relaxation_time=1e-4
    ):
        """Number of Prony terms for the given density per decade of time."""
        shortest, longest = self._get_relaxation_time_range(
            times, 
            min_relaxation_time
        )
        
        return int(
            np.ceil(terms_per_decade * np.log10(longest / shortest))
        ) + 1
    
    
    def _get_prony_relaxation_times(
        self, 
        times,
        number_of_relaxation_times, 
        min_relaxation_time=1e-4
    ):
        """Logarithmically spaced Prony relaxation times."""
        shortest, longest = self._get_relaxation_time_range(
            times,
            min_relaxation_time
        )
        return np.logspace(
            np.log10(shortest), 
            np.log10(longest), 
            number_of_relaxation_times
        )
    def _compute_convolution_integral(
        self,
        times,
        depth_pow_3_2,
        relaxation_time,
    ):
        """Compute the exponential convolution integral."""
        
        if np.allclose(depth_pow_3_2, depth_pow_3_2[0]):
            return (
                depth_pow_3_2[0]
                * np.exp(-(times - times[0]) / relaxation_time)
            )
        
        # Piecewise-linear indentation history
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
        
        
    def _compute_force_matrix(
        self,
        times,
        depth_pow_3_2,
        relaxation_times,
        force_scale,
    ):  
        """Compute the normalised force matrix, so that
            force_matrix @ [1, g_1, ..., g_n]  =  F(t) / force_scale
        
        Column 0 multiplies the fixed coefficient 1 (instantaneous response);
        the remaining columns multiply the Prony coefficients g_i.
        """
        convolution_integrals = [
            self._compute_convolution_integral(
                times,
                depth_pow_3_2,
                relaxation_time
            )
            for relaxation_time in relaxation_times
        ]
        
        return (
            self.contact_coefficient
            * self.instantaneous_shear_modulus
            / force_scale
            * np.column_stack(
                [
                    depth_pow_3_2,
                    *[
                        integral - depth_pow_3_2
                        for integral in convolution_integrals
                    ],
                ]
            )
        )
    
    
    def _compute_regularization_matrix(self, number_of_terms):
        """Compute the second-order regularisation matrix."""
        regularization_matrix = np.zeros(
            (number_of_terms - 2, number_of_terms + 1)
        )
        
        for i in range(number_of_terms - 2):
            regularization_matrix[i, i + 1:i + 4] = [1, -2, 1]
        
        return regularization_matrix
    
    def _objective_function(
        self,
        parameters,
        force_matrix,
        normalized_forces,
        regularization_matrix,
        regularization_parameter,
    ):
        """Minimize ||Ax - b||² + λ²||Dx||²."""
        
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
    
    
    def _objective_gradient(
        self,
        parameters,
        force_matrix,
        normalized_forces,
        regularization_matrix,
        regularization_parameter,
    ):
        """Gradient of ||Ax - b||² + λ²||Dx||² with respect to x."""
        if len(parameters)< 2:
            raise ValueError(
                "At least one Prony coefficient is required."
            )
        residual = force_matrix @ parameters - normalized_forces
        
        gradient = (
            2.0 * force_matrix.T @ residual
            + 2.0 * regularization_parameter ** 2
            * regularization_matrix.T @ (regularization_matrix @ parameters)
        )
        
        return gradient / len(force_matrix)
        
        
    def _compute_numerical_coefficients(
        self,
        force_matrix,
        normalized_forces,
        regularization_matrix,
        regularization_parameter=0.0,
    ):
        """Solve the constrained regularized least-squares problem."""
        
        number_of_parameters = force_matrix.shape[1]
        
        initial_guess = np.zeros(number_of_parameters)
        initial_guess[0] = 1.0 # Normalized instantaneous shear modulus 
        initial_guess[1:] = 0.5 / (number_of_parameters - 1)
        
        jacobian_constraint_instantaneous_modulus = np.zeros(
            number_of_parameters
        )
        jacobian_constraint_instantaneous_modulus[0] = 1.0
        
        jacobian_constraint_prony_coefficients = -np.ones(number_of_parameters)
        jacobian_constraint_prony_coefficients[0] = 0
        
        result = minimize(
            self._objective_function,
            initial_guess,
            args=(
                force_matrix,
                normalized_forces,
                regularization_matrix,
                regularization_parameter,
            ),
            jac=self._objective_gradient,
            method="SLSQP",
            bounds=[(0.0, None)] * number_of_parameters,
            constraints=[
                {
                    "type": "eq",
                    "fun": lambda x: x[0] - 1.0,
                    "jac": lambda x: jacobian_constraint_instantaneous_modulus
                },
                {
                    "type": "ineq",
                    "fun": lambda x: 1.0 - np.sum(x[1:]),
                    "jac": lambda x: jacobian_constraint_prony_coefficients
                },
            ],
            options={"maxiter": 1000, "ftol": 1e-15},
        )
            
        if not result.success:
            raise RuntimeError(result.message)
        
        return result.x
    
    
    def compute_prony_shear_modulus(self):
        """Computes the analytical Prony Series curve over the timeline."""
        times, _, _ = self._get_indentation_history("prony")
        if times is None:
            raise ValueError("No indentation time history is available.")
        
        if self.instantaneous_shear_modulus is None:
            raise ValueError("The instantaneous shear modulus is not defined.")
        
        if self.prony_coefficients is None or self.relaxation_times is None:
            raise ValueError(
                "Prony coefficients and relaxation times must be provided."
            )
        
        if len(self.relaxation_times) != len(self.prony_coefficients):
            raise ValueError(
                "The number of relaxation times must match "
                "the number of Prony coefficients."
            )
        
        prony_shear_modulus = np.full_like(
            times, self.instantaneous_shear_modulus, dtype=float
        )
        
        for tau_i, g_i in zip(self.relaxation_times, self.prony_coefficients):
            prony_shear_modulus -= (
                self.instantaneous_shear_modulus
                * g_i
                * (1 - np.exp(-times / tau_i))
            )
        
        return prony_shear_modulus
    
    
    def compute_closed_form_shear_modulus(self):
        """Compute the shear modulus using the Fadil et al. closed-form
        solution."""
        _, depths, forces = (
            self._get_indentation_history("closed")
        )
        
        indentation_depth_constant = depths[0]
        
        contact_radius_constant = np.sqrt(
            self.indenter_radius * indentation_depth_constant
        )
        
        closed_form_shear_modulus = (
            (3 / 8)
            * forces
            * (1 - self.poisson_ratio)
            * self.indenter_radius
            / contact_radius_constant**3
        )
        return closed_form_shear_modulus
    
    
    def compute_numerical_shear_modulus(self):
        """Compute the shear modulus using the Fadil et al. numerical
        scheme."""
        times, depths, forces = (
            self._get_indentation_history("numerical")
        )
        
        depth_pow_3_2 = depths ** (3 / 2)
        force_scale = max(forces)
        normalized_forces = forces / force_scale
        
        number_of_coefficients = self._get_number_of_coefficients(
            times,
            2,
            1e-4,
        )
        
        fitted_relaxation_times = self._get_prony_relaxation_times(
            times, 
            number_of_coefficients
        )
        
        force_matrix = self._compute_force_matrix(
            times,
            depth_pow_3_2,
            fitted_relaxation_times,
            force_scale,
        )
        
        regularization_matrix = self._compute_regularization_matrix(
            number_of_coefficients
        )
        
        fitted_coefficients = self._compute_numerical_coefficients(
            force_matrix,
            normalized_forces,
            regularization_matrix,
            0.0,
        )
        
        fitted_prony_coefficients = fitted_coefficients[1:]
        
        numerical_shear_modulus = np.full_like(
            times, self.instantaneous_shear_modulus, dtype=float
        )
        
        for tau_i, g_i in zip(
            fitted_relaxation_times, 
            fitted_prony_coefficients
        ):
            numerical_shear_modulus -= (
                self.instantaneous_shear_modulus
                * g_i
                * (1 - np.exp(-times / tau_i))
            )
        return numerical_shear_modulus
    
    
    def plot_shear_moduli(self, xlog=True):
        """Plot Prony, closed-form, and numerical shear moduli."""
        
        if self.closed_form_shear_modulus is None:
            print("Error: No data compiled to plot.")
            return
        
        # Complete time history
        times_prony, _, _ = self._get_indentation_history("prony")
        prony_shear_modulus = self.prony_shear_modulus
        
        times_closed, _, _ = self._get_indentation_history("closed")
        closed_form_shear_modulus = self.closed_form_shear_modulus
        
        times_numerical, _, _ = self._get_indentation_history("numerical")
        numerical_shear_modulus = self.numerical_shear_modulus
        
        # Remove initial point for logarithmic axis
        if xlog:
            times_prony = times_prony[1:]
            prony_shear_modulus = prony_shear_modulus[1:]
            
            times_closed = times_closed[1:]
            closed_form_shear_modulus = closed_form_shear_modulus[1:]
            
            times_numerical = times_numerical[1:]
            numerical_shear_modulus = numerical_shear_modulus[1:]
            
        # Evenly spaced markers in the displayed coordinate system
        spacing = np.geomspace if xlog else np.linspace
        
        marker_targets_closed = spacing(
            times_closed[0],
            times_closed[-1],
            5,
        )
        marker_indices_closed = np.unique(
            np.abs(
                times_closed[:, None] - marker_targets_closed
            ).argmin(axis=0)
        )
        
        marker_targets_numerical = spacing(
            times_numerical[0],
            times_numerical[-1],
            5,
        )
        marker_indices_numerical = np.unique(
            np.abs(
                times_numerical[:, None] - marker_targets_numerical
            ).argmin(axis=0)
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
            fig, ax = plt.subplots(figsize=(6.47, 4))
            
            ax.plot(
                times_prony,
                prony_shear_modulus,
                color="forestgreen",
                linewidth=1.0,
                zorder=2,
                label="Reference Prony series",
            )
            
            ax.plot(
                times_closed,
                closed_form_shear_modulus,
                color="blue",
                linestyle="--",
                linewidth=1.0,
                marker="o",
                markersize=4.2,
                mfc="none",
                mew=1.0,
                markevery=list(marker_indices_closed),
                zorder=4,
                label="Fadil et al., Equation (16)",
            )
            
            ax.plot(
                times_numerical,
                numerical_shear_modulus,
                color="red",
                linestyle="--",
                linewidth=1.0,
                marker="o",
                markersize=4.2,
                mfc="none",
                mew=1.0,
                markevery=list(marker_indices_numerical),
                zorder=4,
                label="Fadil et al., Equation (13)",
            )
            
            if xlog:
                ax.set_xscale("log")
            
            ax.set_xlabel("Time (s)")
            ax.set_ylabel("Shear Modulus (MPa)")
            ax.tick_params(top=True, right=True, width=1.0)
            ax.set_axisbelow(True)
            ax.legend(frameon=True, fancybox=False, edgecolor="black")
            
            fig.tight_layout()
            plt.show()


if __name__ == "__main__":
    indenter_radius = 15.875
    poisson_ratio = 0.48
    youngs_modulus = 476.0142
    
    prony_coefficients = np.array(
        [
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
            
        ]
    )
    
    relaxation_times = np.array(
    [
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
        400.0000,
    ])
    viscoelastic_indentation = ViscoelasticIndentation(
        indenter_radius=indenter_radius, 
        poisson_ratio=poisson_ratio,
        youngs_modulus=youngs_modulus,
        prony_coefficients=prony_coefficients,
        relaxation_times=relaxation_times
    )
    viscoelastic_indentation.plot_shear_moduli()
    viscoelastic_indentation.plot_shear_moduli(False)
    