import numpy as np


class LinearElasticIndentation:

    def __init__(
        self,
        radius,
        depth,
        youngs_modulus,
        poisson_ratio,
        youngs_modulus_indenter=None,
        poisson_ratio_indenter=0.0,
    ):
        self.radius = radius
        self.depth = depth
        self.youngs_modulus = youngs_modulus
        self.poisson_ratio = poisson_ratio

        if youngs_modulus_indenter is None:
            self.reduced_youngs_modulus = youngs_modulus / (
                1 - poisson_ratio**2
            )
        else:
            self.reduced_youngs_modulus = (
                (1 - poisson_ratio**2) / youngs_modulus
                + (1 - poisson_ratio_indenter**2)
                / youngs_modulus_indenter
            ) ** (-1)

    @property
    def contact_radius(self):
        if self.depth <= 0:
            return 0.0
        return (self.radius * self.depth) ** 0.5

    def get_contact_area(self):
        return np.pi * self.radius * self.depth

    def compute_load(self):
        return (
            4
            / 3
            * self.reduced_youngs_modulus
            * self.radius**0.5
            * self.depth**(3 / 2)
        )

    def compute_contact_pressure(self, radial_distance=0.0):
        if radial_distance > self.contact_radius:
            print(
                f"Distance cannot exceed the contact radius distance "
                f"{self.contact_radius}"
            )
            return 0.0

        contact_pressure_0 = (
            2 * self.reduced_youngs_modulus / np.pi
        ) * (self.depth / self.radius) ** 0.5
        contact_pressure = (
            contact_pressure_0
            * (1 - (radial_distance / self.contact_radius) ** 2) ** 0.5
        )
        return contact_pressure

    def compute_radial_tensile_stress(self, radial_distance=None):
        if radial_distance is None:
            radial_distance = self.contact_radius

        contact_stress_0 = self.compute_contact_pressure(0.0)

        max_tensile_stress = (
            contact_stress_0 * (1 - 2 * self.poisson_ratio) / 3
        )

        if radial_distance < self.contact_radius:
            # Inside the footprint, the radial stress turns compressive.
            return 0.0

        tensile_stress = (
            max_tensile_stress * (self.contact_radius / radial_distance) ** 2
        )
        return tensile_stress


if __name__ == "__main__":
    indenter_radius = 15.875
    depth = 0.3
    youngs_modulus = 476.0142
    poisson_ratio = 0.48

    elastic_indentation = LinearElasticIndentation(
        indenter_radius, depth, youngs_modulus, poisson_ratio
    )
    contact_radius = elastic_indentation.contact_radius
    contact_area = elastic_indentation.get_contact_area()

    load = elastic_indentation.compute_load()
    contact_pressure = elastic_indentation.compute_contact_pressure(0)
    tensile_stress = elastic_indentation.compute_radial_tensile_stress()

    print(f"Contact radius: {contact_radius:.4f}")
    print(f"Contact Area: {contact_area:.4f}")
    print(f"Reaction Force: {load:.4f}")
    print(f"Contact Pressure: {contact_pressure:.4f}")
    print(f"Peak Tensile Stress: {tensile_stress:.4f}")
