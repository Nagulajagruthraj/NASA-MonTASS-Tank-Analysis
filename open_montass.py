"""
OpenMonTASS v1.0 — Open-source Monocoque Tank Analysis & Sizing System

Performs thin-shell stress and buckling analysis for rocket propellant tank
structures (domes, cylinders, cones) using classical membrane theory and
NASA SP-8007 buckling knockdown factors.
"""

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Material definition — stores elastic and physical properties for a material.
# ---------------------------------------------------------------------------
class Material:
    """Holds material mechanical properties used throughout the analysis."""

    def __init__(self, name, E, nu, sigma_y, alpha=0.0, rho=0.0):
        self.name = name          # Material designation (e.g. "2219-T87 Aluminum")
        self.E = E                # Young's modulus (psi)
        self.nu = nu              # Poisson's ratio (dimensionless)
        self.sigma_y = sigma_y    # Yield strength (psi)
        self.alpha = alpha        # Coefficient of thermal expansion (1/°F)
        self.rho = rho            # Density (lb/in³)


# Pre-defined materials commonly used in aerospace tank design.
al_2219 = Material("2219-T87 Aluminum", 10.5e6, 0.33, 52000, 13.1e-6, 0.103)
ti_6al4v = Material("Ti-6Al-4V", 16.0e6, 0.31, 130000, 5.0e-6, 0.160)


# ---------------------------------------------------------------------------
# MonocoqueSection — represents a single shell element (cylinder or cone).
# ---------------------------------------------------------------------------
class MonocoqueSection:
    """Computes membrane stresses and buckling loads for a thin-walled shell."""

    def __init__(self, material, R, t, length=None, alpha_cone=0):
        self.mat = material       # Material object for this section
        self.R = R                # Reference radius (in)
        self.t = t                # Wall thickness (in)
        self.length = length      # Axial length of the section (in), optional
        self.alpha = alpha_cone   # Cone half-angle in degrees (0 = cylinder)

    def membrane_stress_pressure(self, P_internal):
        """Return (meridional, hoop) membrane stresses from internal pressure."""
        if self.alpha == 0:
            # Cylinder / sphere: standard thin-wall pressure vessel formulas.
            sigma_mer = P_internal * self.R / (2.0 * self.t)
            sigma_hoop = P_internal * self.R / self.t
        else:
            # Truncated cone: stresses increase by 1/cos(alpha).
            cos_a = np.cos(np.radians(self.alpha))
            sigma_mer = P_internal * self.R / (2.0 * self.t * cos_a)
            sigma_hoop = P_internal * self.R / (self.t * cos_a)
        return sigma_mer, sigma_hoop

    def buckling_critical(self, load_type="axial"):
        """Axial buckling load per unit width using NASA SP-8007 knockdown."""
        if load_type == "axial":
            # Knockdown factor gamma per NASA SP-8007 Rev 2, Eq. 15 (isotropic cylinder).
            gamma = 1.0 - 0.901 * (1.0 - np.exp(-(self.R / self.t) ** 0.16))
            N_cr = (
                gamma
                * (np.pi ** 2 * self.mat.E * self.t ** 2.5)
                / (3.0 * (1.0 - self.mat.nu ** 2) ** 0.75 * np.sqrt(self.R))
            )
            return N_cr
        return None


# ---------------------------------------------------------------------------
# OpenMonTASS — main analysis engine that collects sections and runs checks.
# ---------------------------------------------------------------------------
class OpenMonTASS:
    """Assembles tank sections, runs stress & buckling analysis, and reports."""

    def __init__(self):
        self.sections = []            # List of section dictionaries
        self.results = pd.DataFrame() # Analysis output table

    def add_dome(self, mat, R, t, shape="spherical"):
        """Add a dome (forward or aft) to the tank model."""
        self.sections.append(
            {"type": "dome", "shape": shape, "R": R, "t": t, "mat": mat}
        )

    def add_cone(self, mat, R1, R2, length, t, alpha):
        """Add a conical transition section between two radii."""
        self.sections.append(
            {
                "type": "cone",
                "R1": R1,
                "R2": R2,
                "L": length,
                "t": t,
                "alpha": alpha,
                "mat": mat,
                "R": R1,  # reference radius for analysis
            }
        )

    def analyze(self, P_int=0, P_ext=0, Nx=0, My=0, temp=70):
        """Run membrane stress and buckling analysis on every section.

        Parameters
        ----------
        P_int : float  — internal pressure (psi)
        P_ext : float  — external pressure (psi), reserved for future use
        Nx    : float  — axial line load (lb/in), reserved for future use
        My    : float  — bending moment (lb·in), reserved for future use
        temp  : float  — operating temperature (°F), reserved for future use

        Returns
        -------
        pandas.DataFrame with one row per section showing stresses, margins,
        and buckling loads.
        """
        FACTOR_OF_SAFETY = 1.1  # Design factor of safety applied to yield strength.
        rows = []
        for sec in self.sections:
            radius = sec["R"]
            thickness = sec.get("t", 0)
            cone_angle = sec.get("alpha", 0)

            # Minimum thickness required by pressure (with the design FoS).
            if P_int > 0:
                t_req = (P_int * radius) / (2.0 * sec["mat"].sigma_y / FACTOR_OF_SAFETY)
            else:
                t_req = thickness

            # Build a MonocoqueSection and compute membrane stresses.
            shell = MonocoqueSection(
                sec["mat"], radius, thickness, alpha_cone=cone_angle
            )
            sigma_m, sigma_h = shell.membrane_stress_pressure(P_int)

            # Yield margin of safety using von-Mises equivalent stress.
            von_mises = np.sqrt(sigma_m ** 2 + sigma_h ** 2 - sigma_m * sigma_h)
            MS_yield = (sec["mat"].sigma_y / von_mises) - 1.0

            # Axial buckling critical load per unit circumferential width.
            buckling = shell.buckling_critical()

            rows.append(
                {
                    "Component": sec["type"].capitalize(),
                    "Radius (in)": radius,
                    "Thickness (in)": thickness,
                    "t_required (in)": round(t_req, 4),
                    "Sigma_mer (psi)": round(sigma_m, 3),
                    "Sigma_hoop (psi)": round(sigma_h, 3),
                    "Margin_Yield": round(MS_yield, 3),
                    "Buckling_Load (lb/in)": round(buckling, 3) if buckling else None,
                }
            )

        self.results = pd.DataFrame(rows)
        return self.results


# ========================== DEMO ==========================================
if __name__ == "__main__":
    # 1. Create an analysis session.
    tass = OpenMonTASS()

    # 2. Define tank geometry — LH2-class tank with a dome and a cone.
    tass.add_dome(al_2219, R=60, t=0.12)                          # forward dome
    tass.add_cone(al_2219, R1=60, R2=45, length=120, t=0.10, alpha=15)  # transition cone

    # 3. Run the analysis at 25 psi internal pressure.
    results = tass.analyze(P_int=25, P_ext=0, Nx=50000, temp=70)

    # 4. Print results to the console.
    print("=== OpenMonTASS v1.0 Analysis ===")
    print(results.to_string(index=False))

    # 5. Export to Excel (same column layout as the original MonTASS spreadsheet).
    results.to_excel("OpenMonTASS_Analysis.xlsx", index=False)
    print("\n✅ Excel report saved: OpenMonTASS_Analysis.xlsx")
