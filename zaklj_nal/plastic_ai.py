"""Newton-Raphson plastic iteration for an Euler-Bernoulli cantilever beam.

Model assumptions:
- Cantilever beam with tip force load.
- Only bending normal stress sigma_xx is considered.
- Response is reduced to a single scalar DOF: tip deflection w.
- Bilinear elastoplastic material with isotropic hardening at the critical section
  (fixed end, outer fiber), integrated via 1D return mapping.

This is a compact teaching/example implementation of a global-local iteration:
- global Newton loop solves equilibrium R = P_ext - P_int(w) = 0,
- local constitutive update provides sigma_xx and consistent tangent modulus.
"""

# AI generated code !!!

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Tuple
import matplotlib.pyplot as plt


@dataclass
class MaterialState:
    """Committed internal variables for 1D plasticity."""

    eps_p: float = 0.0  # plastic strain
    alpha: float = 0.0  # equivalent plastic strain (for isotropic hardening)


@dataclass
class BeamResult:
    """Single converged load step output."""

    step: int
    p_ext: float
    w_tip: float
    sigma_xx: float
    iterations: int


def update_stress_1d_bilinear(
    eps: float,
    state_n: MaterialState,
    e_modulus: float,
    sigma_y0: float,
    hardening_modulus: float,
) -> Tuple[float, float, MaterialState]:
    """Return-mapping update for uniaxial bilinear elastoplasticity.

    Args:
        eps: total strain at the current iterate.
        state_n: committed state from previous converged load step.
        e_modulus: Young's modulus E.
        sigma_y0: initial yield stress.
        hardening_modulus: isotropic hardening modulus H.

    Returns:
        sigma: updated stress.
        et: algorithmic tangent modulus d(sigma)/d(eps).
        state_trial: updated internal variables (trial for this global iteration).
    """

    sigma_trial = e_modulus * (eps - state_n.eps_p)
    f_trial = abs(sigma_trial) - (sigma_y0 + hardening_modulus * state_n.alpha)

    if f_trial <= 0.0:
        # Elastic step
        return sigma_trial, e_modulus, MaterialState(state_n.eps_p, state_n.alpha)

    # Plastic correction
    dgamma = f_trial / (e_modulus + hardening_modulus)
    sign_trial = 1.0 if sigma_trial >= 0.0 else -1.0

    eps_p_new = state_n.eps_p + dgamma * sign_trial
    alpha_new = state_n.alpha + dgamma
    sigma = sign_trial * (sigma_y0 + hardening_modulus * alpha_new)

    # Consistent tangent for 1D bilinear isotropic hardening
    et = (e_modulus * hardening_modulus) / (e_modulus + hardening_modulus)

    return sigma, et, MaterialState(eps_p_new, alpha_new)


def solve_cantilever_tip_load(
    length: float,
    inertia: float,
    c_fiber: float,
    e_modulus: float,
    sigma_y0: float,
    hardening_modulus: float,
    p_max: float,
    n_steps: int = 50,
    tol: float = 1e-10,
    max_iter: int = 40,
) -> List[BeamResult]:
    """Incremental-iterative Newton solver for a tip-loaded cantilever beam."""

    # Kinematic/stress mappings for an EB cantilever under tip force
    # sigma_xx(root, outer fiber) = (3 E c / L^2) * w_tip in elasticity.
    strain_per_displacement = 3.0 * c_fiber / (length**2)

    # P = sigma * I / (L c)
    force_per_stress = inertia / (length * c_fiber)

    w_n = 0.0
    state_n = MaterialState()
    results: List[BeamResult] = []

    for step in range(1, n_steps + 1):
        p_ext = p_max * step / n_steps
        w = w_n

        converged = False
        state_trial = MaterialState(state_n.eps_p, state_n.alpha)
        sigma = 0.0

        for it in range(1, max_iter + 1):
            eps = strain_per_displacement * w
            sigma, et, state_trial = update_stress_1d_bilinear(
                eps, state_n, e_modulus, sigma_y0, hardening_modulus
            )

            p_int = force_per_stress * sigma
            residual = p_ext - p_int

            # Global tangent dP_int/dw
            k_tangent = force_per_stress * et * strain_per_displacement

            if abs(residual) < tol:
                converged = True
                results.append(
                    BeamResult(
                        step=step,
                        p_ext=p_ext,
                        w_tip=w,
                        sigma_xx=sigma,
                        iterations=it,
                    )
                )
                w_n = w
                state_n = state_trial
                break

            if abs(k_tangent) < 1e-20:
                raise RuntimeError(
                    f"Step {step}: near-zero tangent stiffness encountered in Newton loop."
                )

            dw = residual / k_tangent
            w += dw

        if not converged:
            raise RuntimeError(f"Step {step}: Newton did not converge in {max_iter} iterations.")

    return results


def _demo() -> None:
    """Run a small monotonic loading example and print selected steps."""

    length = 2.0
    width = 0.05
    height = 0.10

    inertia = width * height**3 / 12.0
    c_fiber = height / 2.0

    e_modulus = 210e9
    sigma_y0 = 250e6
    hardening_modulus = 2.0e9

    p_yield = sigma_y0 * inertia / (length * c_fiber)
    p_max = 1.6 * p_yield

    out = solve_cantilever_tip_load(
        length=length,
        inertia=inertia,
        c_fiber=c_fiber,
        e_modulus=e_modulus,
        sigma_y0=sigma_y0,
        hardening_modulus=hardening_modulus,
        p_max=p_max,
        n_steps=60,
    )

    print("step | P_ext [N] | w_tip [m] | sigma_xx [MPa] | iters")
    print("-" * 58)
    sample_ids = [0, 14, 29, 44, 59]
    for idx in sample_ids:
        r = out[idx]
        print(
            f"{r.step:>4d} | {r.p_ext:>9.2f} | {r.w_tip:>9.6e} |"
            f" {r.sigma_xx / 1e6:>13.3f} | {r.iterations:>5d}"
        )
    
    plt.plot(
        [r.w_tip for r in out],
        [r.sigma_xx / 1e6 for r in out],
        marker="o",
        linestyle="-",
        color="b",
    )
    plt.xlabel("Tip deflection w [m]")
    plt.ylabel("Bending stress sigma_xx [MPa]")
    plt.title("Bending stress vs. tip deflection for bilinear elastoplastic cantilever")
    plt.grid(True)
    plt.show()


if __name__ == "__main__":
    _demo()
