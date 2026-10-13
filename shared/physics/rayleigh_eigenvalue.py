"""
rayleigh_eigenvalue.py — Rayleigh eigenvalue problem for the tanh shear layer.

Solves (U - c)(φ'' - k²φ) = U''φ for the most-unstable c = c_r + i·c_i,
returning σ_nd = k · c_i · δ / U0 (non-dimensional growth rate per Michalke 1964).

Two independent methods are provided for cross-validation:
  rayleigh_shooting  — shooting method with complex RK4 integration
  rayleigh_cheb      — Chebyshev pseudospectral collocation (generalised EVP)

Reference: Michalke, A. (1964). J. Fluid Mech. 19(4):543-556.
           Drazin, P.G., Reid, W.H. (1981). Hydrodynamic Stability. Cambridge. §4.2.
           Orszag, S.A. (1971). J. Fluid Mech. 50:689-703.  [Chebyshev collocation]
"""

import numpy as np
from scipy import linalg, optimize


# ---------------------------------------------------------------------------
# Base-flow helpers for tanh profile
# ---------------------------------------------------------------------------

def _tanh_profile(y, y_c, delta):
    """U(y) = tanh((y - y_c) / delta)"""
    eta = (y - y_c) / delta
    return np.tanh(eta)


def _tanh_d2U(y, y_c, delta):
    """U''(y) = -2 tanh(η) sech²(η) / δ²"""
    eta = (y - y_c) / delta
    t = np.tanh(eta)
    s2 = 1.0 / np.cosh(eta) ** 2
    return -2.0 * t * s2 / delta ** 2


# ---------------------------------------------------------------------------
# Method A: Chebyshev pseudospectral collocation  (Orszag 1971)
# ---------------------------------------------------------------------------

def _cheb_diff_matrix(N):
    """
    Return the N+1 × N+1 Chebyshev differentiation matrix D on z ∈ [-1, 1].
    Nodes: z_j = cos(j π / N), j = 0 … N  (Chebyshev–Gauss–Lobatto).

    Reference: Trefethen (2000) Spectral Methods in MATLAB, Program 6.
    """
    if N == 0:
        return np.array([[0.0]]), np.array([1.0])
    z = np.cos(np.pi * np.arange(N + 1) / N)
    c = np.ones(N + 1)
    c[0] = 2.0
    c[N] = 2.0
    c *= (-1) ** np.arange(N + 1)
    Z = np.tile(z, (N + 1, 1))
    dZ = Z - Z.T
    D = np.outer(c, 1.0 / c) / (dZ + np.eye(N + 1))
    D -= np.diag(D.sum(axis=1))
    return D, z


def rayleigh_cheb(k, Ly, y_c, delta, N=80):
    """
    Solve the Rayleigh eigenvalue problem via Chebyshev pseudospectral collocation.

    Parameters
    ----------
    k     : float  — streamwise wavenumber (physical units, k = km * 2π / Lx)
    Ly    : float  — domain height
    y_c   : float  — shear-layer centre (= Ly/2 for symmetric domain)
    delta : float  — shear-layer momentum thickness
    N     : int    — number of Chebyshev intervals (N+1 nodes incl. boundaries)

    Returns
    -------
    sigma_nd : float — non-dimensional growth rate k·c_i·δ/U0 (Michalke convention)
    c        : complex — phase speed (c_r + i·c_i); c_i > 0 → unstable
    phi      : ndarray (N-1,) complex — eigenfunction on interior nodes
    y_int    : ndarray (N-1,) — y-coordinates of interior nodes
    """
    # Build Chebyshev differentiation matrix on z ∈ [-1, 1]
    D_z, z = _cheb_diff_matrix(N)
    # Map z → y: y = Ly/2 * (z + 1)
    # dy/dz = Ly/2, so D_y = (2/Ly)*D_z, D2_y = (2/Ly)^2 * D_z @ D_z
    scale = 2.0 / Ly
    D1 = scale * D_z
    D2 = scale ** 2 * (D_z @ D_z)

    # y-coordinates at all N+1 nodes
    y_all = 0.5 * Ly * (z + 1.0)

    # Base-flow arrays at all nodes
    U_all = _tanh_profile(y_all, y_c, delta)
    d2U_all = _tanh_d2U(y_all, y_c, delta)

    # Interior indices (exclude boundary rows j=0 and j=N where φ=0)
    interior = slice(1, N)
    n = N - 1  # number of interior points

    D2_int = D2[interior, interior]
    U_int = U_all[interior]
    d2U_int = d2U_all[interior]

    # Rayleigh equation: (U - c)(φ'' - k²φ) = U''φ
    # Rearranged:  [diag(U)·(D2 - k²I) - diag(U'')] φ = c · (D2 - k²I) φ
    # Generalised EVP: A φ = c B φ
    L = D2_int - k ** 2 * np.eye(n)
    A = np.diag(U_int) @ L - np.diag(d2U_int)
    B = L.copy()

    # Handle boundary columns: interior D2 rows already exclude BC; but the
    # off-diagonal D2 coupling to boundary points needs zeroing (BC: φ_bdy = 0)
    # Because we sliced both rows AND columns to [1:N], boundary coupling is
    # already excluded by the column slice — no further action needed.

    # Solve generalised EVP
    eigenvalues, eigenvectors = linalg.eig(A, B)

    # Filter: keep only eigenvalues with Im(c) > 0 (unstable) and |Re(c)| < 1.5
    # Threshold 1e-8 retains near-cutoff modes (Im(c) can be as small as ~0.04).
    mask = (np.imag(eigenvalues) > 1e-8) & (np.abs(np.real(eigenvalues)) < 1.5)
    if not np.any(mask):
        return 0.0, 0.0 + 0.0j, np.zeros(n), y_all[interior]

    eigs_filtered = eigenvalues[mask]
    vecs_filtered = eigenvectors[:, mask]

    # Most-unstable eigenvalue = largest Im(c)
    idx = np.argmax(np.imag(eigs_filtered))
    c = eigs_filtered[idx]
    phi = vecs_filtered[:, idx]

    sigma_nd = k * np.imag(c) * delta  # = Im(c) * k * δ (U0 = 1)
    return float(np.real(sigma_nd)), c, phi, y_all[interior]


# ---------------------------------------------------------------------------
# Method B: Shooting method with complex RK4 integration
# ---------------------------------------------------------------------------

def _rayleigh_rhs(y, phi_vec, k, y_c, delta, c):
    """
    RHS of the Rayleigh ODE written as a first-order system.

    State: [φ, φ'] (complex)
    φ'' = k²φ + U''(y)/(U(y) - c) * φ
    """
    phi, dphi = phi_vec
    U = _tanh_profile(y, y_c, delta)
    d2U = _tanh_d2U(y, y_c, delta)
    denom = U - c
    phi_pp = k ** 2 * phi + d2U / denom * phi
    return np.array([dphi, phi_pp], dtype=complex)


def _rk4_step_complex(y, phi_vec, dy, k, y_c, delta, c):
    """Single RK4 step in y for the complex Rayleigh ODE."""
    k1 = _rayleigh_rhs(y,            phi_vec,          k, y_c, delta, c)
    k2 = _rayleigh_rhs(y + 0.5 * dy, phi_vec + 0.5 * dy * k1, k, y_c, delta, c)
    k3 = _rayleigh_rhs(y + 0.5 * dy, phi_vec + 0.5 * dy * k2, k, y_c, delta, c)
    k4 = _rayleigh_rhs(y + dy,       phi_vec + dy * k3,         k, y_c, delta, c)
    return phi_vec + (dy / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)


def _shoot_from_boundary(k, y_c, delta, c, y_start, y_end, n_steps):
    """
    Integrate the Rayleigh ODE from y_start to y_end using RK4.
    Initial conditions are set by the far-field evanescent behaviour:
      φ ~ exp(-k·|y - y_c|)  far from the critical layer.
    """
    dy = (y_end - y_start) / n_steps
    # Far-field decaying IC: φ ~ exp(-k*|y - y_c|) far from shear layer centre.
    # At y_start (far from y_c): φ'(y_start) = +k*φ if integrating toward y_c from below,
    #                                           = -k*φ if integrating toward y_c from above.
    # Sign: d/dy exp(-k*(y_c - y)) = +k*exp(-k*(y_c-y)) → +k*φ (y < y_c, shooting upward)
    #        d/dy exp(-k*(y - y_c)) = -k*exp(-k*(y-y_c)) → -k*φ (y > y_c, shooting downward)
    phi0 = np.exp(-k * abs(y_start - y_c))
    dphi0 = +np.sign(y_end - y_start) * k * phi0
    state = np.array([phi0, dphi0], dtype=complex)

    y = float(y_start)
    for _ in range(n_steps):
        state = _rk4_step_complex(y, state, dy, k, y_c, delta, c)
        y += dy
    return state


def _wronskian(phi_L, phi_R):
    """
    Wronskian W(φ_L, φ_R) = φ_L * φ'_R - φ'_L * φ_R.
    Zero ↔ linearly dependent solutions ↔ eigenvalue.
    """
    return phi_L[0] * phi_R[1] - phi_L[1] * phi_R[0]


def rayleigh_shooting(k, Ly, y_c, delta, n_steps=400, c_i_lo=1e-4, c_i_hi=0.75,
                      n_bracket=80, tol=1e-8, maxiter=60):
    """
    Solve the Rayleigh eigenvalue problem via the shooting method.

    For the symmetric tanh profile, c_r = 0 by symmetry; we search for c_i > 0.
    Shoots from both boundaries to the midpoint y_c and matches via the Wronskian.

    Parameters
    ----------
    k        : float — streamwise wavenumber
    Ly       : float — domain height
    y_c      : float — shear-layer centre
    delta    : float — shear-layer thickness
    n_steps  : int   — RK4 steps per half-domain per evaluation
    c_i_lo   : float — lower bound for c_i search
    c_i_hi   : float — upper bound for c_i search
    n_bracket: int   — number of points in the initial sign-change scan
    tol      : float — Wronskian bisection convergence tolerance
    maxiter  : int   — maximum bisection iterations

    Returns
    -------
    sigma_nd : float — non-dimensional growth rate k·c_i·δ/U0
    c        : complex — eigenvalue (c_r = 0 for tanh, c_i > 0)
    """
    y_mid = y_c  # matching point at shear layer centre

    def wronskian_abs(c_i):
        """Absolute value of the Wronskian at y_mid for c = 0 + i*c_i.
        Has a clear global minimum at the eigenvalue, unlike Im(W) which is noisy."""
        c = 0.0 + 1j * c_i
        phi_L = _shoot_from_boundary(k, y_c, delta, c, 0.0, y_mid, n_steps)
        phi_R = _shoot_from_boundary(k, y_c, delta, c, Ly,  y_mid, n_steps)
        return abs(_wronskian(phi_L, phi_R))

    # 1. Coarse scan: find the global minimum of |W(c_i)|
    c_i_vals = np.linspace(c_i_lo, c_i_hi, n_bracket)
    W_abs = np.array([wronskian_abs(ci) for ci in c_i_vals])
    i_min = np.argmin(W_abs)

    # Check whether there is a genuine minimum (not just a boundary edge)
    if i_min == 0 or i_min == len(c_i_vals) - 1:
        return 0.0, 0.0 + 0.0j

    # 2. Refine with Brent's method in a bracket around the coarse minimum
    bracket_lo = c_i_vals[max(0, i_min - 2)]
    bracket_hi = c_i_vals[min(len(c_i_vals) - 1, i_min + 2)]
    result = optimize.minimize_scalar(wronskian_abs, bounds=(bracket_lo, bracket_hi),
                                       method='bounded',
                                       options={'xatol': tol, 'maxiter': maxiter})

    c_i_star = float(result.x)
    c_star = 0.0 + 1j * c_i_star
    sigma_nd = k * c_i_star * delta
    return float(sigma_nd), c_star


# ---------------------------------------------------------------------------
# Comparison utility: both methods across Michalke Table 1 wavenumbers
# ---------------------------------------------------------------------------

def compare_methods(Ly=None, delta=None, N_cheb=80, n_steps_shoot=400):
    """
    Run both methods across the k*delta grid from Michalke (1964) Table 1
    and return a comparison table.

    Returns dict with keys 'k_delta', 'sigma_michalke', 'sigma_cheb', 'sigma_shoot'.
    """
    # Domain parameters matching the oracle test
    if Ly is None:
        Ly = 2.0 * np.pi
    if delta is None:
        delta = 0.05 * Ly  # = 0.1π ≈ 0.3142

    y_c = Ly / 2.0

    # Michalke (1964) Table 1 values (k*delta → sigma_nd)
    k_delta_theory = np.array([0.1, 0.2, 0.3, 0.4, 0.4446, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    sigma_michalke  = np.array([0.064, 0.121, 0.162, 0.184, 0.1897, 0.186, 0.171, 0.148, 0.116, 0.079, 0.038])

    # Physical k from k*delta
    k_vals = k_delta_theory / delta

    results = {
        'k_delta':        k_delta_theory,
        'sigma_michalke': sigma_michalke,
        'sigma_cheb':     np.zeros_like(k_delta_theory),
        'sigma_shoot':    np.zeros_like(k_delta_theory),
    }

    for i, (kd, k) in enumerate(zip(k_delta_theory, k_vals)):
        sigma_c, _, _, _ = rayleigh_cheb(k, Ly, y_c, delta, N=N_cheb)
        sigma_s, _       = rayleigh_shooting(k, Ly, y_c, delta, n_steps=n_steps_shoot)
        results['sigma_cheb'][i]  = sigma_c
        results['sigma_shoot'][i] = sigma_s

    return results


def print_comparison(results):
    """Print a formatted comparison table."""
    kd  = results['k_delta']
    sm  = results['sigma_michalke']
    sc  = results['sigma_cheb']
    ss  = results['sigma_shoot']

    print(f"{'k*delta':>8}  {'Michalke':>10}  {'Chebyshev':>10}  {'Shooting':>10}  "
          f"{'err_cheb%':>9}  {'err_shoot%':>10}")
    print("-" * 70)
    for i in range(len(kd)):
        ec = abs(sc[i] - sm[i]) / sm[i] * 100 if sm[i] > 0 else float('nan')
        es = abs(ss[i] - sm[i]) / sm[i] * 100 if sm[i] > 0 else float('nan')
        print(f"{kd[i]:8.4f}  {sm[i]:10.4f}  {sc[i]:10.4f}  {ss[i]:10.4f}  "
              f"{ec:9.2f}  {es:10.2f}")

    # Summary
    valid = sm > 0
    rmse_cheb  = np.sqrt(np.mean((sc[valid] - sm[valid]) ** 2 / sm[valid] ** 2)) * 100
    rmse_shoot = np.sqrt(np.mean((ss[valid] - sm[valid]) ** 2 / sm[valid] ** 2)) * 100
    print(f"\nRMS relative error:  Chebyshev={rmse_cheb:.2f}%   Shooting={rmse_shoot:.2f}%")
    winner = "Chebyshev" if rmse_cheb <= rmse_shoot else "Shooting"
    print(f"Winner: {winner}")
    return winner


# ---------------------------------------------------------------------------
# Convenience: vorticity eigenfunction on a uniform y-grid for simulation seeding
# ---------------------------------------------------------------------------

def rayleigh_cheb_vort_seed(k, Ly, y_c, delta, y_uniform, N=100):
    """
    Return the vorticity eigenfunction ω_eig(y) = k²φ - φ'' normalised to max 1,
    evaluated at the supplied uniform y-grid, using the Chebyshev solution.

    Parameters
    ----------
    k         : float — streamwise wavenumber
    Ly        : float — domain height
    y_c       : float — shear-layer centre
    delta     : float — shear-layer thickness
    y_uniform : ndarray (ny,) — uniform y-coordinates of the simulation grid
    N         : int — Chebyshev collocation points (default 100)

    Returns
    -------
    vort_eig : ndarray (ny,) real — vorticity eigenfunction, max-normalised
    sigma_nd : float — non-dimensional growth rate from the Chebyshev solve
    """
    # --- Chebyshev solve on N+1 nodes ---
    D_z, z = _cheb_diff_matrix(N)
    scale = 2.0 / Ly
    D1 = scale * D_z
    D2 = scale ** 2 * (D_z @ D_z)
    y_all = 0.5 * Ly * (z + 1.0)  # Chebyshev nodes on [0, Ly]

    U_all  = _tanh_profile(y_all, y_c, delta)
    d2U_all = _tanh_d2U(y_all, y_c, delta)

    interior = slice(1, N)
    n = N - 1
    D2_int = D2[interior, interior]
    U_int  = U_all[interior]
    d2U_int = d2U_all[interior]

    L = D2_int - k ** 2 * np.eye(n)
    A = np.diag(U_int) @ L - np.diag(d2U_int)
    B = L.copy()

    eigenvalues, eigenvectors = linalg.eig(A, B)
    mask = (np.imag(eigenvalues) > 1e-8) & (np.abs(np.real(eigenvalues)) < 1.5)
    if not np.any(mask):
        return np.zeros_like(y_uniform), 0.0

    eigs_f = eigenvalues[mask]
    vecs_f = eigenvectors[:, mask]
    idx = np.argmax(np.imag(eigs_f))
    c = eigs_f[idx]
    phi_int = vecs_f[:, idx]
    sigma_nd = float(np.real(k * np.imag(c) * delta))

    # Build full eigenfunction on N+1 nodes (BC: φ=0 at boundaries)
    phi_full = np.zeros(N + 1, dtype=complex)
    phi_full[interior] = phi_int

    # Compute vorticity eigenfunction ω = k²φ - φ'' using full D2
    phi_pp = D2 @ phi_full  # shape (N+1,)
    vort_cheb = k ** 2 * phi_full - phi_pp  # vorticity mode on Chebyshev nodes

    # Normalize: make real by dividing out the phase at the peak
    i_peak = np.argmax(np.abs(vort_cheb))
    phase = vort_cheb[i_peak] / abs(vort_cheb[i_peak])
    vort_real = np.real(vort_cheb / phase)

    # Interpolate from Chebyshev nodes (sorted ascending in y) to uniform grid
    # y_all is from Ly→0 (Chebyshev nodes are in reverse order), so we flip
    y_sorted = y_all[::-1]
    v_sorted = vort_real[::-1]
    vort_uniform = np.interp(y_uniform, y_sorted, v_sorted)

    # Normalise to max amplitude 1
    amp_max = np.max(np.abs(vort_uniform))
    if amp_max > 0:
        vort_uniform /= amp_max

    return vort_uniform, sigma_nd


if __name__ == "__main__":
    print("Rayleigh eigenvalue comparison — tanh shear layer vs Michalke (1964) Table 1")
    print("=" * 70)
    results = compare_methods()
    winner = print_comparison(results)
