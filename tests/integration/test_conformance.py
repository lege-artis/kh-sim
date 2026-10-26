# test_conformance.py -- request-honouring and robustness invariants, all 5 backends
# Task: KH-018 extension (2026-09-26; section 4 and the forcing checks 2026-10-03)
# Finding: development record, finding 2026-09-16 (v0.5)
#
# WHY THIS MODULE EXISTS
#   The reference parameters in conftest.REFERENCE_PARAMS (64x32, 100 steps,
#   Re 1000, ...) are exactly the Pascal backend's built-in defaults. A backend
#   that ignores its request body entirely therefore passes test_simulate.py:
#   it computes the defaults, which ARE the reference. That is how a stale Pascal
#   image that dropped every parameter stayed green from March to September 2026.
#
#   Every request below is deliberately OFF the defaults, and each parameter is
#   checked for effect one at a time. Thresholds were calibrated against the
#   canonical scheme on 2026-09-26 (see SENSITIVITY_FLOOR); they are not guesses.
#
#   Section 4 pins request contract KH-VAL-1 (finding v0.5): requests that no
#   backend can compute are refused by every backend as a 400 naming the field.
#   Before it, a JSON array body made the C++ and Fortran backends compute their
#   defaults and return 200 - the parameter-drop defect once more.
#
#   Section 5 pins contract KH-VAL-2 (development record, modelling rules v0.1,
#   ruling 2026-10-07): inputs that are computable but describe no state of the
#   model (K3) or no field on the grid (K4) are refused as 400 naming the field;
#   nu = |velocity_shear| / Re, so a reversed shear is the mirrored flow and
#   velocity_shear = 0 an inviscid run (D); every 200 carries "warnings" (K5).
#   Before it, velocity_shear = -1 ran with negative viscosity and was reported
#   as "dt beyond the stability bound", and velocity_shear = 0 was a 500 on
#   Pascal and a ZeroDivisionError in the oracle.
#
#   Section 6 pins forcing over HTTP (2026-10-10, ruling 8 and the shared CP-7F-B/D
#   harness): every backend computes forcing as the oracle does, refuses forced
#   modes it cannot represent (K4) and a band without a mode (K3), and continues a
#   run from t_start with the forcing phase intact. Before it the Fortran server
#   refused every forcing request, and a run split into two requests restarted the
#   rotating forcing at t = 0 (2.5e-4 away from the unsplit run in vorticity).
#
#   Section 5 also pins K7, step 2 (2026-10-24): a key the contract does not define
#   is refused as 400 naming it, before any value is read. Until then such a key was
#   logged and ignored; the census of the live logs (2026-10-09, 2026-10-24) found
#   none but the nsteps this module sent on purpose.
#
#   Section 7 pins the request size (2026-10-26): a request the contract allows is read
#   whole - 512^2 with initial_omega is a 5.4 MB body. The Rust backend read at most
#   2 MiB (axum's default) and answered 413 in plain text; the CP-7F-B2 run found it.

import json
import math

import pytest
import requests

from conftest import BACKEND_NAMES

TIMEOUT = 30

# Off every backend's defaults. Grid is a power of two, so it is valid for the
# radix-2 Pascal FFT as well as for general-length FFT libraries.
BASE = {
    "grid_nx": 32, "grid_ny": 16,
    "domain_lx": 1.0, "domain_ly": 0.5,
    "dt": 0.0008, "steps": 17,
    "reynolds_number": 800.0, "velocity_shear": 1.2,
    "perturbation_amplitude": 0.02, "perturbation_mode": 2,
}

# One-at-a-time perturbations. Measured relative change of kinetic_energy on the
# canonical scheme (Python reference == post-5c Pascal == current C++ source):
#   steps 7.1e-3, dt 7.5e-3, reynolds_number 1.5e-2, velocity_shear 4.4e-1,
#   perturbation_amplitude 2.9e-3, perturbation_mode 4.5e-6,
#   domain_lx 2.0e-6, domain_ly 8.6e-3.
PHYSICS_PERTURBATIONS = {
    "steps": 25,
    "dt": 0.0012,
    "reynolds_number": 400.0,
    "velocity_shear": 0.9,
    "perturbation_amplitude": 0.05,
    "perturbation_mode": 3,
    "domain_lx": 1.5,
    "domain_ly": 0.75,
}

# Three orders below the smallest real effect (2.0e-6, domain_lx) and six above
# floating-point rounding (~1e-15). A parameter that moves the result by less
# than this was not honoured.
SENSITIVITY_FLOOR = 1e-9

# Inputs a backend may compute or refuse, but must never be damaged by. A
# non-power-of-two grid is valid for general FFTs and invalid for radix-2; both
# answers are acceptable, a corrupted process is not.
HOSTILE_CASES = {
    "grid_nx_zero":           json.dumps(dict(BASE, grid_nx=0)),
    "grid_nx_negative":       json.dumps(dict(BASE, grid_nx=-8)),
    "grid_nx_three":          json.dumps(dict(BASE, grid_nx=3)),
    "grid_nx_48_non_pow2":    json.dumps(dict(BASE, grid_nx=48)),
    "grid_ny_24_non_pow2":    json.dumps(dict(BASE, grid_ny=24)),
    "domain_lx_zero":         json.dumps(dict(BASE, domain_lx=0.0)),
    "reynolds_number_zero":   json.dumps(dict(BASE, reynolds_number=0.0)),
    "body_not_json":          "this is not json",
    "body_json_array":        "[1, 2, 3]",
}

PROBE = {"grid_nx": 32, "grid_ny": 16, "steps": 10}

# Requests no backend can compute, with the field the 400 must name. Pinned by
# contract KH-VAL-1; every backend implements it (finding v0.5).
INVALID_CASES = {
    "grid_nx_zero":               (dict(BASE, grid_nx=0), "grid_nx"),
    "grid_nx_negative":           (dict(BASE, grid_nx=-8), "grid_nx"),
    "grid_ny_above_axis_cap":     (dict(BASE, grid_ny=131072), "grid_ny"),
    "grid_nx_not_a_number":       (dict(BASE, grid_nx="32"), "grid_nx"),
    "domain_lx_zero":             (dict(BASE, domain_lx=0.0), "domain_lx"),
    "domain_ly_zero":             (dict(BASE, domain_ly=0.0), "domain_ly"),
    "reynolds_number_zero":       (dict(BASE, reynolds_number=0.0), "reynolds_number"),
    "initial_omega_wrong_length": (dict(BASE, initial_omega=[0.0] * 10), "initial_omega"),
    "initial_omega_not_an_array": (dict(BASE, initial_omega="zeros"), "initial_omega"),
    "both_ic_keys":               (dict(BASE, initial_omega=[0.0] * 512, init_omega=[0.0] * 512),
                                   "initial_omega"),
    "unknown_forcing_type":       (dict(BASE, forcing_type="bogus", forcing_amplitude=1.0),
                                   "forcing_type"),
    "body_not_json":              ("this is not json", "body"),
    "body_json_array":            ("[1, 2, 3]", "body"),
    "body_empty":                 ("", "body"),
}

# A dt far beyond the stability bound: the run blows up. A backend must say so
# (500), or return numbers that are all finite - never 200 with nulls.
UNSTABLE = dict(BASE, dt=1.0, steps=5)

# Forcing on a square domain (band forcing requires domain_lx == domain_ly).
# Measured relative change of kinetic_energy against the unforced run on the
# canonical scheme (C++ == Rust, 2026-10-03): mode 1.6e-7, band 7.5e-5,
# forcing_omega0 6.9e-5 (against the unrotated band), linear_drag 1.4e-2.
SQUARE = dict(BASE, domain_ly=1.0)
BAND = dict(SQUARE, forcing_type="band", forcing_amplitude=0.5, forcing_kmin=1.0, forcing_kmax=3.0)
FORCING_CASES = {
    # case: (reference request, perturbed request, field a refusing backend named before
    # 2026-10-10; since then every backend computes forcing and none may refuse)
    "forcing_mode":   (SQUARE, dict(SQUARE, forcing_type="mode", forcing_amplitude=0.5,
                                    forcing_kx=2, forcing_ky=1), "forcing_type"),
    "forcing_band":   (SQUARE, BAND, "forcing_type"),
    "forcing_omega0": (BAND, dict(BAND, forcing_omega0=5.0), "forcing_type"),
    "linear_drag":    (SQUARE, dict(SQUARE, linear_drag=0.5), "linear_drag"),
}


# -- helpers -------------------------------------------------------------------

def _post(url: str, payload, timeout: int = TIMEOUT) -> requests.Response:
    data = payload if isinstance(payload, str) else json.dumps(payload)
    return requests.post(f"{url}/simulate", data=data,
                         headers={"Content-Type": "application/json"}, timeout=timeout)


# Backends that answered at least once in this session. A backend that was
# reachable and later is not has been damaged by something this module sent;
# that must FAIL, never skip. (First version of this module skipped here, and a
# wedged Pascal backend reported as "1 failed, 1 skipped" - the one outcome the
# module exists to catch, reported as a skip.)
_SEEN_ALIVE: dict = {}


def _require(backend_name: str, backend_urls: dict) -> str:
    url = backend_urls[backend_name]
    try:
        requests.get(f"{url}/health", timeout=5)
    except requests.exceptions.RequestException as exc:
        if _SEEN_ALIVE.get(backend_name):
            pytest.fail(f"{backend_name} answered earlier in this session and is now unreachable "
                        f"({type(exc).__name__}) - damaged by a preceding request")
        pytest.skip(f"{backend_name} not reachable at {url}")
    _SEEN_ALIVE[backend_name] = True
    return url


def _ok_body(backend_name: str, resp: requests.Response, what: str) -> dict:
    assert resp.status_code == 200, (
        f"{backend_name}: {what} returned {resp.status_code}: {resp.text[:200]}")
    return resp.json()


def _ke(body: dict) -> float:
    return body["diagnostics"]["kinetic_energy"]


# -- 1. The response states the configuration it computed ----------------------

@pytest.mark.parametrize("backend_name", BACKEND_NAMES)
def test_conformance_echo_of_non_default_request(backend_name: str, backend_urls: dict) -> None:
    """grid, steps_completed and t_final describe THIS request, not the defaults."""
    url = _require(backend_name, backend_urls)
    body = _ok_body(backend_name, _post(url, BASE), "POST /simulate")
    problems = []
    if body.get("grid_nx") != BASE["grid_nx"]:
        problems.append(f"grid_nx {body.get('grid_nx')} != {BASE['grid_nx']}")
    if body.get("grid_ny") != BASE["grid_ny"]:
        problems.append(f"grid_ny {body.get('grid_ny')} != {BASE['grid_ny']}")
    if body.get("steps_completed") != BASE["steps"]:
        problems.append(f"steps_completed {body.get('steps_completed')} != {BASE['steps']}")
    t_final = body.get("t_final")
    expected_t = BASE["steps"] * BASE["dt"]
    if not isinstance(t_final, (int, float)) or abs(t_final - expected_t) > 1e-9:
        problems.append(f"t_final {t_final} != {expected_t}")
    cells = BASE["grid_nx"] * BASE["grid_ny"]
    for f in ("u_velocity", "v_velocity", "vorticity", "pressure"):
        if len(body.get(f, [])) != cells:
            problems.append(f"{f} length {len(body.get(f, []))} != {cells}")
    assert not problems, f"{backend_name}: response does not describe the request: {problems}"


# -- 2. Every parameter is honoured --------------------------------------------

@pytest.mark.parametrize("param", sorted(PHYSICS_PERTURBATIONS))
@pytest.mark.parametrize("backend_name", BACKEND_NAMES)
def test_conformance_parameter_is_honoured(backend_name: str, param: str, backend_urls: dict) -> None:
    """Changing one physics parameter changes kinetic_energy by more than rounding."""
    url = _require(backend_name, backend_urls)
    ke0 = _ke(_ok_body(backend_name, _post(url, BASE), "baseline"))
    ke1 = _ke(_ok_body(backend_name, _post(url, dict(BASE, **{param: PHYSICS_PERTURBATIONS[param]})),
                       f"{param} perturbed"))
    rel = abs(ke1 - ke0) / abs(ke0)
    assert rel > SENSITIVITY_FLOOR, (
        f"{backend_name}: {param} {BASE[param]} -> {PHYSICS_PERTURBATIONS[param]} moved "
        f"kinetic_energy by {rel:.2e} (floor {SENSITIVITY_FLOOR:.0e}); the parameter was not honoured")


@pytest.mark.parametrize("axis", ["grid_nx", "grid_ny"])
@pytest.mark.parametrize("backend_name", BACKEND_NAMES)
def test_conformance_grid_is_honoured(backend_name: str, axis: str, backend_urls: dict) -> None:
    """A grid change is detected by field length: at these grids a resolved flow's
    kinetic_energy may move by as little as 1e-15 (grid_nx 32->64 measured)."""
    url = _require(backend_name, backend_urls)
    req = dict(BASE, **{axis: BASE[axis] * 2})
    body = _ok_body(backend_name, _post(url, req), f"{axis} doubled")
    cells = req["grid_nx"] * req["grid_ny"]
    assert len(body.get("u_velocity", [])) == cells, (
        f"{backend_name}: {axis} doubled but u_velocity has {len(body.get('u_velocity', []))} "
        f"entries, expected {cells}")


@pytest.mark.parametrize("backend_name", BACKEND_NAMES)
def test_conformance_initial_omega_is_honoured(backend_name: str, backend_urls: dict) -> None:
    """initial_omega (the canonical key, shared/physics/kh_physics.py) replaces the
    built-in IC. Zero vorticity stays zero, so kinetic_energy must be exactly 0;
    a backend that ignores the key returns the tanh-shear value (~0.18)."""
    url = _require(backend_name, backend_urls)
    cells = BASE["grid_nx"] * BASE["grid_ny"]
    # perturbation_amplitude 0: initial_omega replaces the analytic field, and a
    # non-zero amplitude sent with it is refused (K4b, ruling 2026-10-09).
    body = _ok_body(backend_name, _post(url, dict(BASE, perturbation_amplitude=0.0,
                                                  initial_omega=[0.0] * cells)),
                    "initial_omega = zeros")
    ke = _ke(body)
    assert abs(ke) < 1e-20, (
        f"{backend_name}: initial_omega of zeros gave kinetic_energy {ke!r}; "
        f"the canonical key 'initial_omega' was ignored")


# -- 3. Hostile input: an honest answer, and no damage --------------------------

@pytest.mark.parametrize("case", sorted(HOSTILE_CASES))
@pytest.mark.parametrize("backend_name", BACKEND_NAMES)
def test_conformance_hostile_input_error_envelope(backend_name: str, case: str, backend_urls: dict) -> None:
    """200, or an error body in the shape openapi.yaml declares:
    400 -> {error, field}; 500 -> {error, ...}. Any other status is a defect."""
    url = _require(backend_name, backend_urls)
    resp = _post(url, HOSTILE_CASES[case], timeout=15)
    assert resp.status_code in (200, 400, 500), (
        f"{backend_name}/{case}: status {resp.status_code}: {resp.text[:200]}")
    if resp.status_code == 200:
        ke = _ke(resp.json())
        assert isinstance(ke, (int, float)) and math.isfinite(ke), (
            f"{backend_name}/{case}: 200 with kinetic_energy {ke!r}")
        return
    try:
        body = resp.json()
    except ValueError:
        pytest.fail(f"{backend_name}/{case}: {resp.status_code} body is not JSON: {resp.text[:200]}")
    assert isinstance(body.get("error"), str) and body["error"], (
        f"{backend_name}/{case}: {resp.status_code} without an 'error' message: {body}")
    if resp.status_code == 400:
        assert isinstance(body.get("field"), str) and body["field"], (
            f"{backend_name}/{case}: 400 without a 'field' naming the rejected input: {body}")


@pytest.mark.parametrize("case", ["grid_nx_48_non_pow2", "grid_ny_24_non_pow2"])
@pytest.mark.parametrize("backend_name", BACKEND_NAMES)
def test_conformance_unsupported_grid_is_a_request_error(backend_name: str, case: str, backend_urls: dict) -> None:
    """A grid length the backend's FFT cannot transform is the caller's error.
    200 (a general-length FFT computes it) or 400 (a radix-2 FFT refuses it) are
    both correct; 500 says the backend tried and failed, which is wrong. This is
    the lock for an HTTP-side grid check: without one, a radix-2 backend's
    internal guard still stops the grid, but as a 500."""
    url = _require(backend_name, backend_urls)
    resp = _post(url, HOSTILE_CASES[case], timeout=15)
    assert resp.status_code in (200, 400), (
        f"{backend_name}/{case}: status {resp.status_code} - an unsupported grid must be "
        f"computed (200) or refused as a request error (400): {resp.text[:200]}")
    if resp.status_code == 200:
        req = json.loads(HOSTILE_CASES[case])
        cells = req["grid_nx"] * req["grid_ny"]
        assert len(resp.json().get("u_velocity", [])) == cells, (
            f"{backend_name}/{case}: computed, but u_velocity has the wrong length")


@pytest.mark.parametrize("case", sorted(HOSTILE_CASES))
@pytest.mark.parametrize("backend_name", BACKEND_NAMES)
def test_conformance_hostile_input_does_no_damage(backend_name: str, case: str, backend_urls: dict) -> None:
    """After any hostile request the backend still answers /health, and a probe
    simulation returns the same result it returned before. Until 2026-09-26 a
    single grid_nx=48 left the Pascal process alive but answering nothing."""
    url = _require(backend_name, backend_urls)
    before = _ke(_ok_body(backend_name, _post(url, PROBE), "probe before"))
    again = _ke(_ok_body(backend_name, _post(url, PROBE), "probe repeated"))
    deterministic = before == again

    try:
        _post(url, HOSTILE_CASES[case], timeout=15)
    except requests.exceptions.RequestException:
        pass  # the verdict is on what follows, not on this reply

    try:
        health = requests.get(f"{url}/health", timeout=5).status_code
    except requests.exceptions.RequestException as exc:
        pytest.fail(f"{backend_name}/{case}: /health unreachable after hostile input ({type(exc).__name__})")
    assert health == 200, f"{backend_name}/{case}: /health returned {health} after hostile input"

    try:
        after = _ke(_ok_body(backend_name, _post(url, PROBE, timeout=15), "probe after"))
    except requests.exceptions.RequestException as exc:
        pytest.fail(f"{backend_name}/{case}: probe simulation failed after hostile input ({type(exc).__name__})")
    if deterministic:
        assert after == before, (
            f"{backend_name}/{case}: probe changed from {before!r} to {after!r} after hostile input")
    else:
        rel = abs(after - before) / abs(before)
        assert rel < 1e-12, (
            f"{backend_name}/{case}: probe drifted by {rel:.2e} after hostile input "
            f"(backend is not bit-reproducible back-to-back; tolerance 1e-12 applied)")


# -- 4. Contract KH-VAL-1: what every backend refuses, and how -------------------

@pytest.mark.parametrize("case", sorted(INVALID_CASES))
@pytest.mark.parametrize("backend_name", BACKEND_NAMES)
def test_conformance_invalid_request_is_400(backend_name: str, case: str, backend_urls: dict) -> None:
    """A request no backend can compute is the caller's error: 400 with the field
    that is wrong. Computing something else instead (the defaults, the built-in
    initial condition, an unforced flow) is the parameter-drop defect."""
    url = _require(backend_name, backend_urls)
    payload, field = INVALID_CASES[case]
    resp = _post(url, payload, timeout=15)
    assert resp.status_code == 400, (
        f"{backend_name}/{case}: expected 400, got {resp.status_code}: {resp.text[:200]}")
    try:
        body = resp.json()
    except ValueError:
        pytest.fail(f"{backend_name}/{case}: 400 body is not JSON: {resp.text[:200]}")
    assert isinstance(body.get("error"), str) and body["error"], (
        f"{backend_name}/{case}: 400 without an 'error' message: {body}")
    assert body.get("field") == field, (
        f"{backend_name}/{case}: field {body.get('field')!r}, expected {field!r}: {body}")


@pytest.mark.parametrize("backend_name", BACKEND_NAMES)
def test_conformance_unstable_run_is_reported(backend_name: str, backend_urls: dict) -> None:
    """A run that blows up is the backend's failure (500 naming the backend), or a
    200 whose numbers are all finite. A 200 carrying nulls hides the failure."""
    url = _require(backend_name, backend_urls)
    resp = _post(url, UNSTABLE)
    assert resp.status_code in (200, 500), (
        f"{backend_name}: unstable run returned {resp.status_code}: {resp.text[:200]}")
    body = resp.json()
    if resp.status_code == 500:
        assert isinstance(body.get("error"), str) and body["error"], body
        assert body.get("backend") == backend_name, (
            f"{backend_name}: a 500 must name the backend: {body}")
        return
    numbers = list(body["diagnostics"].values())
    for f in ("u_velocity", "v_velocity", "vorticity", "pressure"):
        numbers.extend(body.get(f, []))
    bad = [x for x in numbers if not isinstance(x, (int, float)) or not math.isfinite(x)]
    assert not bad, f"{backend_name}: 200 with {len(bad)} non-finite values (first {bad[0]!r})"


@pytest.mark.parametrize("case", sorted(FORCING_CASES))
@pytest.mark.parametrize("backend_name", BACKEND_NAMES)
def test_conformance_forcing_is_honoured(backend_name: str, case: str, backend_urls: dict) -> None:
    """Every backend applies forcing. Until 2026-10-03 the Fortran server dropped
    every forcing parameter and returned the unforced flow; from then until
    2026-10-10 it refused them (400); since ruling 8 it computes them."""
    url = _require(backend_name, backend_urls)
    reference, perturbed, _field = FORCING_CASES[case]
    resp = _post(url, perturbed)
    ke1 = _ke(_ok_body(backend_name, resp, f"{case} request"))
    ke0 = _ke(_ok_body(backend_name, _post(url, reference), f"{case} reference"))
    rel = abs(ke1 - ke0) / abs(ke0)
    assert rel > SENSITIVITY_FLOOR, (
        f"{backend_name}: {case} moved kinetic_energy by {rel:.2e} "
        f"(floor {SENSITIVITY_FLOOR:.0e}); the forcing parameter was accepted and not applied")


# -- 5. Contract KH-VAL-2: the modelling rules (ruling 2026-10-07) ---------------

import pathlib  # noqa: E402
import sys  # noqa: E402

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2] / "shared" / "physics"))
import kh_physics  # noqa: E402

# K3 outside the model and K4 not representable: 400 naming the field.
MODEL_INVALID_CASES = {
    "domain_lx_negative":          (dict(BASE, domain_lx=-1.0), "domain_lx"),
    "domain_ly_negative":          (dict(BASE, domain_ly=-0.5), "domain_ly"),
    "reynolds_number_negative":    (dict(BASE, reynolds_number=-800.0), "reynolds_number"),
    "dt_negative":                 (dict(BASE, dt=-0.0008), "dt"),
    "steps_negative":              (dict(BASE, steps=-1), "steps"),
    "linear_drag_negative":        (dict(BASE, linear_drag=-0.5), "linear_drag"),
    "perturbation_mode_negative":  (dict(BASE, perturbation_mode=-2), "perturbation_mode"),
    "perturbation_mode_zero":      (dict(BASE, perturbation_mode=0), "perturbation_mode"),
    "perturbation_mode_beyond_2_3": (dict(BASE, perturbation_mode=11), "perturbation_mode"),
    "forcing_band_empty":          (dict(BAND, forcing_kmin=3.0, forcing_kmax=1.0), "forcing_kmin"),
    # K4/K3 forcing (2026-10-10). Grid 32x16: the 2/3-rule band ends at 10 along x and
    # 5 along y; a band may reach 5 (min(nx, ny)). 1.1 <= |k| <= 1.3 holds no integer mode.
    "forcing_kx_beyond_2_3":       (dict(SQUARE, forcing_type="mode", forcing_amplitude=0.5,
                                         forcing_kx=11, forcing_ky=1), "forcing_kx"),
    "forcing_ky_beyond_2_3":       (dict(SQUARE, forcing_type="mode", forcing_amplitude=0.5,
                                         forcing_kx=1, forcing_ky=6), "forcing_ky"),
    "forcing_band_beyond_2_3":     (dict(BAND, forcing_kmax=6.0), "forcing_kmax"),
    "forcing_band_without_mode":   (dict(BAND, forcing_kmin=1.1, forcing_kmax=1.3), "forcing_kmax"),
    # K4b (ruling 2026-10-09): BASE sends perturbation_amplitude 0.02, which initial_omega
    # would silently override.
    "amplitude_with_initial_omega": (dict(BASE, initial_omega=[0.0] * 512), "perturbation_amplitude"),
}
# Per-backend field overrides (none since 2026-10-10: the Fortran server refused
# every forcing request first, naming forcing_type, until it computed forcing).
MODEL_INVALID_FIELD = {}

# Agreement with the oracle on the same request: the deploy check of 2026-10-06
# measured ~1e-13 for valid requests; a rule implemented differently moves the
# kinetic energy by orders of magnitude more.
ORACLE_REL = 1e-10

WARNING_PREFIXES = ("cfl_advective:", "cfl_diffusive:", "inviscid:")


def _oracle_ke(payload: dict) -> float:
    req = kh_physics.SimulationRequest(**payload)
    import warnings as _w
    with _w.catch_warnings():
        _w.simplefilter("ignore")
        return kh_physics.simulate(req).diagnostics.kinetic_energy


def _warnings(backend_name: str, body: dict) -> list:
    w = body.get("warnings")
    assert isinstance(w, list) and all(isinstance(x, str) for x in w), (
        f"{backend_name}: 'warnings' must be a list of strings, got {w!r}")
    assert all(x.startswith(WARNING_PREFIXES) for x in w), (
        f"{backend_name}: unknown warning prefix in {w}")
    return [x.split(":")[0] for x in w]


@pytest.mark.parametrize("case", sorted(MODEL_INVALID_CASES))
@pytest.mark.parametrize("backend_name", BACKEND_NAMES)
def test_conformance_outside_the_model_is_400(backend_name: str, case: str, backend_urls: dict) -> None:
    """K3/K4: computable, but no state of the model (negative viscosity, backward
    diffusion, an empty band) or no field on the grid (a mode that aliases or lies
    outside the 2/3-rule band). K4b: an amplitude the initial field would ignore.
    Refused, naming the field - never computed."""
    url = _require(backend_name, backend_urls)
    payload, field = MODEL_INVALID_CASES[case]
    field = MODEL_INVALID_FIELD.get((backend_name, case), field)
    resp = _post(url, payload, timeout=15)
    assert resp.status_code == 400, (
        f"{backend_name}/{case}: expected 400, got {resp.status_code}: {resp.text[:200]}")
    assert resp.json().get("field") == field, (
        f"{backend_name}/{case}: field {resp.json().get('field')!r}, expected {field!r}: {resp.text[:200]}")


@pytest.mark.parametrize("backend_name", BACKEND_NAMES)
def test_conformance_initial_omega_without_amplitude_is_computed(backend_name: str, backend_urls: dict) -> None:
    """K4b refuses only an amplitude that was SENT: the default amplitude (0.01) is
    not the client's statement, so initial_omega alone is computed."""
    url = _require(backend_name, backend_urls)
    payload = {k: v for k, v in BASE.items() if k != "perturbation_amplitude"}
    payload["initial_omega"] = [0.0] * (BASE["grid_nx"] * BASE["grid_ny"])
    body = _ok_body(backend_name, _post(url, payload), "initial_omega without perturbation_amplitude")
    assert abs(_ke(body)) < 1e-20, f"{backend_name}: zero initial_omega gave kinetic_energy {_ke(body)!r}"


@pytest.mark.parametrize("backend_name", BACKEND_NAMES)
def test_conformance_mode_at_the_band_edge_is_computed(backend_name: str, backend_urls: dict) -> None:
    """The largest resolved mode (2*(grid_nx/2)/3 = 10 at grid_nx 32) is computed,
    and equals the oracle: the K4 bound is not off by one."""
    url = _require(backend_name, backend_urls)
    payload = dict(BASE, perturbation_mode=10)
    body = _ok_body(backend_name, _post(url, payload), "mode 10 at grid_nx 32")
    rel = abs(_ke(body) - _oracle_ke(payload)) / abs(_oracle_ke(payload))
    assert rel < ORACLE_REL, f"{backend_name}: differs from the oracle by {rel:.2e}"


@pytest.mark.parametrize("shear", [0.0, -1.2])
@pytest.mark.parametrize("backend_name", BACKEND_NAMES)
def test_conformance_viscosity_is_abs_shear_over_re(backend_name: str, shear: float,
                                                     backend_urls: dict) -> None:
    """nu = |velocity_shear| / Re (D and K3). velocity_shear = 0 is an inviscid run
    (was: ZeroDivisionError in the oracle, 500 on Pascal); a negative shear is the
    mirrored flow (was: negative viscosity, and a blow-up reported as a dt problem).
    Run at Re 100 for 400 steps, where negative viscosity diverged."""
    url = _require(backend_name, backend_urls)
    payload = dict(BASE, velocity_shear=shear, reynolds_number=100.0, steps=400)
    body = _ok_body(backend_name, _post(url, payload), f"velocity_shear={shear}")
    ref = _oracle_ke(payload)
    rel = abs(_ke(body) - ref) / abs(ref)
    assert rel < ORACLE_REL, f"{backend_name}: velocity_shear={shear} differs from the oracle by {rel:.2e}"
    w = _warnings(backend_name, body)
    assert ("inviscid" in w) == (shear == 0.0), f"{backend_name}: warnings {body['warnings']}"


WARNING_CASES = {
    # case: (request, expected warning kinds, in order)
    "none":          (BASE, []),
    "cfl_advective": (dict(BASE, dt=0.05, steps=2), ["cfl_advective"]),
    "cfl_diffusive": (dict(BASE, reynolds_number=1.0, steps=1), ["cfl_diffusive"]),
    "inviscid":      (dict(BASE, velocity_shear=0.0), ["inviscid"]),
}


@pytest.mark.parametrize("case", sorted(WARNING_CASES))
@pytest.mark.parametrize("backend_name", BACKEND_NAMES)
def test_conformance_warnings(backend_name: str, case: str, backend_urls: dict) -> None:
    """K5: a computed request that may not be resolved says so in "warnings" -
    the same kinds, in the same order, as the oracle; an ordinary request has none."""
    url = _require(backend_name, backend_urls)
    payload, expected = WARNING_CASES[case]
    body = _ok_body(backend_name, _post(url, payload), case)
    assert _warnings(backend_name, body) == expected, f"{backend_name}/{case}: {body['warnings']}"
    req = kh_physics.SimulationRequest(**payload)
    import warnings as _w
    with _w.catch_warnings():
        _w.simplefilter("ignore")
        oracle = [x.split(":")[0] for x in kh_physics.simulate(req).warnings]
    assert oracle == expected, f"oracle/{case}: {oracle}"


# K7 (step 2, 2026-10-24): a key the contract does not define. Before, nsteps was
# ignored and the run took the default steps - a plausible wrong answer (R1).
UNKNOWN_KEY_CASES = {
    # case: (extra keys, field the 400 names, the hint the message must carry or None)
    "foreign_name_nsteps": ({"nsteps": 5}, "nsteps", "steps"),
    "no_foreign_name":     ({"frobnicate": 1}, "frobnicate", None),
}


@pytest.mark.parametrize("case", sorted(UNKNOWN_KEY_CASES))
@pytest.mark.parametrize("backend_name", BACKEND_NAMES)
def test_conformance_unknown_key_is_400(backend_name: str, case: str, backend_urls: dict) -> None:
    """K7: an unknown key is refused, 400 naming it; a foreign name with a contract
    equivalent (nsteps -> steps) says which key was meant."""
    url = _require(backend_name, backend_urls)
    extra, field, hint = UNKNOWN_KEY_CASES[case]
    resp = _post(url, dict(BASE, **extra), timeout=15)
    assert resp.status_code == 400, (
        f"{backend_name}/{case}: expected 400, got {resp.status_code}: {resp.text[:200]}")
    body = resp.json()
    assert body.get("field") == field, f"{backend_name}/{case}: {resp.text[:200]}"
    msg = body.get("error", "")
    if hint:
        assert f"did you mean '{hint}'" in msg, f"{backend_name}/{case}: no hint in {msg!r}"
    else:
        assert "did you mean" not in msg, f"{backend_name}/{case}: unexpected hint in {msg!r}"


@pytest.mark.parametrize("backend_name", BACKEND_NAMES)
def test_conformance_unknown_key_is_checked_first(backend_name: str, backend_urls: dict) -> None:
    """K7 is recognised by the parser, before any value: with two unknown keys and an
    invalid grid_nx, the 400 names the first unknown key in byte order - the same
    answer from every backend."""
    url = _require(backend_name, backend_urls)
    resp = _post(url, dict(BASE, grid_nx=0, zz_unknown=1, aa_unknown=2), timeout=15)
    assert resp.status_code == 400, f"{backend_name}: {resp.status_code} {resp.text[:200]}"
    assert resp.json().get("field") == "aa_unknown", f"{backend_name}: {resp.text[:200]}"


# -- 6. Forcing over HTTP (ruling 8, 2026-10-10) -----------------------------------

# Forced runs equal the oracle on the same request (measured 2026-10-10, all four
# local backends: 1e-16 .. 7e-16 in kinetic energy). The mode case sits on the
# 2/3-rule edge (kx 10 of 32, ky 5 of 16), so the K4 bound is not off by one.
FORCING_ORACLE_CASES = {
    "mode_at_band_edge": dict(SQUARE, forcing_type="mode", forcing_amplitude=0.5,
                              forcing_kx=10, forcing_ky=5),
    "band_static_drag":  dict(BAND, linear_drag=0.3),
    "band_rotating":     dict(BAND, forcing_omega0=5.0, linear_drag=0.3),
    "band_rotating_t_start": dict(BAND, forcing_omega0=5.0, linear_drag=0.3, t_start=0.37),
    # A band holding one shell, |k| = sqrt(26) = 5.0990195136 (modes (1,5), (5,1)):
    # 5.09901951 <= |k| <= 5.09901952. A validator taking sqrt in single precision
    # (5.0990195274) finds no mode and answers 400 while the builder, in double,
    # forces the shell - the Pascal server did until 2026-10-11 (FPC types the
    # literal 1.0 in Sqrt(1.0 * n) as Single).
    "band_single_shell_sqrt26": dict(BAND, forcing_kmin=5.09901951, forcing_kmax=5.09901952),
}


@pytest.mark.parametrize("case", sorted(FORCING_ORACLE_CASES))
@pytest.mark.parametrize("backend_name", BACKEND_NAMES)
def test_conformance_forcing_equals_the_oracle(backend_name: str, case: str, backend_urls: dict) -> None:
    """Forcing is computed as the oracle computes it - mode, static band with drag,
    rotating band, and a rotating band started at t_start (its phase moves the
    result by 3e-4 in kinetic energy); t_final = t_start + steps*dt."""
    url = _require(backend_name, backend_urls)
    payload = FORCING_ORACLE_CASES[case]
    body = _ok_body(backend_name, _post(url, payload), case)
    ref = _oracle_ke(payload)
    rel = abs(_ke(body) - ref) / abs(ref)
    assert rel < ORACLE_REL, f"{backend_name}/{case}: differs from the oracle by {rel:.2e}"
    want_t = payload.get("t_start", 0.0) + payload["steps"] * payload["dt"]
    assert abs(body["t_final"] - want_t) < 1e-9, f"{backend_name}/{case}: t_final {body['t_final']} != {want_t}"


@pytest.mark.parametrize("backend_name", BACKEND_NAMES)
def test_conformance_continued_run_equals_the_whole_run(backend_name: str, backend_urls: dict) -> None:
    """A rotating-band run of 2n steps equals n steps continued by n more steps from
    the returned field with t_start = the returned t_final. Without t_start the
    second request restarts the forcing phase at 0 (measured 2.5e-4 in vorticity)."""
    url = _require(backend_name, backend_urls)
    whole = dict(BAND, forcing_omega0=5.0, linear_drag=0.3, steps=2 * BASE["steps"])
    full = _ok_body(backend_name, _post(url, whole), "whole run")
    half = _ok_body(backend_name, _post(url, dict(whole, steps=BASE["steps"])), "first half")
    rest = {k: v for k, v in whole.items() if k != "perturbation_amplitude"}
    rest.update(steps=BASE["steps"], initial_omega=half["vorticity"], t_start=half["t_final"])
    cont = _ok_body(backend_name, _post(url, rest), "second half")
    peak = max(abs(x) for x in full["vorticity"])
    dev = max(abs(a - b) for a, b in zip(cont["vorticity"], full["vorticity"])) / peak
    assert dev < 1e-12, f"{backend_name}: continued run differs from the whole run by {dev:.2e}"
    assert abs(cont["t_final"] - full["t_final"]) < 1e-9, (
        f"{backend_name}: t_final {cont['t_final']} after continuation, {full['t_final']} whole")


# -- 7. Request size (2026-10-26) ---------------------------------------------------

# A 512^2 request carrying initial_omega is a body of ~5.4 MB. The contract allows
# grids up to 4194304 cells (K8, openapi), but the Rust backend read at most 2 MiB
# (axum's default body limit) and answered 413 in plain text before the contract was
# consulted - the CP-7F-B2 run of 2026-10-25 found it (its first 512^2 request).
LARGE_N = 512


@pytest.mark.parametrize("backend_name", BACKEND_NAMES)
def test_conformance_large_body_is_computed(backend_name: str, backend_urls: dict) -> None:
    """A request at 512^2 with initial_omega (a body above 2 MiB) is computed, as the
    oracle computes it."""
    url = _require(backend_name, backend_urls)
    x = [2.0 * math.pi * i / LARGE_N for i in range(LARGE_N)]
    omega = [math.sin(3.0 * a + 0.1) * math.cos(5.0 * b) + 0.25 * math.cos(7.0 * a - 2.0 * b)
             for a in x for b in x]
    payload = {"grid_nx": LARGE_N, "grid_ny": LARGE_N, "domain_lx": 2.0 * math.pi,
               "domain_ly": 2.0 * math.pi, "dt": 0.001, "steps": 1, "reynolds_number": 1000.0,
               "velocity_shear": 1.0, "perturbation_amplitude": 0.0, "perturbation_mode": 1,
               "initial_omega": omega}
    data = json.dumps(payload)
    assert len(data) > 2 * 1024 * 1024, len(data)
    body = _ok_body(backend_name, _post(url, data, timeout=300), "512^2 with initial_omega")
    ref = _oracle_ke(payload)
    rel = abs(_ke(body) - ref) / abs(ref)
    assert rel < ORACLE_REL, f"{backend_name}: differs from the oracle by {rel:.2e}"
