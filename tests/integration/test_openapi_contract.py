# test_openapi_contract.py -- shared/api/openapi.yaml is what the backends implement
# Task: openapi schema pass, 2026-10-09 (modelling rules s7, open item 2)
#
# WHY THIS MODULE EXISTS
#   Until 2026-10-09 openapi.yaml described fields no backend had (density_ratio,
#   initial_state, final_state, log_id, uptime_s, db_connected), nested the result
#   fields under a final_state no backend sends, and omitted initial_omega, the
#   forcing keys, cfl_peak and the /info fields - and nothing failed, because nothing
#   read it. This module makes the document a tested contract:
#
#   1. the request properties ARE the oracle's SimulationRequest fields (plus the
#      legacy alias init_omega), with the same defaults;
#   2. they ARE the key list every backend knows (KH-VAL-2 K7) - read from the five
#      sources, so a key added to one backend and not to the document fails here;
#   3. every backend's /simulate 200, /health, /info and 400 bodies validate against
#      the schemas;
#   4. a request that omits every key gets the documented defaults: 128x64, 100 steps,
#      t_final 0.1, and the oracle's kinetic energy for those defaults.
#
#   The first run of (4) found the Pascal backend computing 64x32 for an empty request
#   while the oracle and the other four compute 128x64.

import dataclasses
import json
import os
import pathlib
import re
import sys
import warnings

import jsonschema
import pytest
import requests
import yaml

from conftest import BACKEND_NAMES

_KH = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_KH / "shared" / "physics"))
import kh_physics  # noqa: E402

SPEC = yaml.safe_load((_KH / "shared" / "api" / "openapi.yaml").read_text(encoding="utf-8"))
SCHEMAS = SPEC["components"]["schemas"]
REQUEST_PROPS = SCHEMAS["SimulationRequest"]["properties"]
LEGACY_ALIASES = {"init_omega"}
TIMEOUT = 60
ORACLE_REL = 1e-10

# Where each backend keeps its K7 known-key list.
KNOWN_KEY_SOURCES = {
    "cpp":     "backends/cpp/src/kh_http_contract.hpp",
    "fortran": "backends/fortran/src/kh_http_contract.hpp",
    "rust":    "backends/rust/src/contract.rs",
    "scala":   "backends/scala/src/main/scala/khsim/Contract.scala",
    "pascal":  "backends/pascal/src/kh_server.pas",
}
_KNOWN_BLOCK = re.compile(r"\b(?:KNOWN|Known|known)\b[^=(\n]*[=(]\s*(?:Set\()?\s*[\[{(]?(.*?)[\]})]\s*;?\s*\n", re.S)


def _validate(instance, ref: str) -> None:
    """Validate against a component of openapi.yaml (OpenAPI 3.1 = JSON Schema 2020-12)."""
    schema = {"$ref": ref, "components": SPEC["components"]}
    jsonschema.Draft202012Validator(schema).validate(instance)


def _url(backend_name: str, backend_urls: dict) -> str:
    url = backend_urls[backend_name]
    try:
        requests.get(f"{url}/health", timeout=5)
    except requests.exceptions.RequestException:
        if os.environ.get("KH_GOLDEN_REQUIRE_SERVER") == "1":
            pytest.fail(f"{backend_name} not reachable at {url}")
        pytest.skip(f"{backend_name} not reachable at {url}")
    return url


_ORACLE_DEFAULT_KE = {}


def _oracle_default_ke() -> float:
    if "ke" not in _ORACLE_DEFAULT_KE:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            _ORACLE_DEFAULT_KE["ke"] = kh_physics.simulate(kh_physics.SimulationRequest()).diagnostics.kinetic_energy
    return _ORACLE_DEFAULT_KE["ke"]


# -- 1, 2: the document against the code ------------------------------------------

def test_request_properties_are_the_oracle_fields_with_its_defaults() -> None:
    fields = {f.name: f for f in dataclasses.fields(kh_physics.SimulationRequest)}
    assert set(REQUEST_PROPS) - LEGACY_ALIASES == set(fields), (
        f"openapi only: {sorted(set(REQUEST_PROPS) - LEGACY_ALIASES - set(fields))}; "
        f"oracle only: {sorted(set(fields) - set(REQUEST_PROPS))}")
    wrong = {}
    for name, f in fields.items():
        doc = REQUEST_PROPS[name].get("default", None)
        if f.default is None or f.default is dataclasses.MISSING:
            if doc is not None:
                wrong[name] = (doc, None)
        elif doc != f.default:
            wrong[name] = (doc, f.default)
    assert not wrong, f"default in openapi.yaml vs oracle: {wrong}"


@pytest.mark.parametrize("backend", sorted(KNOWN_KEY_SOURCES))
def test_request_properties_are_the_backends_known_keys(backend: str) -> None:
    src = (_KH / KNOWN_KEY_SOURCES[backend]).read_text(encoding="utf-8")
    m = _KNOWN_BLOCK.search(src)
    assert m, f"{backend}: no K7 known-key list found in {KNOWN_KEY_SOURCES[backend]}"
    known = set(re.findall(r"[\"']([a-z_0-9]+)[\"']", m.group(1)))
    assert known == set(REQUEST_PROPS), (
        f"{backend}: openapi only {sorted(set(REQUEST_PROPS) - known)}, "
        f"backend only {sorted(known - set(REQUEST_PROPS))}")


# -- 3: every response validates ---------------------------------------------------

@pytest.mark.parametrize("backend_name", BACKEND_NAMES)
def test_simulate_200_matches_the_schema(backend_name: str, backend_urls: dict) -> None:
    url = _url(backend_name, backend_urls)
    r = requests.post(f"{url}/simulate", json={"grid_nx": 16, "grid_ny": 8, "steps": 3}, timeout=TIMEOUT)
    assert r.status_code == 200, f"{backend_name}: {r.status_code} {r.text[:200]}"
    body = r.json()
    _validate(body, "#/components/schemas/SimulationResult")
    for f in ("u_velocity", "v_velocity", "vorticity", "pressure"):
        assert len(body[f]) == 16 * 8, f"{backend_name}: {f} has {len(body[f])} values, not grid_nx*grid_ny"


@pytest.mark.parametrize("backend_name", BACKEND_NAMES)
def test_health_matches_the_schema(backend_name: str, backend_urls: dict) -> None:
    url = _url(backend_name, backend_urls)
    _validate(requests.get(f"{url}/health", timeout=TIMEOUT).json(), "#/components/schemas/HealthResponse")


@pytest.mark.parametrize("backend_name", BACKEND_NAMES)
def test_info_matches_the_schema(backend_name: str, backend_urls: dict) -> None:
    url = _url(backend_name, backend_urls)
    _validate(requests.get(f"{url}/info", timeout=TIMEOUT).json(), "#/components/schemas/BackendInfo")


@pytest.mark.parametrize("backend_name", BACKEND_NAMES)
def test_400_matches_the_schema(backend_name: str, backend_urls: dict) -> None:
    url = _url(backend_name, backend_urls)
    r = requests.post(f"{url}/simulate", json={"grid_nx": 16, "grid_ny": 8, "steps": -1}, timeout=TIMEOUT)
    assert r.status_code == 400, f"{backend_name}: steps -1 returned {r.status_code}"
    body = r.json()
    _validate(body, "#/components/responses/BadRequest/content/application~1json/schema")
    assert body["field"] == "steps", f"{backend_name}: field {body['field']!r}"


# -- 4: the documented defaults are the computed ones -------------------------------

@pytest.mark.parametrize("backend_name", BACKEND_NAMES)
def test_an_empty_request_gets_the_documented_defaults(backend_name: str, backend_urls: dict) -> None:
    url = _url(backend_name, backend_urls)
    r = requests.post(f"{url}/simulate", data="{}", headers={"Content-Type": "application/json"},
                      timeout=TIMEOUT)
    assert r.status_code == 200, f"{backend_name}: {r.status_code} {r.text[:200]}"
    body = r.json()
    want = {k: REQUEST_PROPS[k]["default"] for k in ("grid_nx", "grid_ny", "steps")}
    got = {"grid_nx": body["grid_nx"], "grid_ny": body["grid_ny"], "steps": body["steps_completed"]}
    assert got == want, f"{backend_name}: an empty request computed {got}, openapi.yaml says {want}"
    t_final = REQUEST_PROPS["steps"]["default"] * REQUEST_PROPS["dt"]["default"]
    assert abs(body["t_final"] - t_final) < 1e-12, f"{backend_name}: t_final {body['t_final']} != {t_final}"
    ke, ref = body["diagnostics"]["kinetic_energy"], _oracle_default_ke()
    assert abs(ke - ref) <= ORACLE_REL * abs(ref), (
        f"{backend_name}: KE {ke!r} vs oracle defaults {ref!r} (rel {abs(ke - ref) / abs(ref):.2e})")
