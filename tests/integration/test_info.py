# test_info.py — GET /info on all 5 backends
# Task: KH-018  |  Spec: kh-sim/shared/api/openapi.yaml
#
# Each backend must return HTTP 200 with JSON body containing:
#   { "backend": <str>, "language": <str>, "framework": <str>,
#     "fft_library": <str>, "port": <int>, "openapi_spec": <str> }

import pytest
import requests

from conftest import BACKEND_NAMES, BACKEND_META

REQUIRED_INFO_FIELDS = {"backend", "language", "framework", "fft_library", "port", "openapi_spec"}


@pytest.mark.parametrize("backend_name", BACKEND_NAMES)
def test_backend_info_status_code(backend_name: str, backend_urls: dict) -> None:
    """GET /info returns HTTP 200 for each backend."""
    url = f"{backend_urls[backend_name]}/info"
    try:
        resp = requests.get(url, timeout=5)
    except requests.exceptions.ConnectionError:
        pytest.skip(f"{backend_name} not reachable at {url}")
    assert resp.status_code == 200, (
        f"{backend_name}: /info returned {resp.status_code} (expected 200)"
    )


@pytest.mark.parametrize("backend_name", BACKEND_NAMES)
def test_backend_info_required_fields(backend_name: str, backend_urls: dict) -> None:
    """GET /info response has all required fields."""
    url = f"{backend_urls[backend_name]}/info"
    try:
        resp = requests.get(url, timeout=5)
    except requests.exceptions.ConnectionError:
        pytest.skip(f"{backend_name} not reachable at {url}")

    body = resp.json()
    missing = REQUIRED_INFO_FIELDS - set(body.keys())
    assert not missing, f"{backend_name}: /info response missing fields: {missing}"


@pytest.mark.parametrize("backend_name", BACKEND_NAMES)
def test_backend_info_language_matches_manifest(backend_name: str, backend_urls: dict) -> None:
    """info.language matches expected language from BACKEND_META."""
    url = f"{backend_urls[backend_name]}/info"
    try:
        resp = requests.get(url, timeout=5)
    except requests.exceptions.ConnectionError:
        pytest.skip(f"{backend_name} not reachable at {url}")

    body = resp.json()
    expected_lang = BACKEND_META[backend_name]["language"].lower()
    actual_lang = body.get("language", "").lower()
    assert expected_lang in actual_lang, (
        f"{backend_name}: info.language '{body.get('language')}' "
        f"does not contain expected '{BACKEND_META[backend_name]['language']}'"
    )


@pytest.mark.parametrize("backend_name", BACKEND_NAMES)
def test_backend_info_port_matches_manifest(backend_name: str, backend_urls: dict) -> None:
    """info.port matches expected port from BACKEND_META."""
    url = f"{backend_urls[backend_name]}/info"
    try:
        resp = requests.get(url, timeout=5)
    except requests.exceptions.ConnectionError:
        pytest.skip(f"{backend_name} not reachable at {url}")

    body = resp.json()
    expected_port = BACKEND_META[backend_name]["port"]
    assert body.get("port") == expected_port, (
        f"{backend_name}: info.port {body.get('port')} != expected {expected_port}"
    )


@pytest.mark.parametrize("backend_name", BACKEND_NAMES)
def test_backend_info_openapi_spec_field(backend_name: str, backend_urls: dict) -> None:
    """info.openapi_spec is a non-empty string."""
    url = f"{backend_urls[backend_name]}/info"
    try:
        resp = requests.get(url, timeout=5)
    except requests.exceptions.ConnectionError:
        pytest.skip(f"{backend_name} not reachable at {url}")

    body = resp.json()
    spec = body.get("openapi_spec", "")
    # The repository path of the contract, the same in every backend (Pascal said
    # "shared/api/openapi.yaml" until 2026-10-11).
    assert spec == "kh-sim/shared/api/openapi.yaml", (
        f"{backend_name}: info.openapi_spec is {spec!r}, not kh-sim/shared/api/openapi.yaml"
    )


@pytest.mark.parametrize("backend_name", BACKEND_NAMES)
def test_backend_names_itself_once(backend_name: str, backend_urls: dict) -> None:
    """/health, /info and a /simulate result carry the same "backend" value (Pascal
    answered "pascal", "kh-pascal" and "kh-pascal" until 2026-10-11)."""
    url = backend_urls[backend_name]
    try:
        health = requests.get(f"{url}/health", timeout=5).json()
    except requests.exceptions.ConnectionError:
        pytest.skip(f"{backend_name} not reachable at {url}")
    info = requests.get(f"{url}/info", timeout=5).json()
    sim = requests.post(f"{url}/simulate", json={"grid_nx": 8, "grid_ny": 8, "steps": 1}, timeout=30).json()
    names = {"health": health.get("backend"), "info": info.get("backend"), "simulate": sim.get("backend")}
    assert len(set(names.values())) == 1 and all(names.values()), f"{backend_name}: {names}"
