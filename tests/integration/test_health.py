# test_health.py — GET /health on all 6 LDE services
# Task: KH-018  |  Spec: kh-sim/shared/api/openapi.yaml
#
# Each backend must return HTTP 200 with JSON body containing:
#   { "status": <string>, "backend": <string>, "port": <int> }
# Log service returns HTTP 200 with:
#   { "status": <string>, "mongo": <string>, "uptime_s": <float> }
# Log service may return 503 when MongoDB is unavailable (degraded mode).

import pytest
import requests

from conftest import BACKEND_NAMES, BACKEND_META

VALID_STATUS_VALUES = {"ok", "healthy", "up", "running"}


@pytest.mark.parametrize("backend_name", BACKEND_NAMES)
def test_backend_health_status_code(backend_name: str, backend_urls: dict) -> None:
    """GET /health returns HTTP 200 for each backend."""
    url = f"{backend_urls[backend_name]}/health"
    try:
        resp = requests.get(url, timeout=5)
    except requests.exceptions.ConnectionError:
        pytest.skip(f"{backend_name} not reachable at {url}")
    assert resp.status_code == 200, (
        f"{backend_name}: /health returned {resp.status_code} (expected 200)"
    )


@pytest.mark.parametrize("backend_name", BACKEND_NAMES)
def test_backend_health_json_schema(backend_name: str, backend_urls: dict) -> None:
    """GET /health response body has required fields with correct types."""
    url = f"{backend_urls[backend_name]}/health"
    try:
        resp = requests.get(url, timeout=5)
    except requests.exceptions.ConnectionError:
        pytest.skip(f"{backend_name} not reachable at {url}")

    assert resp.status_code == 200
    body = resp.json()

    assert "status" in body, f"{backend_name}: health response missing 'status'"
    assert "backend" in body, f"{backend_name}: health response missing 'backend'"
    assert "port" in body, f"{backend_name}: health response missing 'port'"

    assert isinstance(body["status"], str), f"{backend_name}: 'status' must be a string"
    assert isinstance(body["port"], int), f"{backend_name}: 'port' must be an int"


@pytest.mark.parametrize("backend_name", BACKEND_NAMES)
def test_backend_health_status_value(backend_name: str, backend_urls: dict) -> None:
    """GET /health returns a recognised status value (ok / healthy / up / running)."""
    url = f"{backend_urls[backend_name]}/health"
    try:
        resp = requests.get(url, timeout=5)
    except requests.exceptions.ConnectionError:
        pytest.skip(f"{backend_name} not reachable at {url}")

    body = resp.json()
    status_val = body.get("status", "").lower()
    assert status_val in VALID_STATUS_VALUES, (
        f"{backend_name}: unexpected status value '{status_val}' "
        f"(expected one of {VALID_STATUS_VALUES})"
    )


@pytest.mark.parametrize("backend_name", BACKEND_NAMES)
def test_backend_health_port_matches_manifest(backend_name: str, backend_urls: dict) -> None:
    """Port reported in /health matches the expected service port from manifest."""
    url = f"{backend_urls[backend_name]}/health"
    try:
        resp = requests.get(url, timeout=5)
    except requests.exceptions.ConnectionError:
        pytest.skip(f"{backend_name} not reachable at {url}")

    body = resp.json()
    expected_port = BACKEND_META[backend_name]["port"]
    assert body.get("port") == expected_port, (
        f"{backend_name}: health.port {body.get('port')} != expected {expected_port}"
    )


def test_log_service_health_status_code(log_url: str) -> None:
    """GET /health on log service returns 200 (connected) or 503 (degraded/no MongoDB)."""
    url = f"{log_url}/health"
    try:
        resp = requests.get(url, timeout=5)
    except requests.exceptions.ConnectionError:
        pytest.skip(f"Log service not reachable at {url}")
    assert resp.status_code in (200, 503), (
        f"log-service: /health returned {resp.status_code} (expected 200 or 503)"
    )


def test_log_service_health_json_schema(log_url: str) -> None:
    """GET /health on log service body has required fields."""
    url = f"{log_url}/health"
    try:
        resp = requests.get(url, timeout=5)
    except requests.exceptions.ConnectionError:
        pytest.skip(f"Log service not reachable at {url}")

    body = resp.json()
    assert "status" in body, "log-service: health response missing 'status'"
    assert "mongo" in body, "log-service: health response missing 'mongo'"
    assert "uptime_s" in body, "log-service: health response missing 'uptime_s'"
    assert isinstance(body["uptime_s"], (int, float)), "'uptime_s' must be numeric"
