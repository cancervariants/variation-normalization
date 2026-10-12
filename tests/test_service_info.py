"""Tests for the GA4GH service-info endpoint."""

from fastapi.testclient import TestClient
from ga4gh.vrs import VRS_VERSION

from variation import __version__
from variation.main import app


def test_service_info() -> None:
    """Expose service identity, VRS version, and reference data releases."""
    response = TestClient(app).get("/service-info")

    assert response.status_code == 200
    info = response.json()
    assert info["id"] == "org.cancervariants.variation_normalizer"
    assert info["name"] == "Variation Normalizer"
    assert info["type"] == {
        "group": "org.cancervariants",
        "artifact": "Variation Normalizer API",
        "version": __version__,
    }
    assert info["organization"]["url"] == "https://cancervariants.org"
    assert info["version"] == __version__
    assert info["specMetadata"] == {"vrsVersion": VRS_VERSION}
    assert info["dataVersions"] == {
        "utaSchema": "uta_20241220",
        "seqrepoVersion": "2024-12-20",
    }
    assert info["documentationUrl"].startswith("https://")
