"""Module for testing gnomad_vcf_to_protein works correctly"""
import pytest
from variation.query import QueryHandler
from tests.conftest import assertion_checks


@pytest.fixture(scope="module")
def test_query_handler():
    """Create test fixture for query handler"""
    return QueryHandler()


def test_substitution(test_query_handler, braf_v600e):
    """Test that substitution queries return correct response"""
    resp, w = test_query_handler.gnomad_vcf_to_protein("7-140753336-A-T")
    assertion_checks(resp, braf_v600e, ignore_id=True)
    assert w == []


def test_silent_mutation(test_query_handler, vhl_silent):
    """Test that silent queries return correct response"""
    resp, w = test_query_handler.gnomad_vcf_to_protein("3-10183714-C-C")
    assertion_checks(resp, vhl_silent, ignore_id=True)
    assert w == []


def test_invalid(test_query_handler):
    """Test that invalid queries return correct response"""
    resp, w = test_query_handler.gnomad_vcf_to_protein("BRAF V600E")
    assert resp.variation.type == "Text"
    assert w == ["BRAF V600E is not a supported gnomad vcf query"]

    resp, w = test_query_handler.gnomad_vcf_to_protein("7-140753336-T-G")
    assert resp.variation.type == "Text"
    assert w == ["Unable to get protein variation for 7-140753336-T-G"]
