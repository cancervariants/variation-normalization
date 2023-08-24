"""Module for testing gnomad_vcf_to_protein works correctly"""
import pytest
from ga4gh.vrs import models

from tests.conftest import assertion_checks
from variation.gnomad_vcf_to_protein_variation import dna_to_rna


@pytest.fixture(scope="module")
def test_handler(test_query_handler):
    """Create test fixture for gnomad vcf to protein handler"""
    return test_query_handler.gnomad_vcf_to_protein_handler


@pytest.fixture(scope="module")
def mmel1_l30m():
    """Create test fixture for MMEL1 L30M"""
    params = {
        "id": "ga4gh:VA.JFidTgmjW5Zh3xxfeWJvvKIuGZAqKv7Q",
        "location": {
            "id": "ga4gh:SL.uxCjsPbsmwDjNg4fVwUC2_490CLIzoZE",
            "end": 30,
            "start": 29,
            "sequence": "ga4gh:SQ.iQ8F_pnsiQOLohiV2qh3OWRZiftUt8jZ",
            "type": "SequenceLocation",
        },
        "state": {"sequence": "M", "type": "LiteralSequenceExpression"},
        "type": "Allele",
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def cdk11a_e314del():
    """Create test fixture for CDK11A Glu314del"""
    params = {
        "id": "ga4gh:VA.A3DMXumsfpY3uVp20m9IJO2_Ey9LgC6K",
        "location": {
            "id": "ga4gh:SL.aQOAiKa1RKPcKJH-pLp-LRMIrJ5lTC-h",
            "end": 321,
            "start": 308,
            "sequence": "ga4gh:SQ.N728VSRRMHJ1SrhJgKqJOCaa3l5Z4sqm",
            "type": "SequenceLocation",
        },
        "state": {"sequence": "EEEEEEEEEEEE", "type": "LiteralSequenceExpression"},
        "type": "Allele",
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def protein_insertion2():
    """Create test fixture for LRP8 p.Gln25_Leu26insArg"""
    params = {
        "id": "ga4gh:VA.VC-a3EBdlGAH2sQrf3p4hSNHFSiXneHg",
        "location": {
            "id": "ga4gh:SL.EyERvAv0SM0R44pskUoTDuabVG_VnjYg",
            "end": 25,
            "start": 25,
            "sequence": "ga4gh:SQ.qgIh8--4F6IpxRwX_lVtD2BhepH5B5Ef",
            "type": "SequenceLocation",
        },
        "state": {"sequence": "R", "type": "LiteralSequenceExpression"},
        "type": "Allele",
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def atad3a_loc():
    """Create test fixture for ATAD3A location"""
    return {
        "id": "ga4gh:SL.mnkdeSGvgf_NAzLsmu4mIfjA8HAerZRq",
        "end": 7,
        "start": 6,
        "sequence": "ga4gh:SQ.MHPOY_7fv8V9SktyvaTxulVFSK6XCxM8",
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def atad3a_i7v(atad3a_loc):
    """Create test fixture for ATAD3A Ile7Val"""
    params = {
        "id": "ga4gh:VA.vVSjx5GoG9oP44Ep303BOJxdNOOVro_a",
        "location": atad3a_loc,
        "state": {"sequence": "V", "type": "LiteralSequenceExpression"},
        "type": "Allele",
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def atad3a_i7t(atad3a_loc):
    """Create test fixture for ATAD3A Ile7Thr"""
    params = {
        "id": "ga4gh:VA.rFe3l5iV4EcL1_LyBysFuXQBoB-YdBFd",
        "location": atad3a_loc,
        "state": {"sequence": "T", "type": "LiteralSequenceExpression"},
        "type": "Allele",
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def atad3a_i7m(atad3a_loc):
    """Create test fixture for ATAD3A Ile7Met"""
    params = {
        "id": "ga4gh:VA.PAspG0Rsk1JGb6qLlECEyuf-u6UfyGz8",
        "location": atad3a_loc,
        "state": {"sequence": "M", "type": "LiteralSequenceExpression"},
        "type": "Allele",
    }
    return models.Allele(**params)


@pytest.fixture(scope="session")
def braf_v600l(braf_600loc):
    """Create test fixture for BRAF Val600Leu."""
    params = {
        "id": "ga4gh:VA.7_QXciVINzAB71-Xx7nezubJAlDUu1Nw",
        "location": braf_600loc,
        "state": {"sequence": "L", "type": "LiteralSequenceExpression"},
        "type": "Allele",
    }
    return models.Allele(**params)


@pytest.fixture(scope="session")
def braf_600_reference_agree(braf_600loc):
    """Create test fixture for BRAF Val600=."""
    params = {
        "id": "ga4gh:VA.fLpGu4zNF8TqU3SPCWDB7JpOmHSHK1E0",
        "location": braf_600loc,
        "state": {"sequence": "V", "type": "LiteralSequenceExpression"},
        "type": "Allele",
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def kras_g12d():
    """Fixture for KRAS G12C"""
    return {
        "id": "ga4gh:VA.8Hx7YSnUL2pSa9qz0jGuUCK3dEHYfLYp",
        "type": "Allele",
        "location": {
            "id": "ga4gh:SL.LDqv0rgy82D_1KR2Tn82b7XqUScZcd44",
            "type": "SequenceLocation",
            "sequence": "ga4gh:SQ.fytWhQSNGnA-86vDiQCxTSzgkk_WfQRS",
            "start": 11,
            "end": 12,
        },
        "state": {"type": "LiteralSequenceExpression", "sequence": "D"},
    }


def test_dna_to_rna():
    """Test that dna_to_rna method works correctly."""
    resp = dna_to_rna("GTA")
    assert resp == "CAU"

    resp = dna_to_rna("AAGTGACA")
    assert resp == "UUCACUGU"


@pytest.mark.asyncio
async def test_substitution(
    test_handler,
    braf_v600e,
    braf_v600l,
    braf_600_reference_agree,
    mmel1_l30m,
    atad3a_i7v,
    atad3a_i7t,
    atad3a_i7m,
    kras_g12d,
):
    """Test that substitution queries return correct response"""
    # Reading Frame 1, Negative Strand
    resp = await test_handler.gnomad_vcf_to_protein("7-140753337-C-A")
    assertion_checks(resp, braf_v600l)
    assert resp.warnings == []

    # Reading Frame 2, Negative Strand
    resp = await test_handler.gnomad_vcf_to_protein("7-140753336-A-T")
    assertion_checks(resp, braf_v600e)
    assert resp.warnings == []

    # Reading Frame 3, Negative Strand
    resp = await test_handler.gnomad_vcf_to_protein("7-140753335-C-A")
    assertion_checks(resp, braf_600_reference_agree)
    assert resp.warnings == []

    # Reading Frame 3, Negative Strand
    resp = await test_handler.gnomad_vcf_to_protein("1-2629397-G-T")
    assertion_checks(resp, mmel1_l30m)
    assert resp.warnings == []

    # Reading Frame 1, Positive Strand
    resp = await test_handler.gnomad_vcf_to_protein("1-1512287-A-G")
    assertion_checks(resp, atad3a_i7v)
    assert resp.warnings == []

    # Reading Frame 2, Positive Strand
    resp = await test_handler.gnomad_vcf_to_protein("1-1512288-T-C")
    assertion_checks(resp, atad3a_i7t)
    assert resp.warnings == []

    # Reading Frame 3, Positive Strand
    resp = await test_handler.gnomad_vcf_to_protein("1-1512289-T-G")
    assertion_checks(resp, atad3a_i7m)
    assert resp.warnings == []

    resp = await test_handler.gnomad_vcf_to_protein("12-25245350-C-T")
    assertion_checks(resp, kras_g12d)


@pytest.mark.asyncio
async def test_reference_agree(test_handler, vhl_reference_agree):
    """Test that reference agree queries return correct response"""
    # https://www.ncbi.nlm.nih.gov/clinvar/variation/379039/?new_evidence=true
    resp = await test_handler.gnomad_vcf_to_protein("3-10142030-C-C")
    assertion_checks(resp, vhl_reference_agree)
    assert resp.warnings == []


@pytest.mark.asyncio
async def test_insertion(test_handler, protein_insertion, protein_insertion2):
    """Test that insertion queries return correct response"""
    resp = await test_handler.gnomad_vcf_to_protein("7-55181319-C-CGGGTTG")
    assertion_checks(resp, protein_insertion)
    assert resp.warnings == []

    resp = await test_handler.gnomad_vcf_to_protein("1-53327836-A-ACGC")
    assertion_checks(resp, protein_insertion2)
    assert resp.warnings == []


@pytest.mark.asyncio
async def test_deletion(test_handler, protein_deletion_np_range, cdk11a_e314del):
    """Test that deletion queries return correct response"""
    resp = await test_handler.gnomad_vcf_to_protein("17-39723966-TTGAGGGAAAACACAT-T")
    assertion_checks(resp, protein_deletion_np_range)
    assert resp.warnings == []

    resp = await test_handler.gnomad_vcf_to_protein("1-1708855-TTCC-T")
    assertion_checks(resp, cdk11a_e314del)
    assert resp.warnings == []


@pytest.mark.asyncio
async def test_invalid(test_handler):
    """Test that invalid queries return correct response"""
    resp = await test_handler.gnomad_vcf_to_protein("dummy")
    assert resp.variation is None

    resp = await test_handler.gnomad_vcf_to_protein("BRAF V600E")
    assert resp.variation is None
    assert resp.warnings == ["BRAF V600E is not a supported gnomad vcf query"]

    resp = await test_handler.gnomad_vcf_to_protein("7-140753336-T-G")
    assert resp.variation is None
    assert set(resp.warnings) == {
        "Expected T but found A on NC_000007.14 at position 140753336"
    }

    resp = await test_handler.gnomad_vcf_to_protein("20-2-TC-TG")
    assert resp.variation is None
    assert resp.warnings == ["20-2-TC-TG is not a valid gnomad vcf query"]
