"""Module for testing gnomad_vcf_to_protein works correctly"""
import pytest
from ga4gh.vrs import models

from tests.conftest import assertion_checks


@pytest.fixture(scope="module")
def test_handler(test_query_handler):
    """Create test fixture for gnomad vcf to protein handler"""
    return test_query_handler.gnomad_vcf_to_protein_handler


@pytest.fixture(scope="module")
def mmel1_l30m():
    """Create test fixture for MMEL1 L30M"""
    params = {
        "id": "ga4gh:VA.OqqETz467CITELOZsYDukkab7JaOWiZf",
        "location": {
            "id": "ga4gh:SL.Q7kfcqUWpIyEOgxcgPK1sRfgWPDv7zKA",
            "end": 30,
            "start": 29,
            "sequenceReference": {
                "type": "SequenceReference",
                "refgetAccession": "SQ.iQ8F_pnsiQOLohiV2qh3OWRZiftUt8jZ",
            },
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
        "id": "ga4gh:VA._CVnGazN6KosqrFnDx7kny-rb6yAZWtB",
        "location": {
            "id": "ga4gh:SL.VqI6HuIFmm4XP3ocOTaobGxwqg4m6Ooi",
            "end": 321,
            "start": 308,
            "sequenceReference": {
                "type": "SequenceReference",
                "refgetAccession": "SQ.N728VSRRMHJ1SrhJgKqJOCaa3l5Z4sqm",
            },
            "type": "SequenceLocation",
        },
        "state": {
            "length": 12,
            "repeatSubunitLength": 1,
            "sequence": "EEEEEEEEEEEE",
            "type": "ReferenceLengthExpression",
        },
        "type": "Allele",
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def protein_insertion2():
    """Create test fixture for LRP8 p.Gln25_Leu26insArg (CA860540)"""
    params = {
        "id": "ga4gh:VA.5KWhsli69ac5zyoGf40Owu4CVNKy27So",
        "location": {
            "id": "ga4gh:SL.I4c4NL0g3vBajHe44faZFQtrcqrbA14d",
            "end": 25,
            "start": 25,
            "sequenceReference": {
                "type": "SequenceReference",
                "refgetAccession": "SQ.qgIh8--4F6IpxRwX_lVtD2BhepH5B5Ef",
            },
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
        "id": "ga4gh:SL.xiP3uciIfJy_f44wNKCBvtsb35BC330Q",
        "end": 7,
        "start": 6,
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.MHPOY_7fv8V9SktyvaTxulVFSK6XCxM8",
        },
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def atad3a_i7v(atad3a_loc):
    """Create test fixture for ATAD3A Ile7Val"""
    params = {
        "id": "ga4gh:VA.i_L_bjPfI4XLMIKmVklV6eDLKEl1f7PD",
        "location": atad3a_loc,
        "state": {"sequence": "V", "type": "LiteralSequenceExpression"},
        "type": "Allele",
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def atad3a_i7t(atad3a_loc):
    """Create test fixture for ATAD3A Ile7Thr"""
    params = {
        "id": "ga4gh:VA.C8QO-YAfG66yj7cEwjEhkEfSd-oCSKfc",
        "location": atad3a_loc,
        "state": {"sequence": "T", "type": "LiteralSequenceExpression"},
        "type": "Allele",
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def atad3a_i7m(atad3a_loc):
    """Create test fixture for ATAD3A Ile7Met"""
    params = {
        "id": "ga4gh:VA.Fhmv3GK3bcIJRXOkigS9QNMzAWGW3WGa",
        "location": atad3a_loc,
        "state": {"sequence": "M", "type": "LiteralSequenceExpression"},
        "type": "Allele",
    }
    return models.Allele(**params)


@pytest.fixture(scope="session")
def braf_v600l(braf_600loc):
    """Create test fixture for BRAF Val600Leu."""
    params = {
        "id": "ga4gh:VA.c6f1MPfquVRPZO46wVzCaGaU8QnXoHNN",
        "location": braf_600loc,
        "state": {"sequence": "L", "type": "LiteralSequenceExpression"},
        "type": "Allele",
    }
    return models.Allele(**params)


@pytest.fixture(scope="session")
def braf_600_reference_agree(braf_600loc):
    """Create test fixture for BRAF Val600=."""
    params = {
        "id": "ga4gh:VA.wS6kJNbPkRJDIWg8F4CjOMQ5mcJzD_X4",
        "location": braf_600loc,
        "state": {"sequence": "V", "type": "LiteralSequenceExpression"},
        "type": "Allele",
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def kras_g12d():
    """Fixture for KRAS G12D"""
    params = {
        "id": "ga4gh:VA.CB571ja_KfZM_Hjn9zjjgV1an3tDWRcl",
        "type": "Allele",
        "location": {
            "id": "ga4gh:SL.OndkjmujtyUEZSjjCv0C-gpwnVbRgfj8",
            "type": "SequenceLocation",
            "sequenceReference": {
                "type": "SequenceReference",
                "refgetAccession": "SQ.fytWhQSNGnA-86vDiQCxTSzgkk_WfQRS",
            },
            "start": 11,
            "end": 12,
        },
        "state": {"type": "LiteralSequenceExpression", "sequence": "D"},
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def multi_nuc_sub_pos():
    """Create test fixture for substitution with more than 1 nucleotide change on the
    positive strand (CA16042245)
    """
    params = {
        "id": "ga4gh:VA.q_fhdDFpLT38y6TPU1EMtc2StwRGRVx0",
        "type": "Allele",
        "location": {
            "id": "ga4gh:SL.xcqWKkKU83UBmLq5Q1yrH6IyLvWkMcWk",
            "type": "SequenceLocation",
            "sequenceReference": {
                "type": "SequenceReference",
                "refgetAccession": "SQ.HtNf7YrmmFih3cwRwYMlylPFMAs7-l9B",
            },
            "start": 242,
            "end": 244,
        },
        "state": {"type": "LiteralSequenceExpression", "sequence": "PS"},
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def multi_nuc_sub_neg():
    """Create test fixture for substitution with more than 1 nucleotide change on the
    negative strand (CA1139661942)
    """
    params = {
        "id": "ga4gh:VA.6K950oAyNXfIkPTmSHzX8f8wpNUroBGK",
        "type": "Allele",
        "location": {
            "id": "ga4gh:SL.RiO0HagK846MBCJZK7NVcDtlYPkdaXC4",
            "type": "SequenceLocation",
            "sequenceReference": {
                "type": "SequenceReference",
                "refgetAccession": "SQ.bg8P_l39rOUQVLwsW3Dme-946Od8-3rB",
            },
            "start": 235,
            "end": 236,
        },
        "state": {"type": "LiteralSequenceExpression", "sequence": "G"},
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def delins_pos():
    """Create test fixture for delins on positive strand (CA645561524)"""
    params = {
        "id": "ga4gh:VA.vBQ2TCfRHiG3ud_vqE88BNZEK7Qw28kg",
        "type": "Allele",
        "location": {
            "id": "ga4gh:SL.1btwhKRj0zZQwo-_CalR-WavTB019t-V",
            "type": "SequenceLocation",
            "sequenceReference": {
                "type": "SequenceReference",
                "refgetAccession": "SQ.vyo55F6mA6n2LgN4cagcdRzOuh38V4mE",
            },
            "start": 746,
            "end": 752,
        },
        "state": {"type": "LiteralSequenceExpression", "sequence": "Q"},
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def delins_neg():
    """Create test fixture for the protein consequence of a delins on negative strand (ClinVar ID 1217291)"""
    params = {
        "id": "ga4gh:VA.iSDLORgPGz21BetTQM5grpXyB3tIfZwl",
        "type": "Allele",
        "location": {
            "id": "ga4gh:SL.cmPqZ9vkQlmV_eWEW3MFKFeVtT_n1-yc",
            "type": "SequenceLocation",
            "sequenceReference": {
                "type": "SequenceReference",
                "refgetAccession": "SQ.aDZLb9cs0cYskMKDIK-AXhaevHRA86JS",
            },
            "start": 239,
            "end": 259,
        },
        "state": {"type": "LiteralSequenceExpression", "sequence": "TLTA"},
    }
    return models.Allele(**params)


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
    multi_nuc_sub_pos,
    multi_nuc_sub_neg,
):
    """Test that substitution queries return correct response"""
    # Reading Frame 1, Negative Strand
    resp = await test_handler.gnomad_vcf_to_protein("7-140753337-C-A")
    assertion_checks(resp, braf_v600l)
    assert resp.warnings == []

    # Reading Frame 2, Negative Strand
    resp = await test_handler.gnomad_vcf_to_protein("7-140753336-A-T")
    assertion_checks(resp, braf_v600e)
    assert resp.gene_context
    assert resp.vrs_ref_allele_seq == "V"
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

    # multi nucleotide ref/alt on positive strand (CA16042245)
    resp = await test_handler.gnomad_vcf_to_protein("2-74530927-TGC-CAT")
    assertion_checks(resp, multi_nuc_sub_pos)

    # multi nucleotide ref/alt on negative strand (CA1139661942)
    resp = await test_handler.gnomad_vcf_to_protein("11-47348490-TG-CA")
    assertion_checks(resp, multi_nuc_sub_neg)


@pytest.mark.asyncio
async def test_reference_agree(test_handler, vhl_reference_agree):
    """Test that reference agree queries return correct response"""
    # https://www.ncbi.nlm.nih.gov/clinvar/variation/379039/?new_evidence=true
    resp = await test_handler.gnomad_vcf_to_protein("3-10142030-C-T")
    assertion_checks(resp, vhl_reference_agree)
    assert resp.vrs_ref_allele_seq == "P"
    assert resp.gene_context
    assert resp.warnings == []


@pytest.mark.asyncio
async def test_insertion(test_handler, protein_insertion, protein_insertion2):
    """Test that insertion queries return correct response"""
    # positive strand (CA645561585)
    resp = await test_handler.gnomad_vcf_to_protein("7-55181319-C-CGGGTTA")
    assertion_checks(resp, protein_insertion)
    assert resp.vrs_ref_allele_seq is None
    assert resp.gene_context
    assert resp.warnings == []

    # negative strand (CA860540)
    resp = await test_handler.gnomad_vcf_to_protein("1-53327836-A-AGCC")
    assertion_checks(resp, protein_insertion2)
    assert resp.vrs_ref_allele_seq is None
    assert resp.gene_context
    assert resp.warnings == []


@pytest.mark.asyncio
async def test_deletion(test_handler, protein_deletion_np_range, cdk11a_e314del):
    """Test that deletion queries return correct response"""
    resp = await test_handler.gnomad_vcf_to_protein("17-39723966-TTGAGGGAAAACACAT-T")
    assertion_checks(resp, protein_deletion_np_range)
    assert resp.vrs_ref_allele_seq == "LRENT"
    assert resp.gene_context
    assert resp.warnings == []

    resp = await test_handler.gnomad_vcf_to_protein("1-1708855-TTCC-T")
    assertion_checks(resp, cdk11a_e314del)
    assert resp.gene_context
    assert resp.warnings == []


@pytest.mark.asyncio
async def test_delins(test_handler, delins_pos, delins_neg):
    """Test that delins queries return correct response"""
    # CA645561524, Positive Strand
    resp = await test_handler.gnomad_vcf_to_protein("7-55174776-TTAAGAGAAGCAACATCT-CAA")
    assertion_checks(resp, delins_pos)
    assert resp.vrs_ref_allele_seq == "LREATS"
    assert resp.gene_context

    # ClinVar ID 1217291, Negative Strand
    resp = await test_handler.gnomad_vcf_to_protein(
        "X-153870419-GCTGCCCCTGCAAGGCCACCAGGTGGCTGCTGGAGTTGGTGGGGAAGAGCAGGCGCGG-CTGTCAATGT"
    )
    assertion_checks(resp, delins_neg)
    assert resp.vrs_ref_allele_seq == "PRLLFPTNSSSHLVALQGQP"
    assert resp.gene_context


@pytest.mark.asyncio
async def test_invalid(test_handler):
    """Test that invalid queries return correct response"""
    resp = await test_handler.gnomad_vcf_to_protein("BRAF V600E")
    assert resp.variation is None
    assert resp.vrs_ref_allele_seq is None
    assert resp.gene_context is None
    assert resp.warnings == [
        "BRAF V600E is not a gnomAD VCF-like query (`chr-pos-ref-alt`)"
    ]

    resp = await test_handler.gnomad_vcf_to_protein("7-140753336-T-G")
    assert resp.variation is None
    assert resp.vrs_ref_allele_seq is None
    assert resp.gene_context is None
    assert set(resp.warnings) == {"Unable to get cDNA and protein representation"}

    resp = await test_handler.gnomad_vcf_to_protein("20-2-TC-TG")
    assert resp.variation is None
    assert resp.vrs_ref_allele_seq is None
    assert resp.gene_context is None
    assert resp.warnings == ["20-2-TC-TG is not a valid gnomad vcf query"]
