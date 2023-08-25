"""Module for testing the normalize endpoint."""
import copy
from datetime import datetime

import pytest
from ga4gh.vrs import models

from tests.conftest import assertion_checks
from variation.main import normalize as normalize_get_response
from variation.main import to_vrs as to_vrs_get_response
from variation.schemas.normalize_response_schema import HGVSDupDelModeOption


@pytest.fixture(scope="module")
def test_handler(test_query_handler):
    """Create test fixture for normalize handler"""
    return test_query_handler.normalize_handler


@pytest.fixture(scope="module")
def dis3_p63a():
    """Create DIS3 P63A test fixture."""
    params = {
        "id": "ga4gh:VA.Firb663fuUVIgiHI0HHEHQXcZPG09Kaq",
        "location": {
            "id": "ga4gh:SL.NDImt2aNbeWtYzzVYce-Wxa-1TjNXFVG",
            "end": 63,
            "start": 62,
            "sequence": "ga4gh:SQ.mlWsxfPKINN3o300stAI8oqN5U7P6kEu",
            "type": "SequenceLocation",
        },
        "state": {"sequence": "A", "type": "LiteralSequenceExpression"},
        "type": "Allele",
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def vhl():
    """Create VHL Tyr185Ter fixture."""
    params = {
        "id": "ga4gh:VA.1pMYEAUEgVyehK1Cc64YF2RfGa54CIsH",
        "location": {
            "id": "ga4gh:SL.lTisPm8tmxS0sJmROtX2m5elUAFAv4fF",
            "end": 185,
            "start": 184,
            "sequence": "ga4gh:SQ.z-Oa0pZkJ6GHJHOYM7h5mY_umc0SJzTu",
            "type": "SequenceLocation",
        },
        "state": {"sequence": "*", "type": "LiteralSequenceExpression"},
        "type": "Allele",
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def braf_v600e_nucleotide(braf_nuc_value):
    """Create a test fixture for BRAF V600E MANE select nucleotide hgvs."""
    variation = copy.deepcopy(braf_nuc_value)
    variation["id"] = "ga4gh:VA.X7YPS2mcsfRvPjzVuuBQkspmpv6ip7js"
    return models.Allele(**variation)


@pytest.fixture(scope="module")
def nm_004448_cdna_delins():
    """Create test fixture for NM_004448.4:c.2326_2327delinsCT."""
    params = {
        "id": "ga4gh:VA.xYgK3CrjWLHXZGEmc0YsE69X2Z5Bznw6",
        "location": {
            "id": "ga4gh:SL.fNzNBM8ZjkOIVZUYFT65PiZUnQbY-d3i",
            "end": 2502,
            "start": 2500,
            "sequence": "ga4gh:SQ.y9b4LVMiCXpZxOg9Xt1NwRtssA03MwWM",
            "type": "SequenceLocation",
        },
        "state": {"sequence": "CT", "type": "LiteralSequenceExpression"},
        "type": "Allele",
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def nm_000551():
    """Create test fixture for NM_000551.4:c.615delinsAA."""
    params = {
        "id": "ga4gh:VA.y3yZYOOYlpDKJhhwgl8Lv2Oz9QE1m6wy",
        "location": {
            "id": "ga4gh:SL.ENrWsh_J9cc2S0ZX5e1eeuVLekbZW9vR",
            "end": 685,
            "start": 684,
            "sequence": "ga4gh:SQ.xBKOKptLLDr-k4hTyCetvARn16pDS_rW",
            "type": "SequenceLocation",
        },
        "state": {"sequence": "AA", "type": "LiteralSequenceExpression"},
        "type": "Allele",
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def braf_nuc_value():
    """Create test fixture for BRAF V600E value on c. coordinate."""
    return {
        "location": {
            "id": "ga4gh:SL.aYiJkv-7VIcMYHV3b3n_qpQOs3rxgQ-9",
            "end": 2025,
            "start": 2024,
            "sequence": "ga4gh:SQ.aKMPEJgmlZXt_F6gRY5cUG3THH2n-GUa",
            "type": "SequenceLocation",
        },
        "state": {"sequence": "A", "type": "LiteralSequenceExpression"},
        "type": "Allele",
    }


@pytest.fixture(scope="module")
def cdna_reference_agree(braf_nuc_value):
    """Create test fixture for NM_004333.4:c.1799=."""
    value = copy.deepcopy(braf_nuc_value)
    value["state"]["sequence"] = "T"
    value["id"] = "ga4gh:VA.tvJL4gSZd7_5XZxBxC1VFrw71W4gpJiz"
    return models.Allele(**value)


@pytest.fixture(scope="module")
def protein_delins():
    """Create test fixture for protein delins."""
    params = {
        "id": "ga4gh:VA.yKy4K4yCQrzxUDq9-C0Yv_ZjbIbs2K5u",
        "location": {
            "id": "ga4gh:SL.vdK2LhK3Fk5c7gxM6xs2ivcViXDZbZCL",
            "end": 751,
            "start": 746,
            "sequence": "ga4gh:SQ.vyo55F6mA6n2LgN4cagcdRzOuh38V4mE",
            "type": "SequenceLocation",
        },
        "state": {"sequence": "P", "type": "LiteralSequenceExpression"},
        "type": "Allele",
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def cdna_deletion():
    """Create test fixture for cdna deletion range with deleted
    sequence.
    """
    params = {
        "id": "ga4gh:VA.bM8kEbQCUg0dxyt_gwG3JHKldHohnAio",
        "location": {
            "id": "ga4gh:SL.bry-Y_Z-XEJWbAxdi6Txl4GHSgQesvfT",
            "end": 2453,
            "start": 2437,
            "sequence": "ga4gh:SQ.y9b4LVMiCXpZxOg9Xt1NwRtssA03MwWM",
            "type": "SequenceLocation",
        },
        "state": {
            "length": 1,
            "repeatSubunitLength": 15,
            "sequence": "T",
            "type": "ReferenceLengthExpression",
        },
        "type": "Allele",
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def genomic_deletion():
    """Create test fixture for genomic deletion range with deleted sequence.
    (CA915940709)
    """
    params = {
        "id": "ga4gh:VA.q96pW8eD17vNcJ-OQmJB9WSfMXLVihSQ",
        "location": {
            "id": "ga4gh:SL.ZcUsDwkeZQKbeeEHK-GgqSrpEtVUfZFV",
            "end": 10146528,
            "start": 10146524,
            "sequence": "ga4gh:SQ.Zu7h9AggXxhTaGVsy7h_EZSChSZGcmgX",
            "type": "SequenceLocation",
        },
        "state": {
            "length": "2",
            "repeatSubunitLength": 2,
            "sequence": "CT",
            "type": "ReferenceLengthExpression",
        },
        "type": "Allele",
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def cdna_insertion():
    """Create test fixture for coding DNA insertion."""
    params = {
        "id": "ga4gh:VA.8M4pujzetYpLADKc-NM1hab5mF9BzDhu",
        "location": {
            "id": "ga4gh:SL.QdZya5aHWvF8y71LXksoB5lo5_V3R4M1",
            "end": 2160,
            "start": 2160,
            "sequence": "ga4gh:SQ.7_mlQyDN-uWH0RlxTQFvFEv6ykd2D-xF",
            "type": "SequenceLocation",
        },
        "state": {"sequence": "A", "type": "LiteralSequenceExpression"},
        "type": "Allele",
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def genomic_insertion():
    """Create a gene insertion test fixture."""
    params = {
        "id": "ga4gh:VA.tucTDF6gu-OpgyVfJcnhcCR5jobIV1fL",
        "location": {
            "id": "ga4gh:SL.doaqZtlQYh108Ct7PAbiEC8WH4MeG6Sh",
            "end": 2500,
            "start": 2488,
            "sequence": "ga4gh:SQ.y9b4LVMiCXpZxOg9Xt1NwRtssA03MwWM",
            "type": "SequenceLocation",
        },
        "state": {
            "length": 24,
            "repeatSubunitLength": 12,
            "sequence": "TACGTGATGGCTTACGTGATGGCT",
            "type": "ReferenceLengthExpression",
        },
        "type": "Allele",
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def genomic_substitution():
    """Create a gene insertion test fixture."""
    params = {
        "id": "ga4gh:VA.TM4PBSbll2OQFtnBtRnG3W92WwcS_-_B",
        "location": {
            "id": "ga4gh:SL.zI2UEoIKmwiRPTa7DOF-yAYJGGlSMWNq",
            "end": 2630,
            "start": 2629,
            "sequence": "ga4gh:SQ.d_QsP29RWJi6bac7GOC9cJ9AO7s_HUMN",
            "type": "SequenceLocation",
        },
        "state": {"sequence": "T", "type": "LiteralSequenceExpression"},
        "type": "Allele",
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def gnomad_vcf_genomic_sub_mnv():
    """Create a genomic substitution mnv test fixture for 5-112175770-GGAA-AGAA."""
    params = {
        "id": "ga4gh:VA.ObKZgkiIXlEvpCKRdSRU5newD7-Qy5Fd",
        "location": {
            "id": "ga4gh:SL.b4dx32Llqi4q7OCx8OMDkSiZXpBrxsVp",
            "end": 112840073,
            "start": 112840072,
            "sequence": "ga4gh:SQ.aUiQCzCPZ2d0csHbMSbh2NzInhonSXwI",
            "type": "SequenceLocation",
        },
        "state": {"sequence": "A", "type": "LiteralSequenceExpression"},
        "type": "Allele",
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def genomic_sub_grch38():
    """Create a genomic substitution GRCh38 test fixture."""
    params = {
        "id": "ga4gh:VA.XZelh8iOTaifb4iezhv5iAOhE_37LXLl",
        "location": {
            "id": "ga4gh:SL.Rj_Hi3p8HwNwKfm3Felsuz_DvcwxqhfE",
            "end": 55181378,
            "start": 55181377,
            "sequence": "ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
            "type": "SequenceLocation",
        },
        "state": {"sequence": "T", "type": "LiteralSequenceExpression"},
        "type": "Allele",
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def grch38_braf_genom_reference_agree():
    """Create a genomic reference agree GRCh38 test fixture for BRAF."""
    params = {
        "id": "ga4gh:VA.O4zd1lJS-OwKukrx6q1n1r4IUUavgZqy",
        "location": {
            "id": "ga4gh:SL.j9b9mPXETlWMYDjZ62S9DIXATtt-vY_e",
            "end": 140753336,
            "start": 140753335,
            "sequence": "ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
            "type": "SequenceLocation",
        },
        "state": {"sequence": "A", "type": "LiteralSequenceExpression"},
        "type": "Allele",
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def grch38_genomic_delins1():
    """Create a test fixture for NC_000007.13:g.140453135_140453136delinsAT."""
    params = {
        "id": "ga4gh:VA.F4m9x5Z73pjeVbrPiO44o_gJElLfc4Xt",
        "location": {
            "id": "ga4gh:SL.Xy7o78eUOMPhpI79TDWYy4-a_x83AFC4",
            "end": 140753336,
            "start": 140753334,
            "sequence": "ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
            "type": "SequenceLocation",
        },
        "state": {"sequence": "AT", "type": "LiteralSequenceExpression"},
        "type": "Allele",
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def grch38_genomic_delins2():
    """Create a test fixture for NC_000003.12:g.10149938delinsAA."""
    params = {
        "id": "ga4gh:VA.WWqzHzKMaEZ3g8KXIjf0uhxMEHz06vy6",
        "location": {
            "id": "ga4gh:SL.TvXHMKWMfYtS-a6c3_ywW59WOOMxcO5M",
            "start": 10149937,
            "end": 10149938,
            "sequence": "ga4gh:SQ.Zu7h9AggXxhTaGVsy7h_EZSChSZGcmgX",
            "type": "SequenceLocation",
        },
        "state": {"sequence": "AA", "type": "LiteralSequenceExpression"},
        "type": "Allele",
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def genomic_delins_gene():
    """Create a test fixture for BRAF g.140453135_140453136delinsAT (CA16602419)."""
    params = {
        "id": "ga4gh:VA.AkXi7lrmLzjAAD6jQBHbl9tEfYb9anaV",
        "location": {
            "id": "ga4gh:SL.zE-4ig0O2S3q9pmKdtyM57ONXgiOGAtM",
            "start": 2024,
            "end": 2026,
            "sequence": "ga4gh:SQ.aKMPEJgmlZXt_F6gRY5cUG3THH2n-GUa",
            "type": "SequenceLocation",
        },
        "state": {"sequence": "AT", "type": "LiteralSequenceExpression"},
        "type": "Allele",
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def gnomad_vcf_genomic_delins1():
    """Create a test fixture for 3-37050340-AAAAGCTTTA-GAGGCTTT.

    https://reg.clinicalgenome.org/redmine/projects/registry/genboree_registry/
    allele?hgvsOrDescriptor=NM_000249.3%3Ac.489_498delinsGAGGCTTT
    """
    params = {
        "id": "ga4gh:VA.QgSSPDvaSQXc9SOWH5NxiY6jGyhjMUcA",
        "location": {
            "id": "ga4gh:SL.VF8R5QN_WMl3XKUNTYq_cq5mFfc8yB-J",
            "start": 37008848,
            "end": 37008858,
            "sequence": "ga4gh:SQ.Zu7h9AggXxhTaGVsy7h_EZSChSZGcmgX",
            "type": "SequenceLocation",
        },
        "state": {"sequence": "GAGGCTTT", "type": "LiteralSequenceExpression"},
        "type": "Allele",
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def gnomad_vcf_genomic_delins2():
    """Create a test fixture for 16-68846036-AG-TGAGTTT (CA396459910)"""
    params = {
        "id": "ga4gh:VA.hb5lqeSL-4BGunjhqrL1ePPDPfTCbbIL",
        "location": {
            "id": "ga4gh:SL.zBZ0KyedHbgR5ElXXweOjKDu0bQQ8tVJ",
            "start": 68812132,
            "end": 68812134,
            "sequence": "ga4gh:SQ.yC_0RBj3fgBlvgyAuycbzdubtLxq-rE0",
            "type": "SequenceLocation",
        },
        "state": {"sequence": "TGAGTTT", "type": "LiteralSequenceExpression"},
        "type": "Allele",
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def gnomad_vcf_genomic_delins3():
    """Create a test fixture for X-70350063-AG-AGGCAGCGCATAAAGCGCATTCTCCG"""
    params = {
        "id": "ga4gh:VA.lv2MkUfJ6gGR7bcw8E11khTzsEKo6P8e",
        "location": {
            "id": "ga4gh:SL.iStfHeIOACHNA5MlWX9RbaG-8gjzn6xT",
            "start": 71130213,
            "end": 71130215,
            "sequence": "ga4gh:SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
            "type": "SequenceLocation",
        },
        "state": {
            "length": 26,
            "repeatSubunitLength": 24,
            "sequence": "GGCAGCGCATAAAGCGCATTCTCCGG",
            "type": "ReferenceLengthExpression",
        },
        "type": "Allele",
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def gnomad_vcf_genomic_delins4():
    """Create a test fixture for 1-55509715-AC-A"""
    params = {
        "id": "ga4gh:VA.d5k4EP2zG3gJurtVrKM3I2xRprSCmJZi",
        "location": {
            "id": "ga4gh:SL._-CxWxxGtqykQa4O4D7iZj_PP1SQRN7z",
            "end": 55044045,
            "start": 55044042,
            "sequence": "ga4gh:SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
            "type": "SequenceLocation",
        },
        "state": {
            "length": 2,
            "repeatSubunitLength": 1,
            "sequence": "CC",
            "type": "ReferenceLengthExpression",
        },
        "type": "Allele",
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def gnomad_vcf_genomic_delins5():
    """Create test fixture for 17-7578455-CGCGG-CGCG (CA497925643)"""
    params = {
        "id": "ga4gh:VA.3mnuk0GOSDtlYkaBsz3d5iumacfw9gbf",
        "type": "Allele",
        "location": {
            "id": "ga4gh:SL.jvGY9Lk18DGlztfAkJuIAtQQ9QgyW4lE",
            "type": "SequenceLocation",
            "sequence": "ga4gh:SQ.dLZ15tNO1Ur0IcGjwc3Sdi_0A6Yf4zm7",
            "start": 7675139,
            "end": 7675141,
        },
        "state": {
            "type": "ReferenceLengthExpression",
            "sequence": "G",
            "length": 1,
            "repeatSubunitLength": 1,
        },
    }
    return models.Allele(**params)


@pytest.mark.asyncio
async def test_protein_substitution(test_handler, braf_v600e, dis3_p63a):
    """Test that protein substitutions normalize correctly."""
    resp = await test_handler.normalize("     BRAF      V600E    ")
    assertion_checks(resp, braf_v600e)

    resp = await test_handler.normalize("NP_004324.2:p.Val600Glu")
    assertion_checks(resp, braf_v600e)

    resp = await test_handler.normalize("braf V512E")
    assertion_checks(resp, braf_v600e)

    resp = await test_handler.normalize(" NP_001365404.1:p.Val512Glu  ")
    assertion_checks(resp, braf_v600e)

    resp = await test_handler.normalize("DIS3 P63A")
    assertion_checks(resp, dis3_p63a)


@pytest.mark.asyncio
async def test_polypeptide_truncation(test_handler, vhl):
    """Test that polypeptide truncations normalize correctly."""
    resp = await test_handler.normalize("NP_000542.1:p.Tyr185Ter")
    assertion_checks(resp, vhl)


@pytest.mark.asyncio
async def test_reference_agree(test_handler, vhl_reference_agree):
    """Test that reference agrees normalize correctly."""
    resp = await test_handler.normalize("NP_000542.1:p.Pro61=")
    assertion_checks(resp, vhl_reference_agree)


@pytest.mark.asyncio
async def test_cdna_and_genomic_substitution(
    test_handler,
    braf_v600e_nucleotide,
    genomic_substitution,
    genomic_sub_grch38,
    braf_v600e_genomic_sub,
    gnomad_vcf_genomic_sub_mnv,
):
    """Test that cdna and genomic substitutions normalize correctly."""
    resp = await test_handler.normalize("NM_004333.4:c.1799T>A")
    assertion_checks(resp, braf_v600e_nucleotide)

    # MANE transcript
    resp = await test_handler.normalize("ENST00000288602.10:c.1799T>A")
    assertion_checks(resp, braf_v600e_nucleotide)

    resp = await test_handler.normalize("BRAF V600E c.1799T>A")
    assertion_checks(resp, braf_v600e_nucleotide)

    resp = await test_handler.normalize("BRAF V600E (c.1799T>A)")
    assertion_checks(resp, braf_v600e_nucleotide)

    resp = await test_handler.normalize("BRAF c.1799T>A")
    assertion_checks(resp, braf_v600e_nucleotide)

    resp = await test_handler.normalize("NC_000007.13:g.140453136A>T")
    assertion_checks(resp, braf_v600e_genomic_sub)

    resp = await test_handler.normalize("7-140453136-A-T")  # 37
    assertion_checks(resp, braf_v600e_genomic_sub)

    resp = await test_handler.normalize("7-140753336-A-T")  # 38
    assertion_checks(resp, braf_v600e_genomic_sub)

    resp = await test_handler.normalize("BRAF V600E (g.140453136A>T)")
    assertion_checks(resp, braf_v600e_nucleotide)

    resp = await test_handler.normalize("BRAF g.140453136A>T")
    assertion_checks(resp, braf_v600e_nucleotide)

    # More than 1 gene (EGFR and EGFR-AS1)
    resp = await test_handler.normalize("NC_000007.13:g.55249071C>T")
    assertion_checks(resp, genomic_sub_grch38)

    resp = await test_handler.normalize("EGFR g.55249071C>T")
    assertion_checks(resp, genomic_substitution)

    # MNV genomic substitution (CA009580)
    q = "5-112175770-GGAA-AGAA"
    resp = await test_handler.normalize(q)
    assertion_checks(resp, gnomad_vcf_genomic_sub_mnv)


@pytest.mark.asyncio
async def test_cdna_reference_agree(test_handler, cdna_reference_agree):
    """Test that cdna Reference Agree normalizes correctly."""
    resp = await test_handler.normalize("NM_004333.4:c.1799= ")
    assertion_checks(resp, cdna_reference_agree)

    resp = await test_handler.normalize("ENST00000288602.11:c.1799=")
    assertion_checks(resp, cdna_reference_agree)

    resp = await test_handler.normalize("BRAF    c.1799=")
    assertion_checks(resp, cdna_reference_agree)

    resp = await test_handler.normalize("  BRAF  V600E  c.1799=  ")
    assertion_checks(resp, cdna_reference_agree)


@pytest.mark.asyncio
async def test_genomic_reference_agree(
    test_handler, cdna_reference_agree, grch38_braf_genom_reference_agree
):
    """Test that genomic reference agree normalizes correctly."""
    resp = await test_handler.normalize("NC_000007.13:g.140453136=")
    assertion_checks(
        resp,
        grch38_braf_genom_reference_agree,
    )

    resp = await test_handler.normalize("7-140453136-A-A")
    assertion_checks(resp, grch38_braf_genom_reference_agree)

    resp = await test_handler.normalize("7-140753336-A-A")
    assertion_checks(resp, grch38_braf_genom_reference_agree)

    q = "7-140753336-ACT-ACT"
    resp = await test_handler.normalize(q)
    assertion_checks(resp, grch38_braf_genom_reference_agree)

    resp = await test_handler.normalize("BRAF g.140453136=")
    assertion_checks(resp, cdna_reference_agree)


@pytest.mark.asyncio
async def test_cdna_delins(test_handler, nm_004448_cdna_delins, nm_000551):
    """Test that cdna DelIns normalizes correctly."""
    resp = await test_handler.normalize("    NM_004448.4:c.2326_2327delinsCT    ")
    assertion_checks(
        resp,
        nm_004448_cdna_delins,
    )

    resp = await test_handler.normalize("NM_000551.3:c.615delinsAA")
    assertion_checks(resp, nm_000551)


@pytest.mark.asyncio
async def test_genomic_delins(
    test_handler,
    grch38_genomic_delins1,
    grch38_genomic_delins2,
    genomic_delins_gene,
    gnomad_vcf_genomic_delins1,
    gnomad_vcf_genomic_delins2,
    gnomad_vcf_genomic_delins3,
    gnomad_vcf_genomic_delins4,
    gnomad_vcf_genomic_delins5,
    genomic_del1_lse,
    genomic_del2_lse,
):
    """Test that Genomic DelIns normalizes correctly."""
    resp = await test_handler.normalize("NC_000007.13:g.140453135_140453136delinsAT")
    assertion_checks(resp, grch38_genomic_delins1)

    resp = await test_handler.normalize("NC_000003.12:g.10149938delinsAA")
    assertion_checks(resp, grch38_genomic_delins2)

    q = "3-10149938-C-AA"
    resp = await test_handler.normalize(q)
    assertion_checks(resp, grch38_genomic_delins2)

    q = "BRAF g.140453135_140453136delinsAT"
    resp = await test_handler.normalize(q)
    assertion_checks(resp, genomic_delins_gene)

    # https://reg.clinicalgenome.org/redmine/projects/registry/genboree_registry/
    # allele?hgvsOrDescriptor=NM_000249.3%3Ac.489_498delinsGAGGCTTT
    q = "3-37050340-AAAAGCTTTA-GAGGCTTT"
    resp = await test_handler.normalize(q)
    assertion_checks(resp, gnomad_vcf_genomic_delins1)

    q = "16-68846036-AG-TGAGTTT"
    resp = await test_handler.normalize(q)
    assertion_checks(resp, gnomad_vcf_genomic_delins2)

    # NC_000023.10:g.70350063_70350064delinsAGGCAGCGCATAAAGCGCATTCTCCG
    # NC_000023.10:g.70350063_70350064insGGCAGCGCATAAAGCGCATTCTCC
    # https://reg.clinicalgenome.org/redmine/projects/registry/genboree_registry/
    # allele?hgvsOrDescriptor=NC_000023.11%3Ag.71130213_71130214insGGCAGCGCATAAAGCGCATTCTCC noqa: E501
    q = "X-70350063-AG-AGGCAGCGCATAAAGCGCATTCTCCG"
    resp = await test_handler.normalize(q)
    assertion_checks(resp, gnomad_vcf_genomic_delins3)

    # CA523275412
    q = "1-55509715-AC-A"
    resp = await test_handler.normalize(q)
    assertion_checks(resp, gnomad_vcf_genomic_delins4)

    # CA497925643
    q = "17-7578455-CGCGG-CGCG"
    resp = await test_handler.normalize(q)
    assertion_checks(resp, gnomad_vcf_genomic_delins5)

    q = "3-10146594-AATGTTGACGGACAGCCTAT-A"
    resp = await test_handler.normalize(q)
    assertion_checks(resp, genomic_del2_lse)

    q = "3-10188278-AATGTTGACGGACAGCCTAT-A"
    resp = await test_handler.normalize(q)
    assertion_checks(resp, genomic_del2_lse)

    q = "3-10149810-CT-C"  # 38
    resp = await test_handler.normalize(q)
    assertion_checks(resp, genomic_del1_lse)

    # gnomad should always return lse even if provided other hgvs dup del mode option
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.COPY_NUMBER_COUNT)
    assertion_checks(resp, genomic_del1_lse)


@pytest.mark.asyncio
async def test_protein_delins(test_handler, protein_delins):
    """Test that Amnio Acid DelIns normalizes correctly."""
    resp = await test_handler.normalize("NP_001333827.1:p.Leu747_Thr751delinsPro")
    assertion_checks(resp, protein_delins)

    resp = await test_handler.normalize("EGFR p.Leu747_Thr751delinsPro")
    assertion_checks(resp, protein_delins)

    resp = await test_handler.normalize("EGFR Leu747_Thr751delinsPro")
    assertion_checks(resp, protein_delins)

    resp = await test_handler.normalize("EGFR L747_T751delinsP")
    assertion_checks(resp, protein_delins)


@pytest.mark.asyncio
async def test_protein_deletion(test_handler, protein_deletion_np_range):
    """Test that Protein Deletion normalizes correctly."""
    resp = await test_handler.normalize("NP_004439.2:p.Leu755_Thr759del")
    assertion_checks(resp, protein_deletion_np_range)

    resp = await test_handler.normalize("ERBB2 p.Leu755_Thr759del")
    assertion_checks(resp, protein_deletion_np_range)

    resp = await test_handler.normalize("ERBB2 Leu755_Thr759del")
    assertion_checks(resp, protein_deletion_np_range)

    resp1 = await test_handler.normalize("EGFR L747_T751del")
    resp2 = await test_handler.normalize("EGFR L747_T751delLREAT")
    assert resp1.variation.id == resp2.variation.id

    # incorrect deleted sequence
    resp = await test_handler.normalize("EGFR L747_T751delLREA")
    assert resp.variation is None


@pytest.mark.asyncio
async def test_cdna_deletion(test_handler, cdna_deletion):
    """Test that cdna deletion normalizes correctly."""
    # https://reg.clinicalgenome.org/redmine/projects/registry/genboree_registry/by_caid?caid=CA645372623  # noqa: E501
    q = "NM_004448.3:c.2264_2278delTGAGGGAAAACACAT"
    resp1 = await test_handler.normalize(q)
    assertion_checks(resp1, cdna_deletion)

    # incorrected deleted sequence
    resp = await test_handler.normalize("NM_004448.3:c.2264_2278delTGAGGGAAAACACTA")
    assert resp.variation is None

    resp2 = await test_handler.normalize("NM_004448.3:c.2264_2278del")
    assert resp1.variation.id == resp2.variation.id

    q = "ERBB2 c.2264_2278delTGAGGGAAAACACAT"
    resp = await test_handler.normalize(q)
    assertion_checks(resp, cdna_deletion)


@pytest.mark.asyncio
async def test_genomic_deletion(test_handler, genomic_deletion):
    """Test that genomic deletion normalizes correctly"""
    # CA915940709
    q = "NC_000003.12:g.10146527_10146528del"
    resp1 = await test_handler.normalize(q)
    assertion_checks(resp1, genomic_deletion)

    resp2 = await test_handler.normalize("NC_000003.12:g.10146527_10146528delCT")
    assert resp2.variation.id == resp1.variation.id

    resp3 = await test_handler.normalize("3-10146526-TCT-T")
    assert resp3.variation.id == resp2.variation.id

    # incorrect deleted sequence
    resp = await test_handler.normalize("NC_000003.12:g.10146527_10146528delCC")
    assert resp.variation is None


@pytest.mark.asyncio
async def test_protein_insertion(test_handler, protein_insertion):
    """Test that protein insertion normalizes correctly."""
    resp = await test_handler.normalize("NP_005219.2:p.Asp770_Asn771insGlyLeu")
    assertion_checks(resp, protein_insertion)

    resp = await test_handler.normalize("EGFR D770_N771insGL")
    assertion_checks(resp, protein_insertion)

    resp = await test_handler.normalize("EGFR p.D770_N771insGL")
    assertion_checks(resp, protein_insertion)

    resp = await test_handler.normalize("EGFR Asp770_Asn771insGlyLeu")
    assertion_checks(resp, protein_insertion)

    resp = await test_handler.normalize("EGFR p.Asp770_Asn771insGlyLeu")
    assertion_checks(resp, protein_insertion)


@pytest.mark.asyncio
async def test_cdna_insertion(test_handler, cdna_insertion):
    """Test that cdna insertion normalizes correctly."""
    resp = await test_handler.normalize("ENST00000331728.9:c.2049_2050insA")
    assertion_checks(resp, cdna_insertion)


@pytest.mark.asyncio
async def test_genomic_insertion(
    test_handler, genomic_insertion, grch38_genomic_insertion_variation
):
    """Test that genomic insertion normalizes correctly."""
    resp = await test_handler.normalize(
        "NC_000017.10:g.37880993_37880994insGCTTACGTGATG"
    )
    assertion_checks(resp, grch38_genomic_insertion_variation)

    resp = await test_handler.normalize("ERBB2 g.37880993_37880994insGCTTACGTGATG")
    assertion_checks(resp, genomic_insertion)

    q = "17-37880993-G-GGCTTACGTGATG"
    resp = await test_handler.normalize(q)
    assertion_checks(resp, grch38_genomic_insertion_variation)


@pytest.mark.asyncio
async def test_amplification(test_handler, braf_amplification, prpf8_amplification):
    """Test that amplification normalizes correctly."""
    q = "BRAF Amplification"
    resp = await test_handler.normalize(q)
    assertion_checks(resp, braf_amplification)

    # Gene with > 1 sequence location
    q = "PRPF8 AMPLIFICATION"
    resp = await test_handler.normalize(q)
    assertion_checks(resp, prpf8_amplification)

    # Gene with no location. This should NOT return a variation
    resp = await test_handler.normalize("IFNR amplification")
    assert resp.variation is None


@pytest.mark.asyncio
async def test_valid_queries(test_handler):
    """Test that valid queries don"t throw exceptions. Used for queries that
    revealed bugs in service.
    """
    assert await test_handler.normalize("CCND1 Y44D")

    resp = await test_handler.normalize("NC_000002.12:g.73448098_73448100delCTC")
    assert resp
    assert resp.variation.state.sequence.root == "CTC"
    assert resp.variation.id == "ga4gh:VA.rr4kgycSnU0VyPKjdFetgcvHk2oGPii1"

    # Test ambiguous IUPAC code N
    for q in [
        "NC_000017.10:g.7572948_7572949insTTTTTTTTTNNNNN",
        "NC_000007.13:g.140453136A>N",
        "NC_000007.13:g.140453135_140453136delinsATN",
        "NM_007294.3:c.2902_2903insTCN",
        "NM_004333.4:c.1799T>N",
        "NM_001289937.1:c.2326_2327delinsCTN",
    ]:
        resp = await test_handler.normalize(q)
        assert resp.variation, q


@pytest.mark.asyncio
async def test_no_matches(test_handler):
    """Test no matches work correctly."""
    queries = [
        "braf",  # no change
        "braf v600e",  # incorrect case
        "braf v600000932092039e",  # invalid pos
        "NP_000213.1:cp.Leu862=",  # cp is invalid
        "NP_000213.1:cp.Leu862",  # cp is invalid
        "BRAF V600E 33",  # not supported query type
        "NP_004324.2:p.Glu600Val",  # not valid ref
        "NP_004324.2:p.Glu600Gal",  # not valid ref
        "NP_004324.2839:p.Glu600Val",  # not valid accession
        "NP_004324.2:t.Glu600Val",  # t is invalid
        "this:c.54G>H",  # not a valid accession
        "NC_000007.13:g.4T<A",  # invalid hgvs format
        "test",
        "131",
        "braf Z600E",  # invalid ref
        "braf E600Z",  # invalid ref
        "Thr790Met",  # no gene/accession
        "p.Tyr365Ter",  # no accession
        "ERBB2 G776delinsVCZ",
        "NP005219.2:p.Glu746_Thr751delinsValAla",
        "NP_005219.2:p.Glu746Thr751delinsValAla",
        "EGFR L747_L474delinsP",
        "NP_005219.2:p.Glu746_Thr751delinssValAla",
        "EGFR delins",
        "NM_004333.4:c.1799_1800delTGinsAT",
        "NM_173851.3(SLC30A8):c.973C>T%20(p.Arg325Trp)",
        "NG_008212.3:g.5426_5445del",  # NG accessions not supported
        "NC_000010.11-87925523-C-G",  # invalid format
        "clinvar:10",
        "    ",
        "",
    ]
    for q in queries:
        resp = await test_handler.normalize(q)
        assert resp.variation is None


@pytest.mark.asyncio
async def test_service_meta():
    """Test that service meta info populates correctly."""
    response = await normalize_get_response("BRAF v600e", "default")
    service_meta = response.service_meta_
    assert service_meta.name == "variation-normalizer"
    assert service_meta.version
    assert isinstance(service_meta.response_datetime, datetime)
    assert (
        service_meta.url == "https://github.com/cancervariants/variation-normalization"
    )

    response = await normalize_get_response("this-wont-normalize", "default")
    service_meta = response.service_meta_
    assert service_meta.name == "variation-normalizer"
    assert service_meta.version
    assert isinstance(service_meta.response_datetime, datetime)
    assert (
        service_meta.url == "https://github.com/cancervariants/variation-normalization"
    )

    response = await to_vrs_get_response("this-wont-normalize")
    service_meta = response.service_meta_
    assert service_meta.name == "variation-normalizer"
    assert service_meta.version
    assert isinstance(service_meta.response_datetime, datetime)
    assert (
        service_meta.url == "https://github.com/cancervariants/variation-normalization"
    )
