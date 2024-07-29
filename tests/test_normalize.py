"""Module for testing the normalize endpoint."""

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
        "location": {
            "end": 63,
            "start": 62,
            "sequenceReference": {
                "type": "SequenceReference",
                "refgetAccession": "SQ.mlWsxfPKINN3o300stAI8oqN5U7P6kEu",
            },
            "sequence": "P",
            "type": "SequenceLocation",
        },
        "state": {"sequence": "A", "type": "LiteralSequenceExpression"},
        "type": "Allele",
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def tp53_g262c():
    """Create TP53 G262C test fixture."""
    params = {
        "location": {
            "start": 261,
            "end": 262,
            "sequenceReference": {
                "type": "SequenceReference",
                "refgetAccession": "SQ.YIlmVwD0rxIqnlvb-8WujHPbR0j3WEGI",
            },
            "type": "SequenceLocation",
        },
        "state": {"sequence": "C", "type": "LiteralSequenceExpression"},
        "type": "Allele",
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def vhl():
    """Create VHL Tyr185Ter fixture."""
    params = {
        "location": {
            "end": 185,
            "start": 184,
            "sequenceReference": {
                "type": "SequenceReference",
                "refgetAccession": "SQ.z-Oa0pZkJ6GHJHOYM7h5mY_umc0SJzTu",
            },
            "sequence": "Y",
            "type": "SequenceLocation",
        },
        "state": {"sequence": "*", "type": "LiteralSequenceExpression"},
        "type": "Allele",
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def nm_004448_cdna_delins():
    """Create test fixture for NM_004448.4:c.2326_2327delinsCT."""
    params = {
        "location": {
            "end": 2502,
            "start": 2500,
            "sequenceReference": {
                "type": "SequenceReference",
                "refgetAccession": "SQ.y9b4LVMiCXpZxOg9Xt1NwRtssA03MwWM",
            },
            "sequence": "GG",
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
        "location": {
            "end": 685,
            "start": 684,
            "sequenceReference": {
                "type": "SequenceReference",
                "refgetAccession": "SQ.xBKOKptLLDr-k4hTyCetvARn16pDS_rW",
            },
            "sequence": "C",
            "type": "SequenceLocation",
        },
        "state": {"sequence": "AA", "type": "LiteralSequenceExpression"},
        "type": "Allele",
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def braf_cdna_seq_loc():
    """Create test fixture for BRAF V600E cDNA representation sequence location"""
    return {
        "end": 2025,
        "start": 2024,
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.aKMPEJgmlZXt_F6gRY5cUG3THH2n-GUa",
        },
        "sequence": "T",
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def braf_v600e_nucleotide(braf_cdna_seq_loc):
    """Create a test fixture for BRAF V600E MANE select nucleotide hgvs."""
    params = {
        "location": braf_cdna_seq_loc,
        "state": {"sequence": "A", "type": "LiteralSequenceExpression"},
        "type": "Allele",
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def cdna_reference_agree(braf_cdna_seq_loc):
    """Create test fixture for NM_004333.4:c.1799=."""
    params = {
        "location": braf_cdna_seq_loc,
        "state": {"sequence": "T", "type": "LiteralSequenceExpression"},
        "type": "Allele",
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def protein_delins():
    """Create test fixture for protein delins."""
    params = {
        "location": {
            "end": 751,
            "start": 746,
            "sequenceReference": {
                "type": "SequenceReference",
                "refgetAccession": "SQ.vyo55F6mA6n2LgN4cagcdRzOuh38V4mE",
            },
            "sequence": "LREAT",
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
        "location": {
            "end": 2453,
            "start": 2437,
            "sequenceReference": {
                "type": "SequenceReference",
                "refgetAccession": "SQ.y9b4LVMiCXpZxOg9Xt1NwRtssA03MwWM",
            },
            "sequence": "TTGAGGGAAAACACAT",
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
        "location": {
            "end": 10146528,
            "start": 10146524,
            "sequenceReference": {
                "type": "SequenceReference",
                "refgetAccession": "SQ.Zu7h9AggXxhTaGVsy7h_EZSChSZGcmgX",
            },
            "sequence": "CTCT",
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
        "location": {
            "end": 2160,
            "start": 2160,
            "sequenceReference": {
                "type": "SequenceReference",
                "refgetAccession": "SQ.7_mlQyDN-uWH0RlxTQFvFEv6ykd2D-xF",
            },
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
        "location": {
            "end": 2500,
            "start": 2488,
            "sequenceReference": {
                "type": "SequenceReference",
                "refgetAccession": "SQ.y9b4LVMiCXpZxOg9Xt1NwRtssA03MwWM",
            },
            "sequence": "TACGTGATGGCT",
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
        "location": {
            "end": 2630,
            "start": 2629,
            "sequenceReference": {
                "type": "SequenceReference",
                "refgetAccession": "SQ.d_QsP29RWJi6bac7GOC9cJ9AO7s_HUMN",
            },
            "sequence": "C",
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
        "location": {
            "end": 112840073,
            "start": 112840072,
            "sequenceReference": {
                "type": "SequenceReference",
                "refgetAccession": "SQ.aUiQCzCPZ2d0csHbMSbh2NzInhonSXwI",
            },
            "sequence": "G",
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
        "location": {
            "end": 55181378,
            "start": 55181377,
            "sequenceReference": {
                "type": "SequenceReference",
                "refgetAccession": "SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
            },
            "sequence": "C",
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
        "location": {
            "end": 140753336,
            "start": 140753335,
            "sequenceReference": {
                "type": "SequenceReference",
                "refgetAccession": "SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
            },
            "sequence": "A",
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
        "location": {
            "end": 140753336,
            "start": 140753334,
            "sequenceReference": {
                "type": "SequenceReference",
                "refgetAccession": "SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
            },
            "sequence": "CA",
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
        "location": {
            "start": 10149937,
            "end": 10149938,
            "sequenceReference": {
                "type": "SequenceReference",
                "refgetAccession": "SQ.Zu7h9AggXxhTaGVsy7h_EZSChSZGcmgX",
            },
            "sequence": "C",
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
        "location": {
            "start": 2024,
            "end": 2026,
            "sequenceReference": {
                "type": "SequenceReference",
                "refgetAccession": "SQ.aKMPEJgmlZXt_F6gRY5cUG3THH2n-GUa",
            },
            "sequence": "TG",
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
        "location": {
            "start": 37008848,
            "end": 37008858,
            "sequenceReference": {
                "type": "SequenceReference",
                "refgetAccession": "SQ.Zu7h9AggXxhTaGVsy7h_EZSChSZGcmgX",
            },
            "sequence": "AAAAGCTTTA",
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
        "location": {
            "start": 68812132,
            "end": 68812134,
            "sequenceReference": {
                "type": "SequenceReference",
                "refgetAccession": "SQ.yC_0RBj3fgBlvgyAuycbzdubtLxq-rE0",
            },
            "sequence": "AG",
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
        "location": {
            "start": 71130213,
            "end": 71130215,
            "sequenceReference": {
                "type": "SequenceReference",
                "refgetAccession": "SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
            },
            "sequence": "GG",
            "type": "SequenceLocation",
        },
        "state": {
            "sequence": "GGCAGCGCATAAAGCGCATTCTCCGG",
            "type": "LiteralSequenceExpression",
        },
        "type": "Allele",
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def gnomad_vcf_genomic_delins4():
    """Create a test fixture for 1-55509715-AC-A"""
    params = {
        "location": {
            "end": 55044045,
            "start": 55044042,
            "sequenceReference": {
                "type": "SequenceReference",
                "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
            },
            "sequence": "CCC",
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
        "type": "Allele",
        "location": {
            "type": "SequenceLocation",
            "sequenceReference": {
                "type": "SequenceReference",
                "refgetAccession": "SQ.dLZ15tNO1Ur0IcGjwc3Sdi_0A6Yf4zm7",
            },
            "start": 7675139,
            "end": 7675141,
            "sequence": "GG",
        },
        "state": {
            "type": "ReferenceLengthExpression",
            "sequence": "G",
            "length": 1,
            "repeatSubunitLength": 1,
        },
    }
    return models.Allele(**params)


@pytest.mark.asyncio()
async def test_protein_substitution(test_handler, braf_v600e, dis3_p63a):
    """Test that protein substitutions normalize correctly."""
    resp = await test_handler.normalize("     BRAF      V600E    ")
    assertion_checks(resp, braf_v600e, check_vrs_id=True)

    resp = await test_handler.normalize("NP_004324.2:p.Val600Glu")
    assertion_checks(resp, braf_v600e)

    resp = await test_handler.normalize("braf V512E")
    assertion_checks(resp, braf_v600e)

    resp = await test_handler.normalize(" NP_001365404.1:p.Val512Glu  ")
    assertion_checks(resp, braf_v600e)

    resp = await test_handler.normalize("DIS3 P63A")
    assertion_checks(resp, dis3_p63a)


@pytest.mark.asyncio()
async def test_polypeptide_truncation(test_handler, vhl):
    """Test that polypeptide truncations normalize correctly."""
    resp = await test_handler.normalize("NP_000542.1:p.Tyr185Ter")
    assertion_checks(resp, vhl)


@pytest.mark.asyncio()
async def test_reference_agree(test_handler, vhl_reference_agree):
    """Test that reference agrees normalize correctly."""
    resp = await test_handler.normalize("NP_000542.1:p.Pro61=")
    assertion_checks(resp, vhl_reference_agree)


@pytest.mark.asyncio()
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
    assertion_checks(resp, braf_v600e_genomic_sub, mane_genes_exts=True)

    resp = await test_handler.normalize("7-140453136-A-T")  # 37
    assertion_checks(resp, braf_v600e_genomic_sub, mane_genes_exts=True)

    resp = await test_handler.normalize("7-140753336-A-T")  # 38
    assertion_checks(resp, braf_v600e_genomic_sub, mane_genes_exts=True)

    resp = await test_handler.normalize("BRAF V600E (g.140453136A>T)")
    assertion_checks(resp, braf_v600e_nucleotide)

    resp = await test_handler.normalize("BRAF g.140453136A>T")
    assertion_checks(resp, braf_v600e_nucleotide)

    # More than 1 gene (EGFR and EGFR-AS1)
    resp = await test_handler.normalize("NC_000007.13:g.55249071C>T")
    assertion_checks(resp, genomic_sub_grch38, mane_genes_exts=True)

    resp = await test_handler.normalize("EGFR g.55249071C>T")
    assertion_checks(resp, genomic_substitution)

    # MNV genomic substitution (CA009580)
    q = "5-112175770-GGAA-AGAA"
    resp = await test_handler.normalize(q)
    assertion_checks(resp, gnomad_vcf_genomic_sub_mnv, mane_genes_exts=True)


@pytest.mark.asyncio()
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


@pytest.mark.asyncio()
async def test_genomic_reference_agree(
    test_handler, cdna_reference_agree, grch38_braf_genom_reference_agree
):
    """Test that genomic reference agree normalizes correctly."""
    resp = await test_handler.normalize("NC_000007.13:g.140453136=")
    assertion_checks(resp, grch38_braf_genom_reference_agree, mane_genes_exts=True)

    resp = await test_handler.normalize("7-140453136-A-A")
    assertion_checks(resp, grch38_braf_genom_reference_agree, mane_genes_exts=True)

    resp = await test_handler.normalize("7-140753336-A-A")
    assertion_checks(resp, grch38_braf_genom_reference_agree, mane_genes_exts=True)

    q = "7-140753336-ACT-ACT"
    resp = await test_handler.normalize(q)
    assertion_checks(resp, grch38_braf_genom_reference_agree, mane_genes_exts=True)

    resp = await test_handler.normalize("BRAF g.140453136=")
    assertion_checks(resp, cdna_reference_agree)


@pytest.mark.asyncio()
async def test_cdna_delins(test_handler, nm_004448_cdna_delins, nm_000551):
    """Test that cdna DelIns normalizes correctly."""
    resp = await test_handler.normalize("    NM_004448.4:c.2326_2327delinsCT    ")
    assertion_checks(
        resp,
        nm_004448_cdna_delins,
    )

    resp = await test_handler.normalize("NM_000551.3:c.615delinsAA")
    assertion_checks(resp, nm_000551)


@pytest.mark.asyncio()
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
    assertion_checks(resp, grch38_genomic_delins1, mane_genes_exts=True)

    resp = await test_handler.normalize("NC_000003.12:g.10149938delinsAA")
    assertion_checks(resp, grch38_genomic_delins2, mane_genes_exts=True)

    q = "3-10149938-C-AA"
    resp = await test_handler.normalize(q)
    assertion_checks(resp, grch38_genomic_delins2, mane_genes_exts=True)

    q = "BRAF g.140453135_140453136delinsAT"
    resp = await test_handler.normalize(q)
    assertion_checks(resp, genomic_delins_gene)

    # https://reg.clinicalgenome.org/redmine/projects/registry/genboree_registry/
    # allele?hgvsOrDescriptor=NM_000249.3%3Ac.489_498delinsGAGGCTTT
    q = "3-37050340-AAAAGCTTTA-GAGGCTTT"
    resp = await test_handler.normalize(q)
    assertion_checks(resp, gnomad_vcf_genomic_delins1, mane_genes_exts=True)

    q = "16-68846036-AG-TGAGTTT"
    resp = await test_handler.normalize(q)
    assertion_checks(resp, gnomad_vcf_genomic_delins2, mane_genes_exts=True)

    # NC_000023.10:g.70350063_70350064delinsAGGCAGCGCATAAAGCGCATTCTCCG
    # NC_000023.10:g.70350063_70350064insGGCAGCGCATAAAGCGCATTCTCC
    # https://reg.clinicalgenome.org/redmine/projects/registry/genboree_registry/
    # allele?hgvsOrDescriptor=NC_000023.11%3Ag.71130213_71130214insGGCAGCGCATAAAGCGCATTCTCC noqa: E501
    q = "X-70350063-AG-AGGCAGCGCATAAAGCGCATTCTCCG"
    resp = await test_handler.normalize(q)
    assertion_checks(resp, gnomad_vcf_genomic_delins3, mane_genes_exts=True)

    # CA523275412
    q = "1-55509715-AC-A"
    resp = await test_handler.normalize(q)
    assertion_checks(resp, gnomad_vcf_genomic_delins4, mane_genes_exts=True)

    # CA497925643
    q = "17-7578455-CGCGG-CGCG"
    resp = await test_handler.normalize(q)
    assertion_checks(resp, gnomad_vcf_genomic_delins5, mane_genes_exts=True)

    q = "3-10146594-AATGTTGACGGACAGCCTAT-A"
    resp = await test_handler.normalize(q)
    assertion_checks(resp, genomic_del2_lse, mane_genes_exts=True)

    q = "3-10188278-AATGTTGACGGACAGCCTAT-A"
    resp = await test_handler.normalize(q)
    assertion_checks(resp, genomic_del2_lse, mane_genes_exts=True)

    q = "3-10149810-CT-C"  # 38
    resp = await test_handler.normalize(q)
    assertion_checks(resp, genomic_del1_lse, mane_genes_exts=True)

    # gnomad should always return lse even if provided other hgvs dup del mode option
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.COPY_NUMBER_COUNT)
    assertion_checks(resp, genomic_del1_lse, mane_genes_exts=True)


@pytest.mark.asyncio()
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


@pytest.mark.asyncio()
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


@pytest.mark.asyncio()
async def test_cdna_deletion(test_handler, cdna_deletion):
    """Test that cdna deletion normalizes correctly."""
    # https://reg.clinicalgenome.org/redmine/projects/registry/genboree_registry/by_caid?caid=CA645372623
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


@pytest.mark.asyncio()
async def test_genomic_deletion(test_handler, genomic_deletion):
    """Test that genomic deletion normalizes correctly"""
    # CA915940709
    q = "NC_000003.12:g.10146527_10146528del"
    resp1 = await test_handler.normalize(q)
    assertion_checks(resp1, genomic_deletion, mane_genes_exts=True)

    resp2 = await test_handler.normalize("NC_000003.12:g.10146527_10146528delCT")
    assert resp2.variation.id == resp1.variation.id

    resp3 = await test_handler.normalize("3-10146526-TCT-T")
    assert resp3.variation.id == resp2.variation.id

    # incorrect deleted sequence
    resp = await test_handler.normalize("NC_000003.12:g.10146527_10146528delCC")
    assert resp.variation is None


@pytest.mark.asyncio()
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


@pytest.mark.asyncio()
async def test_cdna_insertion(test_handler, cdna_insertion):
    """Test that cdna insertion normalizes correctly."""
    resp = await test_handler.normalize("ENST00000331728.9:c.2049_2050insA")
    assertion_checks(resp, cdna_insertion)


@pytest.mark.asyncio()
async def test_genomic_insertion(
    test_handler, genomic_insertion, grch38_genomic_insertion_variation
):
    """Test that genomic insertion normalizes correctly."""
    resp = await test_handler.normalize(
        "NC_000017.10:g.37880993_37880994insGCTTACGTGATG"
    )
    assertion_checks(resp, grch38_genomic_insertion_variation, mane_genes_exts=True)

    resp = await test_handler.normalize("ERBB2 g.37880993_37880994insGCTTACGTGATG")
    assertion_checks(resp, genomic_insertion)

    q = "17-37880993-G-GGCTTACGTGATG"
    resp = await test_handler.normalize(q)
    assertion_checks(resp, grch38_genomic_insertion_variation, mane_genes_exts=True)


@pytest.mark.asyncio()
async def test_amplification(test_handler, braf_amplification, prpf8_amplification):
    """Test that amplification normalizes correctly."""
    q = "BRAF Amplification"
    resp = await test_handler.normalize(q)
    assertion_checks(resp, braf_amplification, check_vrs_id=True)

    # Gene with > 1 sequence location
    q = "PRPF8 AMPLIFICATION"
    resp = await test_handler.normalize(q)
    assertion_checks(resp, prpf8_amplification)

    # Gene with no location. This should NOT return a variation
    resp = await test_handler.normalize("IFNR amplification")
    assert resp.variation is None


@pytest.mark.asyncio()
async def test_valid_queries(test_handler):
    """Test that valid queries don"t throw exceptions. Used for queries that
    revealed bugs in service.
    """
    assert await test_handler.normalize("CCND1 Y44D")

    resp = await test_handler.normalize("NC_000002.12:g.73448098_73448100delCTC")
    assert resp
    assert resp.variation.state.sequence.root == "CTC"
    assert resp.variation.id.startswith("ga4gh:VA.")

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


@pytest.mark.asyncio()
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


@pytest.mark.asyncio()
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
