"""Module for testing gnomad_vcf_to_protein works correctly"""
import pytest
from ga4gh.vrsatile.pydantic.vrsatile_models import VariationDescriptor

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
        "id": "normalize.variation:1-2629397-G-T",
        "type": "VariationDescriptor",
        "variation_id": "ga4gh:VA.9TFIzZ8M_qfiAfZTO1QvZ6toS4I_j4Z9",
        "variation": {
            "_id": "ga4gh:VA.9TFIzZ8M_qfiAfZTO1QvZ6toS4I_j4Z9",
            "location": {
                "_id": "ga4gh:VSL.0h5OViQETz4k5QePhKFiQJE038yy4jpS",
                "interval": {
                    "end": {"value": 30, "type": "Number"},
                    "start": {"value": 29, "type": "Number"},
                    "type": "SequenceInterval",
                },
                "sequence_id": "ga4gh:SQ.iQ8F_pnsiQOLohiV2qh3OWRZiftUt8jZ",
                "type": "SequenceLocation",
            },
            "state": {"sequence": "M", "type": "LiteralSequenceExpression"},
            "type": "Allele",
        },
        "molecule_context": "protein",
        "vrs_ref_allele_seq": "L",
        "gene_context": "hgnc:14668",
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def cdk11a_e314del():
    """Create test fixture for CDK11A Glu314del"""
    params = {
        "id": "normalize.variation:1-1708855-TTCC-T",
        "type": "VariationDescriptor",
        "variation_id": "ga4gh:VA.6GUBUCcortzNbx7Q5DMPzxiFTnvPowhi",
        "variation": {
            "_id": "ga4gh:VA.6GUBUCcortzNbx7Q5DMPzxiFTnvPowhi",
            "location": {
                "_id": "ga4gh:VSL.tnCPFL4wpErBLA6rnIcFrRUGpaO1639_",
                "interval": {
                    "end": {"value": 321, "type": "Number"},
                    "start": {"value": 308, "type": "Number"},
                    "type": "SequenceInterval",
                },
                "sequence_id": "ga4gh:SQ.N728VSRRMHJ1SrhJgKqJOCaa3l5Z4sqm",
                "type": "SequenceLocation",
            },
            "state": {"sequence": "EEEEEEEEEEEE", "type": "LiteralSequenceExpression"},
            "type": "Allele",
        },
        "molecule_context": "protein",
        "vrs_ref_allele_seq": "EEEEEEEEEEEEE",
        "gene_context": "hgnc:1730",
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def protein_insertion2():
    """Create test fixture for LRP8 p.Gln25_Leu26insArg"""
    params = {
        "id": "normalize.variation:1-53327836-A-ACGC",
        "type": "VariationDescriptor",
        "variation_id": "ga4gh:VA.O1MlG0mAJcxpAfd-dRZ_ZXhMdbm3OHHx",
        "variation": {
            "_id": "ga4gh:VA.O1MlG0mAJcxpAfd-dRZ_ZXhMdbm3OHHx",
            "location": {
                "_id": "ga4gh:VSL.5jPV69Go1DH3XWb1bbRpoiz-GeEtIxYv",
                "interval": {
                    "end": {"value": 25, "type": "Number"},
                    "start": {"value": 25, "type": "Number"},
                    "type": "SequenceInterval",
                },
                "sequence_id": "ga4gh:SQ.qgIh8--4F6IpxRwX_lVtD2BhepH5B5Ef",
                "type": "SequenceLocation",
            },
            "state": {"sequence": "R", "type": "LiteralSequenceExpression"},
            "type": "Allele",
        },
        "molecule_context": "protein",
        "gene_context": "hgnc:6700",
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def atad3a_loc():
    """Create test fixture for ATAD3A location"""
    return {
        "_id": "ga4gh:VSL.h7pLvn-VyN4H7GT0vBj6XD5PENOocxOR",
        "interval": {
            "end": {"value": 7, "type": "Number"},
            "start": {"value": 6, "type": "Number"},
            "type": "SequenceInterval",
        },
        "sequence_id": "ga4gh:SQ.MHPOY_7fv8V9SktyvaTxulVFSK6XCxM8",
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def atad3a_i7v(atad3a_loc):
    """Create test fixture for ATAD3A Ile7Val"""
    params = {
        "id": "normalize.variation:1-1512287-A-G",
        "type": "VariationDescriptor",
        "variation_id": "ga4gh:VA.T4StYcC72ctAMe5FGzI9XkyYbLr6VxVk",
        "variation": {
            "_id": "ga4gh:VA.T4StYcC72ctAMe5FGzI9XkyYbLr6VxVk",
            "location": atad3a_loc,
            "state": {"sequence": "V", "type": "LiteralSequenceExpression"},
            "type": "Allele",
        },
        "molecule_context": "protein",
        "vrs_ref_allele_seq": "I",
        "gene_context": "hgnc:25567",
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def atad3a_i7t(atad3a_loc):
    """Create test fixture for ATAD3A Ile7Thr"""
    params = {
        "id": "normalize.variation:1-1512288-T-C",
        "type": "VariationDescriptor",
        "variation_id": "ga4gh:VA.7ohTnL0ULPFaOESzuh_SIRYgp9zGch_Z",
        "variation": {
            "_id": "ga4gh:VA.7ohTnL0ULPFaOESzuh_SIRYgp9zGch_Z",
            "location": atad3a_loc,
            "state": {"sequence": "T", "type": "LiteralSequenceExpression"},
            "type": "Allele",
        },
        "molecule_context": "protein",
        "vrs_ref_allele_seq": "I",
        "gene_context": "hgnc:25567",
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def atad3a_i7m(atad3a_loc):
    """Create test fixture for ATAD3A Ile7Met"""
    params = {
        "id": "normalize.variation:1-1512289-T-G",
        "type": "VariationDescriptor",
        "variation_id": "ga4gh:VA.ua1plOEkzNo23uVhnbSYBSM0jS7to6tN",
        "variation": {
            "_id": "ga4gh:VA.ua1plOEkzNo23uVhnbSYBSM0jS7to6tN",
            "location": atad3a_loc,
            "state": {"sequence": "M", "type": "LiteralSequenceExpression"},
            "type": "Allele",
        },
        "molecule_context": "protein",
        "vrs_ref_allele_seq": "I",
        "gene_context": "hgnc:25567",
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="session")
def braf_v600l(braf_600loc):
    """Create test fixture for BRAF Val600Leu."""
    params = {
        "id": "normalize.variation:7-140753337-C-A",
        "type": "VariationDescriptor",
        "variation_id": "ga4gh:VA.Ktev5asCsmUbHaQG6N-CdSp_g5FyJxLN",
        "variation": {
            "_id": "ga4gh:VA.Ktev5asCsmUbHaQG6N-CdSp_g5FyJxLN",
            "location": braf_600loc,
            "state": {"sequence": "L", "type": "LiteralSequenceExpression"},
            "type": "Allele",
        },
        "molecule_context": "protein",
        "vrs_ref_allele_seq": "V",
        "gene_context": "hgnc:1097",
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="session")
def braf_600_reference_agree(braf_600loc):
    """Create test fixture for BRAF Val600=."""
    params = {
        "id": "normalize.variation:7-140753335-C-A",
        "type": "VariationDescriptor",
        "variation_id": "ga4gh:VA.YUPFxYUZlYnr7Q4nIPiLV5BJznzF3YKl",
        "variation": {
            "_id": "ga4gh:VA.YUPFxYUZlYnr7Q4nIPiLV5BJznzF3YKl",
            "location": braf_600loc,
            "state": {"sequence": "V", "type": "LiteralSequenceExpression"},
            "type": "Allele",
        },
        "molecule_context": "protein",
        "vrs_ref_allele_seq": "V",
        "gene_context": "hgnc:1097",
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def kras_g12d():
    """Fixture for KRAS G12C"""
    return {
        "_id": "ga4gh:VA.NtQTqsdO_Z8G0KpBQ1_z7QsHo_bVN43m",
        "type": "Allele",
        "location": {
            "_id": "ga4gh:VSL.Eiio4mQpHp-kXdQWT_AUHrubE8Q18_br",
            "type": "SequenceLocation",
            "sequence_id": "ga4gh:SQ.fytWhQSNGnA-86vDiQCxTSzgkk_WfQRS",
            "interval": {
                "type": "SequenceInterval",
                "start": {"type": "Number", "value": 11},
                "end": {"type": "Number", "value": 12},
            },
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
    assertion_checks(
        resp.variation_descriptor, braf_v600l, "7-140753337-C-A", ignore_id=True
    )
    assert resp.warnings == []

    # Reading Frame 2, Negative Strand
    resp = await test_handler.gnomad_vcf_to_protein("7-140753336-A-T")
    assertion_checks(
        resp.variation_descriptor, braf_v600e, "7-140753336-A-T", ignore_id=True
    )
    assert resp.warnings == []

    # Reading Frame 3, Negative Strand
    resp = await test_handler.gnomad_vcf_to_protein("7-140753335-C-A")
    assertion_checks(
        resp.variation_descriptor,
        braf_600_reference_agree,
        "7-140753335-C-A",
        ignore_id=True,
    )
    assert resp.warnings == []

    # Reading Frame 3, Negative Strand
    resp = await test_handler.gnomad_vcf_to_protein("1-2629397-G-T")
    assertion_checks(resp.variation_descriptor, mmel1_l30m, "1-2629397-G-T")
    assert resp.warnings == []

    # Reading Frame 1, Positive Strand
    resp = await test_handler.gnomad_vcf_to_protein("1-1512287-A-G")
    assertion_checks(resp.variation_descriptor, atad3a_i7v, "1-1512287-A-G")
    assert resp.warnings == []

    # Reading Frame 2, Positive Strand
    resp = await test_handler.gnomad_vcf_to_protein("1-1512288-T-C")
    assertion_checks(resp.variation_descriptor, atad3a_i7t, "1-1512288-T-C")
    assert resp.warnings == []

    # Reading Frame 3, Positive Strand
    resp = await test_handler.gnomad_vcf_to_protein("1-1512289-T-G")
    assertion_checks(
        resp.variation_descriptor, atad3a_i7m, "1-1512289-T-G", ignore_id=True
    )
    assert resp.warnings == []

    resp = await test_handler.gnomad_vcf_to_protein("12-25245350-C-T")
    assert resp
    variation = resp.variation_descriptor.dict(by_alias=True)["variation"]
    assert variation == kras_g12d


@pytest.mark.asyncio
async def test_reference_agree(test_handler, vhl_reference_agree):
    """Test that reference agree queries return correct response"""
    # https://www.ncbi.nlm.nih.gov/clinvar/variation/379039/?new_evidence=true
    resp = await test_handler.gnomad_vcf_to_protein("3-10142030-C-C")
    assertion_checks(
        resp.variation_descriptor, vhl_reference_agree, "3-10142030-C-C", ignore_id=True
    )
    assert resp.warnings == []


@pytest.mark.asyncio
async def test_insertion(test_handler, protein_insertion, protein_insertion2):
    """Test that insertion queries return correct response"""
    resp = await test_handler.gnomad_vcf_to_protein("7-55181319-C-CGGGTTG")
    assertion_checks(
        resp.variation_descriptor,
        protein_insertion,
        "7-55181319-C-CGGGTTG",
        ignore_id=True,
    )
    assert resp.warnings == []

    resp = await test_handler.gnomad_vcf_to_protein("1-53327836-A-ACGC")
    assertion_checks(resp.variation_descriptor, protein_insertion2, "1-53327836-A-ACGC")
    assert resp.warnings == []


@pytest.mark.asyncio
async def test_deletion(test_handler, protein_deletion_np_range, cdk11a_e314del):
    """Test that deletion queries return correct response"""
    q = "17-39723966-TTGAGGGAAAACACAT-T"
    resp = await test_handler.gnomad_vcf_to_protein(q)
    assertion_checks(
        resp.variation_descriptor, protein_deletion_np_range, q, ignore_id=True
    )
    assert resp.warnings == []

    q = "1-1708855-TTCC-T"
    resp = await test_handler.gnomad_vcf_to_protein(q)
    assertion_checks(resp.variation_descriptor, cdk11a_e314del, q)
    assert resp.warnings == []


@pytest.mark.asyncio
async def test_invalid(test_handler):
    """Test that invalid queries return correct response"""
    resp = await test_handler.gnomad_vcf_to_protein("dummy")
    assert resp.variation_descriptor is None

    resp = await test_handler.gnomad_vcf_to_protein(
        "BRAF V600E", untranslatable_returns_text=True
    )
    assert resp.variation_descriptor.variation.type == "Text"
    assert resp.variation_descriptor.label == "BRAF V600E"
    assert resp.warnings == ["BRAF V600E is not a supported gnomad vcf query"]

    resp = await test_handler.gnomad_vcf_to_protein(
        "7-140753336-T-G", untranslatable_returns_text=True
    )
    assert resp.variation_descriptor.variation.type == "Text"
    assert resp.variation_descriptor.label == "7-140753336-T-G"
    assert set(resp.warnings) == {
        "Expected T but found A on NC_000007.14 at position 140753336"
    }

    resp = await test_handler.gnomad_vcf_to_protein(
        "20-2-TC-TG", untranslatable_returns_text=True
    )
    assert resp.variation_descriptor.variation.type == "Text"
    assert resp.variation_descriptor.label == "20-2-TC-TG"
    assert resp.warnings == ["20-2-TC-TG is not a valid gnomad vcf query"]
