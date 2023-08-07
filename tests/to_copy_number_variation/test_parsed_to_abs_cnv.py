"""Test that parsed_to_cn_var works correctly"""
import pytest
from ga4gh.vrsatile.pydantic.vrs_models import CopyNumberCount

from variation.schemas.service_schema import ClinVarAssembly


@pytest.fixture(scope="module")
def copy_number_gain1():
    """Create test fixture for clinvar copy number gain.
    https://www.ncbi.nlm.nih.gov/clinvar/variation/145208/?new_evidence=true
    """
    variation = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.N6C9rWBjrNuiIhJkPxdPlRKvSGKoFynr",
        "subject": {
            "type": "SequenceLocation",
            "_id": "ga4gh:VSL.JTsxd9PiPZaIPL9Tl3ss78GYYnDeogvf",
            "sequence_id": "ga4gh:SQ.S_KjnFVz-FE7M0W6yoaUDgYxLPc1jyWU",
            "interval": {
                "type": "SequenceInterval",
                "start": {
                    "type": "IndefiniteRange",
                    "value": 143134062,
                    "comparator": "<=",
                },
                "end": {
                    "type": "IndefiniteRange",
                    "value": 143284670,
                    "comparator": ">=",
                },
            },
        },
        "copies": {"type": "Number", "value": 3},
    }
    return CopyNumberCount(**variation)


@pytest.fixture(scope="module")
def copy_number_gain2():
    """Create test fixture for clinvar copy number gain.
    https://www.ncbi.nlm.nih.gov/clinvar/variation/146181/?new_evidence=true
    """
    variation = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.xOEIBXGfoM8TUA2RKNWINRze_hWT1lPP",
        "subject": {
            "type": "SequenceLocation",
            "_id": "ga4gh:VSL.9moblqAMqfEryr9pRUxqZMiOkqbsy5Ml",
            "sequence_id": "ga4gh:SQ.AsXvWL1-2i5U_buw6_niVIxD6zTbAuS6",
            "interval": {
                "type": "SequenceInterval",
                "start": {
                    "type": "IndefiniteRange",
                    "value": 31738808,
                    "comparator": "<=",
                },
                "end": {
                    "type": "IndefiniteRange",
                    "value": 32217725,
                    "comparator": ">=",
                },
            },
        },
        "copies": {"type": "Number", "value": 2},
    }
    return CopyNumberCount(**variation)


@pytest.fixture(scope="module")
def copy_number_loss1():
    """Create test fixture for clinvar copy number loss.
    https://www.ncbi.nlm.nih.gov/clinvar/variation/146181/?new_evidence=true
    """
    variation = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.IhWQwwhYFtjQAG7BejsVYy-KiM4RMyed",
        "subject": {
            "type": "SequenceLocation",
            "_id": "ga4gh:VSL.Szlw1t4YMuaO7lLwFJ-T7fGTcXuhNNKB",
            "sequence_id": "ga4gh:SQ.v7noePfnNpK8ghYXEqZ9NukMXW7YeNsm",
            "interval": {
                "type": "SequenceInterval",
                "start": {
                    "type": "IndefiniteRange",
                    "value": 10491131,
                    "comparator": "<=",
                },
                "end": {
                    "type": "IndefiniteRange",
                    "value": 10535643,
                    "comparator": ">=",
                },
            },
        },
        "copies": {"type": "Number", "value": 1},
    }
    return CopyNumberCount(**variation)


@pytest.fixture(scope="module")
def copy_number_loss2():
    """Create test fixture for clinvar copy number loss.
    https://www.ncbi.nlm.nih.gov/clinvar/variation/148425/?new_evidence=true
    """
    variation = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.sxRwv8l26F1PdcovpJR5HEpgFOY8J95Q",
        "subject": {
            "type": "SequenceLocation",
            "_id": "ga4gh:VSL.Bp-86GeYti1DBmrj_Dtz7qNIMF5ygx5y",
            "sequence_id": "ga4gh:SQ.8_liLu1aycC0tPQPFmUaGXJLDs5SbPZ5",
            "interval": {
                "type": "SequenceInterval",
                "start": {
                    "type": "IndefiniteRange",
                    "value": 10000,
                    "comparator": "<=",
                },
                "end": {
                    "type": "IndefiniteRange",
                    "value": 1223133,
                    "comparator": ">=",
                },
            },
        },
        "copies": {"type": "Number", "value": 0},
    }
    return CopyNumberCount(**variation)


def test_parsed_copy_number_gain(
    test_cnv_handler, copy_number_gain1, copy_number_gain2
):
    """Test that parsed_to_cn_var works for parsed copy number gain queries"""
    # https://www.ncbi.nlm.nih.gov/clinvar/variation/145208/?new_evidence=true
    resp = test_cnv_handler.parsed_to_cn_var(
        143134063, 143284670, 3, assembly=ClinVarAssembly.GRCH37, chr="chr1"
    )
    assert resp.copy_number_count.dict() == copy_number_gain1.dict()
    assert resp.warnings == []

    resp = test_cnv_handler.parsed_to_cn_var(
        143134063, 143284670, 3, assembly=ClinVarAssembly.HG19, chr="chr1"
    )
    assert resp.copy_number_count.dict() == copy_number_gain1.dict()
    assert resp.warnings == []

    resp = test_cnv_handler.parsed_to_cn_var(
        143134063, 143284670, 3, accession="NC_000001.10"
    )
    assert resp.copy_number_count.dict() == copy_number_gain1.dict()
    assert resp.warnings == []

    # https://www.ncbi.nlm.nih.gov/clinvar/variation/146181/?new_evidence=true
    resp = test_cnv_handler.parsed_to_cn_var(
        31738809, 32217725, 2, assembly=ClinVarAssembly.GRCH38, chr="chr15"
    )
    assert resp.copy_number_count.dict() == copy_number_gain2.dict()
    assert resp.warnings == []

    resp = test_cnv_handler.parsed_to_cn_var(
        31738809, 32217725, 2, assembly=ClinVarAssembly.GRCH38, chr="15"
    )
    assert resp.copy_number_count.dict() == copy_number_gain2.dict()
    assert resp.warnings == []

    resp = test_cnv_handler.parsed_to_cn_var(
        31738809, 32217725, 2, assembly=ClinVarAssembly.HG38, chr="chr15"
    )
    assert resp.copy_number_count.dict() == copy_number_gain2.dict()
    assert resp.warnings == []

    resp = test_cnv_handler.parsed_to_cn_var(
        31738809, 32217725, 2, accession="NC_000015.10"
    )
    assert resp.copy_number_count.dict() == copy_number_gain2.dict()
    assert resp.warnings == []


def test_parsed_copy_number_loss(
    test_cnv_handler, copy_number_loss1, copy_number_loss2
):
    """Test that parsed_to_cn_var works for parsed copy number loss queries"""
    # https://www.ncbi.nlm.nih.gov/clinvar/variation/1299222/?new_evidence=true
    resp = test_cnv_handler.parsed_to_cn_var(
        10491132, 10535643, 1, assembly=ClinVarAssembly.GRCH37, chr="chrX"
    )
    assert resp.copy_number_count.dict() == copy_number_loss1.dict()
    assert resp.warnings == []

    resp = test_cnv_handler.parsed_to_cn_var(
        10491132, 10535643, 1, assembly=ClinVarAssembly.HG19, chr="chrX"
    )
    assert resp.copy_number_count.dict() == copy_number_loss1.dict()
    assert resp.warnings == []

    resp = test_cnv_handler.parsed_to_cn_var(
        10491132, 10535643, 1, assembly=ClinVarAssembly.HG19, chr="X"
    )
    assert resp.copy_number_count.dict() == copy_number_loss1.dict()
    assert resp.warnings == []

    resp = test_cnv_handler.parsed_to_cn_var(
        10491132, 10535643, 1, accession="NC_000023.10"
    )
    assert resp.copy_number_count.dict() == copy_number_loss1.dict()
    assert resp.warnings == []

    # https://www.ncbi.nlm.nih.gov/clinvar/variation/148425/?new_evidence=true
    resp = test_cnv_handler.parsed_to_cn_var(
        10001, 1223133, 0, assembly=ClinVarAssembly.GRCH38, chr="chrY"
    )
    assert resp.copy_number_count.dict() == copy_number_loss2.dict()
    assert resp.warnings == []

    resp = test_cnv_handler.parsed_to_cn_var(
        10001, 1223133, 0, assembly=ClinVarAssembly.HG38, chr="chrY"
    )
    assert resp.copy_number_count.dict() == copy_number_loss2.dict()
    assert resp.warnings == []

    resp = test_cnv_handler.parsed_to_cn_var(
        10001, 1223133, 0, accession="NC_000024.10"
    )
    assert resp.copy_number_count.dict() == copy_number_loss2.dict()
    assert resp.warnings == []


def test_invalid(test_cnv_handler):
    """Test invalid queries returns Text variation and warnings"""
    # NCBI36/hg18 assembly
    # https://www.ncbi.nlm.nih.gov/clinvar/variation/443961/?new_evidence=true
    expected_w = ["NCBI36 assembly is not currently supported"]
    resp = test_cnv_handler.parsed_to_cn_var(
        2623228,
        3150942,
        3,
        assembly=ClinVarAssembly.NCBI36,
        chr="chr1",
        untranslatable_returns_text=True,
    )
    assert resp.copy_number_count.type == "Text"
    assert resp.warnings == expected_w

    resp = test_cnv_handler.parsed_to_cn_var(
        2623228,
        3150942,
        3,
        assembly=ClinVarAssembly.HG18,
        chr="chr1",
        untranslatable_returns_text=True,
    )
    assert resp.copy_number_count.type == "Text"
    assert resp.warnings == expected_w

    # Must give both assembly + chr or accession
    expected_w = ["Must provide either `accession` or both `assembly` and `chr`."]
    resp = test_cnv_handler.parsed_to_cn_var(
        31738809,
        32217725,
        2,
        assembly=ClinVarAssembly.HG38,
        untranslatable_returns_text=True,
    )
    assert resp.copy_number_count.type == "Text"
    assert resp.warnings == expected_w

    resp = test_cnv_handler.parsed_to_cn_var(
        31738809, 32217725, 2, chr="chr15", untranslatable_returns_text=True
    )
    assert resp.copy_number_count.type == "Text"
    assert resp.warnings == expected_w

    resp = test_cnv_handler.parsed_to_cn_var(
        31738809, 32217725, 2, untranslatable_returns_text=True
    )
    assert resp.copy_number_count.type == "Text"
    assert resp.warnings == expected_w

    # invalid chr
    resp = test_cnv_handler.parsed_to_cn_var(
        10001, 1223133, 0, assembly=ClinVarAssembly.GRCH38, chr="z"
    )
    assert resp.copy_number_count is None
    assert resp.warnings == [
        "SeqRepo unable to get translated identifiers for GRCh38:z"
    ]

    # invalid assembly
    resp = test_cnv_handler.parsed_to_cn_var(10001, 1223133, 0, assembly="GRCh99")
    assert resp.copy_number_count is None
    assert resp.warnings

    # invalid accession
    resp = test_cnv_handler.parsed_to_cn_var(
        10491132, 10535643, 1, accession="NC_00002310"
    )
    assert resp.copy_number_count is None
    assert resp.warnings == [
        "SeqRepo unable to get translated identifiers for NC_00002310"
    ]

    # Invalid position
    resp = test_cnv_handler.parsed_to_cn_var(
        31738809, 2302991250, 2, accession="NC_000015.10"
    )
    assert resp.copy_number_count is None
    assert resp.warnings == ["Position out of range (2302991250)"]
