"""Test that parsed_to_copy_number works correctly"""
import pytest
from ga4gh.vrsatile.pydantic.vrs_models import (
    CopyNumberCount, CopyNumberChange, CopyChange, VRSTypes
)

from variation.schemas.service_schema import ClinVarAssembly


@pytest.fixture(scope="module")
def cn_gain1():
    """Create test fixture for clinvar copy number gain.
    https://www.ncbi.nlm.nih.gov/clinvar/variation/145208/?new_evidence=true
    """
    variation = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN._IYaKE4CoDa01tkcgOuqPhnYbZ5RuPcj",
        "subject": {
            "type": "SequenceLocation",
            "id": "ga4gh:SL.RIgksXkT_kWCJv3poK4WQ9PK5_YSRBuh",
            "sequence_id": "ga4gh:SQ.S_KjnFVz-FE7M0W6yoaUDgYxLPc1jyWU",
            "start": {
                "type": "IndefiniteRange",
                "value": 143134062,
                "comparator": "<="
            },
            "end": {
                "type": "IndefiniteRange",
                "value": 143284670,
                "comparator": ">="
            }
        },
        "copies": {"type": "Number", "value": 3}
    }
    return CopyNumberCount(**variation)


@pytest.fixture(scope="module")
def cn_gain2():
    """Create test fixture for clinvar copy number gain.
    https://www.ncbi.nlm.nih.gov/clinvar/variation/146181/?new_evidence=true
    """
    variation = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.CrBtCqTLNhQEHVNxYeC1orEMs_jYki8-",
        "subject": {
            "type": "SequenceLocation",
            "id": "ga4gh:SL.xBw2DuIlDmgOswPardadgzXmVdJmrLZF",
            "sequence_id": "ga4gh:SQ.AsXvWL1-2i5U_buw6_niVIxD6zTbAuS6",
            "start": {
                "type": "IndefiniteRange",
                "value": 31738808,
                "comparator": "<="
            },
            "end": {
                "type": "IndefiniteRange",
                "value": 32217725,
                "comparator": ">="
            }
        },
        "copies": {"type": "Number", "value": 2}
    }
    return CopyNumberCount(**variation)


@pytest.fixture(scope="module")
def cn_loss1():
    """Create test fixture for clinvar copy number loss.
    https://www.ncbi.nlm.nih.gov/clinvar/variation/146181/?new_evidence=true
    """
    variation = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.g1agp3no_7j9ekdTmrpmijhcfC7qvrPC",
        "subject": {
            "type": "SequenceLocation",
            "id": "ga4gh:SL.aIbLtN8y8GbKG91n5_sJR3f2dVUrydeD",
            "sequence_id": "ga4gh:SQ.v7noePfnNpK8ghYXEqZ9NukMXW7YeNsm",
            "start": {
                "type": "IndefiniteRange",
                "value": 10491131,
                "comparator": "<="
            },
            "end": {
                "type": "IndefiniteRange",
                "value": 10535643,
                "comparator": ">="
            }
        },
        "copies": {"type": "Number", "value": 1}
    }
    return CopyNumberCount(**variation)


@pytest.fixture(scope="module")
def cn_loss2():
    """Create test fixture for clinvar copy number loss.
    https://www.ncbi.nlm.nih.gov/clinvar/variation/148425/?new_evidence=true
    """
    variation = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.DwoLZiCB408oqbW2gZNTpvME14tVwNHX",
        "subject": {
            "type": "SequenceLocation",
            "id": "ga4gh:SL.n35LAVyENzd6BVllBUUxzI2p5TgJpKAl",
            "sequence_id": "ga4gh:SQ.8_liLu1aycC0tPQPFmUaGXJLDs5SbPZ5",
            "start": {
                "type": "IndefiniteRange",
                "value": 10000,
                "comparator": "<="
            },
            "end": {
                "type": "IndefiniteRange",
                "value": 1223133,
                "comparator": ">="
            }
        },
        "copies": {"type": "Number", "value": 0}
    }
    return CopyNumberCount(**variation)


@pytest.fixture(scope="module")
def cx_numbers():
    """Create test fixture for copy number change using numbers for start and end"""
    variation = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.UirzxujWnAIklYHh4VxSnFglfDROHYv6",
        "subject": {
            "type": "SequenceLocation",
            "id": "ga4gh:SL.x075Sp6tCfGZcpHHmJ1e5oUdAW0CvN0X",
            "sequence_id": "ga4gh:SQ.8_liLu1aycC0tPQPFmUaGXJLDs5SbPZ5",
            "start": {
                "type": "Number",
                "value": 10000
            },
            "end": {
                "type": "Number",
                "value": 1223133
            }
        },
        "copy_change": "efo:0030069"
    }
    return CopyNumberChange(**variation)


@pytest.fixture(scope="module")
def cx_definite_ranges():
    """Create test fixture for copy number change using definite ranges for start and
    end
    """
    variation = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.Z_1wvmOraI1iAxCKHseVhejBGKI27Yw1",
        "subject": {
            "type": "SequenceLocation",
            "id": "ga4gh:SL.VLcZT6AA5puGaCd45BPUiRJBwrGXx9mF",
            "sequence_id": "ga4gh:SQ.8_liLu1aycC0tPQPFmUaGXJLDs5SbPZ5",
            "start": {
                "type": "DefiniteRange",
                "min": 10000,
                "max": 10005
            },
            "end": {
                "type": "DefiniteRange",
                "min": 1223131,
                "max": 1223134
            }
        },
        "copy_change": "efo:0030069"
    }
    return CopyNumberChange(**variation)


@pytest.fixture(scope="module")
def cx_indefinite_ranges():
    """Create test fixture for copy number change using indefinite ranges for start and
    end
    """
    variation = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.AwFAe7hc5-ACPE6jeJpx7qvk25pw9tmK",
        "subject": {
            "type": "SequenceLocation",
            "id": "ga4gh:SL.xKiefTSlYZBCI0621dLUMJib8XmFzzrC",
            "sequence_id": "ga4gh:SQ.8_liLu1aycC0tPQPFmUaGXJLDs5SbPZ5",
            "start": {
                "type": "IndefiniteRange",
                "comparator": "<=",
                "value": 10000
            },
            "end": {
                "type": "IndefiniteRange",
                "comparator": ">=",
                "value": 1223130
            }
        },
        "copy_change": "efo:0030069"
    }
    return CopyNumberChange(**variation)


@pytest.fixture(scope="module")
def cx_number_indefinite():
    """Create test fixture for copy number change using number for start and indefinite
    range for end
    """
    variation = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.CA02DMUPXNWDbj47UPVyrTzcNLNJP2sl",
        "subject": {
            "type": "SequenceLocation",
            "id": "ga4gh:SL.674YRNtetAqnFCg9k_8uB1eBzXGKwVPe",
            "sequence_id": "ga4gh:SQ.8_liLu1aycC0tPQPFmUaGXJLDs5SbPZ5",
            "start": {
                "type": "Number",
                "value": 10000
            },
            "end": {
                "type": "IndefiniteRange",
                "comparator": ">=",
                "value": 1223130
            }
        },
        "copy_change": "efo:0030069"
    }
    return CopyNumberChange(**variation)


def test_parsed_copy_number_gain(test_cnv_handler, cn_gain1, cn_gain2):
    """Test that parsed_to_copy_number works for parsed copy number gain queries"""
    # https://www.ncbi.nlm.nih.gov/clinvar/variation/145208/?new_evidence=true
    resp = test_cnv_handler.parsed_to_copy_number(
        143134063, 143284670, VRSTypes.COPY_NUMBER_COUNT, total_copies=3,
        assembly=ClinVarAssembly.GRCH37, chr="chr1",
        start_pos_type=VRSTypes.INDEFINITE_RANGE, end_pos_type=VRSTypes.INDEFINITE_RANGE
    )
    assert resp.copy_number_count.dict() == cn_gain1.dict()
    assert resp.warnings == []

    resp = test_cnv_handler.parsed_to_copy_number(
        143134063, 143284670, VRSTypes.COPY_NUMBER_COUNT, total_copies=3,
        assembly=ClinVarAssembly.HG19, chr="chr1",
        start_pos_type=VRSTypes.INDEFINITE_RANGE, end_pos_type=VRSTypes.INDEFINITE_RANGE
    )
    assert resp.copy_number_count.dict() == cn_gain1.dict()
    assert resp.warnings == []

    resp = test_cnv_handler.parsed_to_copy_number(
        143134063, 143284670, VRSTypes.COPY_NUMBER_COUNT, total_copies=3,
        accession="NC_000001.10", start_pos_type=VRSTypes.INDEFINITE_RANGE,
        end_pos_type=VRSTypes.INDEFINITE_RANGE
    )
    assert resp.copy_number_count.dict() == cn_gain1.dict()
    assert resp.warnings == []

    # https://www.ncbi.nlm.nih.gov/clinvar/variation/146181/?new_evidence=true
    resp = test_cnv_handler.parsed_to_copy_number(
        31738809, 32217725, VRSTypes.COPY_NUMBER_COUNT, total_copies=2,
        assembly=ClinVarAssembly.GRCH38, chr="chr15",
        start_pos_type=VRSTypes.INDEFINITE_RANGE, end_pos_type=VRSTypes.INDEFINITE_RANGE
    )
    assert resp.copy_number_count.dict() == cn_gain2.dict()
    assert resp.warnings == []

    resp = test_cnv_handler.parsed_to_copy_number(
        31738809, 32217725, VRSTypes.COPY_NUMBER_COUNT, total_copies=2,
        assembly=ClinVarAssembly.GRCH38, chr="15",
        start_pos_type=VRSTypes.INDEFINITE_RANGE, end_pos_type=VRSTypes.INDEFINITE_RANGE
    )
    assert resp.copy_number_count.dict() == cn_gain2.dict()
    assert resp.warnings == []

    resp = test_cnv_handler.parsed_to_copy_number(
        31738809, 32217725, VRSTypes.COPY_NUMBER_COUNT, total_copies=2,
        assembly=ClinVarAssembly.HG38, chr="chr15",
        start_pos_type=VRSTypes.INDEFINITE_RANGE, end_pos_type=VRSTypes.INDEFINITE_RANGE
    )
    assert resp.copy_number_count.dict() == cn_gain2.dict()
    assert resp.warnings == []

    resp = test_cnv_handler.parsed_to_copy_number(
        31738809, 32217725, VRSTypes.COPY_NUMBER_COUNT, total_copies=2,
        accession="NC_000015.10", start_pos_type=VRSTypes.INDEFINITE_RANGE,
        end_pos_type=VRSTypes.INDEFINITE_RANGE
    )
    assert resp.copy_number_count.dict() == cn_gain2.dict()
    assert resp.warnings == []


def test_parsed_copy_number_loss(test_cnv_handler, cn_loss1,
                                 cn_loss2):
    """Test that parsed_to_copy_number works for parsed copy number loss queries"""
    # https://www.ncbi.nlm.nih.gov/clinvar/variation/1299222/?new_evidence=true
    resp = test_cnv_handler.parsed_to_copy_number(
        10491132, 10535643, VRSTypes.COPY_NUMBER_COUNT, total_copies=1,
        assembly=ClinVarAssembly.GRCH37, chr="chrX",
        start_pos_type=VRSTypes.INDEFINITE_RANGE, end_pos_type=VRSTypes.INDEFINITE_RANGE
    )
    assert resp.copy_number_count.dict() == cn_loss1.dict()
    assert resp.warnings == []

    resp = test_cnv_handler.parsed_to_copy_number(
        10491132, 10535643, VRSTypes.COPY_NUMBER_COUNT, total_copies=1,
        assembly=ClinVarAssembly.HG19, chr="chrX",
        start_pos_type=VRSTypes.INDEFINITE_RANGE, end_pos_type=VRSTypes.INDEFINITE_RANGE
    )
    assert resp.copy_number_count.dict() == cn_loss1.dict()
    assert resp.warnings == []

    resp = test_cnv_handler.parsed_to_copy_number(
        10491132, 10535643, VRSTypes.COPY_NUMBER_COUNT, total_copies=1,
        assembly=ClinVarAssembly.HG19, chr="X",
        start_pos_type=VRSTypes.INDEFINITE_RANGE, end_pos_type=VRSTypes.INDEFINITE_RANGE
    )
    assert resp.copy_number_count.dict() == cn_loss1.dict()
    assert resp.warnings == []

    resp = test_cnv_handler.parsed_to_copy_number(
        10491132, 10535643, VRSTypes.COPY_NUMBER_COUNT, total_copies=1,
        accession="NC_000023.10", start_pos_type=VRSTypes.INDEFINITE_RANGE,
        end_pos_type=VRSTypes.INDEFINITE_RANGE
    )
    assert resp.copy_number_count.dict() == cn_loss1.dict()
    assert resp.warnings == []

    # https://www.ncbi.nlm.nih.gov/clinvar/variation/148425/?new_evidence=true
    resp = test_cnv_handler.parsed_to_copy_number(
        10001, 1223133, VRSTypes.COPY_NUMBER_COUNT, total_copies=0,
        assembly=ClinVarAssembly.GRCH38, chr="chrY",
        start_pos_type=VRSTypes.INDEFINITE_RANGE, end_pos_type=VRSTypes.INDEFINITE_RANGE
    )
    assert resp.copy_number_count.dict() == cn_loss2.dict()
    assert resp.warnings == []

    resp = test_cnv_handler.parsed_to_copy_number(
        10001, 1223133, VRSTypes.COPY_NUMBER_COUNT, total_copies=0,
        assembly=ClinVarAssembly.HG38, chr="chrY",
        start_pos_type=VRSTypes.INDEFINITE_RANGE, end_pos_type=VRSTypes.INDEFINITE_RANGE
    )
    assert resp.copy_number_count.dict() == cn_loss2.dict()
    assert resp.warnings == []

    resp = test_cnv_handler.parsed_to_copy_number(
        10001, 1223133, VRSTypes.COPY_NUMBER_COUNT, total_copies=0,
        accession="NC_000024.10", start_pos_type=VRSTypes.INDEFINITE_RANGE,
        end_pos_type=VRSTypes.INDEFINITE_RANGE
    )
    assert resp.copy_number_count.dict() == cn_loss2.dict()
    assert resp.warnings == []


def test_parsed_to_cx_var(
    test_cnv_handler, cx_numbers, cx_definite_ranges, cx_indefinite_ranges,
    cx_number_indefinite
):
    """Test that parsed_to_copy_number works for copy number change"""
    # start and end use number
    resp = test_cnv_handler.parsed_to_copy_number(
        10001, 1223133, VRSTypes.COPY_NUMBER_CHANGE,
        copy_change=CopyChange.COMPLETE_GENOMIC_LOSS, assembly=ClinVarAssembly.GRCH38,
        chr="chrY", start_pos_type=VRSTypes.NUMBER, end_pos_type=VRSTypes.NUMBER
    )
    assert resp.copy_number_change.dict() == cx_numbers.dict()
    assert resp.warnings == []

    # start and end use definite ranges
    resp = test_cnv_handler.parsed_to_copy_number(
        10001, 1223130, VRSTypes.COPY_NUMBER_CHANGE,
        copy_change=CopyChange.COMPLETE_GENOMIC_LOSS, assembly=ClinVarAssembly.GRCH38,
        chr="chrY", start_pos_type=VRSTypes.DEFINITE_RANGE,
        end_pos_type=VRSTypes.DEFINITE_RANGE, start1=10006, end1=1223133
    )
    assert resp.copy_number_change.dict() == cx_definite_ranges.dict()
    assert resp.warnings == []

    # start and end use indefinite ranges
    resp = test_cnv_handler.parsed_to_copy_number(
        10001, 1223130, VRSTypes.COPY_NUMBER_CHANGE,
        copy_change=CopyChange.COMPLETE_GENOMIC_LOSS, assembly=ClinVarAssembly.GRCH38,
        chr="chrY", start_pos_type=VRSTypes.INDEFINITE_RANGE,
        end_pos_type=VRSTypes.INDEFINITE_RANGE
    )
    assert resp.copy_number_change.dict() == cx_indefinite_ranges.dict()
    assert resp.warnings == []

    # start uses number and end use indefinite range
    resp = test_cnv_handler.parsed_to_copy_number(
        10001, 1223130, VRSTypes.COPY_NUMBER_CHANGE,
        copy_change=CopyChange.COMPLETE_GENOMIC_LOSS, assembly=ClinVarAssembly.GRCH38,
        chr="chrY", start_pos_type=VRSTypes.NUMBER,
        end_pos_type=VRSTypes.INDEFINITE_RANGE
    )
    assert resp.copy_number_change.dict() == cx_number_indefinite.dict()
    assert resp.warnings == []


def test_invalid(test_cnv_handler):
    """Test invalid copy number queries returns Text variation and warnings"""
    # Invalid Copy Change
    resp = test_cnv_handler.parsed_to_copy_number(
        10491132, 10535643, VRSTypes.COPY_NUMBER_CHANGE, copy_change="efo:1234",
        accession="NC_000001.10"
    )
    assert resp.copy_number_change is None
    assert "copy_change" in resp.warnings[0]
    assert "value is not a valid enumeration member" in resp.warnings[0]

    # NCBI36/hg18 assembly
    expected_w = ["NCBI36 assembly is not currently supported"]
    resp = test_cnv_handler.parsed_to_copy_number(
        2623228, 3150942, VRSTypes.COPY_NUMBER_CHANGE,
        copy_change=CopyChange.GAIN, assembly=ClinVarAssembly.NCBI36, chr="chr1",
        untranslatable_returns_text=True
    )
    assert resp.copy_number_change.type == "Text"
    assert resp.warnings == expected_w

    resp = test_cnv_handler.parsed_to_copy_number(
        2623228, 3150942, VRSTypes.COPY_NUMBER_CHANGE, copy_change=CopyChange.GAIN,
        assembly=ClinVarAssembly.HG18, chr="chr1", untranslatable_returns_text=True
    )
    assert resp.copy_number_change.type == "Text"
    assert resp.warnings == expected_w

    # Must give both assembly + chr or accession
    expected_w = ["Must provide either `accession` or both `assembly` and `chr`."]
    resp = test_cnv_handler.parsed_to_copy_number(
        31738809, 32217725, VRSTypes.COPY_NUMBER_CHANGE, copy_change=CopyChange.GAIN,
        assembly="hg38", untranslatable_returns_text=True
    )
    assert resp.copy_number_change.type == "Text"
    assert resp.warnings == expected_w

    resp = test_cnv_handler.parsed_to_copy_number(
        31738809, 32217725, VRSTypes.COPY_NUMBER_CHANGE, copy_change=CopyChange.GAIN,
        chr="chr15", untranslatable_returns_text=True
    )
    assert resp.copy_number_change.type == "Text"
    assert resp.warnings == expected_w

    resp = test_cnv_handler.parsed_to_copy_number(
        31738809, 32217725, VRSTypes.COPY_NUMBER_CHANGE, copy_change=CopyChange.GAIN,
        untranslatable_returns_text=True
    )
    assert resp.copy_number_change.type == "Text"
    assert resp.warnings == expected_w

    # invalid chr
    resp = test_cnv_handler.parsed_to_copy_number(
        10001, 1223133, VRSTypes.COPY_NUMBER_CHANGE, copy_change=CopyChange.GAIN,
        assembly=ClinVarAssembly.GRCH38, chr="z"
    )
    assert resp.copy_number_change is None
    assert resp.warnings == \
        ["SeqRepo unable to get translated identifiers for GRCh38:z"]

    # invalid assembly
    resp = test_cnv_handler.parsed_to_copy_number(
        10001, 1223133, VRSTypes.COPY_NUMBER_CHANGE, copy_change=CopyChange.GAIN,
        assembly="GRCh99"
    )
    assert resp.copy_number_change is None
    assert resp.warnings

    # invalid accession
    resp = test_cnv_handler.parsed_to_copy_number(
        10491132, 10535643, VRSTypes.COPY_NUMBER_CHANGE, copy_change=CopyChange.GAIN,
        accession="NC_00002310"
    )
    assert resp.copy_number_change is None
    assert resp.warnings == \
        ["SeqRepo unable to get translated identifiers for NC_00002310"]

    # Invalid position
    resp = test_cnv_handler.parsed_to_copy_number(
        31738809, 2302991250, VRSTypes.COPY_NUMBER_CHANGE, copy_change=CopyChange.GAIN,
        accession="NC_000015.10"
    )
    assert resp.copy_number_change is None
    assert resp.warnings == ["SeqRepo ValueError: Position out of range (2302991249)"]

    # start must be less than end
    resp = test_cnv_handler.parsed_to_copy_number(
        10001, 1223130, VRSTypes.COPY_NUMBER_CHANGE,
        copy_change=CopyChange.COMPLETE_GENOMIC_LOSS, assembly=ClinVarAssembly.GRCH38,
        chr="chrY", start_pos_type=VRSTypes.DEFINITE_RANGE,
        end_pos_type=VRSTypes.DEFINITE_RANGE, start1=1223132, end1=1223133
    )
    assert resp.copy_number_change is None
    assert "end positions must be greater than start" in resp.warnings[0]
