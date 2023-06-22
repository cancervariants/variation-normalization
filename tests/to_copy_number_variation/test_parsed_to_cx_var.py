"""Test that parsed_to_cx_var works correctly"""
import pytest
from ga4gh.vrsatile.pydantic.vrs_models import CopyNumberChange, CopyChange, VRSTypes

from variation.schemas.service_schema import ClinVarAssembly


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


def test_parsed_to_cx_var(
    test_cnv_handler, cx_numbers, cx_definite_ranges, cx_indefinite_ranges,
    cx_number_indefinite
):
    """Test that parsed_to_cx_var works"""
    # start and end use number
    resp = test_cnv_handler.parsed_to_cx_var(
        10001, 1223133, CopyChange.COMPLETE_GENOMIC_LOSS,
        assembly=ClinVarAssembly.GRCH38, chr="chrY", start_pos_type=VRSTypes.NUMBER,
        end_pos_type=VRSTypes.NUMBER
    )
    assert resp.copy_number_change.dict() == cx_numbers.dict()
    assert resp.warnings == []

    # start and end use definite ranges
    resp = test_cnv_handler.parsed_to_cx_var(
        10001, 1223130, CopyChange.COMPLETE_GENOMIC_LOSS,
        assembly=ClinVarAssembly.GRCH38, chr="chrY",
        start_pos_type=VRSTypes.DEFINITE_RANGE, end_pos_type=VRSTypes.DEFINITE_RANGE,
        start1=10006, end1=1223133
    )
    assert resp.copy_number_change.dict() == cx_definite_ranges.dict()
    assert resp.warnings == []

    # start and end use indefinite ranges
    resp = test_cnv_handler.parsed_to_cx_var(
        10001, 1223130, CopyChange.COMPLETE_GENOMIC_LOSS,
        assembly=ClinVarAssembly.GRCH38, chr="chrY",
        start_pos_type=VRSTypes.INDEFINITE_RANGE, end_pos_type=VRSTypes.INDEFINITE_RANGE
    )
    assert resp.copy_number_change.dict() == cx_indefinite_ranges.dict()
    assert resp.warnings == []

    # start uses number and end use indefinite range
    resp = test_cnv_handler.parsed_to_cx_var(
        10001, 1223130, CopyChange.COMPLETE_GENOMIC_LOSS,
        assembly=ClinVarAssembly.GRCH38, chr="chrY",
        start_pos_type=VRSTypes.NUMBER, end_pos_type=VRSTypes.INDEFINITE_RANGE
    )
    assert resp.copy_number_change.dict() == cx_number_indefinite.dict()
    assert resp.warnings == []


def test_invalid(test_cnv_handler):
    """Test invalid queries returns Text variation and warnings"""
    # Invalid Copy Change
    resp = test_cnv_handler.parsed_to_cx_var(
        10491132, 10535643, "efo:1234", accession="NC_000001.10"
    )
    assert resp.copy_number_change is None
    assert "copy_change" in resp.warnings[0]
    assert "value is not a valid enumeration member" in resp.warnings[0]

    # NCBI36/hg18 assembly
    expected_w = ["NCBI36 assembly is not currently supported"]
    resp = test_cnv_handler.parsed_to_cx_var(
        2623228, 3150942, CopyChange.GAIN, assembly=ClinVarAssembly.NCBI36, chr="chr1",
        untranslatable_returns_text=True)
    assert resp.copy_number_change.type == "Text"
    assert resp.warnings == expected_w

    resp = test_cnv_handler.parsed_to_cx_var(
        2623228, 3150942, CopyChange.GAIN, assembly=ClinVarAssembly.HG18, chr="chr1",
        untranslatable_returns_text=True)
    assert resp.copy_number_change.type == "Text"
    assert resp.warnings == expected_w

    # Must give both assembly + chr or accession
    expected_w = ["Must provide either `accession` or both `assembly` and `chr`."]
    resp = test_cnv_handler.parsed_to_cx_var(
        31738809, 32217725, CopyChange.GAIN, assembly="hg38",
        untranslatable_returns_text=True
    )
    assert resp.copy_number_change.type == "Text"
    assert resp.warnings == expected_w

    resp = test_cnv_handler.parsed_to_cx_var(
        31738809, 32217725, CopyChange.GAIN, chr="chr15",
        untranslatable_returns_text=True
    )
    assert resp.copy_number_change.type == "Text"
    assert resp.warnings == expected_w

    resp = test_cnv_handler.parsed_to_cx_var(
        31738809, 32217725, CopyChange.GAIN, untranslatable_returns_text=True)
    assert resp.copy_number_change.type == "Text"
    assert resp.warnings == expected_w

    # invalid chr
    resp = test_cnv_handler.parsed_to_cx_var(
        10001, 1223133, CopyChange.GAIN, assembly=ClinVarAssembly.GRCH38, chr="z")
    assert resp.copy_number_change is None
    assert resp.warnings == \
        ["SeqRepo unable to get translated identifiers for GRCh38:z"]

    # invalid assembly
    resp = test_cnv_handler.parsed_to_cx_var(
        10001, 1223133, CopyChange.GAIN, assembly="GRCh99")
    assert resp.copy_number_change is None
    assert resp.warnings

    # invalid accession
    resp = test_cnv_handler.parsed_to_cx_var(
        10491132, 10535643, CopyChange.GAIN, accession="NC_00002310")
    assert resp.copy_number_change is None
    assert resp.warnings == \
        ["SeqRepo unable to get translated identifiers for NC_00002310"]

    # Invalid position
    resp = test_cnv_handler.parsed_to_cx_var(
        31738809, 2302991250, CopyChange.GAIN, accession="NC_000015.10")
    assert resp.copy_number_change is None
    assert resp.warnings == ["SeqRepo ValueError: Position out of range (2302991249)"]

    # start must be less than end
    resp = test_cnv_handler.parsed_to_cx_var(
        10001, 1223130, CopyChange.COMPLETE_GENOMIC_LOSS,
        assembly=ClinVarAssembly.GRCH38, chr="chrY",
        start_pos_type=VRSTypes.DEFINITE_RANGE, end_pos_type=VRSTypes.DEFINITE_RANGE,
        start1=1223132, end1=1223133
    )
    assert resp.copy_number_change is None
    assert "end positions must be greater than start" in resp.warnings[0]
