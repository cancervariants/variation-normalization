"""Test that parsed_to_copy_number works correctly"""
from copy import deepcopy

import pytest
from ga4gh.vrs import models
from ga4gh.vrsatile.pydantic.vrs_models import (
    CopyNumberCount, CopyNumberChange, CopyChange, VRSTypes, Comparator
)
from pydantic.error_wrappers import ValidationError

from variation.schemas.copy_number_schema import (
    ClinVarAssembly, ParsedToCxVarQuery, ParsedToCnVarQuery
)
from variation.to_copy_number_variation import ToCopyNumberException


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
def cn_gain2_37():
    """Create test fixture for clinvar copy number gain on GRCh37 assembly.
    https://www.ncbi.nlm.nih.gov/clinvar/variation/146181/?new_evidence=true
    """
    variation = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.J2BzX2gg87q75lkCuFV7HylE_uG_cQJN",
        "subject": {
            "type": "SequenceLocation",
            "id": "ga4gh:SL.PZFhnd9S3DxnhQvMUuWd1b3kFqQe937p",
            "sequence_id": "ga4gh:SQ.zIMZb3Ft7RdWa5XYq0PxIlezLY2ccCgt",
            "start": {
                "type": "IndefiniteRange",
                "value": 32031011,
                "comparator": "<="
            },
            "end": {
                "type": "IndefiniteRange",
                "value": 32509926,
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
def cn_definite_number():
    """Create test fixture for copy number count using definite range for start and
    number for end
    """
    variation = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.uzKJQTAw2pYuAGaeHclKATjgPJQHjwmk",
        "subject": {
            "type": "SequenceLocation",
            "id": "ga4gh:SL.wfI7pH2g6czJN8dWwgF1gH9BD4pnJCFG",
            "sequence_id": "ga4gh:SQ.S_KjnFVz-FE7M0W6yoaUDgYxLPc1jyWU",
            "start": {
                "type": "DefiniteRange",
                "min": 143134062,
                "max": 143134064
            },
            "end": {
                "type": "Number",
                "value": 143284670
            }
        },
        "copies": {"type": "Number", "value": 3}
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


def test_get_parsed_ac(test_cnv_handler):
    """Test that _get_parsed_ac works correctly"""
    for assembly in [ClinVarAssembly.GRCH37, ClinVarAssembly.HG19]:
        resp = test_cnv_handler._get_parsed_ac(assembly, "chr7", use_grch38=False)
        assert resp.lifted_over is False
        assert resp.accession == "ga4gh:SQ.IW78mgV5Cqf6M24hy52hPjyyo5tCCd86"

        resp = test_cnv_handler._get_parsed_ac(assembly, "chr7", use_grch38=True)
        assert resp.lifted_over is True
        assert resp.accession == "ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul"

    for assembly in [ClinVarAssembly.GRCH38, ClinVarAssembly.HG38]:
        resp = test_cnv_handler._get_parsed_ac(assembly, "chr7", use_grch38=False)
        assert resp.lifted_over is False
        assert resp.accession == "ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul"

        resp = test_cnv_handler._get_parsed_ac(assembly, "chr7", use_grch38=True)
        assert resp.lifted_over is False
        assert resp.accession == "ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul"

    with pytest.raises(ToCopyNumberException) as e:
        test_cnv_handler._get_parsed_ac(
            ClinVarAssembly.NCBI36, "chr7", use_grch38=False
        )
    assert str(e.value) == "NCBI36 assembly is not currently supported"

    with pytest.raises(ToCopyNumberException) as e:
        test_cnv_handler._get_parsed_ac(
            ClinVarAssembly.HG18, "chr7", use_grch38=False
        )
    assert str(e.value) == "hg18 assembly is not currently supported"


def test_get_parsed_ac_chr(test_cnv_handler):
    """Test that _get_parsed_ac_chr works correctly"""
    resp = test_cnv_handler._get_parsed_ac_chr("NC_000007.13", False)
    assert resp.accession == "ga4gh:SQ.IW78mgV5Cqf6M24hy52hPjyyo5tCCd86"
    assert resp.chromosome == "chr7"
    assert resp.lifted_over is False

    resp = test_cnv_handler._get_parsed_ac_chr("NC_000007.13", True)
    assert resp.accession == "ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul"
    assert resp.chromosome == "chr7"
    assert resp.lifted_over is True

    for do_liftover in [True, False]:
        resp = test_cnv_handler._get_parsed_ac_chr("NC_000007.14", do_liftover)
        assert resp.accession == "ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul"
        assert resp.chromosome == "chr7"
        assert resp.lifted_over is False

    # if genomic ac not provided
    with pytest.raises(ToCopyNumberException) as e:
        test_cnv_handler._get_parsed_ac_chr("NP_000542.1", False)
    assert str(e.value) == "Not a supported genomic accession: NP_000542.1"

    # invalid accession
    with pytest.raises(ToCopyNumberException) as e:
        test_cnv_handler._get_parsed_ac_chr("NC_00000713", False)
    assert str(e.value) == "SeqRepo unable to get translated identifiers for NC_00000713"  # noqa: E501


def test_validate_pos(test_cnv_handler):
    """Test that _validate_ac_pos works correctly"""
    resp = test_cnv_handler._validate_ac_pos("NC_000007.14", 140753336)
    assert resp is None

    # invalid accession
    with pytest.raises(ToCopyNumberException) as e:
        test_cnv_handler._validate_ac_pos("NC_00000714", 140753336)
    assert str(e.value) == "Accession not found in SeqRepo: NC_00000714"

    # invalid pos
    with pytest.raises(ToCopyNumberException) as e:
        test_cnv_handler._validate_ac_pos("NC_000007.14", 159345975)
    assert str(e.value) == "Position (159345975) is not valid on NC_000007.14"

    # invalid pos
    with pytest.raises(ToCopyNumberException) as e:
        test_cnv_handler._validate_ac_pos("NC_000007.14", 99999999999)
    assert str(e.value) == "SeqRepo ValueError: Position out of range (99999999998)"


def test_get_vrs_loc_start_or_end(test_cnv_handler):
    """Test that _get_vrs_loc_start_or_end works correctly"""
    ac = "NC_000007.14"
    pos0 = 140753336
    pos1 = 140753350

    # Number start
    resp = test_cnv_handler._get_vrs_loc_start_or_end(
        ac, pos0, VRSTypes.NUMBER, is_start=True
    )
    assert resp == models.Number(value=140753335, type="Number")

    # Number end
    resp = test_cnv_handler._get_vrs_loc_start_or_end(
        ac, pos0, VRSTypes.NUMBER, is_start=False
    )
    assert resp == models.Number(value=140753336, type="Number")

    # Definite Range start
    resp = test_cnv_handler._get_vrs_loc_start_or_end(
        ac, pos0, VRSTypes.DEFINITE_RANGE, is_start=True, pos1=pos1
    )
    assert resp == models.DefiniteRange(
        min=140753335, max=140753349, type="DefiniteRange"
    )

    # Definite Range end
    resp = test_cnv_handler._get_vrs_loc_start_or_end(
        ac, pos0, VRSTypes.DEFINITE_RANGE, is_start=False, pos1=pos1
    )
    assert resp == models.DefiniteRange(
        min=140753337, max=140753351, type="DefiniteRange"
    )

    # Indefinite Range start
    resp = test_cnv_handler._get_vrs_loc_start_or_end(
        ac, pos0, VRSTypes.INDEFINITE_RANGE, is_start=True,
        comparator=Comparator.LT_OR_EQUAL
    )
    assert resp == models.IndefiniteRange(
        comparator="<=", value=140753335, type="IndefiniteRange"
    )

    # Indefinite Range end
    resp = test_cnv_handler._get_vrs_loc_start_or_end(
        ac, pos0, VRSTypes.INDEFINITE_RANGE, is_start=False,
        comparator=Comparator.GT_OR_EQUAL
    )
    assert resp == models.IndefiniteRange(
        comparator=">=", value=140753336, type="IndefiniteRange"
    )


def test_liftover_pos(test_cnv_handler):
    """Test that _liftover_pos works correctly"""
    resp = test_cnv_handler._liftover_pos("chr7", 140453136, 140453137, None, None)
    assert resp == {
        "start0": 140753336,
        "end0": 140753337,
        "start1": None,
        "end1": None
    }

    resp = test_cnv_handler._liftover_pos(
        "chr7", 140453136, 140453137, 140453138, 140453139
    )
    assert resp == {
        "start0": 140753336,
        "end0": 140753337,
        "start1": 140753338,
        "end1": 140753339
    }

    # invalid pos
    with pytest.raises(ToCopyNumberException) as e:
        test_cnv_handler._liftover_pos(
            "chr7", 159345975, 159345976, None, None
        )
    assert str(e.value) == "Unable to liftover: chr7 with pos 159345975"


def test_parsed_copy_number_gain(test_cnv_handler, cn_gain1, cn_gain2, cn_gain2_37):
    """Test that parsed_to_copy_number works for parsed copy number gain queries"""
    # https://www.ncbi.nlm.nih.gov/clinvar/variation/145208/?new_evidence=true
    rb = ParsedToCnVarQuery(
        start0=143134063, end0=143284670, copies0=3,
        assembly=ClinVarAssembly.GRCH37, chromosome="chr1",
        start_pos_type=VRSTypes.INDEFINITE_RANGE,
        start_pos_comparator=Comparator.LT_OR_EQUAL,
        end_pos_type=VRSTypes.INDEFINITE_RANGE,
        end_pos_comparator=Comparator.GT_OR_EQUAL
    )
    resp = test_cnv_handler.parsed_to_copy_number(rb)
    assert resp.copy_number_count.dict() == cn_gain1.dict()
    assert resp.warnings == []

    rb = ParsedToCnVarQuery(
        start0=143134063, end0=143284670, copies0=3, assembly=ClinVarAssembly.HG19,
        chromosome="chr1", start_pos_type=VRSTypes.INDEFINITE_RANGE,
        start_pos_comparator=Comparator.LT_OR_EQUAL,
        end_pos_type=VRSTypes.INDEFINITE_RANGE,
        end_pos_comparator=Comparator.GT_OR_EQUAL
    )
    resp = test_cnv_handler.parsed_to_copy_number(rb)
    assert resp.copy_number_count.dict() == cn_gain1.dict()
    assert resp.warnings == []

    rb = ParsedToCnVarQuery(
        start0=143134063, end0=143284670, copies0=3, accession="NC_000001.10",
        start_pos_type=VRSTypes.INDEFINITE_RANGE,
        start_pos_comparator=Comparator.LT_OR_EQUAL,
        end_pos_type=VRSTypes.INDEFINITE_RANGE,
        end_pos_comparator=Comparator.GT_OR_EQUAL
    )
    resp = test_cnv_handler.parsed_to_copy_number(rb)
    assert resp.copy_number_count.dict() == cn_gain1.dict()
    assert resp.warnings == []

    # https://www.ncbi.nlm.nih.gov/clinvar/variation/146181/?new_evidence=true
    # 38
    rb = ParsedToCnVarQuery(
        start0=31738809, end0=32217725, copies0=2, assembly=ClinVarAssembly.GRCH38,
        chromosome="chr15", start_pos_type=VRSTypes.INDEFINITE_RANGE,
        start_pos_comparator=Comparator.LT_OR_EQUAL,
        end_pos_type=VRSTypes.INDEFINITE_RANGE,
        end_pos_comparator=Comparator.GT_OR_EQUAL
    )
    resp = test_cnv_handler.parsed_to_copy_number(rb)
    assert resp.copy_number_count.dict() == cn_gain2.dict()
    assert resp.warnings == []

    # 38 with liftover (shouldnt do anything)
    rb = ParsedToCnVarQuery(
        start0=31738809, end0=32217725, copies0=2, assembly=ClinVarAssembly.GRCH38,
        chromosome="chr15", do_liftover=True, start_pos_type=VRSTypes.INDEFINITE_RANGE,
        start_pos_comparator=Comparator.LT_OR_EQUAL,
        end_pos_type=VRSTypes.INDEFINITE_RANGE,
        end_pos_comparator=Comparator.GT_OR_EQUAL
    )
    resp = test_cnv_handler.parsed_to_copy_number(rb)
    assert resp.copy_number_count.dict() == cn_gain2.dict()
    assert resp.warnings == []

    # 38 with liftover (shouldnt do anything)
    rb = ParsedToCnVarQuery(
        start0=31738809, end0=32217725, copies0=2, assembly=ClinVarAssembly.HG38,
        chromosome="chr15", do_liftover=True, start_pos_type=VRSTypes.INDEFINITE_RANGE,
        start_pos_comparator=Comparator.LT_OR_EQUAL,
        end_pos_type=VRSTypes.INDEFINITE_RANGE,
        end_pos_comparator=Comparator.GT_OR_EQUAL
    )
    resp = test_cnv_handler.parsed_to_copy_number(rb)
    assert resp.copy_number_count.dict() == cn_gain2.dict()
    assert resp.warnings == []

    # 38
    rb = ParsedToCnVarQuery(
        start0=31738809, end0=32217725, copies0=2, assembly=ClinVarAssembly.HG38,
        chromosome="chr15", start_pos_type=VRSTypes.INDEFINITE_RANGE,
        start_pos_comparator=Comparator.LT_OR_EQUAL,
        end_pos_type=VRSTypes.INDEFINITE_RANGE,
        end_pos_comparator=Comparator.GT_OR_EQUAL
    )
    resp = test_cnv_handler.parsed_to_copy_number(rb)
    assert resp.copy_number_count.dict() == cn_gain2.dict()
    assert resp.warnings == []

    # 38 accession
    rb = ParsedToCnVarQuery(
        start0=31738809, end0=32217725, copies0=2, accession="NC_000015.10",
        start_pos_type=VRSTypes.INDEFINITE_RANGE,
        start_pos_comparator=Comparator.LT_OR_EQUAL,
        end_pos_type=VRSTypes.INDEFINITE_RANGE,
        end_pos_comparator=Comparator.GT_OR_EQUAL
    )
    resp = test_cnv_handler.parsed_to_copy_number(rb)
    assert resp.copy_number_count.dict() == cn_gain2.dict()
    assert resp.warnings == []

    # 38 accession with liftover (shouldnt do anything)
    rb = ParsedToCnVarQuery(
        start0=31738809, end0=32217725, copies0=2, accession="NC_000015.10",
        start_pos_type=VRSTypes.INDEFINITE_RANGE,
        start_pos_comparator=Comparator.LT_OR_EQUAL,
        end_pos_type=VRSTypes.INDEFINITE_RANGE,
        end_pos_comparator=Comparator.GT_OR_EQUAL
    )
    resp = test_cnv_handler.parsed_to_copy_number(rb)
    assert resp.copy_number_count.dict() == cn_gain2.dict()
    assert resp.warnings == []

    # 37 with liftover
    rb = ParsedToCnVarQuery(
        start0=32031012, end0=32509926, copies0=2, accession="NC_000015.9",
        do_liftover=True, start_pos_type=VRSTypes.INDEFINITE_RANGE,
        start_pos_comparator=Comparator.LT_OR_EQUAL,
        end_pos_type=VRSTypes.INDEFINITE_RANGE,
        end_pos_comparator=Comparator.GT_OR_EQUAL
    )
    resp = test_cnv_handler.parsed_to_copy_number(rb)
    assert resp.copy_number_count.dict() == cn_gain2.dict()
    assert resp.warnings == []

    # 37 chr+accession with liftover
    rb = ParsedToCnVarQuery(
        start0=32031012, end0=32509926, copies0=2, chromosome="chr15",
        assembly=ClinVarAssembly.GRCH37, do_liftover=True,
        start_pos_type=VRSTypes.INDEFINITE_RANGE,
        start_pos_comparator=Comparator.LT_OR_EQUAL,
        end_pos_type=VRSTypes.INDEFINITE_RANGE,
        end_pos_comparator=Comparator.GT_OR_EQUAL
    )
    resp = test_cnv_handler.parsed_to_copy_number(rb)
    assert resp.copy_number_count.dict() == cn_gain2.dict()
    assert resp.warnings == []

    # 37 with no liftover
    rb = ParsedToCnVarQuery(
        start0=32031012, end0=32509926, copies0=2, accession="NC_000015.9",
        start_pos_type=VRSTypes.INDEFINITE_RANGE,
        start_pos_comparator=Comparator.LT_OR_EQUAL,
        end_pos_type=VRSTypes.INDEFINITE_RANGE,
        end_pos_comparator=Comparator.GT_OR_EQUAL
    )
    resp = test_cnv_handler.parsed_to_copy_number(rb)
    assert resp.copy_number_count.dict() == cn_gain2_37.dict()
    assert resp.warnings == []

    # 37 chr+accession with no liftover
    rb = ParsedToCnVarQuery(
        start0=32031012, end0=32509926, copies0=2, chromosome="chr15",
        assembly=ClinVarAssembly.GRCH37, start_pos_type=VRSTypes.INDEFINITE_RANGE,
        start_pos_comparator=Comparator.LT_OR_EQUAL,
        end_pos_type=VRSTypes.INDEFINITE_RANGE,
        end_pos_comparator=Comparator.GT_OR_EQUAL
    )
    resp = test_cnv_handler.parsed_to_copy_number(rb)
    assert resp.copy_number_count.dict() == cn_gain2_37.dict()
    assert resp.warnings == []


def test_parsed_copy_number_loss(test_cnv_handler, cn_loss1,
                                 cn_loss2):
    """Test that parsed_to_copy_number works for parsed copy number loss queries"""
    # https://www.ncbi.nlm.nih.gov/clinvar/variation/1299222/?new_evidence=true
    rb = ParsedToCnVarQuery(
        start0=10491132, end0=10535643, copies0=1, assembly=ClinVarAssembly.GRCH37,
        chromosome="chrX", start_pos_type=VRSTypes.INDEFINITE_RANGE,
        start_pos_comparator=Comparator.LT_OR_EQUAL,
        end_pos_type=VRSTypes.INDEFINITE_RANGE,
        end_pos_comparator=Comparator.GT_OR_EQUAL
    )
    resp = test_cnv_handler.parsed_to_copy_number(rb)
    assert resp.copy_number_count.dict() == cn_loss1.dict()
    assert resp.warnings == []

    rb = ParsedToCnVarQuery(
        start0=10491132, end0=10535643, copies0=1, assembly=ClinVarAssembly.HG19,
        chromosome="chrX", start_pos_type=VRSTypes.INDEFINITE_RANGE,
        start_pos_comparator=Comparator.LT_OR_EQUAL,
        end_pos_type=VRSTypes.INDEFINITE_RANGE,
        end_pos_comparator=Comparator.GT_OR_EQUAL
    )
    resp = test_cnv_handler.parsed_to_copy_number(rb)
    assert resp.copy_number_count.dict() == cn_loss1.dict()
    assert resp.warnings == []

    rb = ParsedToCnVarQuery(
        start0=10491132, end0=10535643, copies0=1, accession="NC_000023.10",
        start_pos_type=VRSTypes.INDEFINITE_RANGE,
        start_pos_comparator=Comparator.LT_OR_EQUAL,
        end_pos_type=VRSTypes.INDEFINITE_RANGE,
        end_pos_comparator=Comparator.GT_OR_EQUAL
    )
    resp = test_cnv_handler.parsed_to_copy_number(rb)
    assert resp.copy_number_count.dict() == cn_loss1.dict()
    assert resp.warnings == []

    # https://www.ncbi.nlm.nih.gov/clinvar/variation/148425/?new_evidence=true
    rb = ParsedToCnVarQuery(
        start0=10001, end0=1223133, copies0=0, assembly=ClinVarAssembly.GRCH38,
        chromosome="chrY", start_pos_type=VRSTypes.INDEFINITE_RANGE,
        start_pos_comparator=Comparator.LT_OR_EQUAL,
        end_pos_type=VRSTypes.INDEFINITE_RANGE,
        end_pos_comparator=Comparator.GT_OR_EQUAL
    )
    resp = test_cnv_handler.parsed_to_copy_number(rb)
    assert resp.copy_number_count.dict() == cn_loss2.dict()
    assert resp.warnings == []

    rb = ParsedToCnVarQuery(
        start0=10001, end0=1223133, copies0=0, assembly=ClinVarAssembly.HG38,
        chromosome="chrY", start_pos_type=VRSTypes.INDEFINITE_RANGE,
        start_pos_comparator=Comparator.LT_OR_EQUAL,
        end_pos_type=VRSTypes.INDEFINITE_RANGE,
        end_pos_comparator=Comparator.GT_OR_EQUAL
    )
    resp = test_cnv_handler.parsed_to_copy_number(rb)
    assert resp.copy_number_count.dict() == cn_loss2.dict()
    assert resp.warnings == []

    rb = ParsedToCnVarQuery(
        start0=10001, end0=1223133, copies0=0, accession="NC_000024.10",
        start_pos_type=VRSTypes.INDEFINITE_RANGE,
        start_pos_comparator=Comparator.LT_OR_EQUAL,
        end_pos_type=VRSTypes.INDEFINITE_RANGE,
        end_pos_comparator=Comparator.GT_OR_EQUAL
    )
    resp = test_cnv_handler.parsed_to_copy_number(rb)
    assert resp.copy_number_count.dict() == cn_loss2.dict()
    assert resp.warnings == []


def test_to_parsed_cn_var(
    test_cnv_handler, cn_definite_number
):
    """Test that parsed_to_copy_number works correctly for copy number count"""
    # start uses definite and end uses number
    rb = ParsedToCnVarQuery(
        start0=143134063, end0=143284670, copies0=3,
        assembly=ClinVarAssembly.GRCH37, chromosome="chr1",
        start_pos_type=VRSTypes.DEFINITE_RANGE, start1=143134065
    )
    resp = test_cnv_handler.parsed_to_copy_number(rb)
    assert resp.copy_number_count.dict() == cn_definite_number.dict()
    assert resp.warnings == []

    # copies is definite range
    rb = ParsedToCnVarQuery(
        start0=143134063, end0=143284670, copies0=3, copies1=5,
        copies_type=VRSTypes.DEFINITE_RANGE,
        assembly=ClinVarAssembly.GRCH37, chromosome="chr1",
        start_pos_type=VRSTypes.DEFINITE_RANGE, start1=143134065
    )
    resp = test_cnv_handler.parsed_to_copy_number(rb)
    expected = deepcopy(cn_definite_number.dict())
    expected["copies"] = {"type": "DefiniteRange", "min": 3, "max": 5}
    expected["id"] = "ga4gh:CN.abS9SjAfusSTVgdX3kIUOGV8q8bmr5uO"
    assert resp.copy_number_count.dict() == expected
    assert resp.warnings == []

    # copies is indefinite range <=
    rb = ParsedToCnVarQuery(
        start0=143134063, end0=143284670, copies0=3,
        copies_comparator=Comparator.LT_OR_EQUAL,
        copies_type=VRSTypes.INDEFINITE_RANGE,
        assembly=ClinVarAssembly.GRCH37, chromosome="chr1",
        start_pos_type=VRSTypes.DEFINITE_RANGE, start1=143134065
    )
    resp = test_cnv_handler.parsed_to_copy_number(rb)
    expected = deepcopy(cn_definite_number.dict())
    expected["copies"] = {"type": "IndefiniteRange", "comparator": "<=", "value": 3}
    expected["id"] = "ga4gh:CN.2aMndMYx2uVCN-XcG-uLimeXkRRHGDCs"
    assert resp.copy_number_count.dict() == expected
    assert resp.warnings == []

    # copies is indefinite range >=
    rb = ParsedToCnVarQuery(
        start0=143134063, end0=143284670, copies0=3,
        copies_comparator=Comparator.GT_OR_EQUAL,
        copies_type=VRSTypes.INDEFINITE_RANGE,
        assembly=ClinVarAssembly.GRCH37, chromosome="chr1",
        start_pos_type=VRSTypes.DEFINITE_RANGE, start1=143134065
    )
    resp = test_cnv_handler.parsed_to_copy_number(rb)
    expected = deepcopy(cn_definite_number.dict())
    expected["copies"] = {"type": "IndefiniteRange", "comparator": ">=", "value": 3}
    expected["id"] = "ga4gh:CN.UDIZLBzw7TiaM7Q7Ani4BPoWyVbstlKu"
    assert resp.copy_number_count.dict() == expected
    assert resp.warnings == []

    # start_pos and end_pos indefinite range
    rb = ParsedToCnVarQuery(
        start0=143134063, end0=143284670, copies0=3,
        assembly=ClinVarAssembly.GRCH37, chromosome="chr1",
        start_pos_type=VRSTypes.INDEFINITE_RANGE,
        end_pos_type=VRSTypes.INDEFINITE_RANGE,
        start_pos_comparator=Comparator.GT_OR_EQUAL,
        end_pos_comparator=Comparator.LT_OR_EQUAL
    )
    resp = test_cnv_handler.parsed_to_copy_number(rb)
    assert resp.copy_number_count.dict() == {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.wHCfVMPVV1AZI6-eURZgAKbNr2srugdL",
        "subject": {
            "type": "SequenceLocation",
            "id": "ga4gh:SL.xHA8jY9-PauRA-3kfu3_RkP-SW4poqiI",
            "sequence_id": "ga4gh:SQ.S_KjnFVz-FE7M0W6yoaUDgYxLPc1jyWU",
            "start": {
                "type": "IndefiniteRange",
                "value": 143134062,
                "comparator": ">="
            },
            "end": {
                "type": "IndefiniteRange",
                "value": 143284670,
                "comparator": "<="
            }
        },
        "copies": {"type": "Number", "value": 3}
    }


def test_parsed_to_cx_var(
    test_cnv_handler, cx_numbers, cx_definite_ranges, cx_indefinite_ranges,
    cx_number_indefinite
):
    """Test that parsed_to_copy_number works for copy number change"""
    # start and end use number
    rb = ParsedToCxVarQuery(
        start0=10001, end0=1223133,
        copy_change=CopyChange.COMPLETE_GENOMIC_LOSS, assembly=ClinVarAssembly.GRCH38,
        chromosome="chrY", start_pos_type=VRSTypes.NUMBER, end_pos_type=VRSTypes.NUMBER
    )
    resp = test_cnv_handler.parsed_to_copy_number(rb)
    assert resp.copy_number_change.dict() == cx_numbers.dict()
    assert resp.warnings == []

    # start and end use definite ranges
    rb = ParsedToCxVarQuery(
        start0=10001, end0=1223130,
        copy_change=CopyChange.COMPLETE_GENOMIC_LOSS, assembly=ClinVarAssembly.GRCH38,
        chromosome="chrY", start_pos_type=VRSTypes.DEFINITE_RANGE,
        end_pos_type=VRSTypes.DEFINITE_RANGE, start1=10006, end1=1223133
    )
    resp = test_cnv_handler.parsed_to_copy_number(rb)
    assert resp.copy_number_change.dict() == cx_definite_ranges.dict()
    assert resp.warnings == []

    # start and end use indefinite ranges
    rb = ParsedToCxVarQuery(
        start0=10001, end0=1223130,
        copy_change=CopyChange.COMPLETE_GENOMIC_LOSS, assembly=ClinVarAssembly.GRCH38,
        chromosome="chrY", start_pos_type=VRSTypes.INDEFINITE_RANGE,
        start_pos_comparator=Comparator.LT_OR_EQUAL,
        end_pos_type=VRSTypes.INDEFINITE_RANGE,
        end_pos_comparator=Comparator.GT_OR_EQUAL
    )
    resp = test_cnv_handler.parsed_to_copy_number(rb)
    assert resp.copy_number_change.dict() == cx_indefinite_ranges.dict()
    assert resp.warnings == []

    # start uses number and end use indefinite range
    rb = ParsedToCxVarQuery(
        start0=10001, end0=1223130,
        copy_change=CopyChange.COMPLETE_GENOMIC_LOSS, assembly=ClinVarAssembly.GRCH38,
        chromosome="chrY", start_pos_type=VRSTypes.NUMBER,
        end_pos_type=VRSTypes.INDEFINITE_RANGE,
        end_pos_comparator=Comparator.GT_OR_EQUAL
    )
    resp = test_cnv_handler.parsed_to_copy_number(rb)
    assert resp.copy_number_change.dict() == cx_number_indefinite.dict()
    assert resp.warnings == []


def test_invalid(test_cnv_handler):
    """Test invalid copy number queries returns Text variation and warnings"""
    # Invalid Copy Change
    with pytest.raises(ValidationError) as e:
        ParsedToCxVarQuery(
            start0=10491132, end0=10535643, copy_change="efo:1234",
            accession="NC_000001.10"
        )
    assert "value is not a valid enumeration member" in str(e.value)

    # NCBI36/hg18 assembly
    rb = ParsedToCxVarQuery(
        start0=2623228, end0=3150942, copy_change=CopyChange.GAIN,
        assembly=ClinVarAssembly.NCBI36, chromosome="chr1",
        untranslatable_returns_text=True
    )
    resp = test_cnv_handler.parsed_to_copy_number(rb)
    assert resp.copy_number_change.type == "Text"
    assert resp.warnings == ["NCBI36 assembly is not currently supported"]

    rb = ParsedToCxVarQuery(
        start0=2623228, end0=3150942, copy_change=CopyChange.GAIN,
        assembly=ClinVarAssembly.HG18, chromosome="chr1",
        untranslatable_returns_text=True
    )
    resp = test_cnv_handler.parsed_to_copy_number(rb)
    assert resp.copy_number_change.type == "Text"
    assert resp.warnings == ["hg18 assembly is not currently supported"]

    # Must give both assembly + chromosome or accession
    ac_assembly_chr_msg = "Must provide either `accession` or both `assembly` and `chromosome`"  # noqa: E501
    with pytest.raises(ValidationError) as e:
        ParsedToCxVarQuery(
            start0=31738809, end0=32217725, copy_change=CopyChange.GAIN,
            assembly="hg38", untranslatable_returns_text=True
        )
    assert ac_assembly_chr_msg in str(e.value)

    # Must give both assembly + chromosome or accession
    with pytest.raises(ValidationError) as e:
        ParsedToCxVarQuery(
            start0=31738809, end0=32217725, copy_change=CopyChange.GAIN,
            chromosome="chr15", untranslatable_returns_text=True
        )
    assert ac_assembly_chr_msg in str(e.value)

    # Must give both assembly + chromosome or accession
    with pytest.raises(ValidationError) as e:
        ParsedToCxVarQuery(
            start0=31738809, end0=32217725, copy_change=CopyChange.GAIN,
            untranslatable_returns_text=True
        )
    assert ac_assembly_chr_msg in str(e.value)

    # invalid chromosome
    with pytest.raises(ValidationError) as e:
        ParsedToCxVarQuery(
            start0=10001, end0=1223133, copy_change=CopyChange.GAIN,
            assembly=ClinVarAssembly.GRCH38, chromosome="z"
        )
    assert "`chromosome`, z, does not match r'^chr(X|Y|([1-9]|1[0-9]|2[0-2]))$'" in str(e.value)  # noqa: E501

    # invalid assembly
    with pytest.raises(ValidationError) as e:
        ParsedToCxVarQuery(
            start0=10001, end0=1223133, copy_change=CopyChange.GAIN, assembly="GRCh99"
        )
    assert "value is not a valid enumeration member" in str(e.value)

    # invalid accession
    rb = ParsedToCxVarQuery(
        start0=10491132, end0=10535643, copy_change=CopyChange.GAIN,
        accession="NC_00002310"
    )
    resp = test_cnv_handler.parsed_to_copy_number(rb)
    assert resp.copy_number_change is None
    assert resp.warnings == \
        ["SeqRepo unable to get translated identifiers for NC_00002310"]

    # Invalid position
    rb = ParsedToCxVarQuery(
        start0=31738809, end0=2302991250, copy_change=CopyChange.GAIN,
        accession="NC_000015.10"
    )
    resp = test_cnv_handler.parsed_to_copy_number(rb)
    assert resp.copy_number_change is None
    assert resp.warnings == ["SeqRepo ValueError: Position out of range (2302991249)"]

    # start must be less than end
    with pytest.raises(ValidationError) as e:
        ParsedToCxVarQuery(
            start0=10001, end0=1223130,
            copy_change=CopyChange.COMPLETE_GENOMIC_LOSS,
            assembly=ClinVarAssembly.GRCH38, chromosome="chrY",
            start_pos_type=VRSTypes.DEFINITE_RANGE,
            end_pos_type=VRSTypes.DEFINITE_RANGE, start1=1223132, end1=1223133
        )
    assert "end positions must be greater than start" in str(e.value)

    # start1 not provided
    with pytest.raises(ValidationError) as e:
        ParsedToCxVarQuery(
            start0=10001, end0=1223130,
            copy_change=CopyChange.COMPLETE_GENOMIC_LOSS,
            assembly=ClinVarAssembly.GRCH38, chromosome="chrY",
            start_pos_type=VRSTypes.DEFINITE_RANGE
        )
    assert "`start1` is required for definite ranges" in str(e.value)

    # copies1 not provided when copies_type is DefiniteRange
    with pytest.raises(ValidationError) as e:
        ParsedToCnVarQuery(
            start0=143134063, end0=143284670, copies0=3,
            copies_type=VRSTypes.DEFINITE_RANGE, assembly=ClinVarAssembly.GRCH37,
            chromosome="chr1", start_pos_type=VRSTypes.INDEFINITE_RANGE,
            start_pos_comparator=Comparator.LT_OR_EQUAL,
            end_pos_type=VRSTypes.INDEFINITE_RANGE,
            end_pos_comparator=Comparator.GT_OR_EQUAL
        )
    assert "`copies1` must be provided for `copies_type == DefiniteRange`" in str(e.value)  # noqa: E501

    # copies_comparator not provided when copies_type is IndefiniteRange
    with pytest.raises(ValidationError) as e:
        ParsedToCnVarQuery(
            start0=143134063, end0=143284670, copies0=3,
            copies_type=VRSTypes.INDEFINITE_RANGE, assembly=ClinVarAssembly.GRCH37,
            chromosome="chr1", start_pos_type=VRSTypes.INDEFINITE_RANGE,
            start_pos_comparator=Comparator.LT_OR_EQUAL,
            end_pos_type=VRSTypes.INDEFINITE_RANGE,
            end_pos_comparator=Comparator.GT_OR_EQUAL
        )
    assert "`copies_comparator` must be provided for `copies_type == IndefiniteRange`" in str(e.value)  # noqa: E501

    # `start_pos_comparator` not provided when start_pos_type is Indefinite Range
    with pytest.raises(ValidationError) as e:
        ParsedToCnVarQuery(
            start0=31738809, end0=32217725, copies0=2, assembly=ClinVarAssembly.GRCH38,
            chromosome="chr15", start_pos_type=VRSTypes.INDEFINITE_RANGE,
            end_pos_type=VRSTypes.NUMBER
        )
    assert "`start_pos_comparator` is required for indefinite ranges" in str(e.value)

    # `end_pos_comparator` not provided when end_pos_type is Indefinite Range
    with pytest.raises(ValidationError) as e:
        ParsedToCnVarQuery(
            start0=31738809, end0=32217725, copies0=2, assembly=ClinVarAssembly.GRCH38,
            chromosome="chr15", start_pos_type=VRSTypes.NUMBER,
            end_pos_type=VRSTypes.INDEFINITE_RANGE
        )
    assert "`end_pos_comparator` is required for indefinite ranges" in str(e.value)
