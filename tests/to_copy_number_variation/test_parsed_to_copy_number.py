"""Test that parsed_to_copy_number works correctly"""

from copy import deepcopy

import pytest
from ga4gh.vrs import models
from pydantic import ValidationError
from tests.conftest import cnv_assertion_checks

from variation.schemas.copy_number_schema import (
    ClinVarAssembly,
    Comparator,
    ParsedPosType,
    ParsedToCnVarQuery,
    ParsedToCxVarQuery,
)
from variation.to_copy_number_variation import ToCopyNumberError


@pytest.fixture(scope="module")
def cn_gain1():
    """Create test fixture for clinvar copy number gain.
    https://www.ncbi.nlm.nih.gov/clinvar/variation/145208/?new_evidence=true
    """
    cn_digest = "pbVk38-x5YGW7yhEtaBnWYjrzcb25L16"
    loc_digest = "6jZXELPqf5JDeN4CpOGde8foTUkHi1jy"
    variation = {
        "type": "CopyNumberCount",
        "id": f"ga4gh:CN.{cn_digest}",
        "digest": cn_digest,
        "location": {
            "type": "SequenceLocation",
            "id": f"ga4gh:SL.{loc_digest}",
            "digest": loc_digest,
            "sequenceReference": {
                "type": "SequenceReference",
                "refgetAccession": "SQ.S_KjnFVz-FE7M0W6yoaUDgYxLPc1jyWU",
            },
            "start": [None, 143134062],
            "end": [143284670, None],
        },
        "copies": 3,
    }
    return models.CopyNumberCount(**variation)


@pytest.fixture(scope="module")
def cn_gain2():
    """Create test fixture for clinvar copy number gain.
    https://www.ncbi.nlm.nih.gov/clinvar/variation/146181/?new_evidence=true
    """
    variation = {
        "type": "CopyNumberCount",
        "location": {
            "type": "SequenceLocation",
            "sequenceReference": {
                "type": "SequenceReference",
                "refgetAccession": "SQ.AsXvWL1-2i5U_buw6_niVIxD6zTbAuS6",
            },
            "start": [None, 31738808],
            "end": [32217725, None],
        },
        "copies": 2,
    }
    return models.CopyNumberCount(**variation)


@pytest.fixture(scope="module")
def cn_gain2_37():
    """Create test fixture for clinvar copy number gain on GRCh37 assembly.
    https://www.ncbi.nlm.nih.gov/clinvar/variation/146181/?new_evidence=true
    """
    variation = {
        "type": "CopyNumberCount",
        "location": {
            "type": "SequenceLocation",
            "sequenceReference": {
                "type": "SequenceReference",
                "refgetAccession": "SQ.zIMZb3Ft7RdWa5XYq0PxIlezLY2ccCgt",
            },
            "start": [None, 32031011],
            "end": [32509926, None],
        },
        "copies": 2,
    }
    return models.CopyNumberCount(**variation)


@pytest.fixture(scope="module")
def cn_loss1():
    """Create test fixture for clinvar copy number loss.
    https://www.ncbi.nlm.nih.gov/clinvar/variation/146181/?new_evidence=true
    """
    variation = {
        "type": "CopyNumberCount",
        "location": {
            "type": "SequenceLocation",
            "sequenceReference": {
                "type": "SequenceReference",
                "refgetAccession": "SQ.v7noePfnNpK8ghYXEqZ9NukMXW7YeNsm",
            },
            "start": [None, 10491131],
            "end": [10535643, None],
        },
        "copies": 1,
    }
    return models.CopyNumberCount(**variation)


@pytest.fixture(scope="module")
def cn_loss2():
    """Create test fixture for clinvar copy number loss.
    https://www.ncbi.nlm.nih.gov/clinvar/variation/148425/?new_evidence=true
    """
    variation = {
        "type": "CopyNumberCount",
        "location": {
            "type": "SequenceLocation",
            "sequenceReference": {
                "type": "SequenceReference",
                "refgetAccession": "SQ.8_liLu1aycC0tPQPFmUaGXJLDs5SbPZ5",
            },
            "start": [None, 10000],
            "end": [1223133, None],
        },
        "copies": 0,
    }
    return models.CopyNumberCount(**variation)


@pytest.fixture(scope="module")
def cn_definite_number():
    """Create test fixture for copy number count using definite range for start and
    number for end
    """
    variation = {
        "type": "CopyNumberCount",
        "location": {
            "type": "SequenceLocation",
            "sequenceReference": {
                "type": "SequenceReference",
                "refgetAccession": "SQ.S_KjnFVz-FE7M0W6yoaUDgYxLPc1jyWU",
            },
            "start": [143134062, 143134064],
            "end": 143284670,
        },
        "copies": 3,
    }
    return models.CopyNumberCount(**variation)


@pytest.fixture(scope="module")
def cx_numbers():
    """Create test fixture for copy number change using numbers for start and end"""
    cx_digest = "5kaJC-7Jj851bfJ6EipsHV413feg1T4T"
    loc_digest = "Iz_azSFTEulx7tCluLgGhE1n0hTLUocb"
    variation = {
        "type": "CopyNumberChange",
        "id": f"ga4gh:CX.{cx_digest}",
        "digest": cx_digest,
        "location": {
            "type": "SequenceLocation",
            "id": f"ga4gh:SL.{loc_digest}",
            "digest": loc_digest,
            "sequenceReference": {
                "type": "SequenceReference",
                "refgetAccession": "SQ.8_liLu1aycC0tPQPFmUaGXJLDs5SbPZ5",
            },
            "start": 10000,
            "end": 1223133,
        },
        "copyChange": "efo:0030069",
    }
    return models.CopyNumberChange(**variation)


@pytest.fixture(scope="module")
def cx_definite_ranges():
    """Create test fixture for copy number change using definite ranges for start and
    end
    """
    variation = {
        "type": "CopyNumberChange",
        "location": {
            "type": "SequenceLocation",
            "sequenceReference": {
                "type": "SequenceReference",
                "refgetAccession": "SQ.8_liLu1aycC0tPQPFmUaGXJLDs5SbPZ5",
            },
            "start": [10000, 10005],
            "end": [1223130, 1223133],
        },
        "copyChange": "efo:0030069",
    }
    return models.CopyNumberChange(**variation)


@pytest.fixture(scope="module")
def cx_indefinite_ranges():
    """Create test fixture for copy number change using indefinite ranges for start and
    end
    """
    variation = {
        "type": "CopyNumberChange",
        "location": {
            "type": "SequenceLocation",
            "sequenceReference": {
                "type": "SequenceReference",
                "refgetAccession": "SQ.8_liLu1aycC0tPQPFmUaGXJLDs5SbPZ5",
            },
            "start": [None, 10000],
            "end": [1223130, None],
        },
        "copyChange": "efo:0030069",
    }
    return models.CopyNumberChange(**variation)


@pytest.fixture(scope="module")
def cx_number_indefinite():
    """Create test fixture for copy number change using number for start and indefinite
    range for end
    """
    variation = {
        "type": "CopyNumberChange",
        "location": {
            "type": "SequenceLocation",
            "sequenceReference": {
                "type": "SequenceReference",
                "refgetAccession": "SQ.8_liLu1aycC0tPQPFmUaGXJLDs5SbPZ5",
            },
            "start": 10000,
            "end": [1223130, None],
        },
        "copyChange": "efo:0030069",
    }
    return models.CopyNumberChange(**variation)


GRCH37_CHR7_VRS_ID = "ga4gh:SQ.IW78mgV5Cqf6M24hy52hPjyyo5tCCd86"
GRCH38_CHR7_VRS_ID = "ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul"


def test_get_parsed_ac(test_cnv_handler):
    """Test that _get_parsed_ac works correctly"""
    for assembly in [ClinVarAssembly.GRCH37, ClinVarAssembly.HG19]:
        resp = test_cnv_handler._get_parsed_ac(assembly, "chr7", use_grch38=False)
        assert resp.lifted_over is False
        assert resp.accession == GRCH37_CHR7_VRS_ID

        resp = test_cnv_handler._get_parsed_ac(assembly, "chr7", use_grch38=True)
        assert resp.lifted_over is True
        assert resp.accession == GRCH38_CHR7_VRS_ID

    for assembly in [ClinVarAssembly.GRCH38, ClinVarAssembly.HG38]:
        resp = test_cnv_handler._get_parsed_ac(assembly, "chr7", use_grch38=False)
        assert resp.lifted_over is False
        assert resp.accession == GRCH38_CHR7_VRS_ID

        resp = test_cnv_handler._get_parsed_ac(assembly, "chr7", use_grch38=True)
        assert resp.lifted_over is False
        assert resp.accession == GRCH38_CHR7_VRS_ID

    with pytest.raises(ToCopyNumberError) as e:
        test_cnv_handler._get_parsed_ac(
            ClinVarAssembly.NCBI36, "chr7", use_grch38=False
        )
    assert str(e.value) == "NCBI36 assembly is not currently supported"

    with pytest.raises(ToCopyNumberError) as e:
        test_cnv_handler._get_parsed_ac(ClinVarAssembly.HG18, "chr7", use_grch38=False)
    assert str(e.value) == "hg18 assembly is not currently supported"


def test_get_parsed_ac_chr(test_cnv_handler):
    """Test that _get_parsed_ac_chr works correctly"""
    resp = test_cnv_handler._get_parsed_ac_chr("NC_000007.13", False)
    assert resp.accession == GRCH37_CHR7_VRS_ID
    assert resp.chromosome == "chr7"
    assert resp.lifted_over is False

    resp = test_cnv_handler._get_parsed_ac_chr("NC_000007.13", True)
    assert resp.accession == GRCH38_CHR7_VRS_ID
    assert resp.chromosome == "chr7"
    assert resp.lifted_over is True

    for do_liftover in [True, False]:
        resp = test_cnv_handler._get_parsed_ac_chr("NC_000007.14", do_liftover)
        assert resp.accession == GRCH38_CHR7_VRS_ID
        assert resp.chromosome == "chr7"
        assert resp.lifted_over is False

    # if genomic ac not provided
    with pytest.raises(ToCopyNumberError) as e:
        test_cnv_handler._get_parsed_ac_chr("NP_000542.1", False)
    assert str(e.value) == "Not a supported genomic accession: NP_000542.1"

    # invalid accession
    with pytest.raises(ToCopyNumberError) as e:
        test_cnv_handler._get_parsed_ac_chr("NC_00000713", False)
    assert (
        str(e.value) == "SeqRepo unable to get translated identifiers for NC_00000713"
    )


def test_validate_pos(test_cnv_handler):
    """Test that _validate_ac_pos works correctly"""
    resp = test_cnv_handler._validate_ac_pos("NC_000007.14", 140753336)
    assert resp is None

    # invalid accession
    with pytest.raises(ToCopyNumberError) as e:
        test_cnv_handler._validate_ac_pos("NC_00000714", 140753336)
    assert str(e.value) == "Accession not found in SeqRepo: NC_00000714"

    # invalid pos
    with pytest.raises(ToCopyNumberError) as e:
        test_cnv_handler._validate_ac_pos("NC_000007.14", 159345975)
    assert str(e.value) == "Position (159345975) is not valid on NC_000007.14"

    # invalid pos
    with pytest.raises(ToCopyNumberError) as e:
        test_cnv_handler._validate_ac_pos("NC_000007.14", 99999999999)
    assert str(e.value) == "SeqRepo ValueError: Position out of range (99999999998)"


def test_get_vrs_loc_start_or_end(test_cnv_handler):
    """Test that _get_vrs_loc_start_or_end works correctly"""
    ac = "NC_000007.14"
    pos0 = 140753336
    pos1 = 140753350

    # Number start
    resp = test_cnv_handler._get_vrs_loc_start_or_end(
        ac, pos0, ParsedPosType.NUMBER, is_start=True
    )
    assert resp == 140753335

    # Number end
    resp = test_cnv_handler._get_vrs_loc_start_or_end(
        ac, pos0, ParsedPosType.NUMBER, is_start=False
    )
    assert resp == 140753336

    # Definite Range start
    resp = test_cnv_handler._get_vrs_loc_start_or_end(
        ac, pos0, ParsedPosType.DEFINITE_RANGE, is_start=True, pos1=pos1
    )
    assert resp == models.Range([140753335, 140753349])

    # Definite Range end
    resp = test_cnv_handler._get_vrs_loc_start_or_end(
        ac, pos0, ParsedPosType.DEFINITE_RANGE, is_start=False, pos1=pos1
    )
    assert resp == models.Range([pos0, pos1])

    # Indefinite Range start
    resp = test_cnv_handler._get_vrs_loc_start_or_end(
        ac,
        pos0,
        ParsedPosType.INDEFINITE_RANGE,
        is_start=True,
        comparator=Comparator.LT_OR_EQUAL,
    )
    assert resp == models.Range([None, 140753335])

    # Indefinite Range end
    resp = test_cnv_handler._get_vrs_loc_start_or_end(
        ac,
        pos0,
        ParsedPosType.INDEFINITE_RANGE,
        is_start=False,
        comparator=Comparator.GT_OR_EQUAL,
    )
    assert resp == models.Range([140753336, None])


def test_liftover_pos(test_cnv_handler):
    """Test that _liftover_pos works correctly"""
    resp = test_cnv_handler._liftover_pos("chr7", 140453136, 140453137, None, None)
    assert resp == {
        "start0": 140753336,
        "end0": 140753337,
        "start1": None,
        "end1": None,
    }

    resp = test_cnv_handler._liftover_pos(
        "chr7", 140453136, 140453137, 140453138, 140453139
    )
    assert resp == {
        "start0": 140753336,
        "end0": 140753337,
        "start1": 140753338,
        "end1": 140753339,
    }

    # invalid pos
    with pytest.raises(ToCopyNumberError) as e:
        test_cnv_handler._liftover_pos("chr7", 159345975, 159345976, None, None)
    assert str(e.value) == "Unable to liftover: chr7 with pos 159345975"


def test_parsed_copy_number_gain(test_cnv_handler, cn_gain1, cn_gain2, cn_gain2_37):
    """Test that parsed_to_copy_number works for parsed copy number gain queries"""
    # https://www.ncbi.nlm.nih.gov/clinvar/variation/145208/?new_evidence=true
    rb = ParsedToCnVarQuery(
        start0=143134063,
        end0=143284670,
        copies0=3,
        assembly=ClinVarAssembly.GRCH37,
        chromosome="chr1",
        start_pos_type=ParsedPosType.INDEFINITE_RANGE,
        start_pos_comparator=Comparator.LT_OR_EQUAL,
        end_pos_type=ParsedPosType.INDEFINITE_RANGE,
        end_pos_comparator=Comparator.GT_OR_EQUAL,
    )
    resp = test_cnv_handler.parsed_to_copy_number(rb)
    cnv_assertion_checks(resp, cn_gain1, check_vrs_id=True)

    rb = ParsedToCnVarQuery(
        start0=143134063,
        end0=143284670,
        copies0=3,
        assembly=ClinVarAssembly.HG19,
        chromosome="chr1",
        start_pos_type=ParsedPosType.INDEFINITE_RANGE,
        start_pos_comparator=Comparator.LT_OR_EQUAL,
        end_pos_type=ParsedPosType.INDEFINITE_RANGE,
        end_pos_comparator=Comparator.GT_OR_EQUAL,
    )
    resp = test_cnv_handler.parsed_to_copy_number(rb)
    cnv_assertion_checks(resp, cn_gain1)

    rb = ParsedToCnVarQuery(
        start0=143134063,
        end0=143284670,
        copies0=3,
        accession="NC_000001.10",
        start_pos_type=ParsedPosType.INDEFINITE_RANGE,
        start_pos_comparator=Comparator.LT_OR_EQUAL,
        end_pos_type=ParsedPosType.INDEFINITE_RANGE,
        end_pos_comparator=Comparator.GT_OR_EQUAL,
    )
    resp = test_cnv_handler.parsed_to_copy_number(rb)
    cnv_assertion_checks(resp, cn_gain1)

    # https://www.ncbi.nlm.nih.gov/clinvar/variation/146181/?new_evidence=true
    # 38
    rb = ParsedToCnVarQuery(
        start0=31738809,
        end0=32217725,
        copies0=2,
        assembly=ClinVarAssembly.GRCH38,
        chromosome="chr15",
        start_pos_type=ParsedPosType.INDEFINITE_RANGE,
        start_pos_comparator=Comparator.LT_OR_EQUAL,
        end_pos_type=ParsedPosType.INDEFINITE_RANGE,
        end_pos_comparator=Comparator.GT_OR_EQUAL,
    )
    resp = test_cnv_handler.parsed_to_copy_number(rb)
    cnv_assertion_checks(resp, cn_gain2)

    # 38 with liftover (shouldnt do anything)
    rb = ParsedToCnVarQuery(
        start0=31738809,
        end0=32217725,
        copies0=2,
        assembly=ClinVarAssembly.GRCH38,
        chromosome="chr15",
        do_liftover=True,
        start_pos_type=ParsedPosType.INDEFINITE_RANGE,
        start_pos_comparator=Comparator.LT_OR_EQUAL,
        end_pos_type=ParsedPosType.INDEFINITE_RANGE,
        end_pos_comparator=Comparator.GT_OR_EQUAL,
    )
    resp = test_cnv_handler.parsed_to_copy_number(rb)
    cnv_assertion_checks(resp, cn_gain2)

    # 38 with liftover (shouldnt do anything)
    rb = ParsedToCnVarQuery(
        start0=31738809,
        end0=32217725,
        copies0=2,
        assembly=ClinVarAssembly.HG38,
        chromosome="chr15",
        do_liftover=True,
        start_pos_type=ParsedPosType.INDEFINITE_RANGE,
        start_pos_comparator=Comparator.LT_OR_EQUAL,
        end_pos_type=ParsedPosType.INDEFINITE_RANGE,
        end_pos_comparator=Comparator.GT_OR_EQUAL,
    )
    resp = test_cnv_handler.parsed_to_copy_number(rb)
    cnv_assertion_checks(resp, cn_gain2)

    # 38
    rb = ParsedToCnVarQuery(
        start0=31738809,
        end0=32217725,
        copies0=2,
        assembly=ClinVarAssembly.HG38,
        chromosome="chr15",
        start_pos_type=ParsedPosType.INDEFINITE_RANGE,
        start_pos_comparator=Comparator.LT_OR_EQUAL,
        end_pos_type=ParsedPosType.INDEFINITE_RANGE,
        end_pos_comparator=Comparator.GT_OR_EQUAL,
    )
    resp = test_cnv_handler.parsed_to_copy_number(rb)
    cnv_assertion_checks(resp, cn_gain2)

    # 38 accession
    rb = ParsedToCnVarQuery(
        start0=31738809,
        end0=32217725,
        copies0=2,
        accession="NC_000015.10",
        start_pos_type=ParsedPosType.INDEFINITE_RANGE,
        start_pos_comparator=Comparator.LT_OR_EQUAL,
        end_pos_type=ParsedPosType.INDEFINITE_RANGE,
        end_pos_comparator=Comparator.GT_OR_EQUAL,
    )
    resp = test_cnv_handler.parsed_to_copy_number(rb)
    cnv_assertion_checks(resp, cn_gain2)

    # 38 accession with liftover (shouldnt do anything)
    rb = ParsedToCnVarQuery(
        start0=31738809,
        end0=32217725,
        copies0=2,
        accession="NC_000015.10",
        start_pos_type=ParsedPosType.INDEFINITE_RANGE,
        start_pos_comparator=Comparator.LT_OR_EQUAL,
        end_pos_type=ParsedPosType.INDEFINITE_RANGE,
        end_pos_comparator=Comparator.GT_OR_EQUAL,
    )
    resp = test_cnv_handler.parsed_to_copy_number(rb)
    cnv_assertion_checks(resp, cn_gain2)

    # 37 with liftover
    rb = ParsedToCnVarQuery(
        start0=32031012,
        end0=32509926,
        copies0=2,
        accession="NC_000015.9",
        do_liftover=True,
        start_pos_type=ParsedPosType.INDEFINITE_RANGE,
        start_pos_comparator=Comparator.LT_OR_EQUAL,
        end_pos_type=ParsedPosType.INDEFINITE_RANGE,
        end_pos_comparator=Comparator.GT_OR_EQUAL,
    )
    resp = test_cnv_handler.parsed_to_copy_number(rb)
    cnv_assertion_checks(resp, cn_gain2)

    # 37 chr+accession with liftover
    rb = ParsedToCnVarQuery(
        start0=32031012,
        end0=32509926,
        copies0=2,
        chromosome="chr15",
        assembly=ClinVarAssembly.GRCH37,
        do_liftover=True,
        start_pos_type=ParsedPosType.INDEFINITE_RANGE,
        start_pos_comparator=Comparator.LT_OR_EQUAL,
        end_pos_type=ParsedPosType.INDEFINITE_RANGE,
        end_pos_comparator=Comparator.GT_OR_EQUAL,
    )
    resp = test_cnv_handler.parsed_to_copy_number(rb)
    cnv_assertion_checks(resp, cn_gain2)

    # 37 with no liftover
    rb = ParsedToCnVarQuery(
        start0=32031012,
        end0=32509926,
        copies0=2,
        accession="NC_000015.9",
        start_pos_type=ParsedPosType.INDEFINITE_RANGE,
        start_pos_comparator=Comparator.LT_OR_EQUAL,
        end_pos_type=ParsedPosType.INDEFINITE_RANGE,
        end_pos_comparator=Comparator.GT_OR_EQUAL,
    )
    resp = test_cnv_handler.parsed_to_copy_number(rb)
    cnv_assertion_checks(resp, cn_gain2_37)

    # 37 chr+accession with no liftover
    rb = ParsedToCnVarQuery(
        start0=32031012,
        end0=32509926,
        copies0=2,
        chromosome="chr15",
        assembly=ClinVarAssembly.GRCH37,
        start_pos_type=ParsedPosType.INDEFINITE_RANGE,
        start_pos_comparator=Comparator.LT_OR_EQUAL,
        end_pos_type=ParsedPosType.INDEFINITE_RANGE,
        end_pos_comparator=Comparator.GT_OR_EQUAL,
    )
    resp = test_cnv_handler.parsed_to_copy_number(rb)
    cnv_assertion_checks(resp, cn_gain2_37)


def test_parsed_copy_number_loss(test_cnv_handler, cn_loss1, cn_loss2):
    """Test that parsed_to_copy_number works for parsed copy number loss queries"""
    # https://www.ncbi.nlm.nih.gov/clinvar/variation/1299222/?new_evidence=true
    rb = ParsedToCnVarQuery(
        start0=10491132,
        end0=10535643,
        copies0=1,
        assembly=ClinVarAssembly.GRCH37,
        chromosome="chrX",
        start_pos_type=ParsedPosType.INDEFINITE_RANGE,
        start_pos_comparator=Comparator.LT_OR_EQUAL,
        end_pos_type=ParsedPosType.INDEFINITE_RANGE,
        end_pos_comparator=Comparator.GT_OR_EQUAL,
    )
    resp = test_cnv_handler.parsed_to_copy_number(rb)
    cnv_assertion_checks(resp, cn_loss1)

    rb = ParsedToCnVarQuery(
        start0=10491132,
        end0=10535643,
        copies0=1,
        assembly=ClinVarAssembly.HG19,
        chromosome="chrX",
        start_pos_type=ParsedPosType.INDEFINITE_RANGE,
        start_pos_comparator=Comparator.LT_OR_EQUAL,
        end_pos_type=ParsedPosType.INDEFINITE_RANGE,
        end_pos_comparator=Comparator.GT_OR_EQUAL,
    )
    resp = test_cnv_handler.parsed_to_copy_number(rb)
    cnv_assertion_checks(resp, cn_loss1)

    rb = ParsedToCnVarQuery(
        start0=10491132,
        end0=10535643,
        copies0=1,
        accession="NC_000023.10",
        start_pos_type=ParsedPosType.INDEFINITE_RANGE,
        start_pos_comparator=Comparator.LT_OR_EQUAL,
        end_pos_type=ParsedPosType.INDEFINITE_RANGE,
        end_pos_comparator=Comparator.GT_OR_EQUAL,
    )
    resp = test_cnv_handler.parsed_to_copy_number(rb)
    cnv_assertion_checks(resp, cn_loss1)

    # https://www.ncbi.nlm.nih.gov/clinvar/variation/148425/?new_evidence=true
    rb = ParsedToCnVarQuery(
        start0=10001,
        end0=1223133,
        copies0=0,
        assembly=ClinVarAssembly.GRCH38,
        chromosome="chrY",
        start_pos_type=ParsedPosType.INDEFINITE_RANGE,
        start_pos_comparator=Comparator.LT_OR_EQUAL,
        end_pos_type=ParsedPosType.INDEFINITE_RANGE,
        end_pos_comparator=Comparator.GT_OR_EQUAL,
    )
    resp = test_cnv_handler.parsed_to_copy_number(rb)
    cnv_assertion_checks(resp, cn_loss2)

    rb = ParsedToCnVarQuery(
        start0=10001,
        end0=1223133,
        copies0=0,
        assembly=ClinVarAssembly.HG38,
        chromosome="chrY",
        start_pos_type=ParsedPosType.INDEFINITE_RANGE,
        start_pos_comparator=Comparator.LT_OR_EQUAL,
        end_pos_type=ParsedPosType.INDEFINITE_RANGE,
        end_pos_comparator=Comparator.GT_OR_EQUAL,
    )
    resp = test_cnv_handler.parsed_to_copy_number(rb)
    cnv_assertion_checks(resp, cn_loss2)

    rb = ParsedToCnVarQuery(
        start0=10001,
        end0=1223133,
        copies0=0,
        accession="NC_000024.10",
        start_pos_type=ParsedPosType.INDEFINITE_RANGE,
        start_pos_comparator=Comparator.LT_OR_EQUAL,
        end_pos_type=ParsedPosType.INDEFINITE_RANGE,
        end_pos_comparator=Comparator.GT_OR_EQUAL,
    )
    resp = test_cnv_handler.parsed_to_copy_number(rb)
    cnv_assertion_checks(resp, cn_loss2)


def test_to_parsed_cn_var(test_cnv_handler, cn_definite_number):
    """Test that parsed_to_copy_number works correctly for copy number count"""
    # start uses definite and end uses number
    rb = ParsedToCnVarQuery(
        start0=143134063,
        end0=143284670,
        copies0=3,
        assembly=ClinVarAssembly.GRCH37,
        chromosome="chr1",
        start_pos_type=ParsedPosType.DEFINITE_RANGE,
        start1=143134065,
    )
    resp = test_cnv_handler.parsed_to_copy_number(rb)
    cnv_assertion_checks(resp, cn_definite_number)

    # copies is definite range
    rb = ParsedToCnVarQuery(
        start0=143134063,
        end0=143284670,
        copies0=3,
        copies1=5,
        copies_type=ParsedPosType.DEFINITE_RANGE,
        assembly=ClinVarAssembly.GRCH37,
        chromosome="chr1",
        start_pos_type=ParsedPosType.DEFINITE_RANGE,
        start1=143134065,
    )
    resp = test_cnv_handler.parsed_to_copy_number(rb)
    expected = deepcopy(cn_definite_number)
    expected.copies = models.Range([3, 5])
    cnv_assertion_checks(resp, expected)

    # copies is indefinite range <=
    rb = ParsedToCnVarQuery(
        start0=143134063,
        end0=143284670,
        copies0=3,
        copies_comparator=Comparator.LT_OR_EQUAL,
        copies_type=ParsedPosType.INDEFINITE_RANGE,
        assembly=ClinVarAssembly.GRCH37,
        chromosome="chr1",
        start_pos_type=ParsedPosType.DEFINITE_RANGE,
        start1=143134065,
    )
    resp = test_cnv_handler.parsed_to_copy_number(rb)
    expected = deepcopy(cn_definite_number)
    expected.copies = models.Range([None, 3])
    cnv_assertion_checks(resp, expected)

    # copies is indefinite range >=
    rb = ParsedToCnVarQuery(
        start0=143134063,
        end0=143284670,
        copies0=3,
        copies_comparator=Comparator.GT_OR_EQUAL,
        copies_type=ParsedPosType.INDEFINITE_RANGE,
        assembly=ClinVarAssembly.GRCH37,
        chromosome="chr1",
        start_pos_type=ParsedPosType.DEFINITE_RANGE,
        start1=143134065,
    )
    resp = test_cnv_handler.parsed_to_copy_number(rb)
    expected = deepcopy(cn_definite_number)
    expected.copies = [3, None]
    cnv_assertion_checks(resp, expected)

    # start_pos and end_pos indefinite range
    rb = ParsedToCnVarQuery(
        start0=143134063,
        end0=143284670,
        copies0=3,
        assembly=ClinVarAssembly.GRCH37,
        chromosome="chr1",
        start_pos_type=ParsedPosType.INDEFINITE_RANGE,
        end_pos_type=ParsedPosType.INDEFINITE_RANGE,
        start_pos_comparator=Comparator.GT_OR_EQUAL,
        end_pos_comparator=Comparator.LT_OR_EQUAL,
    )
    resp = test_cnv_handler.parsed_to_copy_number(rb)
    cnc = resp.copy_number_count.model_dump(exclude_none=True)
    cn_digest = cnc.pop("digest")
    assert cnc.pop("id") == f"ga4gh:CN.{cn_digest}"
    loc_digest = cnc["location"].pop("digest")
    assert cnc["location"].pop("id") == f"ga4gh:SL.{loc_digest}"
    assert cnc == {
        "type": "CopyNumberCount",
        "location": {
            "type": "SequenceLocation",
            "sequenceReference": {
                "type": "SequenceReference",
                "refgetAccession": "SQ.S_KjnFVz-FE7M0W6yoaUDgYxLPc1jyWU",
            },
            "start": [143134062, None],
            "end": [None, 143284670],
        },
        "copies": 3,
    }


def test_parsed_to_cx_var(
    test_cnv_handler,
    cx_numbers,
    cx_definite_ranges,
    cx_indefinite_ranges,
    cx_number_indefinite,
):
    """Test that parsed_to_copy_number works for copy number change"""
    # start and end use number
    rb = ParsedToCxVarQuery(
        start0=10001,
        end0=1223133,
        copy_change=models.CopyChange.EFO_0030069,
        assembly=ClinVarAssembly.GRCH38,
        chromosome="chrY",
        start_pos_type=ParsedPosType.NUMBER,
        end_pos_type=ParsedPosType.NUMBER,
    )
    resp = test_cnv_handler.parsed_to_copy_number(rb)
    cnv_assertion_checks(resp, cx_numbers, check_vrs_id=True)

    # start and end use definite ranges
    rb = ParsedToCxVarQuery(
        start0=10001,
        end0=1223130,
        copy_change=models.CopyChange.EFO_0030069,
        assembly=ClinVarAssembly.GRCH38,
        chromosome="chrY",
        start_pos_type=ParsedPosType.DEFINITE_RANGE,
        end_pos_type=ParsedPosType.DEFINITE_RANGE,
        start1=10006,
        end1=1223133,
    )
    resp = test_cnv_handler.parsed_to_copy_number(rb)
    cnv_assertion_checks(resp, cx_definite_ranges)

    # start and end use indefinite ranges
    rb = ParsedToCxVarQuery(
        start0=10001,
        end0=1223130,
        copy_change=models.CopyChange.EFO_0030069,
        assembly=ClinVarAssembly.GRCH38,
        chromosome="chrY",
        start_pos_type=ParsedPosType.INDEFINITE_RANGE,
        start_pos_comparator=Comparator.LT_OR_EQUAL,
        end_pos_type=ParsedPosType.INDEFINITE_RANGE,
        end_pos_comparator=Comparator.GT_OR_EQUAL,
    )
    resp = test_cnv_handler.parsed_to_copy_number(rb)
    cnv_assertion_checks(resp, cx_indefinite_ranges)

    # start uses number and end use indefinite range
    rb = ParsedToCxVarQuery(
        start0=10001,
        end0=1223130,
        copy_change=models.CopyChange.EFO_0030069,
        assembly=ClinVarAssembly.GRCH38,
        chromosome="chrY",
        start_pos_type=ParsedPosType.NUMBER,
        end_pos_type=ParsedPosType.INDEFINITE_RANGE,
        end_pos_comparator=Comparator.GT_OR_EQUAL,
    )
    resp = test_cnv_handler.parsed_to_copy_number(rb)
    cnv_assertion_checks(resp, cx_number_indefinite)


def test_invalid(test_cnv_handler):
    """Test invalid copy number queries returns no variation and warnings"""
    # Invalid Copy Change
    with pytest.raises(ValidationError) as e:
        ParsedToCxVarQuery(
            start0=10491132,
            end0=10535643,
            copy_change="efo:1234",
            accession="NC_000001.10",
        )
    assert "Input should be 'efo:" in str(e.value)

    # NCBI36/hg18 assembly
    rb = ParsedToCxVarQuery(
        start0=2623228,
        end0=3150942,
        copy_change=models.CopyChange.EFO_0030070,
        assembly=ClinVarAssembly.NCBI36,
        chromosome="chr1",
    )
    resp = test_cnv_handler.parsed_to_copy_number(rb)
    assert resp.copy_number_change is None
    assert resp.warnings == ["NCBI36 assembly is not currently supported"]

    rb = ParsedToCxVarQuery(
        start0=2623228,
        end0=3150942,
        copy_change=models.CopyChange.EFO_0030070,
        assembly=ClinVarAssembly.HG18,
        chromosome="chr1",
    )
    resp = test_cnv_handler.parsed_to_copy_number(rb)
    assert resp.copy_number_change is None
    assert resp.warnings == ["hg18 assembly is not currently supported"]

    # Must give both assembly + chromosome or accession
    ac_assembly_chr_msg = (
        "Must provide either `accession` or both `assembly` and `chromosome`"
    )
    with pytest.raises(ValidationError) as e:
        ParsedToCxVarQuery(
            start0=31738809,
            end0=32217725,
            copy_change=models.CopyChange.EFO_0030070,
            assembly="hg38",
        )
    assert ac_assembly_chr_msg in str(e.value)

    # Must give both assembly + chromosome or accession
    with pytest.raises(ValidationError) as e:
        ParsedToCxVarQuery(
            start0=31738809,
            end0=32217725,
            copy_change=models.CopyChange.EFO_0030070,
            chromosome="chr15",
        )
    assert ac_assembly_chr_msg in str(e.value)

    # Must give both assembly + chromosome or accession
    with pytest.raises(ValidationError) as e:
        ParsedToCxVarQuery(
            start0=31738809,
            end0=32217725,
            copy_change=models.CopyChange.EFO_0030070,
        )
    assert ac_assembly_chr_msg in str(e.value)

    # invalid chromosome
    with pytest.raises(ValidationError) as e:
        ParsedToCxVarQuery(
            start0=10001,
            end0=1223133,
            copy_change=models.CopyChange.EFO_0030070,
            assembly=ClinVarAssembly.GRCH38,
            chromosome="z",
        )
    assert "`chromosome`, z, does not match r'^chr(X|Y|([1-9]|1[0-9]|2[0-2]))$'" in str(
        e.value
    )

    # invalid assembly
    with pytest.raises(ValidationError) as e:
        ParsedToCxVarQuery(
            start0=10001,
            end0=1223133,
            copy_change=models.CopyChange.EFO_0030070,
            assembly="GRCh99",
        )
    assert "Input should be 'GRCh38'," in str(e.value)

    # invalid accession
    rb = ParsedToCxVarQuery(
        start0=10491132,
        end0=10535643,
        copy_change=models.CopyChange.EFO_0030070,
        accession="NC_00002310",
    )
    resp = test_cnv_handler.parsed_to_copy_number(rb)
    assert resp.copy_number_change is None
    assert resp.warnings == [
        "SeqRepo unable to get translated identifiers for NC_00002310"
    ]

    # Invalid position
    rb = ParsedToCxVarQuery(
        start0=31738809,
        end0=2302991250,
        copy_change=models.CopyChange.EFO_0030070,
        accession="NC_000015.10",
    )
    resp = test_cnv_handler.parsed_to_copy_number(rb)
    assert resp.copy_number_change is None
    assert resp.warnings == ["SeqRepo ValueError: Position out of range (2302991249)"]

    # start must be less than end
    with pytest.raises(ValidationError) as e:
        ParsedToCxVarQuery(
            start0=10001,
            end0=1223130,
            copy_change=models.CopyChange.EFO_0030069,
            assembly=ClinVarAssembly.GRCH38,
            chromosome="chrY",
            start_pos_type=ParsedPosType.DEFINITE_RANGE,
            end_pos_type=ParsedPosType.DEFINITE_RANGE,
            start1=1223132,
            end1=1223133,
        )
    assert "end positions must be greater than start" in str(e.value)

    # start1 not provided
    with pytest.raises(ValidationError) as e:
        ParsedToCxVarQuery(
            start0=10001,
            end0=1223130,
            copy_change=models.CopyChange.EFO_0030069,
            assembly=ClinVarAssembly.GRCH38,
            chromosome="chrY",
            start_pos_type=ParsedPosType.DEFINITE_RANGE,
        )
    assert "`start1` is required for definite ranges" in str(e.value)

    # copies1 not provided when copies_type is DefiniteRange
    with pytest.raises(ValidationError) as e:
        ParsedToCnVarQuery(
            start0=143134063,
            end0=143284670,
            copies0=3,
            copies_type=ParsedPosType.DEFINITE_RANGE,
            assembly=ClinVarAssembly.GRCH37,
            chromosome="chr1",
            start_pos_type=ParsedPosType.INDEFINITE_RANGE,
            start_pos_comparator=Comparator.LT_OR_EQUAL,
            end_pos_type=ParsedPosType.INDEFINITE_RANGE,
            end_pos_comparator=Comparator.GT_OR_EQUAL,
        )
    assert (
        "`copies1` must be provided for `copies_type == ParsedPosType.DEFINITE_RANGE`"
        in str(e.value)
    )

    # copies_comparator not provided when copies_type is IndefiniteRange
    with pytest.raises(ValidationError) as e:
        ParsedToCnVarQuery(
            start0=143134063,
            end0=143284670,
            copies0=3,
            copies_type=ParsedPosType.INDEFINITE_RANGE,
            assembly=ClinVarAssembly.GRCH37,
            chromosome="chr1",
            start_pos_type=ParsedPosType.INDEFINITE_RANGE,
            start_pos_comparator=Comparator.LT_OR_EQUAL,
            end_pos_type=ParsedPosType.INDEFINITE_RANGE,
            end_pos_comparator=Comparator.GT_OR_EQUAL,
        )
    assert (
        "`copies_comparator` must be provided for `copies_type == ParsedPosType.INDEFINITE_RANGE`"
        in str(e.value)
    )

    # `start_pos_comparator` not provided when start_pos_type is Indefinite Range
    with pytest.raises(ValidationError) as e:
        ParsedToCnVarQuery(
            start0=31738809,
            end0=32217725,
            copies0=2,
            assembly=ClinVarAssembly.GRCH38,
            chromosome="chr15",
            start_pos_type=ParsedPosType.INDEFINITE_RANGE,
            end_pos_type=ParsedPosType.NUMBER,
        )
    assert "`start_pos_comparator` is required for indefinite ranges" in str(e.value)

    # `end_pos_comparator` not provided when end_pos_type is Indefinite Range
    with pytest.raises(ValidationError) as e:
        ParsedToCnVarQuery(
            start0=31738809,
            end0=32217725,
            copies0=2,
            assembly=ClinVarAssembly.GRCH38,
            chromosome="chr15",
            start_pos_type=ParsedPosType.NUMBER,
            end_pos_type=ParsedPosType.INDEFINITE_RANGE,
        )
    assert "`end_pos_comparator` is required for indefinite ranges" in str(e.value)
