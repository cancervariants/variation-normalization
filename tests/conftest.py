"""Create methods used throughout tests."""

import asyncio
import contextlib
import logging

import pytest
from cool_seq_tool.app import CoolSeqTool
from ga4gh.vrs import models
from gene.database.dynamodb import DynamoDbDatabase
from gene.query import QueryHandler as GeneQueryHandler

from variation.classify import Classify
from variation.query import QueryHandler
from variation.schemas.normalize_response_schema import NormalizeService
from variation.tokenize import Tokenize
from variation.tokenizers import GeneSymbol


def pytest_addoption(parser):
    """Add custom commands to pytest invocation.
    See https://docs.pytest.org/en/7.1.x/reference/reference.html#parser
    """
    parser.addoption(
        "--verbose-logs",
        action="store_true",
        default=False,
        help="show noisy module logs",
    )


def pytest_configure(config):
    """Configure pytest setup."""
    if not config.getoption("--verbose-logs"):
        logging.getLogger("cool_seq_tool").setLevel(logging.INFO)
        logging.getLogger("boto3").setLevel(logging.ERROR)
        logging.getLogger("botocore").setLevel(logging.ERROR)
        logging.getLogger("urllib3.connectionpool").setLevel(logging.ERROR)
        logging.getLogger("hgvs.parser").setLevel(logging.INFO)
        logging.getLogger("biocommons.seqrepo.seqaliasdb.seqaliasdb").setLevel(
            logging.INFO
        )
        logging.getLogger("biocommons.seqrepo.fastadir.fastadir").setLevel(logging.INFO)
        logging.getLogger("asyncio").setLevel(logging.INFO)


@pytest.fixture(scope="session")
def event_loop():
    """Create an instance of the default event loop for each test case."""
    loop = asyncio.get_event_loop_policy().new_event_loop()
    yield loop
    loop.close()


@pytest.fixture(scope="session")
def test_tokenizer():
    """Create test fixture for tokenizer"""
    return Tokenize(GeneSymbol(GeneQueryHandler(DynamoDbDatabase())))


@pytest.fixture(scope="session")
def test_classifier():
    """Create test fixture for classifier"""
    return Classify()


@pytest.fixture(scope="session")
def test_gene_normalizer():
    """Create test fixture for gene normalizer"""
    return GeneQueryHandler(DynamoDbDatabase())


@pytest.fixture(scope="session")
def test_cool_seq_tool():
    """Create test fixture for cool seq tool"""
    return CoolSeqTool()


@pytest.fixture(scope="session")
def val_params(test_cool_seq_tool, test_gene_normalizer):
    """Create test fixture for validator params"""
    return [
        test_cool_seq_tool.seqrepo_access,
        test_cool_seq_tool.transcript_mappings,
        test_cool_seq_tool.uta_db,
        test_gene_normalizer,
        test_cool_seq_tool.liftover,
    ]


@pytest.fixture(scope="session")
def test_query_handler():
    """Build normalize test fixture."""
    return QueryHandler()


@pytest.fixture(scope="session")
def test_cnv_handler(test_query_handler):
    """Create test fixture for copy number variation handler"""
    return test_query_handler.to_copy_number_handler


@pytest.fixture(scope="session")
def braf_ncbi_seq_loc():
    """Create test fixture for BRAF ncbi priority sequence location"""
    digest = "0nPwKHYNnTmJ06G-gSmz8BEhB_NTp-0B"
    return {
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
        },
        "start": 140713327,
        "end": 140924929,
        "id": f"ga4gh:SL.{digest}",
        "digest": digest,
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="session")
def prpf8_ncbi_seq_loc():
    """Create test fixture for PRPF8 ncbi priority sequence location"""
    return {
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.dLZ15tNO1Ur0IcGjwc3Sdi_0A6Yf4zm7",
        },
        "start": 1650628,
        "end": 1684867,
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="session")
def braf_600loc():
    """Create test fixture for BRAF 600 location"""
    digest = "t-3DrWALhgLdXHsupI-e-M00aL3HgK3y"
    return {
        "id": f"ga4gh:SL.{digest}",
        "digest": digest,
        "end": 600,
        "start": 599,
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.cQvw4UsHHRRlogxbWCB8W-mKD4AraM9y",
        },
        "sequence": "V",
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="session")
def braf_v600e(braf_600loc):
    """Create BRAF V600E protein test fixture."""
    digest = "j4XnsLZcdzDIYa5pvvXM7t1wn9OITr0L"
    params = {
        "id": f"ga4gh:VA.{digest}",
        "digest": digest,
        "location": braf_600loc,
        "state": {"sequence": "E", "type": "LiteralSequenceExpression"},
        "type": "Allele",
    }
    return models.Allele(**params)


@pytest.fixture(scope="session")
def vhl_reference_agree():
    """Create NP_000542.1:p.Pro61 fixture."""
    params = {
        "location": {
            "end": 61,
            "start": 60,
            "sequenceReference": {
                "type": "SequenceReference",
                "refgetAccession": "SQ.z-Oa0pZkJ6GHJHOYM7h5mY_umc0SJzTu",
            },
            "sequence": "P",
            "type": "SequenceLocation",
        },
        "state": {"sequence": "P", "type": "LiteralSequenceExpression"},
        "type": "Allele",
    }
    return models.Allele(**params)


@pytest.fixture(scope="session")
def protein_insertion():
    """Create test fixture for NP protein insertion (CA645561585)."""
    params = {
        "location": {
            "end": 770,
            "start": 770,
            "sequenceReference": {
                "type": "SequenceReference",
                "refgetAccession": "SQ.vyo55F6mA6n2LgN4cagcdRzOuh38V4mE",
            },
            "type": "SequenceLocation",
        },
        "state": {"sequence": "GL", "type": "LiteralSequenceExpression"},
        "type": "Allele",
    }
    return models.Allele(**params)


@pytest.fixture(scope="session")
def protein_deletion_np_range():
    """Create test fixture for protein deletion using NP accession and
    range for deletion.
    """
    params = {
        "location": {
            "end": 759,
            "start": 754,
            "sequenceReference": {
                "type": "SequenceReference",
                "refgetAccession": "SQ.AF1UFydIo02-bMplonKSfxlWY2q6ze3m",
            },
            "sequence": "LRENT",
            "type": "SequenceLocation",
        },
        "state": {
            "length": 0,
            "repeatSubunitLength": 5,
            "sequence": "",
            "type": "ReferenceLengthExpression",
        },
        "type": "Allele",
    }
    return models.Allele(**params)


@pytest.fixture(scope="session")
def braf_v600e_genomic_sub():
    """Create test fixture for NC_000007.14:g.140753336A>T"""
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
        "state": {"sequence": "T", "type": "LiteralSequenceExpression"},
        "type": "Allele",
    }
    return models.Allele(**params)


@pytest.fixture(scope="session")
def genomic_dup1_seq_loc_normalized():
    """Create test fixture containing genomic dup1 sequence location normalized"""
    digest = "XyxdODl0lLloyj3EQgIzIOEukl2WiaHw"
    return {
        "id": f"ga4gh:SL.{digest}",
        "digest": digest,
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.Zu7h9AggXxhTaGVsy7h_EZSChSZGcmgX",
        },
        "sequence": "GG",
        "start": 49531260,
        "end": 49531262,
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="session")
def genomic_dup1_seq_loc_not_normalized():
    """Create test fixture containing genomic dup1 sequence location that was
    normalized
    """
    digest = "2vbgFGHGB0QGODwgZNi05fWbROkkjf04"
    return {
        "id": f"ga4gh:SL.{digest}",
        "digest": digest,
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.Zu7h9AggXxhTaGVsy7h_EZSChSZGcmgX",
        },
        "start": 49531261,
        "end": 49531262,
        "sequence": "G",
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="session")
def genomic_dup1_38_cn(genomic_dup1_seq_loc_not_normalized):
    """Create test fixture for copy number count dup1 on GRCh38"""
    digest = "gF1l6Zh6aY3vy_TR7rrat6FTmwiwIukY"
    params = {
        "type": "CopyNumberCount",
        "id": f"ga4gh:CN.{digest}",
        "digest": digest,
        "location": genomic_dup1_seq_loc_not_normalized,
        "copies": 3,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="session")
def genomic_dup2_seq_loc_normalized():
    """Create genomic dup2 sequence location"""
    return {
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
        },
        "start": 33211289,
        "end": 33211293,
        "sequence": "TCTA",
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="session")
def genomic_dup2_38_cn(genomic_dup2_seq_loc_normalized):
    """Create test fixture for copy number count dup2 on GRCh38"""
    params = {
        "type": "CopyNumberCount",
        "location": genomic_dup2_seq_loc_normalized,
        "copies": 3,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="session")
def genomic_del3_dup3_loc_not_normalized():
    """Create genomic del3 dup3 sequence location"""
    return {
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
        },
        "start": [31060226, 31100350],
        "end": [33274278, 33417151],
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="session")
def genomic_dup4_loc():
    """Create genomic dup4 sequence location"""
    return {
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.-A1QmD_MatoqxvgVxBLZTONHz9-c7nQo",
        },
        "start": [None, 30417575],
        "end": [31394018, None],
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="session")
def genomic_dup5_loc():
    """Create genomic dup5 sequence location"""
    return {
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
        },
        "start": [None, 154021811],
        "end": 154092209,
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="session")
def genomic_dup6_loc():
    """Create genomic dup6 sequence location"""
    return {
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
        },
        "start": 154021811,
        "end": [154092209, None],
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="session")
def genomic_del1_seq_loc():
    """Create genomic del1 sequence location"""
    return {
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.Zu7h9AggXxhTaGVsy7h_EZSChSZGcmgX",
        },
        "start": 10149810,
        "end": 10149811,
        "sequence": "T",
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="session")
def genomic_del1_lse(genomic_del1_seq_loc):
    """Create a test fixture for genomic del LSE."""
    params = {
        "type": "Allele",
        "location": genomic_del1_seq_loc,
        "state": {
            "length": 0,
            "repeatSubunitLength": 1,
            "type": "ReferenceLengthExpression",
            "sequence": "",
        },
    }
    return models.Allele(**params)


@pytest.fixture(scope="session")
def genomic_del1_38_cn(genomic_del1_seq_loc):
    """Create test fixture for copy number count del1 on GRCh38"""
    params = {
        "type": "CopyNumberCount",
        "location": genomic_del1_seq_loc,
        "copies": 1,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="session")
def genomic_del2_seq_loc():
    """Create genomic del2 sequence location"""
    return {
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.Zu7h9AggXxhTaGVsy7h_EZSChSZGcmgX",
        },
        "sequence": "ATGTTGACGGACAGCCTAT",
        "start": 10146594,
        "end": 10146613,
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="session")
def genomic_del2_lse(genomic_del2_seq_loc):
    """Create a test fixture for genomic del LSE."""
    params = {
        "type": "Allele",
        "location": genomic_del2_seq_loc,
        "state": {
            "type": "ReferenceLengthExpression",
            "sequence": "",
            "length": 0,
            "repeatSubunitLength": 19,
        },
    }
    return models.Allele(**params)


@pytest.fixture(scope="session")
def genomic_del2_38_cn(genomic_del2_seq_loc):
    """Create test fixture for copy number count del1 on GRCh38"""
    params = {
        "type": "CopyNumberCount",
        "location": genomic_del2_seq_loc,
        "copies": 1,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="session")
def genomic_del4_seq_loc():
    """Create genomic del4 sequence location"""
    return {
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
        },
        "start": [None, 31120495],
        "end": [33339477, None],
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="session")
def genomic_del5_seq_loc():
    """Create genomic del5 sequence location"""
    return {
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
        },
        "start": [None, 18575353],
        "end": 18653629,
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="session")
def genomic_del6_seq_loc():
    """Create genomic del6 sequence location"""
    return {
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.0iKlIQk2oZLoeOG9P1riRU6hvL5Ux8TV",
        },
        "start": 133462763,
        "end": [133464858, None],
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="session")
def grch38_genomic_insertion_seq_loc():
    """Create test fixture for GRCh38 genomic insertion seq location"""
    return {
        "end": 39724743,
        "start": 39724731,
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.dLZ15tNO1Ur0IcGjwc3Sdi_0A6Yf4zm7",
        },
        "sequence": "TACGTGATGGCT",
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="session")
def grch38_genomic_insertion_variation(grch38_genomic_insertion_seq_loc):
    """Create a test fixture for NC_000017.10:g.37880993_37880994insGCTTACGTGATG"""
    params = {
        "location": grch38_genomic_insertion_seq_loc,
        "state": {
            "length": 24,
            "repeatSubunitLength": 12,
            "sequence": "TACGTGATGGCTTACGTGATGGCT",
            "type": "ReferenceLengthExpression",
        },
        "type": "Allele",
    }
    return models.Allele(**params)


@pytest.fixture(scope="session")
def braf_amplification(braf_ncbi_seq_loc):
    """Create test fixture for BRAF Amplification"""
    digest = "uPQaLz6KSwXWdsjNUZ5kRn3znBZF5YwV"
    params = {
        "id": f"ga4gh:CX.{digest}",
        "digest": digest,
        "location": braf_ncbi_seq_loc,
        "copyChange": {"primaryCode": "EFO:0030072"},
        "type": "CopyNumberChange",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="session")
def prpf8_amplification(prpf8_ncbi_seq_loc):
    """Create test fixture for PRPF8 Amplification"""
    params = {
        "location": prpf8_ncbi_seq_loc,
        "copyChange": {"primaryCode": "EFO:0030072"},
        "type": "CopyNumberChange",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_del3_dup3_cn_38(genomic_del3_dup3_loc_not_normalized):
    """Create test fixture copy number variation for del/dup 3 on GRCh38"""
    params = {
        "type": "CopyNumberCount",
        "location": genomic_del3_dup3_loc_not_normalized,
        "copies": 2,
    }
    return models.CopyNumberCount(**params)


def _delete_id_and_digest(vrs_obj_dict):
    """Delete id and digest properties from VRS object"""
    with contextlib.suppress(KeyError):
        # Some fixtures have IDs for other tests
        del vrs_obj_dict["id"]

    with contextlib.suppress(KeyError):
        # Some fixtures have digests for other tests
        del vrs_obj_dict["digest"]


def _vrs_id_and_digest_existence_checks(vrs_obj_dict, prefix=None):
    """Check that VRS id and digest exists.

    Does not check actual values.
    `vrs_obj_dict` will be mutated (id and digest fields removed).
    """
    variation_vrs_digest = vrs_obj_dict.pop("digest")
    variation_vrs_id = vrs_obj_dict.pop("id")

    if not prefix:
        prefix = ("ga4gh:VA.", "ga4gh:CX", "ga4gh:CN.")

    assert variation_vrs_id.startswith(prefix)
    assert variation_vrs_id.endswith(variation_vrs_digest)

    location_vrs_digest = vrs_obj_dict["location"].pop("digest")
    location_vrs_id = vrs_obj_dict["location"].pop("id")
    assert location_vrs_id == f"ga4gh:SL.{location_vrs_digest}"


def _mane_gene_ext_checks(actual_vo: dict) -> None:
    """Check mane gene extensions existence

    :param actual_vo: Actual VRS object represented as a dictionary
    """
    extensions = actual_vo.pop("extensions")
    assert len(extensions) == 1

    mane_genes_ext = extensions[0]
    assert mane_genes_ext["name"] == "mane_genes"
    for mane_gene in mane_genes_ext["value"]:
        assert mane_gene["ncbi_gene_id"]
        assert mane_gene["symbol"]


def assertion_checks(
    normalize_response, test_variation, mane_genes_exts=False, check_vrs_id=False
):
    """Check that normalize_response and test_variation are equal."""
    actual = normalize_response.variation.model_dump(exclude_none=True)
    if not check_vrs_id:
        _vrs_id_and_digest_existence_checks(actual)

    # Check MANE genes existence
    if mane_genes_exts:
        _mane_gene_ext_checks(actual)

    expected = test_variation.model_copy().model_dump(exclude_none=True)
    if not check_vrs_id:
        _delete_id_and_digest(expected)
        _delete_id_and_digest(expected["location"])

    assert actual == expected, "variation"


def cnv_assertion_checks(resp, test_fixture, check_vrs_id=False, mane_genes_exts=False):
    """Check that actual response for to copy number matches expected"""

    def _update_expected_mappings(expected_):
        """Modify test fixture copy to include mappable concept object for CX var"""
        expected_["copyChange"]["mappings"] = [
            {
                "relation": "exactMatch",
                "coding": {
                    "system": "https://www.ebi.ac.uk/efo/",
                    "code": expected_["copyChange"]["primaryCode"],
                },
            }
        ]

    expected = test_fixture.model_copy(deep=True).model_dump(exclude_none=True)

    if isinstance(resp, NormalizeService):
        actual = resp.variation
        if isinstance(actual, models.CopyNumberChange):
            _update_expected_mappings(expected)
            prefix = "ga4gh:CX."
        elif isinstance(actual, models.CopyNumberCount):
            prefix = "ga4gh:CN."

        actual = actual.model_dump(exclude_none=True)
    else:
        try:
            cnc = resp.copy_number_count
        except AttributeError:
            _update_expected_mappings(expected)
            actual = resp.copy_number_change.model_dump(exclude_none=True)
            prefix = "ga4gh:CX."
        else:
            actual = cnc.model_dump(exclude_none=True)
            prefix = "ga4gh:CN."

    # Check MANE genes existence
    if mane_genes_exts:
        _mane_gene_ext_checks(actual)

    if not check_vrs_id:
        _vrs_id_and_digest_existence_checks(actual, prefix=prefix)

    if not check_vrs_id:
        _delete_id_and_digest(expected)
        _delete_id_and_digest(expected["location"])

    assert actual == expected
    assert resp.warnings == []
