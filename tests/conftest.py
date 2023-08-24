"""Create methods used throughout tests."""
import asyncio

import pytest
from cool_seq_tool import CoolSeqTool
from ga4gh.vrs import models
from gene.database.dynamodb import DynamoDbDatabase
from gene.query import QueryHandler as GeneQueryHandler

from variation.classify import Classify
from variation.query import QueryHandler
from variation.tokenize import Tokenize
from variation.tokenizers import GeneSymbol


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
    return {
        "sequence": "ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
        "start": 140713327,
        "end": 140924929,
        "id": "ga4gh:SL.I6Hn1A9YViUPq37PgWuSNnL-BJmU6XgF",
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="session")
def prpf8_ncbi_seq_loc():
    """Create test fixture for PRPF8 ncbi priority sequence location"""
    return {
        "sequence": "ga4gh:SQ.dLZ15tNO1Ur0IcGjwc3Sdi_0A6Yf4zm7",
        "start": 1650628,
        "end": 1684867,
        "id": "ga4gh:SL.dMg5mAjhrEeUbAD7jw_7H-oK-Ry58CWW",
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="session")
def braf_600loc():
    """Create test fixture for BRAF 600 location"""
    return {
        "id": "ga4gh:SL.ko4RJfU-2fvZrbCDpo6i-Ljcfi59TcQI",
        "end": 600,
        "start": 599,
        "sequence": "ga4gh:SQ.cQvw4UsHHRRlogxbWCB8W-mKD4AraM9y",
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="session")
def braf_v600e(braf_600loc):
    """Create BRAF V600E protein test fixture."""
    params = {
        "id": "ga4gh:VA.JENfSejJLoraR6JpXzhBzI1iB3aGMjo3",
        "location": braf_600loc,
        "state": {"sequence": "E", "type": "LiteralSequenceExpression"},
        "type": "Allele",
    }
    return models.Allele(**params)


@pytest.fixture(scope="session")
def vhl_reference_agree():
    """Create NP_000542.1:p.Pro61 fixture."""
    params = {
        "id": "ga4gh:VA.G9YTYX6n-6lX1S7is-CyDZIs5J___DoM",
        "location": {
            "id": "ga4gh:SL.88Oa6i2bJBTaElAQNEjLqMQ8HR0hgdKx",
            "end": 61,
            "start": 60,
            "sequence": "ga4gh:SQ.z-Oa0pZkJ6GHJHOYM7h5mY_umc0SJzTu",
            "type": "SequenceLocation",
        },
        "state": {"sequence": "P", "type": "LiteralSequenceExpression"},
        "type": "Allele",
    }
    return models.Allele(**params)


@pytest.fixture(scope="session")
def protein_insertion():
    """Create test fixture for NP protein insertion."""
    params = {
        "id": "ga4gh:VA.SA-ROrMmAkVSSAmT3xuKXkvED8hgppFY",
        "location": {
            "id": "ga4gh:SL.vUZ6MNJsYdKlL-_Fl3Cm7AnubDRItWv8",
            "end": 770,
            "start": 770,
            "sequence": "ga4gh:SQ.vyo55F6mA6n2LgN4cagcdRzOuh38V4mE",
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
        "id": "ga4gh:VA.uV8vtTCRJGA__U16rXYJX8QmIDlyXwcw",
        "location": {
            "id": "ga4gh:SL.xfIpKI77izLOuEFvJewyQ-VncDr6q059",
            "end": 759,
            "start": 754,
            "sequence": "ga4gh:SQ.AF1UFydIo02-bMplonKSfxlWY2q6ze3m",
            "type": "SequenceLocation",
        },
        "state": {"sequence": "", "type": "LiteralSequenceExpression"},
        "type": "Allele",
    }
    return models.Allele(**params)


@pytest.fixture(scope="session")
def braf_v600e_genomic_sub():
    """Create test fixture for NC_000007.14:g.140753336A>T"""
    params = {
        "id": "ga4gh:VA.oM4qtwQ_V49243VNkCTZYwyiXCWC8gIc",
        "location": {
            "id": "ga4gh:SL.j9b9mPXETlWMYDjZ62S9DIXATtt-vY_e",
            "end": 140753336,
            "start": 140753335,
            "sequence": "ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
            "type": "SequenceLocation",
        },
        "state": {"sequence": "T", "type": "LiteralSequenceExpression"},
        "type": "Allele",
    }
    return models.Allele(**params)


@pytest.fixture(scope="session")
def genomic_dup1_seq_loc_normalized():
    """Create test fixture containing genomic dup1 sequence location normalized"""
    return {
        "id": "ga4gh:SL.LC-4oEU03VAX1bJD4pQii8tEYl9UbiBw",
        "sequence": "ga4gh:SQ.Zu7h9AggXxhTaGVsy7h_EZSChSZGcmgX",
        "start": 49531260,
        "end": 49531262,
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="session")
def genomic_dup1_seq_loc_not_normalized():
    """Create test fixture containing genomic dup1 sequence location that was
    normalized
    """
    return {
        "id": "ga4gh:SL.0NdhAbmG-jjXLdGYeKLarI0lHXRCjro-",
        "sequence": "ga4gh:SQ.Zu7h9AggXxhTaGVsy7h_EZSChSZGcmgX",
        "start": 49531261,
        "end": 49531262,
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="session")
def genomic_dup1_38_cn(genomic_dup1_seq_loc_not_normalized):
    """Create test fixture for copy number count dup1 on GRCh38"""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.xh7DIa3NT8rE10PtEyMIcdwJcuM_rhK7",
        "subject": genomic_dup1_seq_loc_not_normalized,
        "copies": 3,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="session")
def genomic_dup2_seq_loc_normalized():
    """Create genomic dup2 sequence location"""
    return {
        "id": "ga4gh:SL.j0IhhSzTdpGAcbpugs-IW2ZwGLCYoYuY",
        "sequence": "ga4gh:SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
        "start": 33211289,
        "end": 33211293,
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="session")
def genomic_dup2_38_cn(genomic_dup2_seq_loc_normalized):
    """Create test fixture for copy number count dup2 on GRCh38"""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.ST8ehnnfAKFM4kXpCBwmNmeBmv36LG9K",
        "subject": genomic_dup2_seq_loc_normalized,
        "copies": 3,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="session")
def genomic_del3_dup3_loc_not_normalized():
    """Create genomic del3 dup3 sequence location"""
    return {
        "id": "ga4gh:SL.rbCgFyEIFtvtHkVoeV5_L-TVhkHJhzFM",
        "sequence": "ga4gh:SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
        "start": [31060226, 31100350],
        "end": [33274279, 33417152],
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="session")
def genomic_dup4_loc():
    """Create genomic dup4 sequence location"""
    return {
        "id": "ga4gh:SL.gWgZIFB8G8OdYeD-FmtCyuX-5zMmDOpc",
        "sequence": "ga4gh:SQ.-A1QmD_MatoqxvgVxBLZTONHz9-c7nQo",
        "start": [None, 30417575],
        "end": [31394018, None],
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="session")
def genomic_dup5_loc():
    """Create genomic dup5 sequence location"""
    return {
        "id": "ga4gh:SL.2b_skZgj2CdVjgzHQfZt2fhx7prwDmnW",
        "sequence": "ga4gh:SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
        "start": [None, 154021811],
        "end": 154092209,
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="session")
def genomic_dup6_loc():
    """Create genomic dup6 sequence location"""
    return {
        "id": "ga4gh:SL.maCy91pBvLaSCNSS7QqfB01sHOhU8tUX",
        "sequence": "ga4gh:SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
        "start": 154021811,
        "end": [154092209, None],
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="session")
def genomic_del1_seq_loc():
    """Create genomic del1 sequence location"""
    return {
        "id": "ga4gh:SL.wfjWzTQ-yFTd3NrEldcre1AXjbwUAQNe",
        "sequence": "ga4gh:SQ.Zu7h9AggXxhTaGVsy7h_EZSChSZGcmgX",
        "start": 10149810,
        "end": 10149811,
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="session")
def genomic_del1_lse(genomic_del1_seq_loc):
    """Create a test fixture for genomic del LSE."""
    params = {
        "type": "Allele",
        "id": "ga4gh:VA.Kt_0JCbABYTiBhQAYLz0RlxTw-HUHPas",
        "location": genomic_del1_seq_loc,
        "state": {"type": "LiteralSequenceExpression", "sequence": ""},
    }
    return models.Allele(**params)


@pytest.fixture(scope="session")
def genomic_del1_38_cn(genomic_del1_seq_loc):
    """Create test fixture for copy number count del1 on GRCh38"""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.8H0qnVpREcE8hhsxnzdwgIjmkhQ3Y7cK",
        "subject": genomic_del1_seq_loc,
        "copies": 1,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="session")
def genomic_del2_seq_loc():
    """Create genomic del2 sequence location"""
    return {
        "id": "ga4gh:SL.eDHB3W9HlloNvpfz--_7vRYo_aBJW5P1",
        "sequence": "ga4gh:SQ.Zu7h9AggXxhTaGVsy7h_EZSChSZGcmgX",
        "start": 10146594,
        "end": 10146613,
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="session")
def genomic_del2_lse(genomic_del2_seq_loc):
    """Create a test fixture for genomic del LSE."""
    params = {
        "type": "Allele",
        "id": "ga4gh:VA.OhVKf_ePMH0ptie44HofAV_H1rCr-7wC",
        "location": genomic_del2_seq_loc,
        "state": {"type": "LiteralSequenceExpression", "sequence": ""},
    }
    return models.Allele(**params)


@pytest.fixture(scope="session")
def genomic_del2_38_cn(genomic_del2_seq_loc):
    """Create test fixture for copy number count del1 on GRCh38"""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.CYPdGzYwlpkh0za_VwINJk6qPJpIRemA",
        "subject": genomic_del2_seq_loc,
        "copies": 1,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="session")
def genomic_del4_seq_loc():
    """Create genomic del4 sequence location"""
    return {
        "id": "ga4gh:SL.4Pq7vqlY1X8gHO5Znl_TsqKXHReiwpGN",
        "sequence": "ga4gh:SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
        "start": [None, 31120495],
        "end": [33339477, None],
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="session")
def genomic_del5_seq_loc():
    """Create genomic del5 sequence location"""
    return {
        "id": "ga4gh:SL.bLI8ap8wB9ZZfpUUUCx-0qiBW-gEDT3E",
        "sequence": "ga4gh:SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
        "start": [None, 18575353],
        "end": 18653629,
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="session")
def genomic_del6_seq_loc():
    """Create genomic del6 sequence location"""
    return {
        "id": "ga4gh:SL.aBqKQ7nMM-kvniKdh2yZ5XIrNAnWfJ5x",
        "sequence": "ga4gh:SQ.0iKlIQk2oZLoeOG9P1riRU6hvL5Ux8TV",
        "start": 133462763,
        "end": [133464858, None],
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="session")
def grch38_genomic_insertion_seq_loc():
    """Create test fixture for GRCh38 genomic insertion seq location"""
    return {
        "id": "ga4gh:SL.0fSp38w0ROjgtrTiU9UFd_wroFCznvhl",
        "end": 39724743,
        "start": 39724731,
        "sequence": "ga4gh:SQ.dLZ15tNO1Ur0IcGjwc3Sdi_0A6Yf4zm7",
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="session")
def grch38_genomic_insertion_variation(grch38_genomic_insertion_seq_loc):
    """Create a test fixture for NC_000017.10:g.37880993_37880994insGCTTACGTGATG"""
    params = {
        "id": "ga4gh:VA.stx0NrRIPFsy8hBXsx2lVyA6Ewhz7ytI",
        "location": grch38_genomic_insertion_seq_loc,
        "state": {
            "sequence": "TACGTGATGGCTTACGTGATGGCT",
            "type": "LiteralSequenceExpression",
        },
        "type": "Allele",
    }
    return models.Allele(**params)


@pytest.fixture(scope="session")
def braf_amplification(braf_ncbi_seq_loc):
    """Create test fixture for BRAF Amplification"""
    params = {
        "id": "ga4gh:CX.qhzBsMy5O_xODtSBLSbfut6NHa-cQtpU",
        "subject": braf_ncbi_seq_loc,
        "copyChange": "efo:0030072",
        "type": "CopyNumberChange",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="session")
def prpf8_amplification(prpf8_ncbi_seq_loc):
    """Create test fixture for PRPF8 Amplification"""
    params = {
        "id": "ga4gh:CX.j4knp4ubcVMrvTOZ1Z1PGCkqbo7ZvFAz",
        "subject": prpf8_ncbi_seq_loc,
        "copyChange": "efo:0030072",
        "type": "CopyNumberChange",
    }
    return models.CopyNumberChange(**params)


def assertion_checks(normalize_response, test_variation):
    """Check that normalize_response and test_variation are equal."""
    actual = normalize_response.variation.model_dump(exclude_none=True)
    expected = test_variation.model_dump(exclude_none=True)
    assert actual == expected, "variation"


def cnv_assertion_checks(resp, test_fixture):
    """Check that actual response for to copy number matches expected"""
    try:
        getattr(resp, "copy_number_count")
    except AttributeError:
        actual = resp.copy_number_change.model_dump(exclude_none=True)
    else:
        actual = resp.copy_number_count.model_dump(exclude_none=True)
    expected = test_fixture.model_dump(exclude_none=True)
    assert actual == expected
    assert resp.warnings == []
