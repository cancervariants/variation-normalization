"""Create methods used throughout tests."""
import asyncio

import pytest
from cool_seq_tool.app import CoolSeqTool
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
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
        },
        "start": 140713327,
        "end": 140924929,
        "id": "ga4gh:SL.uNBZoxhjhohl24VlIut-JxPJAGfJ7EQE",
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
        "id": "ga4gh:SL.1i49iv3wcBq7SaOA14cs1Kz7SR6DkCw1",
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="session")
def braf_600loc():
    """Create test fixture for BRAF 600 location"""
    return {
        "id": "ga4gh:SL.ZA1XNKhCT_7m2UtmnYb8ZYOVS4eplMEK",
        "end": 600,
        "start": 599,
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.cQvw4UsHHRRlogxbWCB8W-mKD4AraM9y",
        },
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="session")
def braf_v600e(braf_600loc):
    """Create BRAF V600E protein test fixture."""
    params = {
        "id": "ga4gh:VA.4XBXAxSAk-WyAu5H0S1-plrk_SCTW1PO",
        "location": braf_600loc,
        "state": {"sequence": "E", "type": "LiteralSequenceExpression"},
        "type": "Allele",
    }
    return models.Allele(**params)


@pytest.fixture(scope="session")
def vhl_reference_agree():
    """Create NP_000542.1:p.Pro61 fixture."""
    params = {
        "id": "ga4gh:VA.RMmwTvhrPVwfMZ6knsf5zMWQn_F1ukYh",
        "location": {
            "id": "ga4gh:SL.8TZYB8Oqqn93q07zrsNhvRW1JjNpaQXc",
            "end": 61,
            "start": 60,
            "sequenceReference": {
                "type": "SequenceReference",
                "refgetAccession": "SQ.z-Oa0pZkJ6GHJHOYM7h5mY_umc0SJzTu",
            },
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
        "id": "ga4gh:VA.AOCCh_BU5wKkdgoDNqkORF_x4GQwWh1T",
        "location": {
            "id": "ga4gh:SL.ciWb1ylkqUxiviU1djijiuYVZcgsnQnV",
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
        "id": "ga4gh:VA.3Rk_RElDfX820edkQOHsTTYRogr0EMEY",
        "location": {
            "id": "ga4gh:SL.kOTzy0aLlw0yqnmf29Zk8wh65zHQwere",
            "end": 759,
            "start": 754,
            "sequenceReference": {
                "type": "SequenceReference",
                "refgetAccession": "SQ.AF1UFydIo02-bMplonKSfxlWY2q6ze3m",
            },
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
        "id": "ga4gh:VA.LX3ooHBAiZdKY4RfTXcliUmkj48mnD_M",
        "location": {
            "id": "ga4gh:SL.XutGzMvqbzN-vnxmPt2MJf7ehxmB0opi",
            "end": 140753336,
            "start": 140753335,
            "sequenceReference": {
                "type": "SequenceReference",
                "refgetAccession": "SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
            },
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
        "id": "ga4gh:SL.f0nAiaxOC3rPToQEYRRhbVBNO6HKutyc",
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.Zu7h9AggXxhTaGVsy7h_EZSChSZGcmgX",
        },
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
        "id": "ga4gh:SL.y4-cVA2VxMCDxb9gV2oFrzC386yrEVqh",
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.Zu7h9AggXxhTaGVsy7h_EZSChSZGcmgX",
        },
        "start": 49531261,
        "end": 49531262,
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="session")
def genomic_dup1_38_cn(genomic_dup1_seq_loc_not_normalized):
    """Create test fixture for copy number count dup1 on GRCh38"""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.07iM14yvZ80N_AiaM7G_V4f1pCkmFYz4",
        "location": genomic_dup1_seq_loc_not_normalized,
        "copies": 3,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="session")
def genomic_dup2_seq_loc_normalized():
    """Create genomic dup2 sequence location"""
    return {
        "id": "ga4gh:SL.rVXa8TXm6WTEw-_Lom6A347Q45SB7CON",
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
        },
        "start": 33211289,
        "end": 33211293,
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="session")
def genomic_dup2_38_cn(genomic_dup2_seq_loc_normalized):
    """Create test fixture for copy number count dup2 on GRCh38"""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.C8WuNCba5AN1RoXK1enXgALlM1Qz6X6i",
        "location": genomic_dup2_seq_loc_normalized,
        "copies": 3,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="session")
def genomic_del3_dup3_loc_not_normalized():
    """Create genomic del3 dup3 sequence location"""
    return {
        "id": "ga4gh:SL.-zCp7JBaKQ0niPDueJkuCgQhRIQ50hKw",
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
        "id": "ga4gh:SL.o8sCaAaW2a2f_HsNBTsHOCnWRvIyru0y",
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
        "id": "ga4gh:SL.O__pyYq_u7R__2NUbI3koxxkeCBL7WXq",
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
        "id": "ga4gh:SL.Ls2wfxI-2V2OdMY5HHttwlSwgbpNf_j2",
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
        "id": "ga4gh:SL.zMba5wGtQWQmdFd70yEqMYszGoRaYX25",
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.Zu7h9AggXxhTaGVsy7h_EZSChSZGcmgX",
        },
        "start": 10149810,
        "end": 10149811,
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="session")
def genomic_del1_lse(genomic_del1_seq_loc):
    """Create a test fixture for genomic del LSE."""
    params = {
        "type": "Allele",
        "id": "ga4gh:VA.gztc0BFS6p5V1_QVnEYIJ6DwzZQeDCd2",
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
        "id": "ga4gh:CN.wRj3ZKNriLtPDVj0VlPaTCQfklj2ocGU",
        "location": genomic_del1_seq_loc,
        "copies": 1,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="session")
def genomic_del2_seq_loc():
    """Create genomic del2 sequence location"""
    return {
        "id": "ga4gh:SL.usVkXRvjfX0cEXLvP87Oi8eJJGyizjQF",
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.Zu7h9AggXxhTaGVsy7h_EZSChSZGcmgX",
        },
        "start": 10146594,
        "end": 10146613,
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="session")
def genomic_del2_lse(genomic_del2_seq_loc):
    """Create a test fixture for genomic del LSE."""
    params = {
        "type": "Allele",
        "id": "ga4gh:VA.9NmH0sRYerurt-CE6WlF9UaxZiujByIE",
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
        "id": "ga4gh:CN.i7HRf9gge1HJKzazgvtinosa0bE3gHJu",
        "location": genomic_del2_seq_loc,
        "copies": 1,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="session")
def genomic_del4_seq_loc():
    """Create genomic del4 sequence location"""
    return {
        "id": "ga4gh:SL.bWbNmdT__ptImBwTAIYdyNfazhwvEtXD",
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
        "id": "ga4gh:SL.WDxMzftZLrwp2eQJrlasKuY4ns99wG0v",
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
        "id": "ga4gh:SL.TKIwU5OzGgOWIpnzAHfkCLB7vrKupKhD",
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
        "id": "ga4gh:SL.oVzSkGhh3QJ0FAgihm-kNr9CJbF_7Ln2",
        "end": 39724743,
        "start": 39724731,
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.dLZ15tNO1Ur0IcGjwc3Sdi_0A6Yf4zm7",
        },
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="session")
def grch38_genomic_insertion_variation(grch38_genomic_insertion_seq_loc):
    """Create a test fixture for NC_000017.10:g.37880993_37880994insGCTTACGTGATG"""
    params = {
        "id": "ga4gh:VA.eorsJMgis9uDdRRPVd3srYofdPaM_xn2",
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
    params = {
        "id": "ga4gh:CX.89PECTeQjhhXnNW9yg24DheWOQMgmKk2",
        "location": braf_ncbi_seq_loc,
        "copyChange": "efo:0030072",
        "type": "CopyNumberChange",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="session")
def prpf8_amplification(prpf8_ncbi_seq_loc):
    """Create test fixture for PRPF8 Amplification"""
    params = {
        "id": "ga4gh:CX.KH_rYvTqg5Hq0ysqbrh8JR20oeLYa7bk",
        "location": prpf8_ncbi_seq_loc,
        "copyChange": "efo:0030072",
        "type": "CopyNumberChange",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_del3_dup3_cn_38(genomic_del3_dup3_loc_not_normalized):
    """Create test fixture copy number variation for del/dup 3 on GRCh38"""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.rPsK0krAHgmXhDZEw4fqymR0iDQa3UCJ",
        "location": genomic_del3_dup3_loc_not_normalized,
        "copies": 2,
    }
    return models.CopyNumberCount(**params)


def assertion_checks(normalize_response, test_variation):
    """Check that normalize_response and test_variation are equal."""
    actual = normalize_response.variation.model_dump(exclude_none=True)
    expected = test_variation.model_dump(exclude_none=True)
    assert actual == expected, "variation"


def cnv_assertion_checks(resp, test_fixture):
    """Check that actual response for to copy number matches expected"""
    try:
        resp.copy_number_count  # noqa: B018
    except AttributeError:
        actual = resp.copy_number_change.model_dump(exclude_none=True)
    else:
        actual = resp.copy_number_count.model_dump(exclude_none=True)
    expected = test_fixture.model_dump(exclude_none=True)
    assert actual == expected
    assert resp.warnings == []
