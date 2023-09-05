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
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
        },
        "start": 140713327,
        "end": 140924929,
        "id": "ga4gh:SL.rXzVqqlchBvUef98MNQA77FvwSJgiOf5",
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
        "id": "ga4gh:SL._ygg7PVx6ENto_L_BnopQBlcCaWMWQdB",
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="session")
def braf_600loc():
    """Create test fixture for BRAF 600 location"""
    return {
        "id": "ga4gh:SL.pqBCJoUGtDNLKib8F5ke3cyqqrkZzfud",
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
        "id": "ga4gh:VA.tJ7kARp1YlxgcUesm7DXSj_SDXjXg3-u",
        "location": braf_600loc,
        "state": {"sequence": "E", "type": "LiteralSequenceExpression"},
        "type": "Allele",
    }
    return models.Allele(**params)


@pytest.fixture(scope="session")
def vhl_reference_agree():
    """Create NP_000542.1:p.Pro61 fixture."""
    params = {
        "id": "ga4gh:VA.u107jl5_bx268nr2E6VJjmRflj4wQMOG",
        "location": {
            "id": "ga4gh:SL.OiKH3aqBpAjG2vMytW1tNwslBdBkThMr",
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
    """Create test fixture for NP protein insertion."""
    params = {
        "id": "ga4gh:VA.1izVhuyd_6XaNxPIPMwwIB0KeqYrCTUD",
        "location": {
            "id": "ga4gh:SL.qOc3iBPVW7NK7p9yvW1Q6jlZkOTI9Avl",
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
        "id": "ga4gh:VA.poElW1WkhGhPlhSWbvqKHvfQIFc-lhGd",
        "location": {
            "id": "ga4gh:SL.r7VKhO-lvK9ylJ1vDtC9L3gZ2Z9Lyb2f",
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
        "id": "ga4gh:VA.hHoHQ5GTC0b4PlR_xX-TBBxsLECyMfaf",
        "location": {
            "id": "ga4gh:SL.q4OIseiFvw2R0noLIpOQ743VsqyDumZ4",
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
        "id": "ga4gh:SL.ymeWuZvaNt5QpFeq3vC15S_3syWfQVjK",
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
        "id": "ga4gh:SL.L7dW0LicGFRwVQzj30aYac_vdy36Idvk",
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
        "id": "ga4gh:CN.ZelwbEDATHfNcECE25OspdQw9zRW0c4G",
        "subject": genomic_dup1_seq_loc_not_normalized,
        "copies": 3,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="session")
def genomic_dup2_seq_loc_normalized():
    """Create genomic dup2 sequence location"""
    return {
        "id": "ga4gh:SL.RPFDrsfaF3NZDGL-5ExeaDxera4Q8jg6",
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
        "id": "ga4gh:CN.sD5oAN3PXrIWMNtGfTcL3AYAvkFz9e9z",
        "subject": genomic_dup2_seq_loc_normalized,
        "copies": 3,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="session")
def genomic_del3_dup3_loc_not_normalized():
    """Create genomic del3 dup3 sequence location"""
    return {
        "id": "ga4gh:SL.gABz8yV6pFIY5kMY6O8BqvdNcfXt8ulQ",
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
        },
        "start": [31060226, 31100350],
        "end": [33274279, 33417152],
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="session")
def genomic_dup4_loc():
    """Create genomic dup4 sequence location"""
    return {
        "id": "ga4gh:SL.3VRKy0qXlEvL2PZSqLBHchE5TIeJYHhx",
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
        "id": "ga4gh:SL.bF_LzXeTaYlX9EafOpu0l27rkoNZpgU1",
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
        "id": "ga4gh:SL.VW7QImr8yUPvNsphyDKnt5t0ubMRBMXL",
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
        "id": "ga4gh:SL.zSDVB3Qlbw49ezgC0rXrbcTZKFdNZkb6",
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
        "id": "ga4gh:VA.5nh4ErMaI1FPhhoqFD60G8aa4pFEk6Ua",
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
        "id": "ga4gh:CN.ow4omS2-7jShTHhTHTrR3ujTYRtbqmwg",
        "subject": genomic_del1_seq_loc,
        "copies": 1,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="session")
def genomic_del2_seq_loc():
    """Create genomic del2 sequence location"""
    return {
        "id": "ga4gh:SL.ni7iQg4vZPp4RWkpbHjnUMF9QguEC0xp",
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
        "id": "ga4gh:VA.9P3LU5sMg5oQgIop5k0aeiV_og2AtfDW",
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
        "id": "ga4gh:CN.VGwBpEiItF-2J1dW3E8jatsINn9batgA",
        "subject": genomic_del2_seq_loc,
        "copies": 1,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="session")
def genomic_del4_seq_loc():
    """Create genomic del4 sequence location"""
    return {
        "id": "ga4gh:SL.kUqkUmhq2dK4uIjUkKwdCjs2ku_5Ey4K",
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
        "id": "ga4gh:SL.K4NgbWh598IWUZHphvZSZBiRfNY2qYgE",
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
        "id": "ga4gh:SL.i_9j5WJdwQDdPgFOMOL6ikpNlLuBv48s",
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
        "id": "ga4gh:SL.QU6-AyDQcobwTUN7qODvSzhxqIbEsAN0",
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
        "id": "ga4gh:VA.v7qk_YE_i7-H-8w1FZ-cSMpI57luGb7G",
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
        "id": "ga4gh:CX.o8ffM6gQYxFAEY4X0u2a65FqaYCapXmM",
        "subject": braf_ncbi_seq_loc,
        "copyChange": "efo:0030072",
        "type": "CopyNumberChange",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="session")
def prpf8_amplification(prpf8_ncbi_seq_loc):
    """Create test fixture for PRPF8 Amplification"""
    params = {
        "id": "ga4gh:CX.olwM3A0aY-mkYRFhSl0P07_Gtclsau4m",
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
