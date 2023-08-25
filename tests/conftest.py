"""Create methods used throughout tests."""
import asyncio
from copy import deepcopy

import pytest
from cool_seq_tool import CoolSeqTool
from ga4gh.vrsatile.pydantic.vrsatile_models import VariationDescriptor
from gene.database.dynamodb import DynamoDbDatabase
from gene.query import QueryHandler as GeneQueryHandler

from variation.classify import Classify
from variation.query import QueryHandler
from variation.tokenize import Tokenize
from variation.tokenizers import GeneSymbol


@pytest.fixture(scope="session")
def event_loop(request):
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
        "sequence_id": "ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
        "start": {"type": "Number", "value": 140713327},
        "end": {"type": "Number", "value": 140924929},
        "id": "ga4gh:SL.po-AExwyqkstDx3JWYn6plIlxn5eojv4",
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="session")
def prpf8_ncbi_seq_loc():
    """Create test fixture for PRPF8 ncbi priority sequence location"""
    return {
        "sequence_id": "ga4gh:SQ.dLZ15tNO1Ur0IcGjwc3Sdi_0A6Yf4zm7",
        "start": {"type": "Number", "value": 1650628},
        "end": {"type": "Number", "value": 1684867},
        "id": "ga4gh:SL.tiIyPIF9TqW6xrAkH1T0Z_oKIGQf-G3A",
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="session")
def braf_600loc():
    """Create test fixture for BRAF 600 location"""
    return {
        "id": "ga4gh:SL.xfBTztcmMstx8jrrdgPiE_BUoLHLFMMS",
        "end": {"value": 600, "type": "Number"},
        "start": {"value": 599, "type": "Number"},
        "sequence_id": "ga4gh:SQ.cQvw4UsHHRRlogxbWCB8W-mKD4AraM9y",
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="session")
def braf_v600e(braf_600loc):
    """Create BRAF V600E protein test fixture."""
    params = {
        "id": "normalize.variation:BRAF%20V600E",
        "type": "VariationDescriptor",
        "variation": {
            "id": "ga4gh:VA.h313H4CQh6pogbbSJ3H5pI1cPoh9YMm_",
            "location": braf_600loc,
            "state": {"sequence": "E", "type": "LiteralSequenceExpression"},
            "type": "Allele",
        },
        "molecule_context": "protein",
        "vrs_ref_allele_seq": "V",
        "gene_context": "hgnc:1097",
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="session")
def vhl_reference_agree():
    """Create NP_000542.1:p.Pro61 fixture."""
    params = {
        "id": "normalize.variation:NP_000542.1%3Ap.Pro61%3D",
        "type": "VariationDescriptor",
        "variation": {
            "id": "ga4gh:VA.hFvpLOEfU4Qtxd1_pdSx-_XnIJng9Xnb",
            "location": {
                "id": "ga4gh:SL.Kx99ER-oHCNP3RwwatOYn9IN5LRRxiy-",
                "end": {"value": 61, "type": "Number"},
                "start": {"value": 60, "type": "Number"},
                "sequence_id": "ga4gh:SQ.z-Oa0pZkJ6GHJHOYM7h5mY_umc0SJzTu",
                "type": "SequenceLocation",
            },
            "state": {"sequence": "P", "type": "LiteralSequenceExpression"},
            "type": "Allele",
        },
        "molecule_context": "protein",
        "vrs_ref_allele_seq": "P",
        "gene_context": "hgnc:12687",
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="session")
def protein_insertion():
    """Create test fixture for NP protein insertion."""
    params = {
        "id": "normalize.variation:NP_005219.2%3Ap.Asp770_Asn771insGlyLeu",
        "type": "VariationDescriptor",
        "variation": {
            "id": "ga4gh:VA.Daydg17safvIKs_ENTrvKpfoSooKgImP",
            "location": {
                "id": "ga4gh:SL.ozw2OUd_hkRUcAo4zUM_jH40Wlbd_lb0",
                "end": {"value": 770, "type": "Number"},
                "start": {"value": 770, "type": "Number"},
                "sequence_id": "ga4gh:SQ.vyo55F6mA6n2LgN4cagcdRzOuh38V4mE",
                "type": "SequenceLocation",
            },
            "state": {"sequence": "GL", "type": "LiteralSequenceExpression"},
            "type": "Allele",
        },
        "molecule_context": "protein",
        "gene_context": "hgnc:3236",
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="session")
def protein_deletion_np_range():
    """Create test fixture for protein deletion using NP accession and
    range for deletion.
    """
    params = {
        "id": "normalize.variation:NP_004439.2%3Ap.Leu755_Thr759del",
        "type": "VariationDescriptor",
        "variation": {
            "id": "ga4gh:VA.gb8REjJatcpIgb1B3LKiIDI4DhJd70Bk",
            "location": {
                "id": "ga4gh:SL.Ca3e2urICwddGClFXppmmMeGr3Zqg-i8",
                "end": {"value": 759, "type": "Number"},
                "start": {"value": 754, "type": "Number"},
                "sequence_id": "ga4gh:SQ.AF1UFydIo02-bMplonKSfxlWY2q6ze3m",
                "type": "SequenceLocation",
            },
            "state": {"sequence": "", "type": "LiteralSequenceExpression"},
            "type": "Allele",
        },
        "molecule_context": "protein",
        "vrs_ref_allele_seq": "LRENT",
        "gene_context": "hgnc:3430",
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="session")
def braf_v600e_genomic_sub():
    """Create test fixture for NC_000007.14:g.140753336A>T"""
    return {
        "id": "ga4gh:VA.3xFHF399HGbG1JUf5uwcj3oWVKZJ70oX",
        "location": {
            "id": "ga4gh:SL.WBLxdkoypnRME6b8tJtlOWqZKU1ruqY1",
            "end": {"value": 140753336, "type": "Number"},
            "start": {"value": 140753335, "type": "Number"},
            "sequence_id": "ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
            "type": "SequenceLocation",
        },
        "state": {"sequence": "T", "type": "LiteralSequenceExpression"},
        "type": "Allele",
    }


@pytest.fixture(scope="session")
def genomic_dup1_seq_loc_normalized():
    """Create test fixture containing genomic dup1 sequence location normalized"""
    return {
        "id": "ga4gh:SL.KefUQwlqEBGtzoNO-MzOozx7_H1uP-fD",
        "sequence_id": "ga4gh:SQ.Zu7h9AggXxhTaGVsy7h_EZSChSZGcmgX",
        "start": {"value": 49531260, "type": "Number"},
        "end": {"value": 49531262, "type": "Number"},
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="session")
def genomic_dup1_seq_loc_not_normalized():
    """Create test fixture containing genomic dup1 sequence location that was
    normalized
    """
    return {
        "id": "ga4gh:SL.ZDjdLp0iGarNKJcflacojz3Yy37JYQ7T",
        "sequence_id": "ga4gh:SQ.Zu7h9AggXxhTaGVsy7h_EZSChSZGcmgX",
        "start": {"value": 49531261, "type": "Number"},
        "end": {"value": 49531262, "type": "Number"},
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="session")
def genomic_dup1_38_cn(genomic_dup1_seq_loc_not_normalized):
    """Create test fixture for copy number count dup1 on GRCh38"""
    return {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.Zo9E9mi0DK7a93686jPbmeWcmQ_zmSgK",
        "subject": genomic_dup1_seq_loc_not_normalized,
        "copies": {"type": "Number", "value": 3},
    }


@pytest.fixture(scope="session")
def genomic_dup2_seq_loc_normalized():
    """Create genomic dup2 sequence location"""
    return {
        "id": "ga4gh:SL.qT3x2Jj3SS-cTC9PSuIS7rICY39SeENI",
        "sequence_id": "ga4gh:SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
        "start": {"value": 33211289, "type": "Number"},
        "end": {"value": 33211293, "type": "Number"},
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="session")
def genomic_dup2_38_cn(genomic_dup2_seq_loc_normalized):
    """Create test fixture for copy number count dup2 on GRCh38"""
    return {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN._zvzpBI2zyVZ6AMMrmJcrh83EvAaiYCA",
        "subject": genomic_dup2_seq_loc_normalized,
        "copies": {"type": "Number", "value": 3},
    }


@pytest.fixture(scope="session")
def genomic_del3_dup3_loc_not_normalized():
    """Create genomic del3 dup3 sequence location"""
    return {
        "id": "ga4gh:SL.RANaZSqxDM1hoJeXZAfSXErh6XqPrijo",
        "sequence_id": "ga4gh:SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
        "start": {"min": 31060226, "max": 31100350, "type": "DefiniteRange"},
        "end": {"min": 33274279, "max": 33417152, "type": "DefiniteRange"},
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="session")
def genomic_dup4_loc():
    """Create genomic dup4 sequence location"""
    return {
        "id": "ga4gh:SL.uD8efGXIXdiMNyHX4MogVF0jA28jIWb4",
        "sequence_id": "ga4gh:SQ.-A1QmD_MatoqxvgVxBLZTONHz9-c7nQo",
        "start": {"value": 30417575, "comparator": "<=", "type": "IndefiniteRange"},
        "end": {"value": 31394018, "comparator": ">=", "type": "IndefiniteRange"},
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="session")
def genomic_dup5_loc():
    """Create genomic dup5 sequence location"""
    return {
        "id": "ga4gh:SL.GzmuP1MBA9qILR8fVFhp4BdUEcaLwKaR",
        "sequence_id": "ga4gh:SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
        "start": {"value": 154021811, "comparator": "<=", "type": "IndefiniteRange"},
        "end": {"value": 154092209, "type": "Number"},
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="session")
def genomic_dup6_loc():
    """Create genomic dup6 sequence location"""
    return {
        "id": "ga4gh:SL.8j0dwTvx7zKHVk2JCDt1eqoEt-o993hg",
        "sequence_id": "ga4gh:SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
        "start": {"value": 154021811, "type": "Number"},
        "end": {"value": 154092209, "comparator": ">=", "type": "IndefiniteRange"},
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="session")
def genomic_del1_seq_loc():
    """Create genomic del1 sequence location"""
    return {
        "id": "ga4gh:SL.UI2JMyUJaYVNqYd1MI-sAxPsxnolalcw",
        "sequence_id": "ga4gh:SQ.Zu7h9AggXxhTaGVsy7h_EZSChSZGcmgX",
        "start": {"value": 10149810, "type": "Number"},
        "end": {"value": 10149811, "type": "Number"},
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="session")
def genomic_del1():
    """Create test fixture containing params for genomic del VD."""
    params = {
        "id": "normalize.variation:NC_000003.12%3Ag.10149811del",
        "type": "VariationDescriptor",
        "variation": dict(),
        "molecule_context": "genomic",
        "vrs_ref_allele_seq": "T",
    }
    return params


@pytest.fixture(scope="session")
def genomic_del1_lse(genomic_del1, genomic_del1_seq_loc):
    """Create a test fixture for genomic del LSE."""
    params = deepcopy(genomic_del1)
    params["variation"] = {
        "type": "Allele",
        "id": "ga4gh:VA.FVRzUwKTV-A-8zvxPUyREBR9mCunjIPO",
        "location": genomic_del1_seq_loc,
        "state": {"type": "LiteralSequenceExpression", "sequence": ""},
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="session")
def genomic_del1_38_cn(genomic_del1_seq_loc):
    """Create test fixture for copy number count del1 on GRCh38"""
    return {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.4SoY4bWbd7hDMJmIPr2vVjCMWetgtunL",
        "subject": genomic_del1_seq_loc,
        "copies": {"type": "Number", "value": 1},
    }


@pytest.fixture(scope="session")
def genomic_del2_seq_loc():
    """Create genomic del2 sequence location"""
    return {
        "id": "ga4gh:SL.rTvq-fnraKFAvs9aMgduEsi4Z7RTaxO5",
        "sequence_id": "ga4gh:SQ.Zu7h9AggXxhTaGVsy7h_EZSChSZGcmgX",
        "start": {"value": 10146594, "type": "Number"},
        "end": {"value": 10146613, "type": "Number"},
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="session")
def genomic_del2():
    """Create test fixture containing params for genomic del VD."""
    params = {
        "id": "normalize.variation:NC_000003.12%3Ag.10146595_10146613del",
        "type": "VariationDescriptor",
        "variation": dict(),
        "molecule_context": "genomic",
        "vrs_ref_allele_seq": "ATGTTGACGGACAGCCTAT",
    }
    return params


@pytest.fixture(scope="session")
def genomic_del2_lse(genomic_del2, genomic_del2_seq_loc):
    """Create a test fixture for genomic del LSE."""
    genomic_del2["variation"] = {
        "type": "Allele",
        "id": "ga4gh:VA.UgJSDSWAaJFwhRm56LA0Gez47_PYqv0k",
        "location": genomic_del2_seq_loc,
        "state": {"type": "LiteralSequenceExpression", "sequence": ""},
    }
    return VariationDescriptor(**genomic_del2)


@pytest.fixture(scope="session")
def genomic_del2_38_cn(genomic_del2_seq_loc):
    """Create test fixture for copy number count del1 on GRCh38"""
    return {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.-4obfHvLoblo593QiHeJkCORSBT3SbNV",
        "subject": genomic_del2_seq_loc,
        "copies": {"type": "Number", "value": 1},
    }


@pytest.fixture(scope="session")
def genomic_del4_seq_loc():
    """Create genomic del4 sequence location"""
    return {
        "id": "ga4gh:SL.dRc1d9ymsXhbb439OQE830RBELZ4aMXi",
        "sequence_id": "ga4gh:SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
        "start": {"value": 31120495, "comparator": "<=", "type": "IndefiniteRange"},
        "end": {"value": 33339477, "comparator": ">=", "type": "IndefiniteRange"},
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="session")
def genomic_del5_seq_loc():
    """Create genomic del5 sequence location"""
    return {
        "id": "ga4gh:SL.t3PI8S74FRjr39sAp8l3SiJrbcRGRaMx",
        "sequence_id": "ga4gh:SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
        "start": {"value": 18575353, "comparator": "<=", "type": "IndefiniteRange"},
        "end": {"value": 18653629, "type": "Number"},
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="session")
def genomic_del6_seq_loc():
    """Create genomic del6 sequence location"""
    return {
        "id": "ga4gh:SL.sawJf1SKy79nxruXFK4lb4pt2g6AU1Fy",
        "sequence_id": "ga4gh:SQ.0iKlIQk2oZLoeOG9P1riRU6hvL5Ux8TV",
        "start": {"value": 133462763, "type": "Number"},
        "end": {"value": 133464858, "comparator": ">=", "type": "IndefiniteRange"},
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="session")
def grch38_genomic_insertion_seq_loc():
    """Create test fixture for GRCh38 genomic insertion seq location"""
    return {
        "id": "ga4gh:SL.IQbKwy5bXuJrZxrcbeazxmYMCJmEjgsW",
        "end": {"value": 39724743, "type": "Number"},
        "start": {"value": 39724731, "type": "Number"},
        "sequence_id": "ga4gh:SQ.dLZ15tNO1Ur0IcGjwc3Sdi_0A6Yf4zm7",
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="session")
def grch38_genomic_insertion_variation(grch38_genomic_insertion_seq_loc):
    """Create a test fixture for NC_000017.10:g.37880993_37880994insGCTTACGTGATG"""
    return {
        "id": "ga4gh:VA.Zx2gPVfZkaplqGNrv1rzPvj3rjx0soSd",
        "location": grch38_genomic_insertion_seq_loc,
        "state": {
            "sequence": "TACGTGATGGCTTACGTGATGGCT",
            "type": "LiteralSequenceExpression",
        },
        "type": "Allele",
    }


@pytest.fixture(scope="session")
def braf_amplification(braf_ncbi_seq_loc):
    """Create test fixture for BRAF Amplification"""
    params = {
        "id": "normalize.variation:BRAF%20Amplification",
        "type": "VariationDescriptor",
        "variation": {
            "id": "ga4gh:CX.1RJp1zW60x2t4Exc4965_a3CvYFtsL4q",
            "subject": braf_ncbi_seq_loc,
            "copy_change": "efo:0030072",
            "type": "CopyNumberChange",
        },
        "molecule_context": "genomic",
        "gene_context": "hgnc:1097",
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="session")
def prpf8_amplification(prpf8_ncbi_seq_loc):
    """Create test fixture for PRPF8 Amplification"""
    params = {
        "id": "normalize.variation:PRPF8%20AMPLIFICATION",
        "type": "VariationDescriptor",
        "variation": {
            "id": "ga4gh:CX.juBTsOcCmUPHUAMu-s-Bu0oZhj1VktTL",
            "subject": prpf8_ncbi_seq_loc,
            "copy_change": "efo:0030072",
            "type": "CopyNumberChange",
        },
        "molecule_context": "genomic",
        "gene_context": "hgnc:17340",
    }
    return VariationDescriptor(**params)


def assertion_checks(normalize_response, test_variation, label, ignore_id=False):
    """Check that normalize_response and test_variation are equal."""
    if not ignore_id:
        assert normalize_response.id == test_variation.id, "id"
    assert normalize_response.label == label
    assert normalize_response.type == test_variation.type, "type"
    if test_variation.variation.type != "Text":
        if test_variation.variation.id:
            assert (
                normalize_response.variation.id == test_variation.variation.id
            ), "variation.id"
        assert normalize_response.variation == test_variation.variation, "variation"
    else:
        if not ignore_id:
            assert normalize_response.variation.id == test_variation.variation.id
        assert normalize_response.variation.type == test_variation.variation.type
        assert (
            normalize_response.variation.definition
            == test_variation.variation.definition
        )
    assert (
        normalize_response.molecule_context == test_variation.molecule_context
    ), "molecule_context"
    assert (
        normalize_response.vrs_ref_allele_seq == test_variation.vrs_ref_allele_seq
    ), "vrs_ref_allele_seq"

    # Don't need to check actual response from gene normalizer
    # Just check that if it's given in the test fixture, that it's given in the actual
    resp_gene_context = normalize_response.gene_context
    if test_variation.gene_context:
        assert resp_gene_context
    else:
        assert not resp_gene_context
