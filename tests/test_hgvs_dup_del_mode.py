"""Module for testing HGVS Dup Del mode."""
import pytest
from ga4gh.vrs import models

from tests.conftest import assertion_checks
from variation.schemas.normalize_response_schema import HGVSDupDelModeOption


@pytest.fixture(scope="module")
def test_handler(test_query_handler):
    """Create test fixture for normalize handler"""
    return test_query_handler.normalize_handler


@pytest.fixture(scope="module")
def genomic_dup1_lse(genomic_dup1_seq_loc_normalized):
    """Create a test fixture for genomic dup LSE."""
    params = {
        "type": "Allele",
        "id": "ga4gh:VA.J4tgSfdlqvfFtIFW2QY_ux7RKzFco2pd",
        "location": genomic_dup1_seq_loc_normalized,
        "state": {"type": "LiteralSequenceExpression", "sequence": "GGG"},
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def genomic_dup1_cx(genomic_dup1_seq_loc_not_normalized):
    """Create a test fixture for genomic dup copy number change."""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.xyN7zF_CrP_CwIvLIUxlsVUWoIqHNT1Y",
        "subject": genomic_dup1_seq_loc_not_normalized,
        "copy_change": "efo:0030072",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_dup1_rse(genomic_dup1_seq_loc_normalized):
    """Create a test fixture for genomic dup RSE."""
    params = {
        "type": "Allele",
        "id": "ga4gh:VA.C32Z4YjYHzsEZwzbW85O-jz_CBWl1Blu",
        "location": genomic_dup1_seq_loc_normalized,
        "state": {
            "type": "RepeatedSequenceExpression",
            "seq_expr": {
                "type": "DerivedSequenceExpression",
                "location": genomic_dup1_seq_loc_normalized,
                "reverse_complement": False,
            },
            "count": {"type": "Number", "value": 2},
        },
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def genomic_dup1_free_text_seq_loc_normalized():
    """Create genomic dup1 free text sequence subject"""
    return {
        "id": "ga4gh:SL.CnJQs8qYCvlsG8NM2XFEU5IjvjQaCVx1",
        "sequence_id": "ga4gh:SQ.tpvbnWsfEGqip8gJQZnWJAF8-bWDUDKd",
        "start": {"value": 1032, "type": "Number"},
        "end": {"value": 1034, "type": "Number"},
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_dup1_free_text_seq_loc_not_normalized():
    """Create genomic dup1 free text sequence subject"""
    return {
        "id": "ga4gh:SL.5mOwf9ZXKAMLC4i97DeYlq4CwHT1AJVR",
        "sequence_id": "ga4gh:SQ.tpvbnWsfEGqip8gJQZnWJAF8-bWDUDKd",
        "start": {"value": 1033, "type": "Number"},
        "end": {"value": 1034, "type": "Number"},
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_dup1_free_text_lse(genomic_dup1_free_text_seq_loc_normalized):
    """Create a test fixture for genomic dup LSE."""
    params = {
        "type": "Allele",
        "id": "ga4gh:VA.qkuCWe2_OQeWHXjHUa8NQpDitjhYpNhm",
        "location": genomic_dup1_free_text_seq_loc_normalized,
        "state": {"type": "LiteralSequenceExpression", "sequence": "GGG"},
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def genomic_dup1_free_text_cn(genomic_dup1_free_text_seq_loc_not_normalized):
    """Create a test fixture for genomic dup copy number count."""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.klEnlnqSMqfGoqsbK-s5HOeIi8ZVDwdH",
        "subject": genomic_dup1_free_text_seq_loc_not_normalized,
        "copies": {"type": "Number", "value": 3},
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_dup1_free_text_rse(genomic_dup1_free_text_seq_loc_normalized):
    """Create a test fixture for genomic dup RSE."""
    params = {
        "type": "Allele",
        "id": "ga4gh:VA.sY_rNabMagcE7AQW2gAhOII3kvArGjAx",
        "location": genomic_dup1_free_text_seq_loc_normalized,
        "state": {
            "type": "RepeatedSequenceExpression",
            "seq_expr": {
                "type": "DerivedSequenceExpression",
                "location": genomic_dup1_free_text_seq_loc_normalized,
                "reverse_complement": False,
            },
            "count": {"type": "Number", "value": 2},
        },
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def genomic_dup2_lse(genomic_dup2_seq_loc_normalized):
    """Create a test fixture for genomic dup LSE."""
    params = {
        "type": "Allele",
        "id": "ga4gh:VA.hUv2sNAok2K4UTiruVtvVpEHs98SGa8t",
        "location": genomic_dup2_seq_loc_normalized,
        "state": {"type": "LiteralSequenceExpression", "sequence": "TCTATCTA"},
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def genomic_dup2_cx(genomic_dup2_seq_loc_normalized):
    """Create a test fixture for genomic dup copy number change."""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.gMKznAuUtQo4niOGcJOdkgJbieFp8pV9",
        "subject": genomic_dup2_seq_loc_normalized,
        "copy_change": "efo:0030070",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_dup2_rse(genomic_dup2_seq_loc_normalized):
    """Create a test fixture for genomic dup RSE."""
    params = {
        "type": "Allele",
        "id": "ga4gh:VA.BmEpOvDuJjnXg1nU20kSAhWmUHaOt-Ru",
        "location": genomic_dup2_seq_loc_normalized,
        "state": {
            "type": "RepeatedSequenceExpression",
            "seq_expr": {
                "type": "DerivedSequenceExpression",
                "location": genomic_dup2_seq_loc_normalized,
                "reverse_complement": False,
            },
            "count": {"type": "Number", "value": 2},
        },
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def seq_loc_gt_100_bp():
    """Create seq loc for positions 33211290, 33211490 on NC_000023.11"""
    return {
        "id": "ga4gh:SL.Fj378g_oxkQoX7XMiOihnl8ELeNJ6p4_",
        "sequence_id": "ga4gh:SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
        "start": {"value": 33211289, "type": "Number"},
        "end": {"value": 33211490, "type": "Number"},
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_dup2_rse2(seq_loc_gt_100_bp):
    """Create a test fixture for genomic dup RSE where bp > 100."""
    params = {
        "type": "Allele",
        "id": "ga4gh:VA.rnWFuE-BJxh3HMllyS43v9C9MnoOJitS",
        "location": seq_loc_gt_100_bp,
        "state": {
            "type": "RepeatedSequenceExpression",
            "seq_expr": {
                "type": "DerivedSequenceExpression",
                "location": seq_loc_gt_100_bp,
                "reverse_complement": False,
            },
            "count": {"type": "Number", "value": 2},
        },
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def genomic_dup2_free_text_seq_loc():
    """Create genomic dup2 free text sequence subject"""
    return {
        "id": "ga4gh:SL.Ij3eMZP3euBYYAAcxBgayh0miq6wGvcN",
        "sequence_id": "ga4gh:SQ.1DeZLYHMnd-smp3GDlpRxETb9_0AokO7",
        "start": {"value": 256, "type": "Number"},
        "end": {"value": 260, "type": "Number"},
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_dup2_free_text_default(genomic_dup2_free_text_seq_loc):
    """Create a test fixture for genomic dup default and LSE."""
    params = {
        "type": "Allele",
        "id": "ga4gh:VA.ofKjBdAtvrfIAng0zLYiScDbfzIMfZSm",
        "location": genomic_dup2_free_text_seq_loc,
        "state": {"type": "LiteralSequenceExpression", "sequence": "TAGATAGA"},
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def genomic_dup2_free_text_cn(genomic_dup2_free_text_seq_loc):
    """Create a test fixture for genomic dup copy number count."""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.gNWVM7yYOrhA5rzRs3zrj3stQT3Phvyn",
        "subject": genomic_dup2_free_text_seq_loc,
        "copies": {"type": "Number", "value": 3},
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_dup2_free_text_rse(genomic_dup2_free_text_seq_loc):
    """Create a test fixture for genomic dup RSE."""
    params = {
        "type": "Allele",
        "id": "ga4gh:VA.Ys6Td_en3Nc72EMtr6fUEGp2WxDI-CNs",
        "location": genomic_dup2_free_text_seq_loc,
        "state": {
            "type": "RepeatedSequenceExpression",
            "seq_expr": {
                "type": "DerivedSequenceExpression",
                "location": genomic_dup2_free_text_seq_loc,
                "reverse_complement": False,
            },
            "count": {"type": "Number", "value": 2},
        },
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def genomic_dup3_cn(genomic_del3_dup3_loc_not_normalized):
    """Create a test fixture for genomic dup copy number count."""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.a2-YlA01Q-y2zAOgn4r-cP9hiIEGaaqq",
        "subject": genomic_del3_dup3_loc_not_normalized,
        "copies": {"type": "Number", "value": 2},
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_dup3_cx(genomic_del3_dup3_loc_not_normalized):
    """Create a test fixture for genomic dup copy number change."""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.cmQ2h4FXbNHHTZHJ8N0sOXIE01RfyBxb",
        "subject": genomic_del3_dup3_loc_not_normalized,
        "copy_change": "efo:0030070",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_dup3_free_text_subject():
    """Create test fixture for genomic dup3 free text subject"""
    return {
        "id": "ga4gh:SL.tdjgxR-PHAmEUyQJHQSPo1CI0Vgn19Gg",
        "sequence_id": "ga4gh:SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
        "start": {"min": 31147273, "max": 31147277, "type": "DefiniteRange"},
        "end": {"min": 31182738, "max": 31182740, "type": "DefiniteRange"},
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_dup3_free_text_cx(genomic_dup3_free_text_subject):
    """Create a test fixture for genomic dup copy number change."""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.C0njAM6WcOua46XljTLLCVQkNVIChQU0",
        "subject": genomic_dup3_free_text_subject,
        "copy_change": "efo:0030070",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_dup3_free_text_cn(genomic_dup3_free_text, genomic_dup3_free_text_subject):
    """Create a test fixture for genomic dup copy number count."""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.yfsv4o1sfFpFUxKTS1de9ODUN8fSWMQh",
        "subject": genomic_dup3_free_text_subject,
        "copies": {"type": "Number", "value": 4},
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_dup4_cx(genomic_dup4_loc):
    """Create a test fixture for genomic dup copy number change."""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.b7MlHXZcBEe4M81yZviR-Fujo4OKrcw3",
        "subject": genomic_dup4_loc,
        "copy_change": "efo:0030070",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_dup4_cn(genomic_dup4_loc):
    """Create a test fixture for genomic dup copy number count."""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.Fvt2VSCSLNROTjZd20vKZUzdMYxz-cfZ",
        "subject": genomic_dup4_loc,
        "copies": {"type": "Number", "value": 3},
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_dup4_free_text_subject():
    """Create test fixture for genomic dup4 free text subject"""
    return {
        "id": "ga4gh:SL.C177kUIMvRZlQySAXkDN9K4pg5HWqbo0",
        "sequence_id": "ga4gh:SQ.dLZ15tNO1Ur0IcGjwc3Sdi_0A6Yf4zm7",
        "start": {"value": 1674441, "comparator": "<=", "type": "IndefiniteRange"},
        "end": {"value": 1684571, "comparator": ">=", "type": "IndefiniteRange"},
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_dup4_free_text_cx(genomic_dup4_free_text_subject):
    """Create a test fixture for genomic dup copy number change."""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.vhV4v1BK5hFQXHX9qysZFpIYhrgI5-Nf",
        "subject": genomic_dup4_free_text_subject,
        "copy_change": "efo:0030070",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_dup4_free_text_cn(genomic_dup4_free_text_subject):
    """Create a test fixture for genomic dup copy number count."""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.oFSw_z1nfk0Yr1mCIi3uKIF-KGC-x8lQ",
        "subject": genomic_dup4_free_text_subject,
        "copies": {"type": "Number", "value": 3},
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_dup5_cx(genomic_dup5_loc):
    """Create a test fixture for genomic dup5 copy number change."""
    params = {
        {
            "type": "CopyNumberChange",
            "id": "ga4gh:CX.rBUl0BJiar7-NX0pU23NEln2KAifk03p",
            "subject": genomic_dup5_loc,
            "copy_change": "efo:0030070",
        }
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_dup5_cn(genomic_dup5_loc):
    """Create a test fixture for genomic dup5 copy number count."""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.uhQfl3d-UonMZDnK6FAASUtY8FTaHplR",
        "subject": genomic_dup5_loc,
        "copies": {"type": "Number", "value": 3},
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_dup6_cx(genomic_dup6_loc):
    """Create a test fixture for genomic dup copy number change."""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.IbkRumjnx-rxHbWONhJKgl-9zixCnArU",
        "subject": genomic_dup6_loc,
        "copy_change": "efo:0030070",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_dup6_cn(genomic_dup6_loc):
    """Create a test fixture for genomic dup copy number count."""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.Tuvj2C1xYVUF0UL3CSB1GZjCRC6K_eSt",
        "subject": genomic_dup6_loc,
        "copies": {"type": "Number", "value": 2},
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_del1_lse(genomic_del1_seq_loc):
    """Create a test fixture for genomic del LSE."""
    params = {
        "type": "Allele",
        "id": "ga4gh:VA.FVRzUwKTV-A-8zvxPUyREBR9mCunjIPO",
        "location": genomic_del1_seq_loc,
        "state": {"type": "LiteralSequenceExpression", "sequence": ""},
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def genomic_del1_cx(genomic_del1_seq_loc):
    """Create a test fixture for genomic del copy number change."""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.p8k6YurGbJALwJRjvKfq156CRUetenB6",
        "subject": genomic_del1_seq_loc,
        "copy_change": "efo:0030064",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_del1_rse(genomic_del1_seq_loc):
    """Create a test fixture for genomic del RSE."""
    params = {
        "type": "Allele",
        "id": "ga4gh:VA.Mubi52bBVfOHkfemgLXVg1vtl6WLfyxe",
        "location": genomic_del1_seq_loc,
        "state": {
            "type": "RepeatedSequenceExpression",
            "seq_expr": {
                "type": "DerivedSequenceExpression",
                "location": genomic_del1_seq_loc,
                "reverse_complement": False,
            },
            "count": {"type": "Number", "value": 0},
        },
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def genomic_del1_free_text_seq_loc():
    """Create genomic del1 free text sequence subject"""
    return {
        "id": "ga4gh:SL.C8wsPU7c4uq-YG88CXZzEldP888a1FMm",
        "sequence_id": "ga4gh:SQ.xBKOKptLLDr-k4hTyCetvARn16pDS_rW",
        "start": {"value": 557, "type": "Number"},
        "end": {"value": 558, "type": "Number"},
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_del1_free_text_lse(genomic_del1_free_text_seq_loc):
    """Create a test fixture for genomic del LSE."""
    params = {
        "type": "Allele",
        "id": "ga4gh:VA.VbssoYULETSOoFI0FwjGypt2YinbTEOc",
        "location": genomic_del1_free_text_seq_loc,
        "state": {"type": "LiteralSequenceExpression", "sequence": ""},
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def genomic_del1_free_text_cn(genomic_del1_free_text_seq_loc):
    """Create a test fixture for genomic del copy number count."""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.ftCNLk6wg-0szErbIRcEZLhF6aUXd-3D",
        "subject": genomic_del1_free_text_seq_loc,
        "copies": {"type": "Number", "value": 1},
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_del1_free_text_rse(genomic_del1_free_text_seq_loc):
    """Create a test fixture for genomic del RSE."""
    params = {
        "type": "Allele",
        "id": "ga4gh:VA.FMSY98r0s6L1GF5YIELVx0bmmpFmJ4n-",
        "location": genomic_del1_free_text_seq_loc,
        "state": {
            "type": "RepeatedSequenceExpression",
            "seq_expr": {
                "type": "DerivedSequenceExpression",
                "location": genomic_del1_free_text_seq_loc,
                "reverse_complement": False,
            },
            "count": {"type": "Number", "value": 0},
        },
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def genomic_del2_lse(genomic_del2_seq_loc):
    """Create a test fixture for genomic del LSE."""
    params = {
        "type": "Allele",
        "id": "ga4gh:VA.UgJSDSWAaJFwhRm56LA0Gez47_PYqv0k",
        "location": genomic_del2_seq_loc,
        "state": {"type": "LiteralSequenceExpression", "sequence": ""},
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def genomic_del2_lse2(seq_loc_gt_100_bp):
    """Create a test fixture for genomic del LSE where bp > 100."""
    params = {
        "type": "Allele",
        "id": "ga4gh:VA.7YNhEgDjlx_wc3Ci0Dt2iNAtIn9En0BL",
        "location": seq_loc_gt_100_bp,
        "state": {
            "type": "LiteralSequenceExpression",
            "sequence": "",
        },
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def genomic_del2_cx(genomic_del2_seq_loc):
    """Create a test fixture for genomic del copy number change."""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.sfbYfDVrJCocGMvSivZAppmAnZ21Sp3-",
        "subject": genomic_del2_seq_loc,
        "copy_change": "efo:0030069",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_del2_rse(genomic_del2_seq_loc):
    """Create a test fixture for genomic del RSE."""
    params = {
        "type": "Allele",
        "id": "ga4gh:VA._QLDzH5kHqog6-RaKOLg36EmEtWP_qE7",
        "location": genomic_del2_seq_loc,
        "state": {
            "type": "RepeatedSequenceExpression",
            "seq_expr": {
                "type": "DerivedSequenceExpression",
                "location": genomic_del2_seq_loc,
                "reverse_complement": False,
            },
            "count": {"type": "Number", "value": 0},
        },
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def genomic_del2_free_text_seq_loc():
    """Create genomic del2 free text sequence subject"""
    return {
        "id": "ga4gh:SL.0ietJcxUtWPRKr_9XtqQoH3cN4XfCrDM",
        "sequence_id": "ga4gh:SQ.xBKOKptLLDr-k4hTyCetvARn16pDS_rW",
        "start": {"value": 491, "type": "Number"},
        "end": {"value": 510, "type": "Number"},
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_del2_free_text_default(genomic_del2_free_text_seq_loc):
    """Create a test fixture for genomic del default and LSE."""
    params = {
        "type": "Allele",
        "id": "ga4gh:VA.6wTYBh0btGq6SlXDu4V7iEK9UrehXS-6",
        "location": genomic_del2_free_text_seq_loc,
        "state": {"type": "LiteralSequenceExpression", "sequence": ""},
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def genomic_del2_free_text_cnv(genomic_del2_free_text_seq_loc):
    """Create a test fixture for genomic del CNV."""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.4gVfUDeIuWIwLhBU34UrjZWoKbWdG3B7",
        "subject": genomic_del2_free_text_seq_loc,
        "copies": {"type": "Number", "value": 1},
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_del2_free_text_rse(genomic_del2_free_text_seq_loc):
    """Create a test fixture for genomic del RSE."""
    params = {
        "type": "Allele",
        "id": "ga4gh:VA.bXvkvrOgc7R9KG2MM--H7dRJEDv-CEVa",
        "location": genomic_del2_free_text_seq_loc,
        "state": {
            "type": "RepeatedSequenceExpression",
            "seq_expr": {
                "type": "DerivedSequenceExpression",
                "location": genomic_del2_free_text_seq_loc,
                "reverse_complement": False,
            },
            "count": {"type": "Number", "value": 0},
        },
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def genomic_del3_cx(genomic_del3_dup3_loc_not_normalized):
    """Create a test fixture for genomic del copy number change."""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.XOE1DDGAyMEUeknUAXNHFkt7FKsQBUfh",
        "subject": genomic_del3_dup3_loc_not_normalized,
        "copy_change": "efo:0030067",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_del3_cn(genomic_del3_dup3_loc_not_normalized):
    """Create a test fixture for genomic del copy number count."""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.a2-YlA01Q-y2zAOgn4r-cP9hiIEGaaqq",
        "subject": genomic_del3_dup3_loc_not_normalized,
        "copies": {"type": "Number", "value": 2},
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_del3_free_text_subject():
    """Create test fixture for genomic del3 free text subject"""
    return {
        "id": "ga4gh:SL.rqRKvcOMBF9hB5DMhRZ50hu7pnDlCMpi",
        "sequence_id": "ga4gh:SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
        "start": {"min": 68839264, "max": 68839267, "type": "DefiniteRange"},
        "end": {"min": 68841121, "max": 68841126, "type": "DefiniteRange"},
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_del3_free_text_cx(genomic_del3_free_text_subject):
    """Create a test fixture for genomic del copy number change."""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.k7s-BiYoxXGlBaGXpqcUdimdsUygemat",
        "subject": genomic_del3_free_text_subject,
        "copy_change": "efo:0030067",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_del3_free_text_cn(genomic_del3_free_text_subject):
    """Create a test fixture for genomic del copy number count."""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.qWnA8lNUMVz8vNKikQ6dWFjdV5E8FJ3U",
        "subject": genomic_del3_free_text_subject,
        "copies": {"type": "Number", "value": 2},
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_del4_cx(genomic_del4_seq_loc):
    """Create a test fixture for genomic del copy number change."""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.fYi0dG6Q8kACkyY-ICBzzvslv-ONWrPF",
        "subject": genomic_del4_seq_loc,
        "copy_change": "efo:0030067",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_del4_cn(genomic_del4_seq_loc):
    """Create a test fixture for genomic del copy number count."""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.vU4Q8TzrKFU3aQvSn2RiRS3ikh58gw_3",
        "subject": genomic_del4_seq_loc,
        "copies": {"type": "Number", "value": 1},
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_del4_free_text_subject():
    """Create test fixture for genomic del4 free text subject"""
    return {
        "id": "ga4gh:SL.TxgMX9W0YT_v1ix8uLviaLGdcqMd6WRg",
        "sequence_id": "ga4gh:SQ.pnAqCRBrTsUoBghSD1yp_jXWSmlbdh4g",
        "start": {"value": 227022027, "comparator": "<=", "type": "IndefiniteRange"},
        "end": {"value": 227025830, "comparator": ">=", "type": "IndefiniteRange"},
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_del4_free_text_cx(genomic_del4_free_text_subject):
    """Create a test fixture for genomic del copy number change."""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.q5LCBF00ZpNiutyADSXuiM0-u1--cHSA",
        "subject": genomic_del4_free_text_subject,
        "copy_change": "efo:0030067",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_del4_free_text_cn(genomic_del4_free_text_subject):
    """Create a test fixture for genomic del copy number count."""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.FF3yHxoyJB0JUt9FEeYuKx5vVOUe-w39",
        "subject": genomic_del4_free_text_subject,
        "copies": {"type": "Number", "value": 1},
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_uncertain_del_2():
    """Create a genomic uncertain deletion on chr 2 test fixture."""
    params = {
        "id": "ga4gh:CX.0pYLRyzBmiHSiFa_IQqrz8H6iNEYqobh",
        "subject": {
            "id": "ga4gh:SL.gUeB872FGVaphqoSAfI0gz4KXJvpZKL_",
            "sequence_id": "ga4gh:SQ.pnAqCRBrTsUoBghSD1yp_jXWSmlbdh4g",
            "start": {
                "value": 110104899,
                "comparator": "<=",
                "type": "IndefiniteRange",
            },
            "end": {
                "value": 110207160,
                "comparator": ">=",
                "type": "IndefiniteRange",
            },
            "type": "SequenceLocation",
        },
        "copy_change": "efo:0030067",
        "type": "CopyNumberChange",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_uncertain_del_y():
    """Create a genomic uncertain deletion on chr Y test fixture."""
    params = {
        "id": "ga4gh:CX.vcj0GssihvRHGln6RNrWpFnbid70T2Bf",
        "subject": {
            "id": "ga4gh:SL.ykRzA8IFueiCG7oznnN4teL2nXXBshHV",
            "sequence_id": "ga4gh:SQ.8_liLu1aycC0tPQPFmUaGXJLDs5SbPZ5",
            "start": {
                "value": 14076801,
                "comparator": "<=",
                "type": "IndefiniteRange",
            },
            "end": {
                "value": 57165209,
                "comparator": ">=",
                "type": "IndefiniteRange",
            },
            "type": "SequenceLocation",
        },
        "copy_change": "efo:0030067",
        "type": "CopyNumberChange",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_del5_cn_var(genomic_del5_seq_loc):
    """Create genomic del5 copy number count"""
    return {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.ZKXhoSguNUA8MBU9SpXmS52FmnZMhgy3",
        "subject": genomic_del5_seq_loc,
        "copies": {"type": "Number", "value": 3},
    }


@pytest.fixture(scope="module")
def genomic_del5_cx_var(genomic_del5_seq_loc):
    """Create genomic del5 copy number change"""
    return {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.0TqknNhM5e5i-5tAKxGXv4yDJlaQgDRR",
        "subject": genomic_del5_seq_loc,
        "copy_change": "efo:0030067",
    }


@pytest.fixture(scope="module")
def genomic_del6_cx_var(genomic_del6_seq_loc):
    """Create genomic del6 copy number change"""
    return {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.kzLS3Xs3qhMarky22So-4JW_mZAIhVW2",
        "subject": genomic_del6_seq_loc,
        "copy_change": "efo:0030067",
    }


@pytest.fixture(scope="module")
def genomic_del6_cn_var(genomic_del6_seq_loc):
    """Create genomic del6 copy number count"""
    return {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.lfI0nj_cOjRwul_R_4tohZK_vqiZ-LSz",
        "subject": genomic_del6_seq_loc,
        "copies": {"type": "Number", "value": 1},
    }


@pytest.mark.asyncio
async def assert_no_variation(query_list, test_handler):
    """Make assertion checks for invalid queries"""
    for q in query_list:
        resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
        assert resp.variation is None, q


@pytest.mark.asyncio
async def test_genomic_dup1(
    test_handler,
    genomic_dup1_lse,
    genomic_dup1_38_cn,
    genomic_dup1_cx,
    genomic_dup1_rse,
    genomic_dup1_free_text_lse,
    genomic_dup1_free_text_cn,
    genomic_dup1_free_text_rse,
):
    """Test that genomic duplication works correctly."""
    # https://reg.clinicalgenome.org/redmine/projects/registry/genboree_registry/allele?hgvsOrDescriptor=NC_000003.12%3Ag.49531262dup
    q = "NC_000003.12:g.49531262dup"  # 38
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp, genomic_dup1_lse)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp, genomic_dup1_lse)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
    )
    assertion_checks(resp, genomic_dup1_38_cn)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE, copy_change="efo:0030072"
    )
    assertion_checks(resp, genomic_dup1_cx)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.REPEATED_SEQ_EXPR)
    assertion_checks(resp, genomic_dup1_rse)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.LITERAL_SEQ_EXPR)
    assertion_checks(resp, genomic_dup1_lse)

    q = "NC_000003.11:g.49568695dup"  # 37
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp, genomic_dup1_lse)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
    )
    assertion_checks(resp, genomic_dup1_38_cn)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE, copy_change="efo:0030072"
    )
    assertion_checks(resp, genomic_dup1_cx)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.REPEATED_SEQ_EXPR)
    assertion_checks(resp, genomic_dup1_rse)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.LITERAL_SEQ_EXPR)
    assertion_checks(resp, genomic_dup1_lse)

    # Free Text
    for q in ["DAG1 g.49568695dup", "DAG1 g.49531262dup"]:  # 37  # 38
        resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
        assertion_checks(resp, genomic_dup1_free_text_lse)

        resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
        assertion_checks(resp, genomic_dup1_free_text_lse)

        resp = await test_handler.normalize(
            q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
        )
        assertion_checks(resp, genomic_dup1_free_text_cn)

        resp = await test_handler.normalize(q, HGVSDupDelModeOption.REPEATED_SEQ_EXPR)
        assertion_checks(resp, genomic_dup1_free_text_rse)

        resp = await test_handler.normalize(q, HGVSDupDelModeOption.LITERAL_SEQ_EXPR)
        assertion_checks(resp, genomic_dup1_free_text_lse)

    # Invalid
    invalid_queries = [
        "NC_000007.13:g.159138670dup",
        "NC_000007.14:g.159345976dup",
        "BRAF g.140219337dup",
        "BRAF g.141024929dup",
    ]
    await assert_no_variation(invalid_queries, test_handler)


@pytest.mark.asyncio
async def test_genomic_dup2(
    test_handler,
    genomic_dup2_lse,
    genomic_dup2_38_cn,
    genomic_dup2_cx,
    genomic_dup2_rse,
    genomic_dup2_free_text_default,
    genomic_dup2_free_text_cn,
    genomic_dup2_free_text_rse,
    genomic_dup2_rse2,
):
    """Test that genomic duplication works correctly."""
    # https://reg.clinicalgenome.org/redmine/projects/registry/genboree_registry/allele?hgvsOrDescriptor=NM_004006.2%3Ac.20_23dup
    q = "NC_000023.11:g.33211290_33211293dup"  # 38
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp, genomic_dup2_lse)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
    )
    assertion_checks(resp, genomic_dup2_38_cn)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE)
    assertion_checks(resp, genomic_dup2_cx)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.REPEATED_SEQ_EXPR)
    assertion_checks(resp, genomic_dup2_rse)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.LITERAL_SEQ_EXPR)
    assertion_checks(resp, genomic_dup2_lse)

    q = "NC_000023.10:g.33229407_33229410dup"  # 37
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp, genomic_dup2_lse)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
    )
    assertion_checks(resp, genomic_dup2_38_cn)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE)
    assertion_checks(resp, genomic_dup2_cx)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.REPEATED_SEQ_EXPR)
    assertion_checks(resp, genomic_dup2_rse)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.LITERAL_SEQ_EXPR)
    assertion_checks(resp, genomic_dup2_lse)

    # Free text
    for q in ["DMD g.33211290_33211293dup", "DMD g.33229407_33229410dup"]:  # 37  # 38
        resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
        assertion_checks(resp, genomic_dup2_free_text_default)

        resp = await test_handler.normalize(
            q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
        )
        assertion_checks(resp, genomic_dup2_free_text_cn)

        resp = await test_handler.normalize(q, HGVSDupDelModeOption.REPEATED_SEQ_EXPR)
        assertion_checks(resp, genomic_dup2_free_text_rse)

        resp = await test_handler.normalize(q, HGVSDupDelModeOption.LITERAL_SEQ_EXPR)
        assertion_checks(resp, genomic_dup2_free_text_default)

    # Greater than 100 bps -> rse
    q = "NC_000023.11:g.33211290_33211490dup"
    resp = await test_handler.normalize(q)
    assertion_checks(resp, genomic_dup2_rse2)

    # Invalid
    invalid_queries = [
        "NC_000007.13:g.140413127_159138670dup",
        "NC_000007.14:g.140413127_159345976dup",
        "BRAF g.140219337_140924929dup",
        "BRAF g.140719326_141024929dup",
    ]
    await assert_no_variation(invalid_queries, test_handler)


@pytest.mark.asyncio
async def test_genomic_dup3(
    test_handler,
    genomic_dup3_cx,
    genomic_dup3_cn,
    genomic_dup3_rse_lse,
    genomic_dup3_free_text_cn,
    genomic_dup3_free_text_cx,
    genomic_dup3_free_text_rse_lse,
):
    """Test that genomic duplication works correctly."""
    q = "NC_000023.11:g.(31060227_31100351)_(33274278_33417151)dup"  # 38
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp, genomic_dup3_cx)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=1
    )
    assertion_checks(resp, genomic_dup3_cn)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE, copy_change="efo:0030070"
    )
    assertion_checks(resp, genomic_dup3_cx)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.REPEATED_SEQ_EXPR)
    assertion_checks(resp, genomic_dup3_rse_lse)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.LITERAL_SEQ_EXPR)
    assertion_checks(resp, genomic_dup3_rse_lse)

    q = "NC_000023.10:g.(31078344_31118468)_(33292395_33435268)dup"  # 37
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp, genomic_dup3_cx)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=1
    )
    assertion_checks(resp, genomic_dup3_cn)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE)
    assertion_checks(resp, genomic_dup3_cx)

    genomic_dup3_rse_lse.variation.definition = q
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.REPEATED_SEQ_EXPR)
    assertion_checks(resp, genomic_dup3_rse_lse)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.LITERAL_SEQ_EXPR)
    assertion_checks(resp, genomic_dup3_rse_lse)

    # Free Text
    for q in ["DMD g.(31147274_31147278)_(31182737_31182739)dup"]:  # 38
        resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
        assertion_checks(resp, genomic_dup3_free_text_cx)

        resp = await test_handler.normalize(
            q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=3
        )
        assertion_checks(resp, genomic_dup3_free_text_cn)

        genomic_dup3_rse_lse.variation.definition = q
        resp = await test_handler.normalize(q, HGVSDupDelModeOption.REPEATED_SEQ_EXPR)
        assertion_checks(resp, genomic_dup3_free_text_rse_lse)

        resp = await test_handler.normalize(q, HGVSDupDelModeOption.LITERAL_SEQ_EXPR)
        assertion_checks(resp, genomic_dup3_free_text_rse_lse)

    # Invalid
    invalid_queries = [
        "NC_000023.10:g.(31119221_31119227)_(31119300_155270562)dup",
        "NC_000023.11:g.(31119221_31119227)_(31119300_156040899)dup",
        "DMD g.(31060227_31100351)_(33274278_33417151)dup",
    ]
    await assert_no_variation(invalid_queries, test_handler)


@pytest.mark.asyncio
async def test_genomic_dup4(
    test_handler,
    genomic_dup4_cn,
    genomic_dup4_cx,
    genomic_dup4_rse_lse,
    genomic_dup4_free_text_cn,
    genomic_dup4_free_text_cx,
    genomic_dup4_free_text_rse_lse,
):
    """Test that genomic duplication works correctly."""
    q = "NC_000020.11:g.(?_30417576)_(31394018_?)dup"  # 38
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp, genomic_dup4_cx)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
    )
    assertion_checks(resp, genomic_dup4_cn)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE)
    assertion_checks(resp, genomic_dup4_cx)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.REPEATED_SEQ_EXPR)
    assertion_checks(resp, genomic_dup4_rse_lse)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.LITERAL_SEQ_EXPR)
    assertion_checks(resp, genomic_dup4_rse_lse)

    q = "NC_000020.10:g.(?_29652252)_(29981821_?)dup"  # 37
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp, genomic_dup4_cx)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
    )
    assertion_checks(resp, genomic_dup4_cn)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE)
    assertion_checks(resp, genomic_dup4_cx)

    genomic_dup4_rse_lse.variation.definition = q
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.REPEATED_SEQ_EXPR)
    assertion_checks(resp, genomic_dup4_rse_lse)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.LITERAL_SEQ_EXPR)
    assertion_checks(resp, genomic_dup4_rse_lse)

    # Free Text
    for q in [
        "PRPF8 g.(?_1577736)_(1587865_?)dup",  # 37
        "PRPF8 g.(?_1674442)_(1684571_?)dup",  # 38
    ]:
        resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
        assertion_checks(resp, genomic_dup4_free_text_cx)

        resp = await test_handler.normalize(
            q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
        )
        assertion_checks(resp, genomic_dup4_free_text_cn)

        genomic_dup4_rse_lse.variation.definition = q
        resp = await test_handler.normalize(q, HGVSDupDelModeOption.REPEATED_SEQ_EXPR)
        genomic_dup4_free_text_rse_lse.variation.definition = q
        assertion_checks(resp, genomic_dup4_free_text_rse_lse)

        resp = await test_handler.normalize(q, HGVSDupDelModeOption.LITERAL_SEQ_EXPR)
        assertion_checks(resp, genomic_dup4_free_text_rse_lse)

    # Invalid
    invalid_queries = [
        "NC_000020.10:g.(?_29652252)_(63025530_?)dup",
        "NC_000020.11:g.(?_29652252)_(64444169_?)dup",
        "PRPF8 g.(?_1650628)_(1684571_?)dup",
    ]
    await assert_no_variation(invalid_queries, test_handler)


@pytest.mark.asyncio
async def test_genomic_dup5(
    test_handler,
    genomic_dup5_cn,
    genomic_dup5_cx,
    genomic_dup5_rse_lse,
    genomic_dup5_free_text_rse_lse,
):
    """Test that genomic duplication works correctly."""
    q = "NC_000023.11:g.(?_154021812)_154092209dup"  # 38
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp, genomic_dup5_cx)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
    )
    assertion_checks(resp, genomic_dup5_cn)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE)
    assertion_checks(resp, genomic_dup5_cx)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.REPEATED_SEQ_EXPR)
    assertion_checks(resp, genomic_dup5_rse_lse)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.LITERAL_SEQ_EXPR)
    assertion_checks(resp, genomic_dup5_rse_lse)

    q = "NC_000023.10:g.(?_153287263)_153357667dup"  # 37
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp, genomic_dup5_cx)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
    )
    assertion_checks(resp, genomic_dup5_cn)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE)
    assertion_checks(resp, genomic_dup5_cx)

    genomic_dup5_rse_lse.variation.definition = q
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.REPEATED_SEQ_EXPR)
    assertion_checks(resp, genomic_dup5_rse_lse)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.LITERAL_SEQ_EXPR)
    assertion_checks(resp, genomic_dup5_rse_lse)

    # Free Text
    for q in [
        "MECP2 g.(?_153287263)_153357667dup",  # 37
        "MECP2 g.(?_154021812)_154092209dup",  # 38
    ]:
        resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
        assertion_checks(resp, genomic_dup5_cx)

        resp = await test_handler.normalize(
            q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
        )
        assertion_checks(resp, genomic_dup5_cn)

        genomic_dup5_free_text_rse_lse.variation.definition = q
        resp = await test_handler.normalize(q, HGVSDupDelModeOption.REPEATED_SEQ_EXPR)
        assertion_checks(resp, genomic_dup5_free_text_rse_lse)

        resp = await test_handler.normalize(q, HGVSDupDelModeOption.LITERAL_SEQ_EXPR)
        assertion_checks(resp, genomic_dup5_free_text_rse_lse)

    # Invalid
    for q in [
        "NC_000023.10:g.(?_153287263)_155270561dup",
        "NC_000023.11:g.(?_154021812)_156040896dup",
        "MECP2 g.(?_154021812)_154097733dup"  # 37
        "MECP2 g.(?_154021572)_154092209dup",  # 38
    ]:
        resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
        assert resp.variation is None, q


@pytest.mark.asyncio
async def test_genomic_dup6(
    test_handler,
    genomic_dup6_cn,
    genomic_dup6_cx,
    genomic_dup6_rse_lse,
    genomic_dup6_free_text_rse_lse,
):
    """Test that genomic duplication works correctly."""
    q = "NC_000023.11:g.154021812_(154092209_?)dup"  # 38
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp, genomic_dup6_cx)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=1
    )
    assertion_checks(resp, genomic_dup6_cn)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE)
    assertion_checks(resp, genomic_dup6_cx)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.REPEATED_SEQ_EXPR)
    assertion_checks(resp, genomic_dup6_rse_lse)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.LITERAL_SEQ_EXPR)
    assertion_checks(resp, genomic_dup6_rse_lse)

    q = "NC_000023.10:g.153287263_(153357667_?)dup"  # 37
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp, genomic_dup6_cx)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=1
    )
    assertion_checks(resp, genomic_dup6_cn)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE)
    assertion_checks(resp, genomic_dup6_cx)

    genomic_dup6_rse_lse.variation.definition = q
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.REPEATED_SEQ_EXPR)
    assertion_checks(resp, genomic_dup6_rse_lse)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.LITERAL_SEQ_EXPR)
    assertion_checks(resp, genomic_dup6_rse_lse)

    # Free Text
    for q in [
        "MECP2 g.153287263_(153357667_?)dup",  # 37
        "MECP2 g.154021812_(154092209_?)dup",  # 38
    ]:
        resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
        assertion_checks(resp, genomic_dup6_cx)

        resp = await test_handler.normalize(
            q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=1
        )
        assertion_checks(resp, genomic_dup6_cn)

        genomic_dup6_free_text_rse_lse.variation.definition = q
        resp = await test_handler.normalize(q, HGVSDupDelModeOption.REPEATED_SEQ_EXPR)
        assertion_checks(resp, genomic_dup6_free_text_rse_lse)

        resp = await test_handler.normalize(q, HGVSDupDelModeOption.LITERAL_SEQ_EXPR)
        assertion_checks(resp, genomic_dup6_free_text_rse_lse)

    # Invalid
    for q in [
        "NC_000023.10:g.153287263_(155270561_?)dup",
        "NC_000023.11:g.154021812_(156040896_?)dup",
        "MECP2 g.154021812_(154097733_?)dup"  # 37
        "MECP2 g.154021572_(154092209_?)dup",  # 38
    ]:
        resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
        assert resp.variation is None, q


@pytest.mark.asyncio
async def test_genomic_del1(
    test_handler,
    genomic_del1_lse,
    genomic_del1_38_cn,
    genomic_del1_cx,
    genomic_del1_rse,
    genomic_del1_free_text_lse,
    genomic_del1_free_text_cn,
    genomic_del1_free_text_rse,
):
    """Test that genomic deletion works correctly."""
    q = "NC_000003.12:g.10149811del"  # 38
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp, genomic_del1_lse)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
    )
    assertion_checks(resp, genomic_del1_38_cn)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE, copy_change="efo:0030064"
    )
    assertion_checks(resp, genomic_del1_cx)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.REPEATED_SEQ_EXPR)
    assertion_checks(resp, genomic_del1_rse)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.LITERAL_SEQ_EXPR)
    assertion_checks(resp, genomic_del1_lse)

    q = "NC_000003.11:g.10191495del"  # 37
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp, genomic_del1_lse)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
    )
    assertion_checks(resp, genomic_del1_38_cn)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE, copy_change="efo:0030064"
    )
    assertion_checks(resp, genomic_del1_cx)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.REPEATED_SEQ_EXPR)
    assertion_checks(resp, genomic_del1_rse)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.LITERAL_SEQ_EXPR)
    assertion_checks(resp, genomic_del1_lse)

    # Free text
    for q in ["VHL g.10191495del", "VHL g.10149811del"]:  # 37  # 38
        resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
        assertion_checks(resp, genomic_del1_free_text_lse)

        resp = await test_handler.normalize(
            q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
        )
        assertion_checks(resp, genomic_del1_free_text_cn)

        resp = await test_handler.normalize(q, HGVSDupDelModeOption.REPEATED_SEQ_EXPR)
        assertion_checks(resp, genomic_del1_free_text_rse)

        resp = await test_handler.normalize(q, HGVSDupDelModeOption.LITERAL_SEQ_EXPR)
        assertion_checks(resp, genomic_del1_free_text_lse)

    # Invalid
    invalid_queries = [
        "NC_000003.11:g.198022431del",
        "NC_000003.12:g.198295567del",
        "BRAF g.140413127del",
        "BRAF g.141024929del",
    ]
    await assert_no_variation(invalid_queries, test_handler)


@pytest.mark.asyncio
async def test_genomic_del2(
    test_handler,
    genomic_del2_lse,
    genomic_del2_38_cn,
    genomic_del2_cx,
    genomic_del2_rse,
    genomic_del2_free_text_default,
    genomic_del2_free_text_cnv,
    genomic_del2_free_text_rse,
    genomic_del2_lse2,
):
    """Test that genomic deletion works correctly."""
    q = "NC_000003.12:g.10146595_10146613del"  # 38
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp, genomic_del2_lse)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
    )
    assertion_checks(resp, genomic_del2_38_cn)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE, copy_change="efo:0030069"
    )
    assertion_checks(resp, genomic_del2_cx)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.REPEATED_SEQ_EXPR)
    assertion_checks(resp, genomic_del2_rse)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.LITERAL_SEQ_EXPR)
    assertion_checks(resp, genomic_del2_lse)

    q = "NC_000003.11:g.10188279_10188297del"  # 37
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp, genomic_del2_lse)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
    )
    assertion_checks(resp, genomic_del2_38_cn)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE, copy_change="efo:0030069"
    )
    assertion_checks(resp, genomic_del2_cx)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.REPEATED_SEQ_EXPR)
    assertion_checks(resp, genomic_del2_rse)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.LITERAL_SEQ_EXPR)
    assertion_checks(resp, genomic_del2_lse)

    # Free text
    for q in ["VHL g.10188279_10188297del", "VHL g.10146595_10146613del"]:  # 37  # 38
        resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
        assertion_checks(resp, genomic_del2_free_text_default)

        resp = await test_handler.normalize(
            q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
        )
        assertion_checks(resp, genomic_del2_free_text_cnv)

        resp = await test_handler.normalize(q, HGVSDupDelModeOption.REPEATED_SEQ_EXPR)
        assertion_checks(resp, genomic_del2_free_text_rse)

        resp = await test_handler.normalize(q, HGVSDupDelModeOption.LITERAL_SEQ_EXPR)
        assertion_checks(resp, genomic_del2_free_text_default)

    # Check that del > 100 bps returns LSE
    q = "NC_000023.11:g.33211290_33211490del"
    resp = await test_handler.normalize(q)
    assertion_checks(resp, genomic_del2_lse2)

    # gnomad vcf
    q = "3-10146594-AATGTTGACGGACAGCCTAT-A"
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp, genomic_del2_lse)

    q = "3-10188278-AATGTTGACGGACAGCCTAT-A"
    resp = await test_handler.normalize(q)
    assertion_checks(resp, genomic_del2_lse)

    # Invalid
    invalid_queries = [
        "NC_000003.12:g.10146595_198295580del",
        "NC_000003.11:g.198022435_198022437del",
        "BRAF g.140413127_140419136del",
        "BRAF g.140719326_141024929del",
    ]
    await assert_no_variation(invalid_queries, test_handler)


@pytest.mark.asyncio
async def test_genomic_del3(
    test_handler,
    genomic_del3_cn,
    genomic_del3_cx,
    genomic_del3_rse_lse,
    genomic_del3_free_text_cn,
    genomic_del3_free_text_cx,
    genomic_del3_free_text_rse_lse,
):
    """Test that genomic deletion works correctly."""
    q = "NC_000023.11:g.(31060227_31100351)_(33274278_33417151)del"  # 38
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp, genomic_del3_cx)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=3
    )
    assertion_checks(resp, genomic_del3_cn)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE)
    assertion_checks(resp, genomic_del3_cx)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.REPEATED_SEQ_EXPR)
    assertion_checks(resp, genomic_del3_rse_lse)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.LITERAL_SEQ_EXPR)
    assertion_checks(resp, genomic_del3_rse_lse)

    q = "NC_000023.10:g.(31078344_31118468)_(33292395_33435268)del"  # 37
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp, genomic_del3_cx)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=3
    )
    assertion_checks(resp, genomic_del3_cn)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE)
    assertion_checks(resp, genomic_del3_cx)

    genomic_del3_rse_lse.variation.definition = q
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.REPEATED_SEQ_EXPR)
    assertion_checks(resp, genomic_del3_rse_lse)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.LITERAL_SEQ_EXPR)
    assertion_checks(resp, genomic_del3_rse_lse)

    # Free Text
    for q in [
        "EFNB1 g.(68059108_68059111)_(68060963_68060968)del",  # 37
        "EFNB1 g.(68839265_68839268)_(68841120_68841125)del",  # 38
    ]:
        resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
        assertion_checks(resp, genomic_del3_free_text_cx)

        resp = await test_handler.normalize(
            q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=3
        )
        assertion_checks(resp, genomic_del3_free_text_cn)

        genomic_del3_free_text_rse_lse.variation.definition = q
        resp = await test_handler.normalize(q, HGVSDupDelModeOption.REPEATED_SEQ_EXPR)
        assertion_checks(resp, genomic_del3_free_text_rse_lse)

        resp = await test_handler.normalize(q, HGVSDupDelModeOption.LITERAL_SEQ_EXPR)
        assertion_checks(resp, genomic_del3_free_text_rse_lse)

    # Invalid
    invalid_queries = [
        "NC_000023.11:g.(156040880_156040883)_(156040896_156040899)del",
        "NC_000023.10:g.(155270550_155270555)_(155270560_155270562)del",
        "EFNB1 g.(68048863_68048870)_(68842150_68842152)del",  # 37
        "EFNB1 g.(68829022_68829030)_(68842150_68842161)del",  # 38
    ]
    await assert_no_variation(invalid_queries, test_handler)


@pytest.mark.asyncio
async def test_genomic_del4(
    test_handler,
    genomic_del4_cn,
    genomic_del4_cx,
    genomic_del4_rse_lse,
    genomic_uncertain_del_2,
    genomic_uncertain_del_y,
    genomic_del4_free_text_cn,
    genomic_del4_free_text_rse_lse,
    genomic_del4_free_text_cx,
):
    """Test that genomic deletion works correctly."""
    q = "NC_000023.11:g.(?_31120496)_(33339477_?)del"  # 38
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp, genomic_del4_cx)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
    )
    assertion_checks(resp, genomic_del4_cn)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE)
    assertion_checks(resp, genomic_del4_cx)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.REPEATED_SEQ_EXPR)
    assertion_checks(resp, genomic_del4_rse_lse)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.LITERAL_SEQ_EXPR)
    assertion_checks(resp, genomic_del4_rse_lse)

    q = "NC_000023.10:g.(?_31138613)_(33357594_?)del"  # 37
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp, genomic_del4_cx)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
    )
    assertion_checks(resp, genomic_del4_cn)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE)
    assertion_checks(resp, genomic_del4_cx)

    genomic_del4_rse_lse.variation.definition = q
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.REPEATED_SEQ_EXPR)
    assertion_checks(resp, genomic_del4_rse_lse)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.LITERAL_SEQ_EXPR)
    assertion_checks(resp, genomic_del4_rse_lse)

    q = "NC_000002.12:g.(?_110104900)_(110207160_?)del"
    resp = await test_handler.normalize(q)
    assertion_checks(resp, genomic_uncertain_del_2)

    q = "NC_000024.10:g.(?_14076802)_(57165209_?)del"
    resp = await test_handler.normalize(q)
    assertion_checks(resp, genomic_uncertain_del_y)

    # Free Text
    for q in ["COL4A4 g.(?_227022028)_(227025830_?)del"]:  # 38
        resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
        assertion_checks(resp, genomic_del4_free_text_cx)

        resp = await test_handler.normalize(
            q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
        )
        assertion_checks(resp, genomic_del4_free_text_cn)

        resp = await test_handler.normalize(q, HGVSDupDelModeOption.REPEATED_SEQ_EXPR)
        assertion_checks(resp, genomic_del4_free_text_rse_lse)

        resp = await test_handler.normalize(q, HGVSDupDelModeOption.LITERAL_SEQ_EXPR)
        assertion_checks(resp, genomic_del4_free_text_rse_lse)

    # Invalid
    invalid_queries = [
        "NC_000023.11:g.(?_156040899)_(156040900_?)del",
        "NC_000024.10:g.(?_155270565)_(155270568_?)del",
        "COL4A4 g.(?_227002710)_(227003710_?)del",
        "COL4A4 g.(?_227867430)_(228029276_?)del",
    ]
    await assert_no_variation(invalid_queries, test_handler)


@pytest.mark.asyncio
async def test_genomic_del5(
    test_handler,
    genomic_del5_cn_var,
    genomic_del5_cx_var,
    genomic_del5_rse_lse,
    genomic_del5_free_text_cn,
    genomic_del5_free_text_cx,
    genomic_del5_free_text_rse_lse,
):
    """Test that genomic deletion works correctly."""
    q = "NC_000023.11:g.(?_18575354)_18653629del"  # 38
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp, genomic_del5_cx_var)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=4
    )
    assertion_checks(resp, genomic_del5_cn_var)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE)
    assertion_checks(resp, genomic_del5_cx_var)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.REPEATED_SEQ_EXPR)
    assertion_checks(resp, genomic_del5_rse_lse)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.LITERAL_SEQ_EXPR)
    assertion_checks(resp, genomic_del5_rse_lse)

    q = "NC_000023.10:g.(?_18593474)_18671749del"  # 37
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp, genomic_del5_cx_var)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=4
    )
    assertion_checks(resp, genomic_del5_cn_var)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE)
    assertion_checks(resp, genomic_del5_cx_var)

    genomic_del5_rse_lse.variation.definition = q
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.REPEATED_SEQ_EXPR)
    assertion_checks(resp, genomic_del5_rse_lse)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.LITERAL_SEQ_EXPR)
    assertion_checks(resp, genomic_del5_rse_lse)

    # Free text
    for q in ["CDKL5 g.(?_18575354)_18653629del"]:
        resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
        assertion_checks(resp, genomic_del5_free_text_cx)

        resp = await test_handler.normalize(
            q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=4
        )
        assertion_checks(resp, genomic_del5_free_text_cn)

        genomic_del5_free_text_rse_lse.variation.definition = q
        resp = await test_handler.normalize(q, HGVSDupDelModeOption.REPEATED_SEQ_EXPR)
        assertion_checks(resp, genomic_del5_free_text_rse_lse)

        resp = await test_handler.normalize(q, HGVSDupDelModeOption.LITERAL_SEQ_EXPR)
        assertion_checks(resp, genomic_del5_free_text_rse_lse)

    # Invalid
    invalid_queries = [
        "NC_000023.10:g.(?_155270550)_155270570del",
        "NC_000023.11:g.(?_18593474)_18671749del"
        "CDKL5  g.(?_18443702)_18671700del",  # 37
        "CDKL5  g.(?_18425585)_18653631del",  # 38
        "CDKL5  g.(?_18425582)_18653500del",  # 38
    ]
    await assert_no_variation(invalid_queries, test_handler)


@pytest.mark.asyncio
async def test_genomic_del6(
    test_handler,
    genomic_del6_cn_var,
    genomic_del6_cx_var,
    genomic_del6_rse_lse,
    genomic_del6_free_text_cn,
    genomic_del6_free_text_cx,
    genomic_del6_free_text_rse_lse,
):
    """Test that genomic deletion works correctly."""
    q = "NC_000006.12:g.133462764_(133464858_?)del"  # 38
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp, genomic_del6_cx_var)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
    )
    assertion_checks(resp, genomic_del6_cn_var)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE)
    assertion_checks(resp, genomic_del6_cx_var)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.REPEATED_SEQ_EXPR)
    assertion_checks(resp, genomic_del6_rse_lse)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.LITERAL_SEQ_EXPR)
    assertion_checks(resp, genomic_del6_rse_lse)

    q = "NC_000006.11:g.133783902_(133785996_?)del"  # 37
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp, genomic_del6_cx_var)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
    )
    assertion_checks(resp, genomic_del6_cn_var)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE)
    assertion_checks(resp, genomic_del6_cx_var)

    genomic_del6_rse_lse.variation.definition = q
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.REPEATED_SEQ_EXPR)
    assertion_checks(resp, genomic_del6_rse_lse)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.LITERAL_SEQ_EXPR)
    assertion_checks(resp, genomic_del6_rse_lse)

    # Free text
    for q in ["EYA4 g.133462764_(133464858_?)del"]:  # 38
        resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
        assertion_checks(resp, genomic_del6_free_text_cx)

        resp = await test_handler.normalize(
            q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
        )
        assertion_checks(resp, genomic_del6_free_text_cn)

        genomic_del6_rse_lse.variation.definition = q
        resp = await test_handler.normalize(q, HGVSDupDelModeOption.REPEATED_SEQ_EXPR)
        assertion_checks(resp, genomic_del6_free_text_rse_lse)

        resp = await test_handler.normalize(q, HGVSDupDelModeOption.LITERAL_SEQ_EXPR)
        assertion_checks(resp, genomic_del6_free_text_rse_lse)

    # Invalid
    invalid_queries = [
        "NC_000006.11:g.171115069_(171115080_?)del",
        "NC_000006.12:g.170805981_(170805989_?)del"
        "EYA4 g.133561700_(133853270_?)del",  # 37
        "EYA4 g.133561651_(133561708_?)del",  # 37
        "EYA4 g.133240513_(133240600_?)del",  # 38
        "EYA4 g.133240515_(133532130_?)del",  # 38
    ]
    await assert_no_variation(invalid_queries, test_handler)


@pytest.mark.asyncio
async def test_parameters(test_handler):
    """Check that valid and invalid parameters work as intended."""
    resp = await test_handler.normalize("7-140453136-A-T")
    assert resp.variation
    assert resp.warnings == []

    q = "NC_000003.12:g.49531262dup"
    resp = await test_handler.normalize(q)
    assert resp.variation
    assert resp.warnings == []

    resp = await test_handler.normalize(q, hgvs_dup_del_mode=None)
    assert resp.variation
    assert resp.warnings == []

    resp = await test_handler.normalize(
        q, hgvs_dup_del_mode=HGVSDupDelModeOption.COPY_NUMBER_COUNT
    )
    assert resp.variation is None
    assert resp.warnings == ["copy_number_count mode requires `baseline_copies`"]
