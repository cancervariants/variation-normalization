"""Module for testing HGVS Dup Del mode."""
from copy import deepcopy

import pytest
from ga4gh.vrsatile.pydantic.vrsatile_models import VariationDescriptor

from tests.conftest import assertion_checks
from variation.schemas.normalize_response_schema import HGVSDupDelModeOption


@pytest.fixture(scope="module")
def test_handler(test_query_handler):
    """Create test fixture for normalize handler"""
    return test_query_handler.normalize_handler


@pytest.fixture(scope="module")
def genomic_dup1_normalized():
    """Create test fixture containing params for genomic dup VD."""
    params = {
        "id": "normalize.variation:NC_000003.12%3Ag.49531262dup",
        "type": "VariationDescriptor",
        "variation": dict(),
        "molecule_context": "genomic",
        "vrs_ref_allele_seq": "GG",
    }
    return params


@pytest.fixture(scope="module")
def genomic_dup1_not_normalized():
    """Create test fixture containing params for genomic dup not normalized VD."""
    params = {
        "id": "normalize.variation:NC_000003.12%3Ag.49531262dup",
        "type": "VariationDescriptor",
        "variation": dict(),
        "molecule_context": "genomic",
        "vrs_ref_allele_seq": "G",
    }
    return params


@pytest.fixture(scope="module")
def genomic_dup1_lse(genomic_dup1_normalized, genomic_dup1_seq_loc_normalized):
    """Create a test fixture for genomic dup LSE."""
    params = deepcopy(genomic_dup1_normalized)
    params["variation"] = {
        "type": "Allele",
        "id": "ga4gh:VA.J4tgSfdlqvfFtIFW2QY_ux7RKzFco2pd",
        "location": genomic_dup1_seq_loc_normalized,
        "state": {"type": "LiteralSequenceExpression", "sequence": "GGG"},
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_dup1_cn(genomic_dup1_not_normalized, genomic_dup1_38_cn):
    """Create a test fixture for genomic dup copy number count."""
    params = deepcopy(genomic_dup1_not_normalized)
    params["variation"] = genomic_dup1_38_cn
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_dup1_cx(genomic_dup1_not_normalized, genomic_dup1_seq_loc_not_normalized):
    """Create a test fixture for genomic dup copy number change."""
    params = deepcopy(genomic_dup1_not_normalized)
    params["variation"] = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.xyN7zF_CrP_CwIvLIUxlsVUWoIqHNT1Y",
        "subject": genomic_dup1_seq_loc_not_normalized,
        "copy_change": "efo:0030072",
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_dup1_rse(genomic_dup1_normalized, genomic_dup1_seq_loc_normalized):
    """Create a test fixture for genomic dup RSE."""
    params = deepcopy(genomic_dup1_normalized)
    params["variation"] = {
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
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_dup1_free_text_normalized():
    """Create test fixture containing params for genomic dup VD."""
    params = {
        "id": "normalize.variation:DAG1%20g.49568695dup",
        "type": "VariationDescriptor",
        "variation": dict(),
        "molecule_context": "transcript",
        "vrs_ref_allele_seq": "GG",
        "gene_context": "hgnc:2666",
    }
    return params


@pytest.fixture(scope="module")
def genomic_dup1_free_text_not_normalized():
    """Create test fixture containing params for genomic dup VD."""
    params = {
        "id": "normalize.variation:DAG1%20g.49568695dup",
        "type": "VariationDescriptor",
        "variation": dict(),
        "molecule_context": "transcript",
        "vrs_ref_allele_seq": "G",
        "gene_context": "hgnc:2666",
    }
    return params


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
def genomic_dup1_free_text_lse(
    genomic_dup1_free_text_normalized, genomic_dup1_free_text_seq_loc_normalized
):
    """Create a test fixture for genomic dup LSE."""
    params = deepcopy(genomic_dup1_free_text_normalized)
    params["variation"] = {
        "type": "Allele",
        "id": "ga4gh:VA.qkuCWe2_OQeWHXjHUa8NQpDitjhYpNhm",
        "location": genomic_dup1_free_text_seq_loc_normalized,
        "state": {"type": "LiteralSequenceExpression", "sequence": "GGG"},
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_dup1_free_text_cn(
    genomic_dup1_free_text_not_normalized, genomic_dup1_free_text_seq_loc_not_normalized
):
    """Create a test fixture for genomic dup copy number count."""
    params = deepcopy(genomic_dup1_free_text_not_normalized)
    params["variation"] = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.klEnlnqSMqfGoqsbK-s5HOeIi8ZVDwdH",
        "subject": genomic_dup1_free_text_seq_loc_not_normalized,
        "copies": {"type": "Number", "value": 3},
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_dup1_free_text_rse(
    genomic_dup1_free_text_normalized, genomic_dup1_free_text_seq_loc_normalized
):
    """Create a test fixture for genomic dup RSE."""
    params = deepcopy(genomic_dup1_free_text_normalized)
    params["variation"] = {
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
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_dup2_normalized():
    """Create test fixture containing params for genomic dup VD."""
    params = {
        "id": "normalize.variation:NC_000023.11%3Ag.33211290_33211293dup",
        "type": "VariationDescriptor",
        "variation": dict(),
        "molecule_context": "genomic",
        "vrs_ref_allele_seq": "TCTA",
    }
    return params


@pytest.fixture(scope="module")
def genomic_dup2_not_normalized():
    """Create test fixture containing params for genomic dup VD."""
    params = {
        "id": "normalize.variation:NC_000023.11%3Ag.33211290_33211293dup",
        "type": "VariationDescriptor",
        "variation": dict(),
        "molecule_context": "genomic",
        "vrs_ref_allele_seq": "TCTA",
    }
    return params


@pytest.fixture(scope="module")
def genomic_dup2_lse(genomic_dup2_normalized, genomic_dup2_seq_loc_normalized):
    """Create a test fixture for genomic dup LSE."""
    params = deepcopy(genomic_dup2_normalized)
    params["variation"] = {
        "type": "Allele",
        "id": "ga4gh:VA.hUv2sNAok2K4UTiruVtvVpEHs98SGa8t",
        "location": genomic_dup2_seq_loc_normalized,
        "state": {"type": "LiteralSequenceExpression", "sequence": "TCTATCTA"},
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_dup2_cn(genomic_dup2_not_normalized, genomic_dup2_38_cn):
    """Create a test fixture for genomic dup copy number count."""
    params = deepcopy(genomic_dup2_not_normalized)
    params["variation"] = genomic_dup2_38_cn
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_dup2_cx(genomic_dup2_normalized, genomic_dup2_seq_loc_normalized):
    """Create a test fixture for genomic dup copy number change."""
    params = deepcopy(genomic_dup2_normalized)
    params["variation"] = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.gMKznAuUtQo4niOGcJOdkgJbieFp8pV9",
        "subject": genomic_dup2_seq_loc_normalized,
        "copy_change": "efo:0030070",
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_dup2_rse(genomic_dup2_normalized, genomic_dup2_seq_loc_normalized):
    """Create a test fixture for genomic dup RSE."""
    params = deepcopy(genomic_dup2_normalized)
    params["variation"] = {
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
    return VariationDescriptor(**params)


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
def vrs_ref_allele_seq_gt_100_bp():
    """Create vrs_ref_allele_seq for positions 33211290, 33211490 on NC_000023.11"""
    return "TCTACTTCTTCCCACCAAAGCATTTTGAAAAGTGTATATCAAGGCAGCGATAAAAAAAACCTGGTAAAAGTTCTTCAAACTTTATTGCTCCAGTAGGCTTAAAAACAATGAGAAACCAACAAACTTCAGCAGCTTTAAAAAAAGTAACACTTCAGTTTTTCCTATTCGTTTTTCTCCGAAGGTAATTGCCTCCCAGATCTG"  # noqa: E501


@pytest.fixture(scope="module")
def genomic_dup2_rse2(seq_loc_gt_100_bp, vrs_ref_allele_seq_gt_100_bp):
    """Create a test fixture for genomic dup RSE where bp > 100."""
    params = {
        "id": "normalize.variation:NC_000023.11%3Ag.33211290_33211490dup",
        "type": "VariationDescriptor",
        "variation": {
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
        },
        "molecule_context": "genomic",
        "vrs_ref_allele_seq": vrs_ref_allele_seq_gt_100_bp,
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_dup2_free_text():
    """Create test fixture containing params for genomic dup VD."""
    params = {
        "id": "normalize.variation:TSC2%20g.2137939_2137949dup",
        "type": "VariationDescriptor",
        "variation": dict(),
        "molecule_context": "transcript",
        "vrs_ref_allele_seq": "TAGA",
        "gene_context": "hgnc:2928",
    }
    return params


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
def genomic_dup2_free_text_default(
    genomic_dup2_free_text, genomic_dup2_free_text_seq_loc
):
    """Create a test fixture for genomic dup default and LSE."""
    params = deepcopy(genomic_dup2_free_text)
    params["variation"] = {
        "type": "Allele",
        "id": "ga4gh:VA.ofKjBdAtvrfIAng0zLYiScDbfzIMfZSm",
        "location": genomic_dup2_free_text_seq_loc,
        "state": {"type": "LiteralSequenceExpression", "sequence": "TAGATAGA"},
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_dup2_free_text_cn(genomic_dup2_free_text, genomic_dup2_free_text_seq_loc):
    """Create a test fixture for genomic dup copy number count."""
    params = deepcopy(genomic_dup2_free_text)
    params["variation"] = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.gNWVM7yYOrhA5rzRs3zrj3stQT3Phvyn",
        "subject": genomic_dup2_free_text_seq_loc,
        "copies": {"type": "Number", "value": 3},
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_dup2_free_text_rse(genomic_dup2_free_text, genomic_dup2_free_text_seq_loc):
    """Create a test fixture for genomic dup RSE."""
    params = deepcopy(genomic_dup2_free_text)
    params["variation"] = {
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
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_dup3():
    """Create test fixture containing params for genomic dup VD."""
    params = {
        "id": "normalize.variation:NC_000023.11%3Ag.%2831060227_31100351%29_%2833274278_33417151%29dup",  # noqa: E501
        "type": "VariationDescriptor",
        "variation": dict(),
        "molecule_context": "genomic",
        "vrs_ref_allele_seq": None,
    }
    return params


@pytest.fixture(scope="module")
def genomic_dup3_cn(genomic_dup3, genomic_del3_dup3_loc_not_normalized):
    """Create a test fixture for genomic dup copy number count."""
    params = deepcopy(genomic_dup3)
    params["variation"] = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.a2-YlA01Q-y2zAOgn4r-cP9hiIEGaaqq",
        "subject": genomic_del3_dup3_loc_not_normalized,
        "copies": {"type": "Number", "value": 2},
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_dup3_cx(genomic_dup3, genomic_del3_dup3_loc_not_normalized):
    """Create a test fixture for genomic dup copy number change."""
    params = deepcopy(genomic_dup3)
    params["variation"] = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.cmQ2h4FXbNHHTZHJ8N0sOXIE01RfyBxb",
        "subject": genomic_del3_dup3_loc_not_normalized,
        "copy_change": "efo:0030070",
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_dup3_rse_lse(genomic_dup3):
    """Create test fixture for genomic dup rse and lse."""
    params = {
        "id": genomic_dup3["id"],
        "type": genomic_dup3["type"],
        "variation": {
            "id": "ga4gh:VT.15sKDgSyoCPOgfrFHvSea-fHVeu7huVT",
            "type": "Text",
            "definition": "NC_000023.11:g.(31060227_31100351)_(33274278_33417151)dup",
        },
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_dup3_free_text():
    """Create test fixture containing params for genomic dup VD."""
    params = {
        "id": "normalize.variation:DMD%20g.%2831147274_31147278%29_%2831182737_31182739%29dup",  # noqa: E501
        "type": "VariationDescriptor",
        "variation": dict(),
        "molecule_context": "genomic",
        "vrs_ref_allele_seq": None,
        "gene_context": "hgnc:2928",
    }
    return params


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
def genomic_dup3_free_text_cx(genomic_dup3_free_text, genomic_dup3_free_text_subject):
    """Create a test fixture for genomic dup copy number change."""
    params = deepcopy(genomic_dup3_free_text)
    params["variation"] = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.C0njAM6WcOua46XljTLLCVQkNVIChQU0",
        "subject": genomic_dup3_free_text_subject,
        "copy_change": "efo:0030070",
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_dup3_free_text_cn(genomic_dup3_free_text, genomic_dup3_free_text_subject):
    """Create a test fixture for genomic dup copy number count."""
    params = deepcopy(genomic_dup3_free_text)
    params["variation"] = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.yfsv4o1sfFpFUxKTS1de9ODUN8fSWMQh",
        "subject": genomic_dup3_free_text_subject,
        "copies": {"type": "Number", "value": 4},
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_dup3_free_text_rse_lse(genomic_dup3_free_text):
    """Create test fixture for genomic dup rse and lse."""
    params = {
        "id": genomic_dup3_free_text["id"],
        "type": genomic_dup3_free_text["type"],
        "variation": {
            "id": "ga4gh:VT.F0AX-RkMN4U8KLkIE68ECGU83Y-ICWXh",
            "type": "Text",
            "definition": "DMD g.(31147274_31147278)_(31182737_31182739)dup",
        },
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_dup4():
    """Create test fixture containing params for genomic dup VD."""
    params = {
        "id": "normalize.variation:NC_000020.11%3Ag.%28%3F_30417576%29_%2831394018_%3F%29dup",  # noqa: E501
        "type": "VariationDescriptor",
        "variation": dict(),
        "molecule_context": "genomic",
        "vrs_ref_allele_seq": None,
    }
    return params


@pytest.fixture(scope="module")
def genomic_dup4_cx(genomic_dup4, genomic_dup4_loc):
    """Create a test fixture for genomic dup copy number change."""
    params = deepcopy(genomic_dup4)
    params["variation"] = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.b7MlHXZcBEe4M81yZviR-Fujo4OKrcw3",
        "subject": genomic_dup4_loc,
        "copy_change": "efo:0030070",
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_dup4_cn(genomic_dup4, genomic_dup4_loc):
    """Create a test fixture for genomic dup copy number count."""
    params = deepcopy(genomic_dup4)
    params["variation"] = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.Fvt2VSCSLNROTjZd20vKZUzdMYxz-cfZ",
        "subject": genomic_dup4_loc,
        "copies": {"type": "Number", "value": 3},
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_dup4_rse_lse(genomic_dup4):
    """Create test fixture for genomic dup rse and lse."""
    params = {
        "id": genomic_dup4["id"],
        "type": genomic_dup4["type"],
        "variation": {
            "id": "ga4gh:VT.Pga4IH82qga2iZAodjxYw9OXhB4Xa2g8",
            "type": "Text",
            "definition": "NC_000020.11:g.(?_30417576)_(31394018_?)dup",
        },
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_dup4_free_text():
    """Create test fixture containing params for genomic dup VD."""
    params = {
        "id": "normalize.variation:PRF8%20g.%28%3F_1577736%29_%281587865_%3F%29",
        "type": "VariationDescriptor",
        "variation": dict(),
        "molecule_context": "genomic",
        "vrs_ref_allele_seq": None,
        "gene_context": "hgnc:17340",
    }
    return params


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
def genomic_dup4_free_text_cx(genomic_dup4_free_text, genomic_dup4_free_text_subject):
    """Create a test fixture for genomic dup copy number change."""
    params = deepcopy(genomic_dup4_free_text)
    params["variation"] = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.vhV4v1BK5hFQXHX9qysZFpIYhrgI5-Nf",
        "subject": genomic_dup4_free_text_subject,
        "copy_change": "efo:0030070",
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_dup4_free_text_cn(genomic_dup4_free_text, genomic_dup4_free_text_subject):
    """Create a test fixture for genomic dup copy number count."""
    params = deepcopy(genomic_dup4_free_text)
    params["variation"] = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.oFSw_z1nfk0Yr1mCIi3uKIF-KGC-x8lQ",
        "subject": genomic_dup4_free_text_subject,
        "copies": {"type": "Number", "value": 3},
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_dup4_free_text_rse_lse(genomic_dup4_free_text):
    """Create test fixture for genomic dup rse and lse."""
    params = {
        "id": genomic_dup4_free_text["id"],
        "type": genomic_dup4_free_text["type"],
        "variation": {
            "id": "ga4gh:VT.Pga4IH82qga2iZAodjxYw9OXhB4Xa2g8",
            "type": "Text",
            "definition": "NC_000020.11:g.(?_30417576)_(31394018_?)dup",
        },
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_dup5():
    """Create test fixture containing params for genomic dup VD."""
    params = {
        "id": "normalize.variation:NC_000023.11%3Ag.%28%3F_154021812%29_154092209dup",
        "type": "VariationDescriptor",
        "variation": dict(),
        "molecule_context": "genomic",
        "vrs_ref_allele_seq": None,
    }
    return params


def genomic_dup5_cn_var(params, genomic_dup5_loc):
    """Create genomic dup5 copy number count"""
    params["variation"] = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.uhQfl3d-UonMZDnK6FAASUtY8FTaHplR",
        "subject": genomic_dup5_loc,
        "copies": {"type": "Number", "value": 3},
    }


def genomic_dup5_cx_var(params, genomic_dup5_loc):
    """Create genomic dup4 copy number change"""
    params["variation"] = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.rBUl0BJiar7-NX0pU23NEln2KAifk03p",
        "subject": genomic_dup5_loc,
        "copy_change": "efo:0030070",
    }


@pytest.fixture(scope="module")
def genomic_dup5_cx(genomic_dup5, genomic_dup5_loc):
    """Create a test fixture for genomic dup5 copy number change."""
    genomic_dup5_cx_var(genomic_dup5, genomic_dup5_loc)
    return VariationDescriptor(**genomic_dup5)


@pytest.fixture(scope="module")
def genomic_dup5_cn(genomic_dup5, genomic_dup5_loc):
    """Create a test fixture for genomic dup5 copy number count."""
    genomic_dup5_cn_var(genomic_dup5, genomic_dup5_loc)
    return VariationDescriptor(**genomic_dup5)


@pytest.fixture(scope="module")
def genomic_dup5_rse_lse(genomic_dup5):
    """Create test fixture for genomic dup rse and lse."""
    params = {
        "id": genomic_dup5["id"],
        "type": genomic_dup5["type"],
        "variation": {
            "id": "ga4gh:VT.of16BEeHyU9od62SrjSCQ4LyUtbbGoKi",
            "type": "Text",
            "definition": "NC_000023.11:g.(?_154021812)_154092209dup",
        },
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_dup5_free_text():
    """Create test fixture containing params for genomic dup VD."""
    params = {
        "id": "normalize.variation:MECP2%20g.%28%3F_154021812%29_154092209dup",
        "type": "VariationDescriptor",
        "variation": dict(),
        "molecule_context": "genomic",
        "vrs_ref_allele_seq": None,
        "gene_context": "hgnc:6990",
    }
    return params


@pytest.fixture(scope="module")
def genomic_dup5_free_text_cx(genomic_dup5_free_text, genomic_dup5_loc):
    """Create a test fixture for genomic dup copy number change."""
    genomic_dup5_cx_var(genomic_dup5_free_text, genomic_dup5_loc)
    return VariationDescriptor(**genomic_dup5_free_text)


@pytest.fixture(scope="module")
def genomic_dup5_free_text_cn(genomic_dup5_free_text, genomic_dup5_loc):
    """Create a test fixture for genomic dup copy number count."""
    genomic_dup5_cn_var(genomic_dup5_free_text, genomic_dup5_loc)
    return VariationDescriptor(**genomic_dup5_free_text)


@pytest.fixture(scope="module")
def genomic_dup5_free_text_rse_lse(genomic_dup5_free_text):
    """Create test fixture for genomic dup rse and lse."""
    params = {
        "id": genomic_dup5_free_text["id"],
        "type": genomic_dup5_free_text["type"],
        "variation": {
            "id": "ga4gh:VT.Kw18bSFpQp9xdKg88fqW7zUx4_VXFIiW",
            "type": "Text",
            "definition": "MECP2 g.(?_154021812)_154092209dup",
        },
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_dup6():
    """Create test fixture containing params for genomic dup VD."""
    params = {
        "id": "normalize.variation:NC_000023.11%3Ag.154021812_%28154092209_%3F%29dup",
        "type": "VariationDescriptor",
        "variation": dict(),
        "molecule_context": "genomic",
        "vrs_ref_allele_seq": None,
    }
    return params


def genomic_dup6_cx_var(params, genomic_dup6_loc):
    """Create genomic dup6 copy number change"""
    params["variation"] = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.IbkRumjnx-rxHbWONhJKgl-9zixCnArU",
        "subject": genomic_dup6_loc,
        "copy_change": "efo:0030070",
    }


def genomic_dup6_cn_var(params, genomic_dup6_loc):
    """Create genomic dup6 copy number count"""
    params["variation"] = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.Tuvj2C1xYVUF0UL3CSB1GZjCRC6K_eSt",
        "subject": genomic_dup6_loc,
        "copies": {"type": "Number", "value": 2},
    }


@pytest.fixture(scope="module")
def genomic_dup6_cx(genomic_dup6, genomic_dup6_loc):
    """Create a test fixture for genomic dup copy number change."""
    genomic_dup6_cx_var(genomic_dup6, genomic_dup6_loc)
    return VariationDescriptor(**genomic_dup6)


@pytest.fixture(scope="module")
def genomic_dup6_cn(genomic_dup6, genomic_dup6_loc):
    """Create a test fixture for genomic dup copy number count."""
    genomic_dup6_cn_var(genomic_dup6, genomic_dup6_loc)
    return VariationDescriptor(**genomic_dup6)


@pytest.fixture(scope="module")
def genomic_dup6_rse_lse(genomic_dup6):
    """Create test fixture for genomic dup rse and lse."""
    params = {
        "id": genomic_dup6["id"],
        "type": genomic_dup6["type"],
        "variation": {
            "id": "ga4gh:VT.2k5AWTbGJxvLVT6bUW0pUMq6XGAcEjXW",
            "type": "Text",
            "definition": "NC_000023.11:g.154021812_(154092209_?)dup",
        },
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_dup6_free_text():
    """Create test fixture containing params for genomic dup VD."""
    params = {
        "id": "normalize.variation:MECP2%20g.154021812_%28154092209_%3F%29dup",
        "type": "VariationDescriptor",
        "variation": dict(),
        "molecule_context": "genomic",
        "vrs_ref_allele_seq": None,
        "gene_context": "hgnc:6990",
    }
    return params


@pytest.fixture(scope="module")
def genomic_dup6_free_text_cx(genomic_dup6_free_text, genomic_dup6_loc):
    """Create a test fixture for genomic dup copy number change."""
    genomic_dup6_cx_var(genomic_dup6_free_text, genomic_dup6_loc)
    return VariationDescriptor(**genomic_dup6_free_text)


@pytest.fixture(scope="module")
def genomic_dup6_free_text_cn(genomic_dup6_free_text, genomic_dup6_loc):
    """Create a test fixture for genomic dup copy number count."""
    genomic_dup6_cn_var(genomic_dup6_free_text, genomic_dup6_loc)
    return VariationDescriptor(**genomic_dup6_free_text)


@pytest.fixture(scope="module")
def genomic_dup6_free_text_rse_lse(genomic_dup6_free_text):
    """Create test fixture for genomic dup rse and lse."""
    params = {
        "id": genomic_dup6_free_text["id"],
        "type": genomic_dup6_free_text["type"],
        "variation": {
            "id": "ga4gh:VT.LbAqiLmJs1t9-FgEKD0-KDKwzvM3AAlz",
            "type": "Text",
            "definition": "MECP2 g.154021812_(154092209_?)dup",
        },
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
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


@pytest.fixture(scope="module")
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


@pytest.fixture(scope="module")
def genomic_del1_cn(genomic_del1, genomic_del1_38_cn):
    """Create a test fixture for genomic del copy number count."""
    genomic_del1["variation"] = genomic_del1_38_cn
    return VariationDescriptor(**genomic_del1)


@pytest.fixture(scope="module")
def genomic_del1_cx(genomic_del1, genomic_del1_seq_loc):
    """Create a test fixture for genomic del copy number change."""
    params = deepcopy(genomic_del1)
    params["variation"] = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.p8k6YurGbJALwJRjvKfq156CRUetenB6",
        "subject": genomic_del1_seq_loc,
        "copy_change": "efo:0030064",
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_del1_rse(genomic_del1, genomic_del1_seq_loc):
    """Create a test fixture for genomic del RSE."""
    params = deepcopy(genomic_del1)
    params["variation"] = {
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
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_del1_free_text():
    """Create test fixture containing params for genomic del VD."""
    params = {
        "id": "normalize.variation:VHL%20g.10191495del",
        "type": "VariationDescriptor",
        "variation": dict(),
        "molecule_context": "transcript",
        "vrs_ref_allele_seq": "T",
        "gene_context": "hgnc:12687",
    }
    return params


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
def genomic_del1_free_text_lse(genomic_del1_free_text, genomic_del1_free_text_seq_loc):
    """Create a test fixture for genomic del LSE."""
    params = deepcopy(genomic_del1_free_text)
    params["variation"] = {
        "type": "Allele",
        "id": "ga4gh:VA.VbssoYULETSOoFI0FwjGypt2YinbTEOc",
        "location": genomic_del1_free_text_seq_loc,
        "state": {"type": "LiteralSequenceExpression", "sequence": ""},
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_del1_free_text_cn(genomic_del1_free_text, genomic_del1_free_text_seq_loc):
    """Create a test fixture for genomic del copy number count."""
    params = deepcopy(genomic_del1_free_text)
    params["variation"] = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.ftCNLk6wg-0szErbIRcEZLhF6aUXd-3D",
        "subject": genomic_del1_free_text_seq_loc,
        "copies": {"type": "Number", "value": 1},
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_del1_free_text_rse(genomic_del1_free_text, genomic_del1_free_text_seq_loc):
    """Create a test fixture for genomic del RSE."""
    params = deepcopy(genomic_del1_free_text)
    params["variation"] = {
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
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
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


@pytest.fixture(scope="module")
def genomic_del2_lse(genomic_del2, genomic_del2_seq_loc):
    """Create a test fixture for genomic del LSE."""
    params = deepcopy(genomic_del2)
    params["variation"] = {
        "type": "Allele",
        "id": "ga4gh:VA.UgJSDSWAaJFwhRm56LA0Gez47_PYqv0k",
        "location": genomic_del2_seq_loc,
        "state": {"type": "LiteralSequenceExpression", "sequence": ""},
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_del2_lse2(seq_loc_gt_100_bp, vrs_ref_allele_seq_gt_100_bp):
    """Create a test fixture for genomic del LSE where bp > 100."""
    params = {
        "id": "normalize.variation:NC_000023.11%3Ag.33211290_33211490del",
        "type": "VariationDescriptor",
        "variation": {
            "type": "Allele",
            "id": "ga4gh:VA.7YNhEgDjlx_wc3Ci0Dt2iNAtIn9En0BL",
            "location": seq_loc_gt_100_bp,
            "state": {
                "type": "LiteralSequenceExpression",
                "sequence": "",
            },
        },
        "molecule_context": "genomic",
        "vrs_ref_allele_seq": vrs_ref_allele_seq_gt_100_bp,
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_del2_cn(genomic_del2, genomic_del2_38_cn):
    """Create a test fixture for genomic del copy number count."""
    params = deepcopy(genomic_del2)
    params["variation"] = genomic_del2_38_cn
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_del2_cx(genomic_del2, genomic_del2_seq_loc):
    """Create a test fixture for genomic del copy number change."""
    params = deepcopy(genomic_del2)
    params["variation"] = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.sfbYfDVrJCocGMvSivZAppmAnZ21Sp3-",
        "subject": genomic_del2_seq_loc,
        "copy_change": "efo:0030069",
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_del2_rse(genomic_del2, genomic_del2_seq_loc):
    """Create a test fixture for genomic del RSE."""
    params = deepcopy(genomic_del2)
    params["variation"] = {
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
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_del2_free_text():
    """Create test fixture containing params for genomic del VD."""
    params = {
        "id": "normalize.variation:VHL%20g.10188279_10188297del",
        "type": "VariationDescriptor",
        "variation": dict(),
        "molecule_context": "transcript",
        "vrs_ref_allele_seq": "ATGTTGACGGACAGCCTAT",
        "gene_context": "hgnc:12687",
    }
    return params


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
def genomic_del2_free_text_default(
    genomic_del2_free_text, genomic_del2_free_text_seq_loc
):
    """Create a test fixture for genomic del default and LSE."""
    params = deepcopy(genomic_del2_free_text)
    params["variation"] = {
        "type": "Allele",
        "id": "ga4gh:VA.6wTYBh0btGq6SlXDu4V7iEK9UrehXS-6",
        "location": genomic_del2_free_text_seq_loc,
        "state": {"type": "LiteralSequenceExpression", "sequence": ""},
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_del2_free_text_cnv(genomic_del2_free_text, genomic_del2_free_text_seq_loc):
    """Create a test fixture for genomic del CNV."""
    params = deepcopy(genomic_del2_free_text)
    params["variation"] = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.4gVfUDeIuWIwLhBU34UrjZWoKbWdG3B7",
        "subject": genomic_del2_free_text_seq_loc,
        "copies": {"type": "Number", "value": 1},
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_del2_free_text_rse(genomic_del2_free_text, genomic_del2_free_text_seq_loc):
    """Create a test fixture for genomic del RSE."""
    params = deepcopy(genomic_del2_free_text)
    params["variation"] = {
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
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_del3():
    """Create test fixture containing params for genomic del VD."""
    params = {
        "id": "normalize.variation:NC_000023.11%3Ag.%2831060227_31100351%29_%2833274278_33417151%29del",  # noqa: E501
        "type": "VariationDescriptor",
        "variation": dict(),
        "molecule_context": "genomic",
        "vrs_ref_allele_seq": None,
    }
    return params


@pytest.fixture(scope="module")
def genomic_del3_cx(genomic_del3, genomic_del3_dup3_loc_not_normalized):
    """Create a test fixture for genomic del copy number change."""
    params = deepcopy(genomic_del3)
    params["variation"] = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.XOE1DDGAyMEUeknUAXNHFkt7FKsQBUfh",
        "subject": genomic_del3_dup3_loc_not_normalized,
        "copy_change": "efo:0030067",
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_del3_cn(genomic_del3, genomic_del3_dup3_loc_not_normalized):
    """Create a test fixture for genomic del copy number count."""
    params = deepcopy(genomic_del3)
    params["variation"] = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.a2-YlA01Q-y2zAOgn4r-cP9hiIEGaaqq",
        "subject": genomic_del3_dup3_loc_not_normalized,
        "copies": {"type": "Number", "value": 2},
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_del3_rse_lse(genomic_del3):
    """Create test fixture for genomic del rse and lse."""
    params = {
        "id": genomic_del3["id"],
        "type": genomic_del3["type"],
        "variation": {
            "id": "ga4gh:VT.tmA3mpMy9HKUweaB8aYsq6uuejEx9iK7",
            "type": "Text",
            "definition": "NC_000023.11:g.(31060227_31100351)_(33274278_33417151)del",
        },
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_del3_free_text():
    """Create test fixture containing params for genomic del VD."""
    params = {
        "id": "normalize.variation:EFNB1%20g.%2868839265_68839268%29_%2868841120_68841125%29del",  # noqa: E501
        "type": "VariationDescriptor",
        "variation": dict(),
        "molecule_context": "genomic",
        "vrs_ref_allele_seq": None,
        "gene_context": "hgnc:3226",
    }
    return params


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
def genomic_del3_free_text_cx(genomic_del3_free_text, genomic_del3_free_text_subject):
    """Create a test fixture for genomic del copy number change."""
    params = deepcopy(genomic_del3_free_text)
    params["variation"] = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.k7s-BiYoxXGlBaGXpqcUdimdsUygemat",
        "subject": genomic_del3_free_text_subject,
        "copy_change": "efo:0030067",
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_del3_free_text_cn(genomic_del3_free_text, genomic_del3_free_text_subject):
    """Create a test fixture for genomic del copy number count."""
    params = deepcopy(genomic_del3_free_text)
    params["variation"] = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.qWnA8lNUMVz8vNKikQ6dWFjdV5E8FJ3U",
        "subject": genomic_del3_free_text_subject,
        "copies": {"type": "Number", "value": 2},
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_del3_free_text_rse_lse(genomic_del3_free_text):
    """Create test fixture for genomic del rse and lse."""
    params = {
        "id": genomic_del3_free_text["id"],
        "type": genomic_del3_free_text["type"],
        "variation": {
            "id": "ga4gh:VT.9mGg0U_Z7NZCFV3jrLdGxSQU03g7z3Z1",
            "type": "Text",
            "definition": "EFNB1 g.(68839265_68839268)_(68841120_68841125)del",
        },
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_del4():
    """Create test fixture containing params for genomic del VD."""
    params = {
        "id": "normalize.variation:NC_000023.11%3Ag.%28%3F_31120496%29_%2833339477_%3F%29del",  # noqa: E501
        "type": "VariationDescriptor",
        "variation": dict(),
        "molecule_context": "genomic",
        "vrs_ref_allele_seq": None,
    }
    return params


@pytest.fixture(scope="module")
def genomic_del4_cx(genomic_del4, genomic_del4_seq_loc):
    """Create a test fixture for genomic del copy number change."""
    params = deepcopy(genomic_del4)
    params["variation"] = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.fYi0dG6Q8kACkyY-ICBzzvslv-ONWrPF",
        "subject": genomic_del4_seq_loc,
        "copy_change": "efo:0030067",
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_del4_cn(genomic_del4, genomic_del4_seq_loc):
    """Create a test fixture for genomic del copy number count."""
    params = deepcopy(genomic_del4)
    params["variation"] = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.vU4Q8TzrKFU3aQvSn2RiRS3ikh58gw_3",
        "subject": genomic_del4_seq_loc,
        "copies": {"type": "Number", "value": 1},
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_del4_rse_lse(genomic_del4):
    """Create test fixture for genomic dup rse and lse."""
    params = {
        "id": genomic_del4["id"],
        "type": genomic_del4["type"],
        "variation": {
            "id": "ga4gh:VT.whBY5P24WVxF1wneDcI8x8btqorJUWXQ",
            "type": "Text",
            "definition": "NC_000023.11:g.(?_31120496)_(33339477_?)del",
        },
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_del4_free_text():
    """Create test fixture containing params for genomic del VD."""
    params = {
        "id": "normalize.variation:COL4A4%20g.%28%3F_227022028%29_%28227025830_%3F%29del",  # noqa: E501
        "type": "VariationDescriptor",
        "variation": dict(),
        "molecule_context": "genomic",
        "vrs_ref_allele_seq": None,
        "gene_context": "hgnc:2206",
    }
    return params


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
def genomic_del4_free_text_cx(genomic_del4_free_text, genomic_del4_free_text_subject):
    """Create a test fixture for genomic del copy number change."""
    params = deepcopy(genomic_del4_free_text)
    params["variation"] = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.q5LCBF00ZpNiutyADSXuiM0-u1--cHSA",
        "subject": genomic_del4_free_text_subject,
        "copy_change": "efo:0030067",
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_del4_free_text_cn(genomic_del4_free_text, genomic_del4_free_text_subject):
    """Create a test fixture for genomic del copy number count."""
    params = deepcopy(genomic_del4_free_text)
    params["variation"] = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.FF3yHxoyJB0JUt9FEeYuKx5vVOUe-w39",
        "subject": genomic_del4_free_text_subject,
        "copies": {"type": "Number", "value": 1},
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_del4_free_text_rse_lse(genomic_del4_free_text):
    """Create test fixture for genomic dup rse and lse."""
    params = {
        "id": genomic_del4_free_text["id"],
        "type": genomic_del4_free_text["type"],
        "variation": {
            "id": "ga4gh:VT.lT0rFYhOGFLA9MYA8ypnCf5q-CkV8dJv",
            "type": "Text",
            "definition": "COL4A4 g.(?_227022028)_(227025830_?)del",
        },
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_uncertain_del_2():
    """Create a genomic uncertain deletion on chr 2 test fixture."""
    params = {
        "id": "normalize.variation:NC_000002.12%3Ag.%28%3F_110104900%29_%28110207160_%3F%29del",  # noqa: E501
        "type": "VariationDescriptor",
        "variation": {
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
        },
        "molecule_context": "genomic",
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_uncertain_del_y():
    """Create a genomic uncertain deletion on chr Y test fixture."""
    params = {
        "id": "normalize.variation:NC_000024.10%3Ag.%28%3F_14076802%29_%2857165209_%3F%29del",  # noqa: E501
        "type": "VariationDescriptor",
        "variation": {
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
        },
        "molecule_context": "genomic",
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_del5():
    """Create test fixture containing params for genomic del VD."""
    params = {
        "id": "normalize.variation:NC_000023.11%3Ag.%28%3F_18575354%29_18653629del",
        "type": "VariationDescriptor",
        "variation": dict(),
        "molecule_context": "genomic",
        "vrs_ref_allele_seq": None,
    }
    return params


def genomic_del5_cn_var(params, genomic_del5_seq_loc):
    """Create genomic del5 copy number count"""
    params["variation"] = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.ZKXhoSguNUA8MBU9SpXmS52FmnZMhgy3",
        "subject": genomic_del5_seq_loc,
        "copies": {"type": "Number", "value": 3},
    }


def genomic_del5_cx_var(params, genomic_del5_seq_loc):
    """Create genomic del5 copy number change"""
    params["variation"] = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.0TqknNhM5e5i-5tAKxGXv4yDJlaQgDRR",
        "subject": genomic_del5_seq_loc,
        "copy_change": "efo:0030067",
    }


@pytest.fixture(scope="module")
def genomic_del5_cx(genomic_del5, genomic_del5_seq_loc):
    """Create a test fixture for genomic del copy number change."""
    genomic_del5_cx_var(genomic_del5, genomic_del5_seq_loc)
    return VariationDescriptor(**genomic_del5)


@pytest.fixture(scope="module")
def genomic_del5_cn(genomic_del5, genomic_del5_seq_loc):
    """Create a test fixture for genomic del copy number count."""
    genomic_del5_cn_var(genomic_del5, genomic_del5_seq_loc)
    return VariationDescriptor(**genomic_del5)


@pytest.fixture(scope="module")
def genomic_del5_rse_lse(genomic_del5):
    """Create test fixture for genomic del rse and lse."""
    params = {
        "id": genomic_del5["id"],
        "type": genomic_del5["type"],
        "variation": {
            "id": "ga4gh:VT.xCLHh3GpCebrP6KDMsWZRdIiW7Sti27H",
            "type": "Text",
            "definition": "NC_000023.11:g.(?_18575354)_18653629del",
        },
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_del5_free_text():
    """Create test fixture containing params for genomic del VD."""
    params = {
        "id": "normalize.variation:CDKL5%20g.%28%3F_18575354%29_18653629del",
        "type": "VariationDescriptor",
        "variation": dict(),
        "molecule_context": "genomic",
        "vrs_ref_allele_seq": None,
        "gene_context": "hgnc:11411",
    }
    return params


@pytest.fixture(scope="module")
def genomic_del5_free_text_cx(genomic_del5_free_text, genomic_del5_seq_loc):
    """Create a test fixture for genomic del copy number change."""
    genomic_del5_cx_var(genomic_del5_free_text, genomic_del5_seq_loc)
    return VariationDescriptor(**genomic_del5_free_text)


@pytest.fixture(scope="module")
def genomic_del5_free_text_cn(genomic_del5_free_text, genomic_del5_seq_loc):
    """Create a test fixture for genomic del copy number count."""
    genomic_del5_cn_var(genomic_del5_free_text, genomic_del5_seq_loc)
    return VariationDescriptor(**genomic_del5_free_text)


@pytest.fixture(scope="module")
def genomic_del5_free_text_rse_lse(genomic_del5_free_text):
    """Create test fixture for genomic del rse and lse."""
    params = {
        "id": genomic_del5_free_text["id"],
        "type": genomic_del5_free_text["type"],
        "variation": {
            "id": "ga4gh:VT.xCLHh3GpCebrP6KDMsWZRdIiW7Sti27H",
            "type": "Text",
            "definition": "NC_000023.11:g.(?_18575354)_18653629del",
        },
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_del6():
    """Create test fixture containing params for genomic del VD."""
    params = {
        "id": "normalize.variation:NC_000006.12%3Ag.133462764_%28133464858_%3F%29del",
        "type": "VariationDescriptor",
        "variation": dict(),
        "molecule_context": "genomic",
        "vrs_ref_allele_seq": None,
    }
    return params


def genomic_del6_cx_var(params, genomic_del6_seq_loc):
    """Create genomic del6 copy number change"""
    params["variation"] = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.kzLS3Xs3qhMarky22So-4JW_mZAIhVW2",
        "subject": genomic_del6_seq_loc,
        "copy_change": "efo:0030067",
    }


def genomic_del6_cn_var(params, genomic_del6_seq_loc):
    """Create genomic del6 copy number count"""
    params["variation"] = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.lfI0nj_cOjRwul_R_4tohZK_vqiZ-LSz",
        "subject": genomic_del6_seq_loc,
        "copies": {"type": "Number", "value": 1},
    }


@pytest.fixture(scope="module")
def genomic_del6_cx(genomic_del6, genomic_del6_seq_loc):
    """Create a test fixture for genomic del copy number change."""
    genomic_del6_cx_var(genomic_del6, genomic_del6_seq_loc)
    return VariationDescriptor(**genomic_del6)


@pytest.fixture(scope="module")
def genomic_del6_cn(genomic_del6, genomic_del6_seq_loc):
    """Create a test fixture for genomic del copy number count."""
    genomic_del6_cn_var(genomic_del6, genomic_del6_seq_loc)
    return VariationDescriptor(**genomic_del6)


@pytest.fixture(scope="module")
def genomic_del6_rse_lse(genomic_del6):
    """Create test fixture for genomic del rse and lse."""
    params = {
        "id": genomic_del6["id"],
        "type": genomic_del6["type"],
        "variation": {
            "id": "ga4gh:VT.Df49jbB-kZ2LSm180uA9wn4TT_p215yX",
            "type": "Text",
            "definition": "NC_000006.12:g.133462764_(133464858_?)del",
        },
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_del6_free_text():
    """Create test fixture containing params for genomic del VD."""
    params = {
        "id": "normalize.variation:EYA4%20g.133462764_%28133464858_%3F%29del",
        "type": "VariationDescriptor",
        "variation": dict(),
        "molecule_context": "genomic",
        "vrs_ref_allele_seq": None,
        "gene_context": "hgnc:3522",
    }
    return params


@pytest.fixture(scope="module")
def genomic_del6_free_text_cx(genomic_del6_free_text, genomic_del6_seq_loc):
    """Create a test fixture for genomic del copy number change."""
    genomic_del6_cx_var(genomic_del6_free_text, genomic_del6_seq_loc)
    return VariationDescriptor(**genomic_del6_free_text)


@pytest.fixture(scope="module")
def genomic_del6_free_text_cn(genomic_del6_free_text, genomic_del6_seq_loc):
    """Create a test fixture for genomic del copy number count."""
    genomic_del6_cn_var(genomic_del6_free_text, genomic_del6_seq_loc)
    return VariationDescriptor(**genomic_del6_free_text)


@pytest.fixture(scope="module")
def genomic_del6_free_text_rse_lse(genomic_del6_free_text):
    """Create test fixture for genomic del rse and lse."""
    params = {
        "id": genomic_del6_free_text["id"],
        "type": genomic_del6_free_text["type"],
        "variation": {
            "id": "ga4gh:VT.a3kXhodtO3tgsdPlEL39Ql4jOuCpOc0s",
            "type": "Text",
            "definition": "EYA4 g.133462764_(133464858_?)del",
        },
    }
    return VariationDescriptor(**params)


@pytest.mark.asyncio
async def assert_text_variation(query_list, test_handler):
    """Make assertion checks for invalid queries"""
    for q in query_list:
        resp = await test_handler.normalize(
            q, HGVSDupDelModeOption.DEFAULT, untranslatable_returns_text=True
        )
        assert resp.variation_descriptor.label == q.strip(), q
        assert resp.variation_descriptor.variation.type == "Text", q


@pytest.mark.asyncio
async def test_genomic_dup1(
    test_handler,
    genomic_dup1_lse,
    genomic_dup1_cn,
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
    assertion_checks(resp.variation_descriptor, genomic_dup1_lse, q)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp.variation_descriptor, genomic_dup1_lse, q)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
    )
    assertion_checks(resp.variation_descriptor, genomic_dup1_cn, q)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE, copy_change="efo:0030072"
    )
    assertion_checks(resp.variation_descriptor, genomic_dup1_cx, q)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.REPEATED_SEQ_EXPR)
    assertion_checks(resp.variation_descriptor, genomic_dup1_rse, q)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.LITERAL_SEQ_EXPR)
    assertion_checks(resp.variation_descriptor, genomic_dup1_lse, q)

    q = "NC_000003.11:g.49568695dup"  # 37
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp.variation_descriptor, genomic_dup1_lse, q, ignore_id=True)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
    )
    assertion_checks(resp.variation_descriptor, genomic_dup1_cn, q, ignore_id=True)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE, copy_change="efo:0030072"
    )
    assertion_checks(resp.variation_descriptor, genomic_dup1_cx, q, ignore_id=True)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.REPEATED_SEQ_EXPR)
    assertion_checks(resp.variation_descriptor, genomic_dup1_rse, q, ignore_id=True)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.LITERAL_SEQ_EXPR)
    assertion_checks(resp.variation_descriptor, genomic_dup1_lse, q, ignore_id=True)

    # Free Text
    for q in ["DAG1 g.49568695dup", "DAG1 g.49531262dup"]:  # 37  # 38
        resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
        assertion_checks(
            resp.variation_descriptor, genomic_dup1_free_text_lse, q, ignore_id=True
        )

        resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
        assertion_checks(
            resp.variation_descriptor, genomic_dup1_free_text_lse, q, ignore_id=True
        )

        resp = await test_handler.normalize(
            q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
        )
        assertion_checks(
            resp.variation_descriptor, genomic_dup1_free_text_cn, q, ignore_id=True
        )

        resp = await test_handler.normalize(q, HGVSDupDelModeOption.REPEATED_SEQ_EXPR)
        assertion_checks(
            resp.variation_descriptor, genomic_dup1_free_text_rse, q, ignore_id=True
        )

        resp = await test_handler.normalize(q, HGVSDupDelModeOption.LITERAL_SEQ_EXPR)
        assertion_checks(
            resp.variation_descriptor, genomic_dup1_free_text_lse, q, ignore_id=True
        )

    # Invalid
    invalid_queries = [
        "NC_000007.13:g.159138670dup",
        "NC_000007.14:g.159345976dup",
        "BRAF g.140219337dup",
        "BRAF g.141024929dup",
    ]
    await assert_text_variation(invalid_queries, test_handler)


@pytest.mark.asyncio
async def test_genomic_dup2(
    test_handler,
    genomic_dup2_lse,
    genomic_dup2_cn,
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
    assertion_checks(resp.variation_descriptor, genomic_dup2_lse, q)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
    )
    assertion_checks(resp.variation_descriptor, genomic_dup2_cn, q)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE)
    assertion_checks(resp.variation_descriptor, genomic_dup2_cx, q)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.REPEATED_SEQ_EXPR)
    assertion_checks(resp.variation_descriptor, genomic_dup2_rse, q)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.LITERAL_SEQ_EXPR)
    assertion_checks(resp.variation_descriptor, genomic_dup2_lse, q)

    q = "NC_000023.10:g.33229407_33229410dup"  # 37
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp.variation_descriptor, genomic_dup2_lse, q, ignore_id=True)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
    )
    assertion_checks(resp.variation_descriptor, genomic_dup2_cn, q, ignore_id=True)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE)
    assertion_checks(resp.variation_descriptor, genomic_dup2_cx, q, ignore_id=True)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.REPEATED_SEQ_EXPR)
    assertion_checks(resp.variation_descriptor, genomic_dup2_rse, q, ignore_id=True)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.LITERAL_SEQ_EXPR)
    assertion_checks(resp.variation_descriptor, genomic_dup2_lse, q, ignore_id=True)

    # Free text
    for q in ["DMD g.33211290_33211293dup", "DMD g.33229407_33229410dup"]:  # 37  # 38
        resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
        assertion_checks(
            resp.variation_descriptor, genomic_dup2_free_text_default, q, ignore_id=True
        )

        resp = await test_handler.normalize(
            q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
        )
        assertion_checks(
            resp.variation_descriptor, genomic_dup2_free_text_cn, q, ignore_id=True
        )

        resp = await test_handler.normalize(q, HGVSDupDelModeOption.REPEATED_SEQ_EXPR)
        assertion_checks(
            resp.variation_descriptor, genomic_dup2_free_text_rse, q, ignore_id=True
        )

        resp = await test_handler.normalize(q, HGVSDupDelModeOption.LITERAL_SEQ_EXPR)
        assertion_checks(
            resp.variation_descriptor, genomic_dup2_free_text_default, q, ignore_id=True
        )

    # Greater than 100 bps -> rse
    q = "NC_000023.11:g.33211290_33211490dup"
    resp = await test_handler.normalize(q)
    assertion_checks(resp.variation_descriptor, genomic_dup2_rse2, q)

    # Invalid
    invalid_queries = [
        "NC_000007.13:g.140413127_159138670dup",
        "NC_000007.14:g.140413127_159345976dup",
        "BRAF g.140219337_140924929dup",
        "BRAF g.140719326_141024929dup",
    ]
    await assert_text_variation(invalid_queries, test_handler)


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
    assertion_checks(resp.variation_descriptor, genomic_dup3_cx, q)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=1
    )
    assertion_checks(resp.variation_descriptor, genomic_dup3_cn, q)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE, copy_change="efo:0030070"
    )
    assertion_checks(resp.variation_descriptor, genomic_dup3_cx, q)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.REPEATED_SEQ_EXPR, untranslatable_returns_text=True
    )
    assertion_checks(resp.variation_descriptor, genomic_dup3_rse_lse, q)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.LITERAL_SEQ_EXPR, untranslatable_returns_text=True
    )
    assertion_checks(resp.variation_descriptor, genomic_dup3_rse_lse, q)

    q = "NC_000023.10:g.(31078344_31118468)_(33292395_33435268)dup"  # 37
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp.variation_descriptor, genomic_dup3_cx, q, ignore_id=True)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=1
    )
    assertion_checks(resp.variation_descriptor, genomic_dup3_cn, q, ignore_id=True)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE)
    assertion_checks(resp.variation_descriptor, genomic_dup3_cx, q, ignore_id=True)

    genomic_dup3_rse_lse.variation.definition = q
    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.REPEATED_SEQ_EXPR, untranslatable_returns_text=True
    )
    assertion_checks(resp.variation_descriptor, genomic_dup3_rse_lse, q, ignore_id=True)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.LITERAL_SEQ_EXPR, untranslatable_returns_text=True
    )
    assertion_checks(resp.variation_descriptor, genomic_dup3_rse_lse, q, ignore_id=True)

    # Free Text
    for q in ["DMD g.(31147274_31147278)_(31182737_31182739)dup"]:  # 38
        resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
        assertion_checks(
            resp.variation_descriptor, genomic_dup3_free_text_cx, q, ignore_id=True
        )

        resp = await test_handler.normalize(
            q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=3
        )
        assertion_checks(
            resp.variation_descriptor, genomic_dup3_free_text_cn, q, ignore_id=True
        )

        genomic_dup3_rse_lse.variation.definition = q
        resp = await test_handler.normalize(
            q, HGVSDupDelModeOption.REPEATED_SEQ_EXPR, untranslatable_returns_text=True
        )
        assertion_checks(
            resp.variation_descriptor, genomic_dup3_free_text_rse_lse, q, ignore_id=True
        )

        resp = await test_handler.normalize(
            q, HGVSDupDelModeOption.LITERAL_SEQ_EXPR, untranslatable_returns_text=True
        )
        assertion_checks(
            resp.variation_descriptor, genomic_dup3_free_text_rse_lse, q, ignore_id=True
        )

    # Invalid
    invalid_queries = [
        "NC_000023.10:g.(31119221_31119227)_(31119300_155270562)dup",
        "NC_000023.11:g.(31119221_31119227)_(31119300_156040899)dup",
        "DMD g.(31060227_31100351)_(33274278_33417151)dup",
    ]
    await assert_text_variation(invalid_queries, test_handler)


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
    assertion_checks(resp.variation_descriptor, genomic_dup4_cx, q)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
    )
    assertion_checks(resp.variation_descriptor, genomic_dup4_cn, q)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE)
    assertion_checks(resp.variation_descriptor, genomic_dup4_cx, q)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.REPEATED_SEQ_EXPR, untranslatable_returns_text=True
    )
    assertion_checks(resp.variation_descriptor, genomic_dup4_rse_lse, q)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.LITERAL_SEQ_EXPR, untranslatable_returns_text=True
    )
    assertion_checks(resp.variation_descriptor, genomic_dup4_rse_lse, q)

    q = "NC_000020.10:g.(?_29652252)_(29981821_?)dup"  # 37
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp.variation_descriptor, genomic_dup4_cx, q, ignore_id=True)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
    )
    assertion_checks(resp.variation_descriptor, genomic_dup4_cn, q, ignore_id=True)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE)
    assertion_checks(resp.variation_descriptor, genomic_dup4_cx, q, ignore_id=True)

    genomic_dup4_rse_lse.variation.definition = q
    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.REPEATED_SEQ_EXPR, untranslatable_returns_text=True
    )
    assertion_checks(resp.variation_descriptor, genomic_dup4_rse_lse, q, ignore_id=True)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.LITERAL_SEQ_EXPR, untranslatable_returns_text=True
    )
    assertion_checks(resp.variation_descriptor, genomic_dup4_rse_lse, q, ignore_id=True)

    # Free Text
    for q in [
        "PRPF8 g.(?_1577736)_(1587865_?)dup",  # 37
        "PRPF8 g.(?_1674442)_(1684571_?)dup",  # 38
    ]:
        resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
        assertion_checks(
            resp.variation_descriptor, genomic_dup4_free_text_cx, q, ignore_id=True
        )

        resp = await test_handler.normalize(
            q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
        )
        assertion_checks(
            resp.variation_descriptor, genomic_dup4_free_text_cn, q, ignore_id=True
        )

        genomic_dup4_rse_lse.variation.definition = q
        resp = await test_handler.normalize(
            q, HGVSDupDelModeOption.REPEATED_SEQ_EXPR, untranslatable_returns_text=True
        )
        genomic_dup4_free_text_rse_lse.variation.definition = q
        assertion_checks(
            resp.variation_descriptor, genomic_dup4_free_text_rse_lse, q, ignore_id=True
        )

        resp = await test_handler.normalize(
            q, HGVSDupDelModeOption.LITERAL_SEQ_EXPR, untranslatable_returns_text=True
        )
        assertion_checks(
            resp.variation_descriptor, genomic_dup4_free_text_rse_lse, q, ignore_id=True
        )

    # Invalid
    invalid_queries = [
        "NC_000020.10:g.(?_29652252)_(63025530_?)dup",
        "NC_000020.11:g.(?_29652252)_(64444169_?)dup",
        "PRPF8 g.(?_1650628)_(1684571_?)dup",
    ]
    await assert_text_variation(invalid_queries, test_handler)


@pytest.mark.asyncio
async def test_genomic_dup5(
    test_handler,
    genomic_dup5_cn,
    genomic_dup5_cx,
    genomic_dup5_rse_lse,
    genomic_dup5_free_text_cn,
    genomic_dup5_free_text_cx,
    genomic_dup5_free_text_rse_lse,
):
    """Test that genomic duplication works correctly."""
    q = "NC_000023.11:g.(?_154021812)_154092209dup"  # 38
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp.variation_descriptor, genomic_dup5_cx, q)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
    )
    assertion_checks(resp.variation_descriptor, genomic_dup5_cn, q)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE)
    assertion_checks(resp.variation_descriptor, genomic_dup5_cx, q)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.REPEATED_SEQ_EXPR, untranslatable_returns_text=True
    )
    assertion_checks(resp.variation_descriptor, genomic_dup5_rse_lse, q)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.LITERAL_SEQ_EXPR, untranslatable_returns_text=True
    )
    assertion_checks(resp.variation_descriptor, genomic_dup5_rse_lse, q)

    q = "NC_000023.10:g.(?_153287263)_153357667dup"  # 37
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp.variation_descriptor, genomic_dup5_cx, q, ignore_id=True)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
    )
    assertion_checks(resp.variation_descriptor, genomic_dup5_cn, q, ignore_id=True)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE)
    assertion_checks(resp.variation_descriptor, genomic_dup5_cx, q, ignore_id=True)

    genomic_dup5_rse_lse.variation.definition = q
    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.REPEATED_SEQ_EXPR, untranslatable_returns_text=True
    )
    assertion_checks(resp.variation_descriptor, genomic_dup5_rse_lse, q, ignore_id=True)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.LITERAL_SEQ_EXPR, untranslatable_returns_text=True
    )
    assertion_checks(resp.variation_descriptor, genomic_dup5_rse_lse, q, ignore_id=True)

    # Free Text
    for q in [
        "MECP2 g.(?_153287263)_153357667dup",  # 37
        "MECP2 g.(?_154021812)_154092209dup",  # 38
    ]:
        resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
        assertion_checks(
            resp.variation_descriptor, genomic_dup5_free_text_cx, q, ignore_id=True
        )

        resp = await test_handler.normalize(
            q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
        )
        assertion_checks(
            resp.variation_descriptor, genomic_dup5_free_text_cn, q, ignore_id=True
        )

        genomic_dup5_free_text_rse_lse.variation.definition = q
        resp = await test_handler.normalize(
            q, HGVSDupDelModeOption.REPEATED_SEQ_EXPR, untranslatable_returns_text=True
        )
        assertion_checks(
            resp.variation_descriptor, genomic_dup5_free_text_rse_lse, q, ignore_id=True
        )

        resp = await test_handler.normalize(
            q, HGVSDupDelModeOption.LITERAL_SEQ_EXPR, untranslatable_returns_text=True
        )
        assertion_checks(
            resp.variation_descriptor, genomic_dup5_free_text_rse_lse, q, ignore_id=True
        )

    # Invalid
    for q in [
        "NC_000023.10:g.(?_153287263)_155270561dup",
        "NC_000023.11:g.(?_154021812)_156040896dup",
        "MECP2 g.(?_154021812)_154097733dup"  # 37
        "MECP2 g.(?_154021572)_154092209dup",  # 38
    ]:
        resp = await test_handler.normalize(
            q, HGVSDupDelModeOption.DEFAULT, untranslatable_returns_text=True
        )
        assert resp.variation_descriptor.variation.type == "Text", q


@pytest.mark.asyncio
async def test_genomic_dup6(
    test_handler,
    genomic_dup6_cn,
    genomic_dup6_cx,
    genomic_dup6_rse_lse,
    genomic_dup6_free_text_cn,
    genomic_dup6_free_text_cx,
    genomic_dup6_free_text_rse_lse,
):
    """Test that genomic duplication works correctly."""
    q = "NC_000023.11:g.154021812_(154092209_?)dup"  # 38
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp.variation_descriptor, genomic_dup6_cx, q)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=1
    )
    assertion_checks(resp.variation_descriptor, genomic_dup6_cn, q)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE)
    assertion_checks(resp.variation_descriptor, genomic_dup6_cx, q)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.REPEATED_SEQ_EXPR, untranslatable_returns_text=True
    )
    assertion_checks(resp.variation_descriptor, genomic_dup6_rse_lse, q)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.LITERAL_SEQ_EXPR, untranslatable_returns_text=True
    )
    assertion_checks(resp.variation_descriptor, genomic_dup6_rse_lse, q)

    q = "NC_000023.10:g.153287263_(153357667_?)dup"  # 37
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp.variation_descriptor, genomic_dup6_cx, q, ignore_id=True)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=1
    )
    assertion_checks(resp.variation_descriptor, genomic_dup6_cn, q, ignore_id=True)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE)
    assertion_checks(resp.variation_descriptor, genomic_dup6_cx, q, ignore_id=True)

    genomic_dup6_rse_lse.variation.definition = q
    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.REPEATED_SEQ_EXPR, untranslatable_returns_text=True
    )
    assertion_checks(resp.variation_descriptor, genomic_dup6_rse_lse, q, ignore_id=True)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.LITERAL_SEQ_EXPR, untranslatable_returns_text=True
    )
    assertion_checks(resp.variation_descriptor, genomic_dup6_rse_lse, q, ignore_id=True)

    # Free Text
    for q in [
        "MECP2 g.153287263_(153357667_?)dup",  # 37
        "MECP2 g.154021812_(154092209_?)dup",  # 38
    ]:
        resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
        assertion_checks(
            resp.variation_descriptor, genomic_dup6_free_text_cx, q, ignore_id=True
        )

        resp = await test_handler.normalize(
            q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=1
        )
        assertion_checks(
            resp.variation_descriptor, genomic_dup6_free_text_cn, q, ignore_id=True
        )

        genomic_dup6_free_text_rse_lse.variation.definition = q
        resp = await test_handler.normalize(
            q, HGVSDupDelModeOption.REPEATED_SEQ_EXPR, untranslatable_returns_text=True
        )
        assertion_checks(
            resp.variation_descriptor, genomic_dup6_free_text_rse_lse, q, ignore_id=True
        )

        resp = await test_handler.normalize(
            q, HGVSDupDelModeOption.LITERAL_SEQ_EXPR, untranslatable_returns_text=True
        )
        assertion_checks(
            resp.variation_descriptor, genomic_dup6_free_text_rse_lse, q, ignore_id=True
        )

    # Invalid
    for q in [
        "NC_000023.10:g.153287263_(155270561_?)dup",
        "NC_000023.11:g.154021812_(156040896_?)dup",
        "MECP2 g.154021812_(154097733_?)dup"  # 37
        "MECP2 g.154021572_(154092209_?)dup",  # 38
    ]:
        resp = await test_handler.normalize(
            q, HGVSDupDelModeOption.DEFAULT, untranslatable_returns_text=True
        )
        assert resp.variation_descriptor.variation.type == "Text", q


@pytest.mark.asyncio
async def test_genomic_del1(
    test_handler,
    genomic_del1_lse,
    genomic_del1_cn,
    genomic_del1_cx,
    genomic_del1_rse,
    genomic_del1_free_text_lse,
    genomic_del1_free_text_cn,
    genomic_del1_free_text_rse,
):
    """Test that genomic deletion works correctly."""
    q = "NC_000003.12:g.10149811del"  # 38
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp.variation_descriptor, genomic_del1_lse, q)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
    )
    assertion_checks(resp.variation_descriptor, genomic_del1_cn, q)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE, copy_change="efo:0030064"
    )
    assertion_checks(resp.variation_descriptor, genomic_del1_cx, q)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.REPEATED_SEQ_EXPR)
    assertion_checks(resp.variation_descriptor, genomic_del1_rse, q)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.LITERAL_SEQ_EXPR)
    assertion_checks(resp.variation_descriptor, genomic_del1_lse, q)

    q = "NC_000003.11:g.10191495del"  # 37
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp.variation_descriptor, genomic_del1_lse, q, ignore_id=True)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
    )
    assertion_checks(resp.variation_descriptor, genomic_del1_cn, q, ignore_id=True)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE, copy_change="efo:0030064"
    )
    assertion_checks(resp.variation_descriptor, genomic_del1_cx, q, ignore_id=True)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.REPEATED_SEQ_EXPR)
    assertion_checks(resp.variation_descriptor, genomic_del1_rse, q, ignore_id=True)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.LITERAL_SEQ_EXPR)
    assertion_checks(resp.variation_descriptor, genomic_del1_lse, q, ignore_id=True)

    # Free text
    for q in ["VHL g.10191495del", "VHL g.10149811del"]:  # 37  # 38
        resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
        assertion_checks(
            resp.variation_descriptor, genomic_del1_free_text_lse, q, ignore_id=True
        )

        resp = await test_handler.normalize(
            q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
        )
        assertion_checks(
            resp.variation_descriptor, genomic_del1_free_text_cn, q, ignore_id=True
        )

        resp = await test_handler.normalize(q, HGVSDupDelModeOption.REPEATED_SEQ_EXPR)
        assertion_checks(
            resp.variation_descriptor, genomic_del1_free_text_rse, q, ignore_id=True
        )

        resp = await test_handler.normalize(q, HGVSDupDelModeOption.LITERAL_SEQ_EXPR)
        assertion_checks(
            resp.variation_descriptor, genomic_del1_free_text_lse, q, ignore_id=True
        )

    # Invalid
    invalid_queries = [
        "NC_000003.11:g.198022431del",
        "NC_000003.12:g.198295567del",
        "BRAF g.140413127del",
        "BRAF g.141024929del",
    ]
    await assert_text_variation(invalid_queries, test_handler)


@pytest.mark.asyncio
async def test_genomic_del2(
    test_handler,
    genomic_del2_lse,
    genomic_del2_cn,
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
    assertion_checks(resp.variation_descriptor, genomic_del2_lse, q)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
    )
    assertion_checks(resp.variation_descriptor, genomic_del2_cn, q)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE, copy_change="efo:0030069"
    )
    assertion_checks(resp.variation_descriptor, genomic_del2_cx, q)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.REPEATED_SEQ_EXPR)
    assertion_checks(resp.variation_descriptor, genomic_del2_rse, q)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.LITERAL_SEQ_EXPR)
    assertion_checks(resp.variation_descriptor, genomic_del2_lse, q)

    q = "NC_000003.11:g.10188279_10188297del"  # 37
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp.variation_descriptor, genomic_del2_lse, q, ignore_id=True)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
    )
    assertion_checks(resp.variation_descriptor, genomic_del2_cn, q, ignore_id=True)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE, copy_change="efo:0030069"
    )
    assertion_checks(resp.variation_descriptor, genomic_del2_cx, q, ignore_id=True)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.REPEATED_SEQ_EXPR)
    assertion_checks(resp.variation_descriptor, genomic_del2_rse, q, ignore_id=True)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.LITERAL_SEQ_EXPR)
    assertion_checks(resp.variation_descriptor, genomic_del2_lse, q, ignore_id=True)

    # Free text
    for q in ["VHL g.10188279_10188297del", "VHL g.10146595_10146613del"]:  # 37  # 38
        resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
        assertion_checks(
            resp.variation_descriptor, genomic_del2_free_text_default, q, ignore_id=True
        )

        resp = await test_handler.normalize(
            q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
        )
        assertion_checks(
            resp.variation_descriptor, genomic_del2_free_text_cnv, q, ignore_id=True
        )

        resp = await test_handler.normalize(q, HGVSDupDelModeOption.REPEATED_SEQ_EXPR)
        assertion_checks(
            resp.variation_descriptor, genomic_del2_free_text_rse, q, ignore_id=True
        )

        resp = await test_handler.normalize(q, HGVSDupDelModeOption.LITERAL_SEQ_EXPR)
        assertion_checks(
            resp.variation_descriptor, genomic_del2_free_text_default, q, ignore_id=True
        )

    # Check that del > 100 bps returns LSE
    q = "NC_000023.11:g.33211290_33211490del"
    resp = await test_handler.normalize(q)
    assertion_checks(resp.variation_descriptor, genomic_del2_lse2, q)

    # gnomad vcf
    q = "3-10146594-AATGTTGACGGACAGCCTAT-A"
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp.variation_descriptor, genomic_del2_lse, q, ignore_id=True)

    q = "3-10188278-AATGTTGACGGACAGCCTAT-A"
    resp = await test_handler.normalize(q)
    assertion_checks(resp.variation_descriptor, genomic_del2_lse, q, ignore_id=True)

    # Invalid
    invalid_queries = [
        "NC_000003.12:g.10146595_198295580del",
        "NC_000003.11:g.198022435_198022437del",
        "BRAF g.140413127_140419136del",
        "BRAF g.140719326_141024929del",
    ]
    await assert_text_variation(invalid_queries, test_handler)


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
    assertion_checks(resp.variation_descriptor, genomic_del3_cx, q)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=3
    )
    assertion_checks(resp.variation_descriptor, genomic_del3_cn, q)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE)
    assertion_checks(resp.variation_descriptor, genomic_del3_cx, q)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.REPEATED_SEQ_EXPR, untranslatable_returns_text=True
    )
    assertion_checks(resp.variation_descriptor, genomic_del3_rse_lse, q)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.LITERAL_SEQ_EXPR, untranslatable_returns_text=True
    )
    assertion_checks(resp.variation_descriptor, genomic_del3_rse_lse, q)

    q = "NC_000023.10:g.(31078344_31118468)_(33292395_33435268)del"  # 37
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp.variation_descriptor, genomic_del3_cx, q, ignore_id=True)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=3
    )
    assertion_checks(resp.variation_descriptor, genomic_del3_cn, q, ignore_id=True)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE)
    assertion_checks(resp.variation_descriptor, genomic_del3_cx, q, ignore_id=True)

    genomic_del3_rse_lse.variation.definition = q
    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.REPEATED_SEQ_EXPR, untranslatable_returns_text=True
    )
    assertion_checks(resp.variation_descriptor, genomic_del3_rse_lse, q, ignore_id=True)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.LITERAL_SEQ_EXPR, untranslatable_returns_text=True
    )
    assertion_checks(resp.variation_descriptor, genomic_del3_rse_lse, q, ignore_id=True)

    # Free Text
    for q in [
        "EFNB1 g.(68059108_68059111)_(68060963_68060968)del",  # 37
        "EFNB1 g.(68839265_68839268)_(68841120_68841125)del",  # 38
    ]:
        resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
        assertion_checks(
            resp.variation_descriptor, genomic_del3_free_text_cx, q, ignore_id=True
        )

        resp = await test_handler.normalize(
            q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=3
        )
        assertion_checks(
            resp.variation_descriptor, genomic_del3_free_text_cn, q, ignore_id=True
        )

        genomic_del3_free_text_rse_lse.variation.definition = q
        resp = await test_handler.normalize(
            q, HGVSDupDelModeOption.REPEATED_SEQ_EXPR, untranslatable_returns_text=True
        )
        assertion_checks(
            resp.variation_descriptor, genomic_del3_free_text_rse_lse, q, ignore_id=True
        )

        resp = await test_handler.normalize(
            q, HGVSDupDelModeOption.LITERAL_SEQ_EXPR, untranslatable_returns_text=True
        )
        assertion_checks(
            resp.variation_descriptor, genomic_del3_free_text_rse_lse, q, ignore_id=True
        )

    # Invalid
    invalid_queries = [
        "NC_000023.11:g.(156040880_156040883)_(156040896_156040899)del",
        "NC_000023.10:g.(155270550_155270555)_(155270560_155270562)del",
        "EFNB1 g.(68048863_68048870)_(68842150_68842152)del",  # 37
        "EFNB1 g.(68829022_68829030)_(68842150_68842161)del",  # 38
    ]
    await assert_text_variation(invalid_queries, test_handler)


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
    assertion_checks(resp.variation_descriptor, genomic_del4_cx, q)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
    )
    assertion_checks(resp.variation_descriptor, genomic_del4_cn, q)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE)
    assertion_checks(resp.variation_descriptor, genomic_del4_cx, q)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.REPEATED_SEQ_EXPR, untranslatable_returns_text=True
    )
    assertion_checks(resp.variation_descriptor, genomic_del4_rse_lse, q)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.LITERAL_SEQ_EXPR, untranslatable_returns_text=True
    )
    assertion_checks(resp.variation_descriptor, genomic_del4_rse_lse, q)

    q = "NC_000023.10:g.(?_31138613)_(33357594_?)del"  # 37
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp.variation_descriptor, genomic_del4_cx, q, ignore_id=True)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
    )
    assertion_checks(resp.variation_descriptor, genomic_del4_cn, q, ignore_id=True)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE)
    assertion_checks(resp.variation_descriptor, genomic_del4_cx, q, ignore_id=True)

    genomic_del4_rse_lse.variation.definition = q
    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.REPEATED_SEQ_EXPR, untranslatable_returns_text=True
    )
    assertion_checks(resp.variation_descriptor, genomic_del4_rse_lse, q, ignore_id=True)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.LITERAL_SEQ_EXPR, untranslatable_returns_text=True
    )
    assertion_checks(resp.variation_descriptor, genomic_del4_rse_lse, q, ignore_id=True)

    q = "NC_000002.12:g.(?_110104900)_(110207160_?)del"
    resp = await test_handler.normalize(q)
    assertion_checks(resp.variation_descriptor, genomic_uncertain_del_2, q)

    q = "NC_000024.10:g.(?_14076802)_(57165209_?)del"
    resp = await test_handler.normalize(q)
    assertion_checks(resp.variation_descriptor, genomic_uncertain_del_y, q)

    # Free Text
    for q in ["COL4A4 g.(?_227022028)_(227025830_?)del"]:  # 38
        resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
        assertion_checks(
            resp.variation_descriptor, genomic_del4_free_text_cx, q, ignore_id=True
        )

        resp = await test_handler.normalize(
            q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
        )
        assertion_checks(
            resp.variation_descriptor, genomic_del4_free_text_cn, q, ignore_id=True
        )

        resp = await test_handler.normalize(
            q, HGVSDupDelModeOption.REPEATED_SEQ_EXPR, untranslatable_returns_text=True
        )
        assertion_checks(
            resp.variation_descriptor, genomic_del4_free_text_rse_lse, q, ignore_id=True
        )

        resp = await test_handler.normalize(
            q, HGVSDupDelModeOption.LITERAL_SEQ_EXPR, untranslatable_returns_text=True
        )
        assertion_checks(
            resp.variation_descriptor, genomic_del4_free_text_rse_lse, q, ignore_id=True
        )

    # Invalid
    invalid_queries = [
        "NC_000023.11:g.(?_156040899)_(156040900_?)del",
        "NC_000024.10:g.(?_155270565)_(155270568_?)del",
        "COL4A4 g.(?_227002710)_(227003710_?)del",
        "COL4A4 g.(?_227867430)_(228029276_?)del",
    ]
    await assert_text_variation(invalid_queries, test_handler)


@pytest.mark.asyncio
async def test_genomic_del5(
    test_handler,
    genomic_del5_cn,
    genomic_del5_cx,
    genomic_del5_rse_lse,
    genomic_del5_free_text_cn,
    genomic_del5_free_text_cx,
    genomic_del5_free_text_rse_lse,
):
    """Test that genomic deletion works correctly."""
    q = "NC_000023.11:g.(?_18575354)_18653629del"  # 38
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp.variation_descriptor, genomic_del5_cx, q)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=4
    )
    assertion_checks(resp.variation_descriptor, genomic_del5_cn, q)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE)
    assertion_checks(resp.variation_descriptor, genomic_del5_cx, q)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.REPEATED_SEQ_EXPR, untranslatable_returns_text=True
    )
    assertion_checks(resp.variation_descriptor, genomic_del5_rse_lse, q)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.LITERAL_SEQ_EXPR, untranslatable_returns_text=True
    )
    assertion_checks(resp.variation_descriptor, genomic_del5_rse_lse, q)

    q = "NC_000023.10:g.(?_18593474)_18671749del"  # 37
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp.variation_descriptor, genomic_del5_cx, q, ignore_id=True)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=4
    )
    assertion_checks(resp.variation_descriptor, genomic_del5_cn, q, ignore_id=True)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE)
    assertion_checks(resp.variation_descriptor, genomic_del5_cx, q, ignore_id=True)

    genomic_del5_rse_lse.variation.definition = q
    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.REPEATED_SEQ_EXPR, untranslatable_returns_text=True
    )
    assertion_checks(resp.variation_descriptor, genomic_del5_rse_lse, q, ignore_id=True)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.LITERAL_SEQ_EXPR, untranslatable_returns_text=True
    )
    assertion_checks(resp.variation_descriptor, genomic_del5_rse_lse, q, ignore_id=True)

    # Free text
    for q in ["CDKL5 g.(?_18575354)_18653629del"]:
        resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
        assertion_checks(
            resp.variation_descriptor, genomic_del5_free_text_cx, q, ignore_id=True
        )

        resp = await test_handler.normalize(
            q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=4
        )
        assertion_checks(
            resp.variation_descriptor, genomic_del5_free_text_cn, q, ignore_id=True
        )

        genomic_del5_free_text_rse_lse.variation.definition = q
        resp = await test_handler.normalize(
            q, HGVSDupDelModeOption.REPEATED_SEQ_EXPR, untranslatable_returns_text=True
        )
        assertion_checks(
            resp.variation_descriptor, genomic_del5_free_text_rse_lse, q, ignore_id=True
        )

        resp = await test_handler.normalize(
            q, HGVSDupDelModeOption.LITERAL_SEQ_EXPR, untranslatable_returns_text=True
        )
        assertion_checks(
            resp.variation_descriptor, genomic_del5_free_text_rse_lse, q, ignore_id=True
        )

    # Invalid
    invalid_queries = [
        "NC_000023.10:g.(?_155270550)_155270570del",
        "NC_000023.11:g.(?_18593474)_18671749del"
        "CDKL5  g.(?_18443702)_18671700del",  # 37
        "CDKL5  g.(?_18425585)_18653631del",  # 38
        "CDKL5  g.(?_18425582)_18653500del",  # 38
    ]
    await assert_text_variation(invalid_queries, test_handler)


@pytest.mark.asyncio
async def test_genomic_del6(
    test_handler,
    genomic_del6_cn,
    genomic_del6_cx,
    genomic_del6_rse_lse,
    genomic_del6_free_text_cn,
    genomic_del6_free_text_cx,
    genomic_del6_free_text_rse_lse,
):
    """Test that genomic deletion works correctly."""
    q = "NC_000006.12:g.133462764_(133464858_?)del"  # 38
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp.variation_descriptor, genomic_del6_cx, q)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
    )
    assertion_checks(resp.variation_descriptor, genomic_del6_cn, q)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE)
    assertion_checks(resp.variation_descriptor, genomic_del6_cx, q)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.REPEATED_SEQ_EXPR, untranslatable_returns_text=True
    )
    assertion_checks(resp.variation_descriptor, genomic_del6_rse_lse, q)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.LITERAL_SEQ_EXPR, untranslatable_returns_text=True
    )
    assertion_checks(resp.variation_descriptor, genomic_del6_rse_lse, q)

    q = "NC_000006.11:g.133783902_(133785996_?)del"  # 37
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp.variation_descriptor, genomic_del6_cx, q, ignore_id=True)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
    )
    assertion_checks(resp.variation_descriptor, genomic_del6_cn, q, ignore_id=True)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE)
    assertion_checks(resp.variation_descriptor, genomic_del6_cx, q, ignore_id=True)

    genomic_del6_rse_lse.variation.definition = q
    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.REPEATED_SEQ_EXPR, untranslatable_returns_text=True
    )
    assertion_checks(resp.variation_descriptor, genomic_del6_rse_lse, q, ignore_id=True)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.LITERAL_SEQ_EXPR, untranslatable_returns_text=True
    )
    assertion_checks(resp.variation_descriptor, genomic_del6_rse_lse, q, ignore_id=True)

    # Free text
    for q in ["EYA4 g.133462764_(133464858_?)del"]:  # 38
        resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
        assertion_checks(
            resp.variation_descriptor, genomic_del6_free_text_cx, q, ignore_id=True
        )

        resp = await test_handler.normalize(
            q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
        )
        assertion_checks(
            resp.variation_descriptor, genomic_del6_free_text_cn, q, ignore_id=True
        )

        genomic_del6_rse_lse.variation.definition = q
        resp = await test_handler.normalize(
            q, HGVSDupDelModeOption.REPEATED_SEQ_EXPR, untranslatable_returns_text=True
        )
        assertion_checks(
            resp.variation_descriptor, genomic_del6_free_text_rse_lse, q, ignore_id=True
        )

        resp = await test_handler.normalize(
            q, HGVSDupDelModeOption.LITERAL_SEQ_EXPR, untranslatable_returns_text=True
        )
        assertion_checks(
            resp.variation_descriptor, genomic_del6_free_text_rse_lse, q, ignore_id=True
        )

    # Invalid
    invalid_queries = [
        "NC_000006.11:g.171115069_(171115080_?)del",
        "NC_000006.12:g.170805981_(170805989_?)del"
        "EYA4 g.133561700_(133853270_?)del",  # 37
        "EYA4 g.133561651_(133561708_?)del",  # 37
        "EYA4 g.133240513_(133240600_?)del",  # 38
        "EYA4 g.133240515_(133532130_?)del",  # 38
    ]
    await assert_text_variation(invalid_queries, test_handler)


@pytest.mark.asyncio
async def test_parameters(test_handler):
    """Check that valid and invalid parameters work as intended."""
    resp = await test_handler.normalize("7-140453136-A-T")
    assert resp.variation_descriptor
    assert resp.warnings == []

    q = "NC_000003.12:g.49531262dup"
    resp = await test_handler.normalize(q)
    assert resp.variation_descriptor
    assert resp.warnings == []

    resp = await test_handler.normalize(q, hgvs_dup_del_mode=None)
    assert resp.variation_descriptor
    assert resp.warnings == []

    resp = await test_handler.normalize(
        q, hgvs_dup_del_mode=HGVSDupDelModeOption.COPY_NUMBER_COUNT
    )
    assert resp.variation_descriptor is None
    assert resp.warnings == ["copy_number_count mode requires `baseline_copies`"]
