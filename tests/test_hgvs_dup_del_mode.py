"""Module for testing HGVS Dup Del mode."""

import pytest
from ga4gh.vrs import models

from tests.conftest import assertion_checks, cnv_assertion_checks
from variation.schemas.normalize_response_schema import HGVSDupDelModeOption


@pytest.fixture(scope="module")
def test_handler(test_query_handler):
    """Create test fixture for normalize handler"""
    return test_query_handler.normalize_handler


@pytest.fixture(scope="module")
def genomic_dup1_lse(genomic_dup1_seq_loc_normalized):
    """Create a test fixture for genomic dup LSE."""
    digest = "vfLfV0PTIdjGBINwgHKFBoVjPSkZ7s5-"
    params = {
        "type": "Allele",
        "id": f"ga4gh:VA.{digest}",
        "digest": digest,
        "location": genomic_dup1_seq_loc_normalized,
        "state": {
            "type": "ReferenceLengthExpression",
            "repeatSubunitLength": 1,
            "length": 3,
            "sequence": "GGG",
        },
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def genomic_dup1_cx(genomic_dup1_seq_loc_not_normalized):
    """Create a test fixture for genomic dup copy number change."""
    digest = "GU4QwgBqWmlwshiv6CCtvynDEHkHBVv2"
    params = {
        "type": "CopyNumberChange",
        "id": f"ga4gh:CX.{digest}",
        "digest": digest,
        "location": genomic_dup1_seq_loc_not_normalized,
        "copyChange": {"primaryCode": "EFO:0030072"},
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_dup1_free_text_seq_loc_normalized():
    """Create genomic dup1 free text sequence location"""
    return {
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.tpvbnWsfEGqip8gJQZnWJAF8-bWDUDKd",
        },
        "start": 1032,
        "end": 1034,
        "sequence": "GG",
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_dup1_free_text_seq_loc_not_normalized():
    """Create genomic dup1 free text sequence location"""
    return {
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.tpvbnWsfEGqip8gJQZnWJAF8-bWDUDKd",
        },
        "start": 1033,
        "end": 1034,
        "sequence": "G",
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_dup1_free_text_lse(genomic_dup1_free_text_seq_loc_normalized):
    """Create a test fixture for genomic dup LSE."""
    params = {
        "type": "Allele",
        "location": genomic_dup1_free_text_seq_loc_normalized,
        "state": {
            "type": "ReferenceLengthExpression",
            "repeatSubunitLength": 1,
            "length": 3,
            "sequence": "GGG",
        },
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def genomic_dup1_free_text_cn(genomic_dup1_free_text_seq_loc_not_normalized):
    """Create a test fixture for genomic dup copy number count."""
    params = {
        "type": "CopyNumberCount",
        "location": genomic_dup1_free_text_seq_loc_not_normalized,
        "copies": 3,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_dup2_lse(genomic_dup2_seq_loc_normalized):
    """Create a test fixture for genomic dup LSE."""
    params = {
        "type": "Allele",
        "location": genomic_dup2_seq_loc_normalized,
        "state": {
            "type": "ReferenceLengthExpression",
            "repeatSubunitLength": 4,
            "length": 8,
            "sequence": "TCTATCTA",
        },
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def genomic_dup2_cx(genomic_dup2_seq_loc_normalized):
    """Create a test fixture for genomic dup copy number change."""
    params = {
        "type": "CopyNumberChange",
        "location": genomic_dup2_seq_loc_normalized,
        "copyChange": {"primaryCode": "EFO:0030070"},
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def seq_loc_gt_100_bp():
    """Create seq loc for positions 33211290, 33211490 on NC_000023.11"""
    return {
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
        },
        "start": 33211289,
        "end": 33211490,
        "sequence": "TCTACTTCTTCCCACCAAAGCATTTTGAAAAGTGTATATCAAGGCAGCGATAAAAAAAACCTGGTAAAAGTTCTTCAAACTTTATTGCTCCAGTAGGCTTAAAAACAATGAGAAACCAACAAACTTCAGCAGCTTTAAAAAAAGTAACACTTCAGTTTTTCCTATTCGTTTTTCTCCGAAGGTAATTGCCTCCCAGATCTG",
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_dup2_rle2(seq_loc_gt_100_bp):
    """Create a test fixture for genomic dup RSE where bp > 100."""
    params = {
        "type": "Allele",
        "location": seq_loc_gt_100_bp,
        "state": {
            "type": "ReferenceLengthExpression",
            "repeatSubunitLength": 201,
            "length": 402,
        },
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def genomic_dup2_free_text_seq_loc():
    """Create genomic dup2 free text sequence location"""
    return {
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.1DeZLYHMnd-smp3GDlpRxETb9_0AokO7",
        },
        "start": 256,
        "end": 260,
        "sequence": "TAGA",
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_dup2_free_text_default(genomic_dup2_free_text_seq_loc):
    """Create a test fixture for genomic dup default and LSE."""
    params = {
        "type": "Allele",
        "location": genomic_dup2_free_text_seq_loc,
        "state": {
            "type": "ReferenceLengthExpression",
            "repeatSubunitLength": 4,
            "length": 8,
            "sequence": "TAGATAGA",
        },
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def genomic_dup2_free_text_cn(genomic_dup2_free_text_seq_loc):
    """Create a test fixture for genomic dup copy number count."""
    params = {
        "type": "CopyNumberCount",
        "location": genomic_dup2_free_text_seq_loc,
        "copies": 3,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_dup3_cx(genomic_del3_dup3_loc_not_normalized):
    """Create a test fixture for genomic dup copy number change."""
    params = {
        "type": "CopyNumberChange",
        "location": genomic_del3_dup3_loc_not_normalized,
        "copyChange": {"primaryCode": "EFO:0030070"},
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_dup3_free_text_subject():
    """Create test fixture for genomic dup3 free text location"""
    return {
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
        },
        "start": [31147273, 31147277],
        "end": [31182737, 31182739],
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_dup3_free_text_cx(genomic_dup3_free_text_subject):
    """Create a test fixture for genomic dup copy number change."""
    params = {
        "type": "CopyNumberChange",
        "location": genomic_dup3_free_text_subject,
        "copyChange": {"primaryCode": "EFO:0030070"},
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_dup3_free_text_cn(genomic_dup3_free_text_subject):
    """Create a test fixture for genomic dup copy number count."""
    params = {
        "type": "CopyNumberCount",
        "location": genomic_dup3_free_text_subject,
        "copies": 4,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_dup4_cx(genomic_dup4_loc):
    """Create a test fixture for genomic dup copy number change."""
    params = {
        "type": "CopyNumberChange",
        "location": genomic_dup4_loc,
        "copyChange": {"primaryCode": "EFO:0030070"},
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_dup4_cn(genomic_dup4_loc):
    """Create a test fixture for genomic dup copy number count."""
    params = {
        "type": "CopyNumberCount",
        "location": genomic_dup4_loc,
        "copies": 3,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_dup4_free_text_subject():
    """Create test fixture for genomic dup4 free text location"""
    return {
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.dLZ15tNO1Ur0IcGjwc3Sdi_0A6Yf4zm7",
        },
        "start": [None, 1674441],
        "end": [1684571, None],
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_dup4_free_text_cx(genomic_dup4_free_text_subject):
    """Create a test fixture for genomic dup copy number change."""
    params = {
        "type": "CopyNumberChange",
        "location": genomic_dup4_free_text_subject,
        "copyChange": {"primaryCode": "EFO:0030070"},
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_dup4_free_text_cn(genomic_dup4_free_text_subject):
    """Create a test fixture for genomic dup copy number count."""
    params = {
        "type": "CopyNumberCount",
        "location": genomic_dup4_free_text_subject,
        "copies": 3,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_dup5_cx(genomic_dup5_loc):
    """Create a test fixture for genomic dup5 copy number change."""
    params = {
        "type": "CopyNumberChange",
        "location": genomic_dup5_loc,
        "copyChange": {"primaryCode": "EFO:0030070"},
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_dup5_cn(genomic_dup5_loc):
    """Create a test fixture for genomic dup5 copy number count."""
    params = {
        "type": "CopyNumberCount",
        "location": genomic_dup5_loc,
        "copies": 3,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_dup6_cx(genomic_dup6_loc):
    """Create a test fixture for genomic dup copy number change."""
    params = {
        "type": "CopyNumberChange",
        "location": genomic_dup6_loc,
        "copyChange": {"primaryCode": "EFO:0030070"},
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_dup6_cn(genomic_dup6_loc):
    """Create a test fixture for genomic dup copy number count."""
    params = {
        "type": "CopyNumberCount",
        "location": genomic_dup6_loc,
        "copies": 2,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_del1_lse(genomic_del1_seq_loc):
    """Create a test fixture for genomic del LSE."""
    params = {
        "type": "Allele",
        "location": genomic_del1_seq_loc,
        "state": {
            "type": "ReferenceLengthExpression",
            "repeatSubunitLength": 1,
            "length": 0,
            "sequence": "",
        },
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def genomic_del1_cx(genomic_del1_seq_loc):
    """Create a test fixture for genomic del copy number change."""
    params = {
        "type": "CopyNumberChange",
        "location": genomic_del1_seq_loc,
        "copyChange": {"primaryCode": "EFO:0030064"},
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_del1_rle(genomic_del1_seq_loc):
    """Create a test fixture for genomic del RSE."""
    params = {
        "type": "Allele",
        "location": genomic_del1_seq_loc,
        "state": {
            "type": "ReferenceLengthExpression",
            "repeatSubunitLength": 1,
            "length": 2,
        },
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def genomic_del1_free_text_seq_loc():
    """Create genomic del1 free text sequence location"""
    return {
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.xBKOKptLLDr-k4hTyCetvARn16pDS_rW",
        },
        "start": 557,
        "end": 558,
        "sequence": "T",
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_del1_free_text_lse(genomic_del1_free_text_seq_loc):
    """Create a test fixture for genomic del LSE."""
    params = {
        "type": "Allele",
        "location": genomic_del1_free_text_seq_loc,
        "state": {
            "type": "ReferenceLengthExpression",
            "repeatSubunitLength": 1,
            "length": 0,
            "sequence": "",
        },
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def genomic_del1_free_text_cn(genomic_del1_free_text_seq_loc):
    """Create a test fixture for genomic del copy number count."""
    params = {
        "type": "CopyNumberCount",
        "location": genomic_del1_free_text_seq_loc,
        "copies": 1,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_del1_free_text_rle(genomic_del1_free_text_seq_loc):
    """Create a test fixture for genomic del RSE."""
    params = {
        "type": "Allele",
        "location": genomic_del1_free_text_seq_loc,
        "state": {
            "type": "ReferenceLengthExpression",
            "repeatSubunitLength": 1,
            "length": 0,
        },
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def genomic_del2_lse(genomic_del2_seq_loc):
    """Create a test fixture for genomic del LSE."""
    params = {
        "type": "Allele",
        "location": genomic_del2_seq_loc,
        "state": {
            "type": "ReferenceLengthExpression",
            "repeatSubunitLength": 19,
            "length": 0,
            "sequence": "",
        },
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def genomic_del2_lse2(seq_loc_gt_100_bp):
    """Create a test fixture for genomic del LSE where bp > 100."""
    params = {
        "type": "Allele",
        "location": seq_loc_gt_100_bp,
        "state": {
            "type": "ReferenceLengthExpression",
            "repeatSubunitLength": 201,
            "length": 0,
            "sequence": "",
        },
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def genomic_del2_cx(genomic_del2_seq_loc):
    """Create a test fixture for genomic del copy number change."""
    params = {
        "type": "CopyNumberChange",
        "location": genomic_del2_seq_loc,
        "copyChange": {"primaryCode": "EFO:0030069"},
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_del2_rle(genomic_del2_seq_loc):
    """Create a test fixture for genomic del RSE."""
    params = {
        "type": "Allele",
        "location": genomic_del2_seq_loc,
        "state": {
            "type": "ReferenceLengthExpression",
            "repeatSubunitLength": 19,
            "length": 0,
        },
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def genomic_del2_free_text_seq_loc():
    """Create genomic del2 free text sequence location"""
    return {
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.xBKOKptLLDr-k4hTyCetvARn16pDS_rW",
        },
        "start": 491,
        "end": 510,
        "sequence": "ATGTTGACGGACAGCCTAT",
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_del2_free_text_default(genomic_del2_free_text_seq_loc):
    """Create a test fixture for genomic del default and LSE."""
    params = {
        "type": "Allele",
        "location": genomic_del2_free_text_seq_loc,
        "state": {
            "type": "ReferenceLengthExpression",
            "repeatSubunitLength": 19,
            "length": 0,
            "sequence": "",
        },
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def genomic_del2_free_text_cnv(genomic_del2_free_text_seq_loc):
    """Create a test fixture for genomic del CNV."""
    params = {
        "type": "CopyNumberCount",
        "location": genomic_del2_free_text_seq_loc,
        "copies": 1,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_del3_cx(genomic_del3_dup3_loc_not_normalized):
    """Create a test fixture for genomic del copy number change."""
    params = {
        "type": "CopyNumberChange",
        "location": genomic_del3_dup3_loc_not_normalized,
        "copyChange": {"primaryCode": "EFO:0030067"},
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_del3_free_text_subject():
    """Create test fixture for genomic del3 free text location"""
    return {
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
        },
        "start": [68839264, 68839267],
        "end": [68841120, 68841125],
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_del3_free_text_cx(genomic_del3_free_text_subject):
    """Create a test fixture for genomic del copy number change."""
    params = {
        "type": "CopyNumberChange",
        "location": genomic_del3_free_text_subject,
        "copyChange": {"primaryCode": "EFO:0030067"},
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_del3_free_text_cn(genomic_del3_free_text_subject):
    """Create a test fixture for genomic del copy number count."""
    params = {
        "type": "CopyNumberCount",
        "location": genomic_del3_free_text_subject,
        "copies": 2,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_del4_cx(genomic_del4_seq_loc):
    """Create a test fixture for genomic del copy number change."""
    params = {
        "type": "CopyNumberChange",
        "location": genomic_del4_seq_loc,
        "copyChange": {"primaryCode": "EFO:0030067"},
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_del4_cn(genomic_del4_seq_loc):
    """Create a test fixture for genomic del copy number count."""
    params = {
        "type": "CopyNumberCount",
        "location": genomic_del4_seq_loc,
        "copies": 1,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_del4_free_text_subject():
    """Create test fixture for genomic del4 free text location"""
    return {
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.pnAqCRBrTsUoBghSD1yp_jXWSmlbdh4g",
        },
        "start": [None, 227022027],
        "end": [227025830, None],
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_del4_free_text_cx(genomic_del4_free_text_subject):
    """Create a test fixture for genomic del copy number change."""
    params = {
        "type": "CopyNumberChange",
        "location": genomic_del4_free_text_subject,
        "copyChange": {"primaryCode": "EFO:0030067"},
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_del4_free_text_cn(genomic_del4_free_text_subject):
    """Create a test fixture for genomic del copy number count."""
    params = {
        "type": "CopyNumberCount",
        "location": genomic_del4_free_text_subject,
        "copies": 1,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_uncertain_del_2():
    """Create a genomic uncertain deletion on chr 2 test fixture."""
    params = {
        "location": {
            "sequenceReference": {
                "type": "SequenceReference",
                "refgetAccession": "SQ.pnAqCRBrTsUoBghSD1yp_jXWSmlbdh4g",
            },
            "start": [None, 110104899],
            "end": [110207160, None],
            "type": "SequenceLocation",
        },
        "copyChange": {"primaryCode": "EFO:0030067"},
        "type": "CopyNumberChange",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_uncertain_del_y():
    """Create a genomic uncertain deletion on chr Y test fixture."""
    params = {
        "location": {
            "sequenceReference": {
                "type": "SequenceReference",
                "refgetAccession": "SQ.8_liLu1aycC0tPQPFmUaGXJLDs5SbPZ5",
            },
            "start": [None, 14076801],
            "end": [57165209, None],
            "type": "SequenceLocation",
        },
        "copyChange": {"primaryCode": "EFO:0030067"},
        "type": "CopyNumberChange",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_del5_cn_var(genomic_del5_seq_loc):
    """Create genomic del5 copy number count"""
    params = {
        "type": "CopyNumberCount",
        "location": genomic_del5_seq_loc,
        "copies": 3,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_del5_cx_var(genomic_del5_seq_loc):
    """Create genomic del5 copy number change"""
    params = {
        "type": "CopyNumberChange",
        "location": genomic_del5_seq_loc,
        "copyChange": {"primaryCode": "EFO:0030067"},
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_del6_cx_var(genomic_del6_seq_loc):
    """Create genomic del6 copy number change"""
    params = {
        "type": "CopyNumberChange",
        "location": genomic_del6_seq_loc,
        "copyChange": {"primaryCode": "EFO:0030067"},
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_del6_cn_var(genomic_del6_seq_loc):
    """Create genomic del6 copy number count"""
    params = {
        "type": "CopyNumberCount",
        "location": genomic_del6_seq_loc,
        "copies": 1,
    }
    return models.CopyNumberCount(**params)


def no_variation_check(resp, q):
    """Check that variation is None in normalize response"""
    assert resp.variation is None, q


@pytest.mark.asyncio()
async def invalid_query_list_checks(query_list, test_handler):
    """Check that invalid queries in query list do not normalize"""
    for q in query_list:
        resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
        no_variation_check(resp, q)


@pytest.mark.asyncio()
async def test_genomic_dup1(
    test_handler,
    genomic_dup1_lse,
    genomic_dup1_38_cn,
    genomic_dup1_cx,
    genomic_dup1_free_text_lse,
    genomic_dup1_free_text_cn,
):
    """Test that genomic duplication works correctly."""
    # https://reg.clinicalgenome.org/redmine/projects/registry/genboree_registry/allele?hgvsOrDescriptor=NC_000003.12%3Ag.49531262dup
    q = "NC_000003.12:g.49531262dup"  # 38
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp, genomic_dup1_lse, check_vrs_id=True, mane_genes_exts=True)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp, genomic_dup1_lse, check_vrs_id=True, mane_genes_exts=True)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
    )
    cnv_assertion_checks(
        resp, genomic_dup1_38_cn, check_vrs_id=True, mane_genes_exts=True
    )

    resp = await test_handler.normalize(
        q,
        HGVSDupDelModeOption.COPY_NUMBER_CHANGE,
        copy_change=models.CopyChange.EFO_0030072,
    )
    cnv_assertion_checks(resp, genomic_dup1_cx, check_vrs_id=True, mane_genes_exts=True)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.ALLELE)
    assertion_checks(resp, genomic_dup1_lse, check_vrs_id=True, mane_genes_exts=True)

    q = "NC_000003.11:g.49568695dup"  # 37
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp, genomic_dup1_lse, check_vrs_id=True, mane_genes_exts=True)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
    )
    cnv_assertion_checks(
        resp, genomic_dup1_38_cn, check_vrs_id=True, mane_genes_exts=True
    )

    resp = await test_handler.normalize(
        q,
        HGVSDupDelModeOption.COPY_NUMBER_CHANGE,
        copy_change=models.CopyChange.EFO_0030072,
    )
    cnv_assertion_checks(resp, genomic_dup1_cx, check_vrs_id=True, mane_genes_exts=True)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.ALLELE)
    assertion_checks(resp, genomic_dup1_lse, check_vrs_id=True, mane_genes_exts=True)

    # Free Text
    for q in ["DAG1 g.49568695dup", "DAG1 g.49531262dup"]:  # 37  # 38
        resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
        assertion_checks(resp, genomic_dup1_free_text_lse, mane_genes_exts=True)

        resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
        assertion_checks(resp, genomic_dup1_free_text_lse, mane_genes_exts=True)

        resp = await test_handler.normalize(
            q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
        )
        cnv_assertion_checks(resp, genomic_dup1_free_text_cn, mane_genes_exts=True)

        resp = await test_handler.normalize(q, HGVSDupDelModeOption.ALLELE)
        assertion_checks(resp, genomic_dup1_free_text_lse, mane_genes_exts=True)

    # Invalid
    invalid_queries = [
        "NC_000007.13:g.159138670dup",
        "NC_000007.14:g.159345976dup",
        "BRAF g.140219337dup",
        "BRAF g.141024929dup",
    ]
    await invalid_query_list_checks(invalid_queries, test_handler)


@pytest.mark.asyncio()
async def test_genomic_dup2(
    test_handler,
    genomic_dup2_lse,
    genomic_dup2_38_cn,
    genomic_dup2_cx,
    genomic_dup2_free_text_default,
    genomic_dup2_free_text_cn,
    genomic_dup2_rle2,
):
    """Test that genomic duplication works correctly."""
    # https://reg.clinicalgenome.org/redmine/projects/registry/genboree_registry/allele?hgvsOrDescriptor=NM_004006.2%3Ac.20_23dup
    q = "NC_000023.11:g.33211290_33211293dup"  # 38
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp, genomic_dup2_lse, mane_genes_exts=True)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
    )
    cnv_assertion_checks(resp, genomic_dup2_38_cn, mane_genes_exts=True)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE)
    cnv_assertion_checks(resp, genomic_dup2_cx, mane_genes_exts=True)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.ALLELE)
    assertion_checks(resp, genomic_dup2_lse, mane_genes_exts=True)

    q = "NC_000023.10:g.33229407_33229410dup"  # 37
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp, genomic_dup2_lse, mane_genes_exts=True)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
    )
    cnv_assertion_checks(resp, genomic_dup2_38_cn, mane_genes_exts=True)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE)
    cnv_assertion_checks(resp, genomic_dup2_cx, mane_genes_exts=True)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.ALLELE)
    assertion_checks(resp, genomic_dup2_lse, mane_genes_exts=True)

    # Free text
    for q in ["DMD g.33211290_33211293dup", "DMD g.33229407_33229410dup"]:  # 37  # 38
        resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
        assertion_checks(resp, genomic_dup2_free_text_default, mane_genes_exts=True)

        resp = await test_handler.normalize(
            q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
        )
        cnv_assertion_checks(resp, genomic_dup2_free_text_cn, mane_genes_exts=True)

        resp = await test_handler.normalize(q, HGVSDupDelModeOption.ALLELE)
        assertion_checks(resp, genomic_dup2_free_text_default, mane_genes_exts=True)

    # Greater than 100 bps -> rse
    q = "NC_000023.11:g.33211290_33211490dup"
    resp = await test_handler.normalize(q)
    assertion_checks(resp, genomic_dup2_rle2, mane_genes_exts=True)

    # Invalid
    invalid_queries = [
        "NC_000007.13:g.140413127_159138670dup",
        "NC_000007.14:g.140413127_159345976dup",
        "BRAF g.140219337_140924929dup",
        "BRAF g.140719326_141024929dup",
    ]
    await invalid_query_list_checks(invalid_queries, test_handler)


@pytest.mark.asyncio()
async def test_genomic_dup3(
    test_handler,
    genomic_dup3_cx,
    genomic_del3_dup3_cn_38,
    genomic_dup3_free_text_cn,
    genomic_dup3_free_text_cx,
):
    """Test that genomic duplication works correctly."""
    q = "NC_000023.11:g.(31060227_31100351)_(33274278_33417151)dup"  # 38
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    cnv_assertion_checks(resp, genomic_dup3_cx)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=1
    )
    cnv_assertion_checks(resp, genomic_del3_dup3_cn_38)

    resp = await test_handler.normalize(
        q,
        HGVSDupDelModeOption.COPY_NUMBER_CHANGE,
        copy_change=models.CopyChange.EFO_0030070,
    )
    cnv_assertion_checks(resp, genomic_dup3_cx)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.ALLELE)
    no_variation_check(resp, q)

    q = "NC_000023.10:g.(31078344_31118468)_(33292395_33435268)dup"  # 37
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    cnv_assertion_checks(resp, genomic_dup3_cx)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=1
    )
    cnv_assertion_checks(resp, genomic_del3_dup3_cn_38)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE)
    cnv_assertion_checks(resp, genomic_dup3_cx)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.ALLELE)
    no_variation_check(resp, q)

    # Free Text
    for q in ["DMD g.(31147274_31147278)_(31182737_31182739)dup"]:  # 38
        resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
        cnv_assertion_checks(resp, genomic_dup3_free_text_cx)

        resp = await test_handler.normalize(
            q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=3
        )
        cnv_assertion_checks(resp, genomic_dup3_free_text_cn)

        resp = await test_handler.normalize(q, HGVSDupDelModeOption.ALLELE)
        no_variation_check(resp, q)

    # Invalid
    invalid_queries = [
        "NC_000023.10:g.(31119221_31119227)_(31119300_155270562)dup",
        "NC_000023.11:g.(31119221_31119227)_(31119300_156040899)dup",
        "DMD g.(31060227_31100351)_(33274278_33417151)dup",
    ]
    await invalid_query_list_checks(invalid_queries, test_handler)


@pytest.mark.asyncio()
async def test_genomic_dup4(
    test_handler,
    genomic_dup4_cn,
    genomic_dup4_cx,
    genomic_dup4_free_text_cn,
    genomic_dup4_free_text_cx,
):
    """Test that genomic duplication works correctly."""
    q = "NC_000020.11:g.(?_30417576)_(31394018_?)dup"  # 38
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    cnv_assertion_checks(resp, genomic_dup4_cx)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
    )
    cnv_assertion_checks(resp, genomic_dup4_cn)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE)
    cnv_assertion_checks(resp, genomic_dup4_cx)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.ALLELE)
    no_variation_check(resp, q)

    q = "NC_000020.10:g.(?_29652252)_(29981821_?)dup"  # 37
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    cnv_assertion_checks(resp, genomic_dup4_cx)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
    )
    cnv_assertion_checks(resp, genomic_dup4_cn)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE)
    cnv_assertion_checks(resp, genomic_dup4_cx)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.ALLELE)
    no_variation_check(resp, q)

    # Free Text
    for q in [
        "PRPF8 g.(?_1577736)_(1587865_?)dup",  # 37
        "PRPF8 g.(?_1674442)_(1684571_?)dup",  # 38
    ]:
        resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
        cnv_assertion_checks(resp, genomic_dup4_free_text_cx)

        resp = await test_handler.normalize(
            q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
        )
        cnv_assertion_checks(resp, genomic_dup4_free_text_cn)

        resp = await test_handler.normalize(q, HGVSDupDelModeOption.ALLELE)
        no_variation_check(resp, q)

    # Invalid
    invalid_queries = [
        "NC_000020.10:g.(?_29652252)_(63025530_?)dup",
        "NC_000020.11:g.(?_29652252)_(64444169_?)dup",
        "PRPF8 g.(?_1650628)_(1684571_?)dup",
    ]
    await invalid_query_list_checks(invalid_queries, test_handler)


@pytest.mark.asyncio()
async def test_genomic_dup5(
    test_handler,
    genomic_dup5_cn,
    genomic_dup5_cx,
):
    """Test that genomic duplication works correctly."""
    q = "NC_000023.11:g.(?_154021812)_154092209dup"  # 38
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    cnv_assertion_checks(resp, genomic_dup5_cx)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
    )
    cnv_assertion_checks(resp, genomic_dup5_cn)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE)
    cnv_assertion_checks(resp, genomic_dup5_cx)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.ALLELE)
    no_variation_check(resp, q)

    q = "NC_000023.10:g.(?_153287263)_153357667dup"  # 37
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    cnv_assertion_checks(resp, genomic_dup5_cx)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
    )
    cnv_assertion_checks(resp, genomic_dup5_cn)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE)
    cnv_assertion_checks(resp, genomic_dup5_cx)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.ALLELE)
    no_variation_check(resp, q)

    # Free Text
    for q in [
        "MECP2 g.(?_153287263)_153357667dup",  # 37
        "MECP2 g.(?_154021812)_154092209dup",  # 38
    ]:
        resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
        cnv_assertion_checks(resp, genomic_dup5_cx)

        resp = await test_handler.normalize(
            q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
        )
        cnv_assertion_checks(resp, genomic_dup5_cn)

        resp = await test_handler.normalize(q, HGVSDupDelModeOption.ALLELE)
        no_variation_check(resp, q)

    # Invalid
    for q in [
        "NC_000023.10:g.(?_153287263)_155270561dup",
        "NC_000023.11:g.(?_154021812)_156040896dup",
        "MECP2 g.(?_154021812)_154097733dup"  # 37
        "MECP2 g.(?_154021572)_154092209dup",  # 38
    ]:
        resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
        assert resp.variation is None, q


@pytest.mark.asyncio()
async def test_genomic_dup6(
    test_handler,
    genomic_dup6_cn,
    genomic_dup6_cx,
):
    """Test that genomic duplication works correctly."""
    q = "NC_000023.11:g.154021812_(154092209_?)dup"  # 38
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    cnv_assertion_checks(resp, genomic_dup6_cx)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=1
    )
    cnv_assertion_checks(resp, genomic_dup6_cn)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE)
    cnv_assertion_checks(resp, genomic_dup6_cx)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.ALLELE)
    no_variation_check(resp, q)

    q = "NC_000023.10:g.153287263_(153357667_?)dup"  # 37
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    cnv_assertion_checks(resp, genomic_dup6_cx)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=1
    )
    cnv_assertion_checks(resp, genomic_dup6_cn)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE)
    cnv_assertion_checks(resp, genomic_dup6_cx)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.ALLELE)
    no_variation_check(resp, q)

    # Free Text
    for q in [
        "MECP2 g.153287263_(153357667_?)dup",  # 37
        "MECP2 g.154021812_(154092209_?)dup",  # 38
    ]:
        resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
        cnv_assertion_checks(resp, genomic_dup6_cx)

        resp = await test_handler.normalize(
            q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=1
        )
        cnv_assertion_checks(resp, genomic_dup6_cn)

        resp = await test_handler.normalize(q, HGVSDupDelModeOption.ALLELE)
        no_variation_check(resp, q)

    # Invalid
    for q in [
        "NC_000023.10:g.153287263_(155270561_?)dup",
        "NC_000023.11:g.154021812_(156040896_?)dup",
        "MECP2 g.154021812_(154097733_?)dup"  # 37
        "MECP2 g.154021572_(154092209_?)dup",  # 38
    ]:
        resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
        assert resp.variation is None, q


@pytest.mark.asyncio()
async def test_genomic_del1(
    test_handler,
    genomic_del1_lse,
    genomic_del1_38_cn,
    genomic_del1_cx,
    genomic_del1_free_text_lse,
    genomic_del1_free_text_cn,
):
    """Test that genomic deletion works correctly."""
    q = "NC_000003.12:g.10149811del"  # 38
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp, genomic_del1_lse, mane_genes_exts=True)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
    )
    cnv_assertion_checks(resp, genomic_del1_38_cn, mane_genes_exts=True)

    resp = await test_handler.normalize(
        q,
        HGVSDupDelModeOption.COPY_NUMBER_CHANGE,
        copy_change=models.CopyChange.EFO_0030064,
    )
    cnv_assertion_checks(resp, genomic_del1_cx, mane_genes_exts=True)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.ALLELE)
    assertion_checks(resp, genomic_del1_lse, mane_genes_exts=True)

    q = "NC_000003.11:g.10191495del"  # 37
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp, genomic_del1_lse, mane_genes_exts=True)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
    )
    cnv_assertion_checks(resp, genomic_del1_38_cn, mane_genes_exts=True)

    resp = await test_handler.normalize(
        q,
        HGVSDupDelModeOption.COPY_NUMBER_CHANGE,
        copy_change=models.CopyChange.EFO_0030064,
    )
    cnv_assertion_checks(resp, genomic_del1_cx, mane_genes_exts=True)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.ALLELE)
    assertion_checks(resp, genomic_del1_lse, mane_genes_exts=True)

    # Free text
    for q in ["VHL g.10191495del", "VHL g.10149811del"]:  # 37  # 38
        resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
        assertion_checks(resp, genomic_del1_free_text_lse, mane_genes_exts=True)

        resp = await test_handler.normalize(
            q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
        )
        cnv_assertion_checks(resp, genomic_del1_free_text_cn, mane_genes_exts=True)

        resp = await test_handler.normalize(q, HGVSDupDelModeOption.ALLELE)
        assertion_checks(resp, genomic_del1_free_text_lse, mane_genes_exts=True)

    # Invalid
    invalid_queries = [
        "NC_000003.11:g.198022431del",
        "NC_000003.12:g.198295567del",
        "BRAF g.140413127del",
        "BRAF g.141024929del",
    ]
    await invalid_query_list_checks(invalid_queries, test_handler)


@pytest.mark.asyncio()
async def test_genomic_del2(
    test_handler,
    genomic_del2_lse,
    genomic_del2_38_cn,
    genomic_del2_cx,
    genomic_del2_free_text_default,
    genomic_del2_free_text_cnv,
    genomic_del2_lse2,
):
    """Test that genomic deletion works correctly."""
    q = "NC_000003.12:g.10146595_10146613del"  # 38
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp, genomic_del2_lse, mane_genes_exts=True)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
    )
    cnv_assertion_checks(resp, genomic_del2_38_cn, mane_genes_exts=True)

    resp = await test_handler.normalize(
        q,
        HGVSDupDelModeOption.COPY_NUMBER_CHANGE,
        copy_change=models.CopyChange.EFO_0030069,
    )
    cnv_assertion_checks(resp, genomic_del2_cx, mane_genes_exts=True)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.ALLELE)
    assertion_checks(resp, genomic_del2_lse, mane_genes_exts=True)

    q = "NC_000003.11:g.10188279_10188297del"  # 37
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp, genomic_del2_lse, mane_genes_exts=True)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
    )
    cnv_assertion_checks(resp, genomic_del2_38_cn, mane_genes_exts=True)

    resp = await test_handler.normalize(
        q,
        HGVSDupDelModeOption.COPY_NUMBER_CHANGE,
        copy_change=models.CopyChange.EFO_0030069,
    )
    cnv_assertion_checks(resp, genomic_del2_cx, mane_genes_exts=True)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.ALLELE)
    assertion_checks(resp, genomic_del2_lse, mane_genes_exts=True)

    # Free text
    for q in ["VHL g.10188279_10188297del", "VHL g.10146595_10146613del"]:  # 37  # 38
        resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
        assertion_checks(resp, genomic_del2_free_text_default, mane_genes_exts=True)

        resp = await test_handler.normalize(
            q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
        )
        cnv_assertion_checks(resp, genomic_del2_free_text_cnv, mane_genes_exts=True)

        resp = await test_handler.normalize(q, HGVSDupDelModeOption.ALLELE)
        assertion_checks(resp, genomic_del2_free_text_default, mane_genes_exts=True)

    # Check that del > 100 bps returns LSE
    q = "NC_000023.11:g.33211290_33211490del"
    resp = await test_handler.normalize(q)
    assertion_checks(resp, genomic_del2_lse2, mane_genes_exts=True)

    # gnomad vcf
    q = "3-10146594-AATGTTGACGGACAGCCTAT-A"
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp, genomic_del2_lse, mane_genes_exts=True)

    q = "3-10188278-AATGTTGACGGACAGCCTAT-A"
    resp = await test_handler.normalize(q)
    assertion_checks(resp, genomic_del2_lse, mane_genes_exts=True)

    # Invalid
    invalid_queries = [
        "NC_000003.12:g.10146595_198295580del",
        "NC_000003.11:g.198022435_198022437del",
        "BRAF g.140413127_140419136del",
        "BRAF g.140719326_141024929del",
    ]
    await invalid_query_list_checks(invalid_queries, test_handler)


@pytest.mark.asyncio()
async def test_genomic_del3(
    test_handler,
    genomic_del3_dup3_cn_38,
    genomic_del3_cx,
    genomic_del3_free_text_cn,
    genomic_del3_free_text_cx,
):
    """Test that genomic deletion works correctly."""
    q = "NC_000023.11:g.(31060227_31100351)_(33274278_33417151)del"  # 38
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    cnv_assertion_checks(resp, genomic_del3_cx)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=3
    )
    cnv_assertion_checks(resp, genomic_del3_dup3_cn_38)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE)
    cnv_assertion_checks(resp, genomic_del3_cx)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.ALLELE)
    no_variation_check(resp, q)

    q = "NC_000023.10:g.(31078344_31118468)_(33292395_33435268)del"  # 37
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    cnv_assertion_checks(resp, genomic_del3_cx)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=3
    )
    cnv_assertion_checks(resp, genomic_del3_dup3_cn_38)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE)
    cnv_assertion_checks(resp, genomic_del3_cx)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.ALLELE)
    no_variation_check(resp, q)

    # Free Text
    for q in [
        "EFNB1 g.(68059108_68059111)_(68060963_68060968)del",  # 37
        "EFNB1 g.(68839265_68839268)_(68841120_68841125)del",  # 38
    ]:
        resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
        cnv_assertion_checks(resp, genomic_del3_free_text_cx)

        resp = await test_handler.normalize(
            q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=3
        )
        cnv_assertion_checks(resp, genomic_del3_free_text_cn)

        resp = await test_handler.normalize(q, HGVSDupDelModeOption.ALLELE)
        no_variation_check(resp, q)

    # Invalid
    invalid_queries = [
        "NC_000023.11:g.(156040880_156040883)_(156040896_156040899)del",
        "NC_000023.10:g.(155270550_155270555)_(155270560_155270562)del",
        "EFNB1 g.(68048863_68048870)_(68842150_68842152)del",  # 37
        "EFNB1 g.(68829022_68829030)_(68842150_68842161)del",  # 38
    ]
    await invalid_query_list_checks(invalid_queries, test_handler)


@pytest.mark.asyncio()
async def test_genomic_del4(
    test_handler,
    genomic_del4_cn,
    genomic_del4_cx,
    genomic_uncertain_del_2,
    genomic_uncertain_del_y,
    genomic_del4_free_text_cn,
    genomic_del4_free_text_cx,
):
    """Test that genomic deletion works correctly."""
    q = "NC_000023.11:g.(?_31120496)_(33339477_?)del"  # 38
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    cnv_assertion_checks(resp, genomic_del4_cx)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
    )
    cnv_assertion_checks(resp, genomic_del4_cn)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE)
    cnv_assertion_checks(resp, genomic_del4_cx)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.ALLELE)
    no_variation_check(resp, q)

    q = "NC_000023.10:g.(?_31138613)_(33357594_?)del"  # 37
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    cnv_assertion_checks(resp, genomic_del4_cx)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
    )
    cnv_assertion_checks(resp, genomic_del4_cn)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE)
    cnv_assertion_checks(resp, genomic_del4_cx)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.ALLELE)
    no_variation_check(resp, q)

    q = "NC_000002.12:g.(?_110104900)_(110207160_?)del"
    resp = await test_handler.normalize(q)
    cnv_assertion_checks(resp, genomic_uncertain_del_2)

    q = "NC_000024.10:g.(?_14076802)_(57165209_?)del"
    resp = await test_handler.normalize(q)
    cnv_assertion_checks(resp, genomic_uncertain_del_y)

    # Free Text
    for q in ["COL4A4 g.(?_227022028)_(227025830_?)del"]:  # 38
        resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
        cnv_assertion_checks(resp, genomic_del4_free_text_cx)

        resp = await test_handler.normalize(
            q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
        )
        cnv_assertion_checks(resp, genomic_del4_free_text_cn)

        resp = await test_handler.normalize(q, HGVSDupDelModeOption.ALLELE)
        no_variation_check(resp, q)

    # Invalid
    invalid_queries = [
        "NC_000023.11:g.(?_156040899)_(156040900_?)del",
        "NC_000024.10:g.(?_155270565)_(155270568_?)del",
        "COL4A4 g.(?_227002710)_(227003710_?)del",
        "COL4A4 g.(?_227867430)_(228029276_?)del",
    ]
    await invalid_query_list_checks(invalid_queries, test_handler)


@pytest.mark.asyncio()
async def test_genomic_del5(
    test_handler,
    genomic_del5_cn_var,
    genomic_del5_cx_var,
):
    """Test that genomic deletion works correctly."""
    q = "NC_000023.11:g.(?_18575354)_18653629del"  # 38
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    cnv_assertion_checks(resp, genomic_del5_cx_var)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=4
    )
    cnv_assertion_checks(resp, genomic_del5_cn_var)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE)
    cnv_assertion_checks(resp, genomic_del5_cx_var)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.ALLELE)
    no_variation_check(resp, q)

    q = "NC_000023.10:g.(?_18593474)_18671749del"  # 37
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    cnv_assertion_checks(resp, genomic_del5_cx_var)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=4
    )
    cnv_assertion_checks(resp, genomic_del5_cn_var)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE)
    cnv_assertion_checks(resp, genomic_del5_cx_var)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.ALLELE)
    no_variation_check(resp, q)

    # Free text
    for q in ["CDKL5 g.(?_18575354)_18653629del"]:
        resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
        cnv_assertion_checks(resp, genomic_del5_cx_var)

        resp = await test_handler.normalize(
            q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=4
        )
        cnv_assertion_checks(resp, genomic_del5_cn_var)

        resp = await test_handler.normalize(q, HGVSDupDelModeOption.ALLELE)
        no_variation_check(resp, q)

    # Invalid
    invalid_queries = [
        "NC_000023.10:g.(?_155270550)_155270570del",
        "NC_000023.11:g.(?_18593474)_18671749del"
        "CDKL5  g.(?_18443702)_18671700del",  # 37
        "CDKL5  g.(?_18425585)_18653631del",  # 38
        "CDKL5  g.(?_18425582)_18653500del",  # 38
    ]
    await invalid_query_list_checks(invalid_queries, test_handler)


@pytest.mark.asyncio()
async def test_genomic_del6(
    test_handler,
    genomic_del6_cn_var,
    genomic_del6_cx_var,
):
    """Test that genomic deletion works correctly."""
    q = "NC_000006.12:g.133462764_(133464858_?)del"  # 38
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    cnv_assertion_checks(resp, genomic_del6_cx_var)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
    )
    cnv_assertion_checks(resp, genomic_del6_cn_var)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE)
    cnv_assertion_checks(resp, genomic_del6_cx_var)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.ALLELE)
    no_variation_check(resp, q)

    q = "NC_000006.11:g.133783902_(133785996_?)del"  # 37
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    cnv_assertion_checks(resp, genomic_del6_cx_var)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
    )
    cnv_assertion_checks(resp, genomic_del6_cn_var)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE)
    cnv_assertion_checks(resp, genomic_del6_cx_var)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.ALLELE)
    no_variation_check(resp, q)

    # Free text
    for q in ["EYA4 g.133462764_(133464858_?)del"]:  # 38
        resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
        cnv_assertion_checks(resp, genomic_del6_cx_var)

        resp = await test_handler.normalize(
            q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
        )
        cnv_assertion_checks(resp, genomic_del6_cn_var)

        resp = await test_handler.normalize(q, HGVSDupDelModeOption.ALLELE)
        no_variation_check(resp, q)

    # Invalid
    invalid_queries = [
        "NC_000006.11:g.171115069_(171115080_?)del",
        "NC_000006.12:g.170805981_(170805989_?)del"
        "EYA4 g.133561700_(133853270_?)del",  # 37
        "EYA4 g.133561651_(133561708_?)del",  # 37
        "EYA4 g.133240513_(133240600_?)del",  # 38
        "EYA4 g.133240515_(133532130_?)del",  # 38
    ]
    await invalid_query_list_checks(invalid_queries, test_handler)


@pytest.mark.asyncio()
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
