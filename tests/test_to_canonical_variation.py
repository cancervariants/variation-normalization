"""Module for testing to_canonical_variation"""
import pytest
from ga4gh.vrsatile.pydantic.vrsatile_models import CanonicalVariation

from variation.schemas.normalize_response_schema import HGVSDupDelMode \
    as HGVSDupDelModeEnum


@pytest.fixture(scope="module")
def test_handler(test_query_handler):
    """Create test fixture for to canonical variation handler"""
    return test_query_handler.to_canonical_handler


@pytest.fixture(scope="module")
def variation1_seq_loc():
    """Create test fixture for variation1 sequence location"""
    return {
        "id": "ga4gh:SL.lMqM3PtNhpC6avJf_XTqvuLtZNDV1XVr",
        "type": "SequenceLocation",
        "sequence_id": "ga4gh:SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT",
        "start": {"type": "Number", "value": 20189346},
        "end": {"type": "Number", "value": 20189349}
    }


@pytest.fixture(scope="module")
def variation1_lse(variation1_seq_loc):
    """Create test fixture for NC_000013.11:20189346:GGG:GG"""
    params = {
        "id": "ga4gh:CAN.tg0nG9q1DMP_J-vcWiaASPW44GMEh47k",
        "type": "CanonicalVariation",
        "canonical_context": {
            "id": "ga4gh:VA.jkAILAe4dK4tQ3y2hz-GHtZRAnbVC__T",
            "type": "Allele",
            "location": variation1_seq_loc,
            "state": {
                "type": "LiteralSequenceExpression",
                "sequence": "GG"
            }
        }
    }
    return CanonicalVariation(**params)


@pytest.fixture(scope="module")
def variation_del_lse():
    """Create test fixture for NC_000013.11:20003096:C:"""
    params = {
        "id": "ga4gh:CAN.6_wBT_bhV-hwjaqDxq3kEs3nyILkF4du",
        "type": "CanonicalVariation",
        "canonical_context": {
            "id": "ga4gh:VA.l55oQYOlWUoYwAxb4trpbqmMNaknTa1U",
            "type": "Allele",
            "location": {
                "id": "ga4gh:SL.aEhVhpZMgXXidD4yA4wQbqxKPAf2gnjq",
                "end": {"value": 20003097, "type": "Number"},
                "start": {"value": 20003096, "type": "Number"},
                "sequence_id": "ga4gh:SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT",
                "type": "SequenceLocation"
            },
            "state": {
                "type": "LiteralSequenceExpression",
                "sequence": ""
            }
        }
    }
    return CanonicalVariation(**params)


@pytest.fixture(scope="module")
def variation1_cn_var(variation1_seq_loc):
    """Create test fixture for variation1 represented as copy number count"""
    params = {
        "id": "ga4gh:CAN.COZ3w7stbXexmpHiSIcVkGDXjRs_uJWn",
        "type": "CanonicalVariation",
        "canonical_context": {
            "id": "ga4gh:CN.BEQZ7WHuRThM0k0pKzZEPe3hLBNjMDWN",
            "type": "CopyNumberCount",
            "subject": variation1_seq_loc,
            "copies": {"type": "Number", "value": 2}
        }
    }
    return CanonicalVariation(**params)


@pytest.fixture(scope="module")
def variation1_cx_var(variation1_seq_loc):
    """Create test fixture for variation1 represented as copy number change"""
    params = {
        "id": "ga4gh:CAN.o0zxV2pHejFzHNNNNefrE4EbxlGPN5sU",
        "type": "CanonicalVariation",
        "canonical_context": {
            "id": "ga4gh:CX.mKsge5NlTOn_azljwHIQ7cvdjmyHNI8E",
            "type": "CopyNumberChange",
            "subject": variation1_seq_loc,
            "copy_change": "efo:0030069"
        }
    }
    return CanonicalVariation(**params)


@pytest.fixture(scope="module")
def variation1_rse(variation1_seq_loc):
    """Create test fixture for variation1 represented as RSE"""
    params = {
        "id": "ga4gh:CAN.oizMWSwBdIddFvs_vA8YqxR7YWupBMrF",
        "type": "CanonicalVariation",
        "canonical_context": {
            "id": "ga4gh:VA.91WFk_XWzbzUEI-SfYZip8r1g7I5wqBo",
            "type": "Allele",
            "location": variation1_seq_loc,
            "state": {
                "type": "RepeatedSequenceExpression",
                "seq_expr": {
                    "type": "DerivedSequenceExpression",
                    "location": variation1_seq_loc,
                    "reverse_complement": False
                },
                "count": {
                    "type": "Number",
                    "value": 0
                }
            }
        }
    }
    return CanonicalVariation(**params)


@pytest.fixture(scope="module")
def variation2(braf_v600e_genomic_sub):
    """Create test fixture for NC_000007.14:140753335:A:T"""
    params = {
        "id": "ga4gh:CAN.dP6z4p7SoGJFmlFQcjOQo2d1mXuo1QiY",
        "type": "CanonicalVariation",
        "canonical_context": braf_v600e_genomic_sub
    }
    return CanonicalVariation(**params)


@pytest.fixture(scope="module")
def variation3_lse(grch38_genomic_insertion_variation):
    """Create test fixture for NC_000017.10:g.37880993_37880994insGCTTACGTGATG"""
    params = {
        "id": "ga4gh:CAN.8Pi46FGQsmKIb-6Q0NYKQh0baDBOMvFF",
        "type": "CanonicalVariation",
        "canonical_context": grch38_genomic_insertion_variation
    }
    return CanonicalVariation(**params)


@pytest.fixture(scope="module")
def variation3_cn_var(grch38_genomic_insertion_seq_loc):
    """Create test fixture for variation3 represented as copy number count"""
    params = {
        "id": "ga4gh:CAN.e0xl0wFhoZamh4rMZ9h5F8-bIvrMk2Jp",
        "type": "CanonicalVariation",
        "canonical_context": {
            "id": "ga4gh:CN.AGpJnGGrAb0EyFNYbQm8Xl3ycImaoR0s",
            "type": "CopyNumberCount",
            "subject": grch38_genomic_insertion_seq_loc,
            "copies": {"type": "Number", "value": 2}
        }
    }
    return CanonicalVariation(**params)


@pytest.fixture(scope="module")
def variation3_cx_var(grch38_genomic_insertion_seq_loc):
    """Create test fixture for variation3 represented as copy number change"""
    params = {
        "id": "ga4gh:CAN.cMEZCpwf-Y1830nNv2VWxee-VdOKl-cg",
        "type": "CanonicalVariation",
        "canonical_context": {
            "id": "ga4gh:CX.upJUYyGZaVMvnwsDGJeAC4Ngdf6XIeiJ",
            "type": "CopyNumberChange",
            "subject": grch38_genomic_insertion_seq_loc,
            "copy_change": "efo:0030072"
        }
    }
    return CanonicalVariation(**params)


@pytest.fixture(scope="module")
def variation3_rse(grch38_genomic_insertion_seq_loc):
    """Create test fixture for variation3 represented as RSE"""
    params = {
        "id": "ga4gh:CAN.Ohalqu5SLmNmwjaE00v26KskGfggXwjq",
        "type": "CanonicalVariation",
        "canonical_context": {
            "id": "ga4gh:VA.j6BnT9kvqTO_BQCjTsOzhcnjiwNlhMHv",
            "type": "Allele",
            "location": grch38_genomic_insertion_seq_loc,
            "state": {
                "type": "RepeatedSequenceExpression",
                "seq_expr": {
                    "type": "DerivedSequenceExpression",
                    "location": grch38_genomic_insertion_seq_loc,
                    "reverse_complement": False
                },
                "count": {
                    "type": "Number",
                    "value": 2
                }
            }
        }
    }
    return CanonicalVariation(**params)


@pytest.fixture(scope="module")
def variation4():
    """Create test fixture for NC_000001.11:g.2229202_2229203insCTC"""
    params = {
        "id": "ga4gh:CAN.azUwFImJO4pTH6rSRx6UItCoxJcTvxln",
        "type": "CanonicalVariation",
        "canonical_context": {
            "id": "ga4gh:VA.M9Ekcss52lqr1IoX3wZeLTxVsrPW1MSq",
            "type": "Allele",
            "location": {
                "id": "ga4gh:SL.UdTlDoFUhChBvTu2PXFg_QHA8QsxAy0m",
                "type": "SequenceLocation",
                "sequence_id": "ga4gh:SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
                "start": {"type": "Number", "value": 2229201},
                "end": {"type": "Number", "value": 2229202}
            },
            "state": {"type": "LiteralSequenceExpression", "sequence": "CCTC"}
        }
    }
    return CanonicalVariation(**params)


@pytest.mark.asyncio
async def test_to_canonical_variation_deletion(
    test_handler, variation1_lse, variation1_cn_var, variation1_cx_var,
    variation1_rse, variation_del_lse
):
    """Test that to_canonical_variation works correctly for deletions"""
    # https://www.ncbi.nlm.nih.gov/clinvar/variation/17014/?new_evidence=true
    q = " NC_000013.11:20189346:GGG:GG "  # 38
    resp = await test_handler.to_canonical_variation(q, fmt="spdi")
    assert resp.canonical_variation == variation1_lse
    assert resp.warnings == []

    q = " NC_000013.10:20763485:GGG:GG "  # 37
    resp = await test_handler.to_canonical_variation(q, fmt="spdi", do_liftover=True)
    assert resp.canonical_variation == variation1_lse
    assert resp.warnings == []

    q = " NC_000013.11:g.20189349del "  # 38
    resp = await test_handler.to_canonical_variation(q, fmt="hgvs")
    assert resp.canonical_variation == variation1_lse
    assert resp.warnings == []

    q = " NC_000013.10:g.20763488del "  # 37
    resp = await test_handler.to_canonical_variation(q, fmt="hgvs", do_liftover=True)
    assert resp.canonical_variation == variation1_lse
    assert resp.warnings == []

    # Deletion
    # https://www.ncbi.nlm.nih.gov/clinvar/variation/1373966/?new_evidence=true
    resp = await test_handler.to_canonical_variation(
        "NC_000013.11:20003096:C:", fmt="spdi")
    assert resp.canonical_variation == variation_del_lse
    assert resp.warnings == []

    resp = await test_handler.to_canonical_variation(
        "NC_000013.11:20003096:1:", fmt="spdi")
    assert resp.canonical_variation == variation_del_lse
    assert resp.warnings == []

    resp = await test_handler.to_canonical_variation(
        "NC_000013.11:g.20003097del", fmt="hgvs")
    assert resp.canonical_variation == variation_del_lse
    assert resp.warnings == []

    # Testing hgvs dup del mode params
    resp = await test_handler.to_canonical_variation(
        q, fmt="hgvs", do_liftover=True, hgvs_dup_del_mode="default")
    assert resp.canonical_variation == variation1_lse
    assert resp.warnings == []

    resp = await test_handler.to_canonical_variation(
        q, fmt="hgvs", do_liftover=True, hgvs_dup_del_mode="copy_number_count",
        baseline_copies=3)
    assert resp.canonical_variation == variation1_cn_var
    assert resp.warnings == []

    resp = await test_handler.to_canonical_variation(
        q, fmt="hgvs", do_liftover=True, hgvs_dup_del_mode="copy_number_change",
        copy_change="efo:0030069")
    assert resp.canonical_variation == variation1_cx_var
    assert resp.warnings == []

    resp = await test_handler.to_canonical_variation(
        q, fmt="hgvs", do_liftover=True, hgvs_dup_del_mode="repeated_seq_expr")
    assert resp.canonical_variation == variation1_rse
    assert resp.warnings == []

    resp = await test_handler.to_canonical_variation(
        q, fmt="hgvs", do_liftover=True, hgvs_dup_del_mode="literal_seq_expr")
    assert resp.canonical_variation == variation1_lse
    assert resp.warnings == []


@pytest.mark.asyncio
async def test_to_canonical_variation_substitution(test_handler, variation2):
    """Test that to_canonical_variation works correctly for substitutions"""
    # https://www.ncbi.nlm.nih.gov/clinvar/variation/13961/
    q = "NC_000007.13:140453135:A:T"  # 37
    resp = await test_handler.to_canonical_variation(q, fmt="spdi", do_liftover=True)
    assert resp.canonical_variation == variation2
    assert resp.warnings == []

    q = "NC_000007.14:140753335:A:T"  # 38
    resp = await test_handler.to_canonical_variation(q, fmt="spdi")
    assert resp.canonical_variation == variation2
    assert resp.warnings == []

    q = "NC_000007.13:g.140453136A>T"  # 37
    resp = await test_handler.to_canonical_variation(q, fmt="hgvs", do_liftover=True)
    assert resp.canonical_variation == variation2
    assert resp.warnings == []

    q = "NC_000007.14:g.140753336A>T"  # 38
    resp = await test_handler.to_canonical_variation(q, fmt="hgvs")
    assert resp.canonical_variation == variation2
    assert resp.warnings == []

    # HGVS Dup Del Mode should not affect
    for mode in [m.value for m in HGVSDupDelModeEnum.__members__.values()]:
        resp = await test_handler.to_canonical_variation(
            q, fmt="hgvs", hgvs_dup_del_mode=mode, baseline_copies=2)
        assert resp.canonical_variation == variation2
        assert resp.warnings == []


@pytest.mark.asyncio
async def test_to_canonical_variation_duplication(
    test_handler, variation3_lse, variation3_cn_var, variation3_cx_var,
    variation3_rse
):
    """Test that to_canonical_variation works correctly for duplications"""
    # https://www.ncbi.nlm.nih.gov/clinvar/variation/44985/
    q = "NC_000017.10:g.37880993_37880994insGCTTACGTGATG"  # 37
    resp = await test_handler.to_canonical_variation(
        q, fmt="hgvs", do_liftover=True)
    assert resp.canonical_variation == variation3_lse
    assert resp.warnings == []

    q = "NC_000017.11:g.39724740_39724741insGCTTACGTGATG"  # 38
    resp = await test_handler.to_canonical_variation(
        q, fmt="hgvs", do_liftover=True)  # even tho it's set it wont liftover
    assert resp.canonical_variation == variation3_lse
    assert resp.warnings == []

    q = "NC_000017.11:39724731:TACGTGATGGCT:TACGTGATGGCTTACGTGATGGCT"  # 38
    resp = await test_handler.to_canonical_variation(q, fmt="spdi")
    assert resp.canonical_variation == variation3_lse
    assert resp.warnings == []

    q = "NC_000017.11:g.39724732_39724743dup"  # 38
    resp = await test_handler.to_canonical_variation(q, fmt="hgvs")
    assert resp.canonical_variation == variation3_lse
    assert resp.warnings == []

    q = "NC_000017.10:g.37880985_37880996dup"  # 37
    resp = await test_handler.to_canonical_variation(
        q, fmt="hgvs", do_liftover=True)
    assert resp.canonical_variation == variation3_lse
    assert resp.warnings == []

    # Testing hgvs dup del mode params
    resp = await test_handler.to_canonical_variation(
        q, fmt="hgvs", do_liftover=True, hgvs_dup_del_mode="default")
    assert resp.canonical_variation == variation3_lse
    assert resp.warnings == []

    resp = await test_handler.to_canonical_variation(
        q, fmt="hgvs", do_liftover=True, hgvs_dup_del_mode="copy_number_count",
        baseline_copies=1)
    assert resp.canonical_variation == variation3_cn_var
    assert resp.warnings == []

    resp = await test_handler.to_canonical_variation(
        q, fmt="hgvs", do_liftover=True, hgvs_dup_del_mode="copy_number_change",
        copy_change="efo:0030072")
    assert resp.canonical_variation == variation3_cx_var
    assert resp.warnings == []

    resp = await test_handler.to_canonical_variation(
        q, fmt="hgvs", do_liftover=True, hgvs_dup_del_mode="repeated_seq_expr")
    assert resp.canonical_variation == variation3_rse
    assert resp.warnings == []

    resp = await test_handler.to_canonical_variation(
        q, fmt="hgvs", do_liftover=True, hgvs_dup_del_mode="literal_seq_expr")
    assert resp.canonical_variation == variation3_lse
    assert resp.warnings == []


@pytest.mark.asyncio
async def test_to_canonical_variation_insertions(test_handler, variation4):
    """Test that to_canonical_variation works correctly for duplications"""
    # https://www.ncbi.nlm.nih.gov/clinvar/variation/440270/
    q = "NC_000001.10:g.2160641_2160642insCTC"  # 37
    resp = await test_handler.to_canonical_variation(
        q, fmt="hgvs", do_liftover=True)
    assert resp.canonical_variation == variation4
    assert resp.warnings == []

    q = "NC_000001.11:g.2229202_2229203insCTC"  # 38
    resp = await test_handler.to_canonical_variation(q, fmt="hgvs")
    assert resp.canonical_variation == variation4
    assert resp.warnings == []

    # HGVS Dup Del Mode should not affect
    for mode in [m.value for m in HGVSDupDelModeEnum.__members__.values()]:
        resp = await test_handler.to_canonical_variation(
            q, fmt="hgvs", hgvs_dup_del_mode=mode, baseline_copies=2)
        assert resp.canonical_variation == variation4
        assert resp.warnings == []

    q = "NC_000001.10:2160640:C:CCTC"  # 37
    resp = await test_handler.to_canonical_variation(
        q, fmt="spdi", do_liftover=True)
    assert resp.canonical_variation == variation4
    assert resp.warnings == []

    q = "NC_000001.11:2229201:C:CCTC"  # 38
    resp = await test_handler.to_canonical_variation(q, fmt="spdi")
    assert resp.canonical_variation == variation4
    assert resp.warnings == []


@pytest.mark.asyncio
async def test_invalid(test_handler):
    """Test that invalid queries return the correct response"""
    resp = await test_handler.to_canonical_variation(
        "NC_000013.11:201845654659346:GGG:GG", fmt="spdi",
        untranslatable_returns_text=True)
    assert resp.canonical_variation.canonical_context.type == "Text"
    assert resp.warnings == ["start out of range (201845654659346)"]

    resp = await test_handler.to_canonical_variation(
        "NC_000013.11:2018459346:GGG:GG", fmt="spdi", untranslatable_returns_text=True)
    assert resp.canonical_variation.canonical_context.type == "Text"
    assert resp.warnings == ["Position, 2018459346, does not exist on NC_000013.11"]

    resp = await test_handler.to_canonical_variation(
        "NC_000013.1:20189346:GGG:GG", fmt="spdi", untranslatable_returns_text=True)
    assert resp.canonical_variation.canonical_context.type == "Text"
    assert resp.warnings == ["vrs-python translator raised error: seqrepo could not "
                             "translate identifier 'refseq:NC_000013.1'"]

    resp = await test_handler.to_canonical_variation(
        "NP_004324.2:p.Val600Glu", fmt="spdi", untranslatable_returns_text=True)
    assert resp.canonical_variation.canonical_context.type == "Text"
    assert resp.warnings == ["NP_004324.2:p.Val600Glu is not a valid SPDI expression"]

    resp = await test_handler.to_canonical_variation(
        "NC_000013.11:20189346:GCG:GG", fmt="spdi", untranslatable_returns_text=True)
    assert resp.canonical_variation.canonical_context.type == "Text"
    assert resp.warnings ==\
           ["Expected to find reference sequence GCG but found GGG on NC_000013.11"]

    resp = await test_handler.to_canonical_variation(
        "NC_000007.14:140753335:A:T", fmt="hgvs")
    assert resp.canonical_variation is None
    assert resp.warnings == \
        ["NC_000007.14:140753335:A:T is not a valid HGVS expression"]

    resp = await test_handler.to_canonical_variation(
        "NC_000007.14:g.140753336464564654A>T", fmt="hgvs")
    assert resp.canonical_variation is None
    assert resp.warnings == ["Unable to find valid result for classifications:"
                             " {'genomic substitution'}"]

    q = "NC_000001.10:2160640:A:ACTC"
    resp = await test_handler.to_canonical_variation(
        q, fmt="spdi", do_liftover=True)
    assert resp.canonical_variation is None
    assert resp.warnings == \
        ["Expected to find reference sequence A but found C on NC_000001.11"]

    q = " NC_000013.11:g.20189349del "  # 38
    resp = await test_handler.to_canonical_variation(
        q, fmt="hgvs", hgvs_dup_del_mode="copy_number_count")
    assert resp.canonical_variation is None
    assert resp.warnings == ["copy_number_count requires `baseline_copies`"]
