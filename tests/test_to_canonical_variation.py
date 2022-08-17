"""Module for testing to_canonical_variation"""
import copy

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
        "_id": "ga4gh:VSL.p4e9kMEY9PrKZ1BbNRuFr6n30DkwXWlX",
        "type": "SequenceLocation",
        "sequence_id": "ga4gh:SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT",
        "interval": {
            "type": "SequenceInterval",
            "start": {
                "type": "Number",
                "value": 20189346
            },
            "end": {
                "type": "Number",
                "value": 20189349
            }
        }
    }


@pytest.fixture(scope="module")
def variation1_lse(variation1_seq_loc):
    """Create test fixture for NC_000013.11:20189346:GGG:GG"""
    params = {
        "_id": "ga4gh:VCC.6FxWtQdkEyVSMIOnvRj0bOEgvHgN3pRh",
        "complement": False,
        "type": "CanonicalVariation",
        "variation": {
            "_id": "ga4gh:VA.KGopzor-bEw8Ot5sAQQ5o5SVx4o7TuLN",
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
def variation1_abs_cnv(variation1_seq_loc):
    """Create test fixture for variation1 represented as absolute cnv"""
    params = {
        "_id": "ga4gh:VCC.SRZiwHEKF9qJ5hCd1hKPDpjBpK29PO_4",
        "complement": False,
        "type": "CanonicalVariation",
        "variation": {
            "_id": "ga4gh:VAC.WkHL0pE5vn2yZfTlB4vo6Qr9gvFgnm2s",
            "type": "AbsoluteCopyNumber",
            "subject": variation1_seq_loc,
            "copies": {"type": "Number", "value": 2}
        }
    }
    return CanonicalVariation(**params)


@pytest.fixture(scope="module")
def variation1_rel_cnv(variation1_seq_loc):
    """Create test fixture for variation1 represented as relative cnv"""
    params = {
        "_id": "ga4gh:VCC.hpdklR_18hDE3YwHTYYFJmbQd__pwBeH",
        "complement": False,
        "type": "CanonicalVariation",
        "variation": {
            "_id": "ga4gh:VRC.d0tA1C6eSo-rkg1wVQ76Bh8XLSwz36se",
            "type": "RelativeCopyNumber",
            "subject": variation1_seq_loc,
            "relative_copy_class": "complete loss"
        }
    }
    return CanonicalVariation(**params)


@pytest.fixture(scope="module")
def variation1_rse(variation1_seq_loc):
    """Create test fixture for variation1 represented as RSE"""
    params = {
        "_id": "ga4gh:VCC.8bwdX-_oSSoT36vvA2Z6CVG_Mfw-jqEO",
        "complement": False,
        "type": "CanonicalVariation",
        "variation": {
            "_id": "ga4gh:VA.0TP_vHDqgjDm2Wme2U6f6upSlYxAqw4l",
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
        "_id": "ga4gh:VCC.jh-D_L4J74__OdCUI8bn6j4cI1RWZIAA",
        "complement": True,
        "type": "CanonicalVariation",
        "variation": braf_v600e_genomic_sub
    }
    return CanonicalVariation(**params)


@pytest.fixture(scope="module")
def variation3_lse(grch38_genomic_insertion_variation):
    """Create test fixture for NC_000017.10:g.37880993_37880994insGCTTACGTGATG"""
    params = {
        "_id": "ga4gh:VCC.Wvq61Rejsn80kUCJkqSXKzCSnRwaRZNm",
        "complement": False,
        "type": "CanonicalVariation",
        "variation": grch38_genomic_insertion_variation
    }
    return CanonicalVariation(**params)


@pytest.fixture(scope="module")
def variation3_abs_cnv(grch38_genomic_insertion_seq_loc):
    """Create test fixture for variation3 represented as absolute cnv"""
    params = {
        "_id": "ga4gh:VCC.DSncpwIaoyP94KOVf_45gcQGJ0_f41Pa",
        "complement": False,
        "type": "CanonicalVariation",
        "variation": {
            "_id": "ga4gh:VAC.c1xp6z9seT4j8Ae4e0gePB7LiPzjykmO",
            "type": "AbsoluteCopyNumber",
            "subject": grch38_genomic_insertion_seq_loc,
            "copies": {"type": "Number", "value": 2}
        }
    }
    return CanonicalVariation(**params)


@pytest.fixture(scope="module")
def variation3_rel_cnv(grch38_genomic_insertion_seq_loc):
    """Create test fixture for variation3 represented as relative cnv"""
    params = {
        "_id": "ga4gh:VCC.fUeN_b2i0iBgad1rizfs1fXOnrZvptc9",
        "complement": False,
        "type": "CanonicalVariation",
        "variation": {
            "_id": "ga4gh:VRC.TJSKDRpI2hW8JbqGsc9yag6KrP1YVhwU",
            "type": "RelativeCopyNumber",
            "subject": grch38_genomic_insertion_seq_loc,
            "relative_copy_class": "high-level gain"
        }
    }
    return CanonicalVariation(**params)


@pytest.fixture(scope="module")
def variation3_rse(grch38_genomic_insertion_seq_loc):
    """Create test fixture for variation3 represented as RSE"""
    params = {
        "_id": "ga4gh:VCC.2Dr1zKcf6_tmFlcQG_1ZgzjFsUMuSVf3",
        "complement": False,
        "type": "CanonicalVariation",
        "variation": {
            "_id": "ga4gh:VA.Psy669nQn5V84LRrECrl4DKLt4V64_5h",
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
        "_id": "ga4gh:VCC.DSkrWcvpWq3U1KkxLDXauJs7BmbOZVBA",
        "complement": False,
        "type": "CanonicalVariation",
        "variation": {
            "_id": "ga4gh:VA.PZ-BIW_H8sMhl65ZU4bgn_99XGirUPpz",
            "type": "Allele",
            "location": {
                "_id": "ga4gh:VSL.ivN4CeEpVq7IPlav6JWMx5cJZ8vO3NKz",
                "type": "SequenceLocation",
                "sequence_id": "ga4gh:SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
                "interval": {
                    "type": "SequenceInterval",
                    "start": {"type": "Number", "value": 2229201},
                    "end": {"type": "Number", "value": 2229202}
                }
            },
            "state": {"type": "LiteralSequenceExpression", "sequence": "CCTC"}
        }
    }
    return CanonicalVariation(**params)


@pytest.mark.asyncio
async def test_to_canonical_variation_deletion(
    test_handler, variation1_lse, variation1_abs_cnv, variation1_rel_cnv,
    variation1_rse
):
    """Test that to_canonical_variation works correctly for deletions"""
    # https://www.ncbi.nlm.nih.gov/clinvar/variation/17014/?new_evidence=true
    q = " NC_000013.11:20189346:GGG:GG "  # 38
    resp = await test_handler.to_canonical_variation(q, fmt="spdi")
    assert resp.canonical_variation == variation1_lse
    assert resp.warnings == []

    q = " NC_000013.10:20763485:GGG:GG "  # 37
    resp = await test_handler.to_canonical_variation(
        q, fmt="spdi", do_liftover=True)
    assert resp.canonical_variation == variation1_lse
    assert resp.warnings == []

    q = " NC_000013.11:g.20189349del "  # 38
    resp = await test_handler.to_canonical_variation(q, fmt="hgvs")
    assert resp.canonical_variation == variation1_lse
    assert resp.warnings == []

    q = " NC_000013.10:g.20763488del "  # 37
    resp = await test_handler.to_canonical_variation(
        q, fmt="hgvs", do_liftover=True)
    assert resp.canonical_variation == variation1_lse
    assert resp.warnings == []

    # Testing hgvs dup del mode params
    resp = await test_handler.to_canonical_variation(
        q, fmt="hgvs", do_liftover=True, hgvs_dup_del_mode="default")
    assert resp.canonical_variation == variation1_lse
    assert resp.warnings == []

    resp = await test_handler.to_canonical_variation(
        q, fmt="hgvs", do_liftover=True, hgvs_dup_del_mode="absolute_cnv",
        baseline_copies=3)
    assert resp.canonical_variation == variation1_abs_cnv
    assert resp.warnings == []

    resp = await test_handler.to_canonical_variation(
        q, fmt="hgvs", do_liftover=True, hgvs_dup_del_mode="relative_cnv",
        relative_copy_class="complete loss")
    assert resp.canonical_variation == variation1_rel_cnv
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
    resp = await test_handler.to_canonical_variation(
        q, fmt="spdi", complement=True, do_liftover=True)
    assert resp.canonical_variation == variation2
    assert resp.warnings == []

    q = "NC_000007.14:140753335:A:T"  # 38
    resp = await test_handler.to_canonical_variation(
        q, fmt="spdi", complement=True)
    assert resp.canonical_variation == variation2
    assert resp.warnings == []

    resp = await test_handler.to_canonical_variation(
        q, fmt="spdi", complement=False)
    cpy_variation2 = copy.deepcopy(variation2)
    cpy_variation2.complement = False
    cpy_variation2.id = "ga4gh:VCC.W0r_NF_ecKXjgvTwcMNkyVS1pB_CXMj9"
    assert resp.canonical_variation == cpy_variation2
    assert resp.warnings == []

    q = "NC_000007.13:g.140453136A>T"  # 37
    resp = await test_handler.to_canonical_variation(
        q, fmt="hgvs", complement=True, do_liftover=True)
    assert resp.canonical_variation == variation2
    assert resp.warnings == []

    q = "NC_000007.14:g.140753336A>T"  # 38
    resp = await test_handler.to_canonical_variation(
        q, fmt="hgvs", complement=True)
    assert resp.canonical_variation == variation2
    assert resp.warnings == []

    resp = await test_handler.to_canonical_variation(
        q, fmt="hgvs", complement=False)
    cpy_variation2 = copy.deepcopy(variation2)
    cpy_variation2.complement = False
    cpy_variation2.id = "ga4gh:VCC.W0r_NF_ecKXjgvTwcMNkyVS1pB_CXMj9"
    assert resp.canonical_variation == cpy_variation2
    assert resp.warnings == []

    # HGVS Dup Del Mode should not affect
    for mode in [m.value for m in HGVSDupDelModeEnum.__members__.values()]:
        resp = await test_handler.to_canonical_variation(
            q, fmt="hgvs", complement=True, hgvs_dup_del_mode=mode,
            baseline_copies=2)
        assert resp.canonical_variation == variation2
        assert resp.warnings == []


@pytest.mark.asyncio
async def test_to_canonical_variation_duplication(
    test_handler, variation3_lse, variation3_abs_cnv, variation3_rel_cnv,
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
        q, fmt="hgvs", do_liftover=True, hgvs_dup_del_mode="absolute_cnv",
        baseline_copies=1)
    assert resp.canonical_variation == variation3_abs_cnv
    assert resp.warnings == []

    resp = await test_handler.to_canonical_variation(
        q, fmt="hgvs", do_liftover=True, hgvs_dup_del_mode="relative_cnv",
        relative_copy_class="high-level gain")
    assert resp.canonical_variation == variation3_rel_cnv
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
    assert resp.canonical_variation.variation.type == "Text"
    assert resp.warnings == ["start out of range (201845654659346)"]

    resp = await test_handler.to_canonical_variation(
        "NC_000013.11:2018459346:GGG:GG", fmt="spdi", untranslatable_returns_text=True)
    assert resp.canonical_variation.variation.type == "Text"
    assert resp.warnings == ["Position, 2018459346, does not exist on NC_000013.11"]

    resp = await test_handler.to_canonical_variation(
        "NC_000013.1:20189346:GGG:GG", fmt="spdi", untranslatable_returns_text=True)
    assert resp.canonical_variation.variation.type == "Text"
    assert resp.warnings == ["vrs-python translator raised error: seqrepo could not "
                             "translate identifier 'refseq:NC_000013.1'"]

    resp = await test_handler.to_canonical_variation(
        "NP_004324.2:p.Val600Glu", fmt="spdi", untranslatable_returns_text=True)
    assert resp.canonical_variation.variation.type == "Text"
    assert resp.warnings == ["NP_004324.2:p.Val600Glu is not a valid SPDI expression"]

    resp = await test_handler.to_canonical_variation(
        "NC_000013.11:20189346:GCG:GG", fmt="spdi", untranslatable_returns_text=True)
    assert resp.canonical_variation.variation.type == "Text"
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
        q, fmt="hgvs", hgvs_dup_del_mode="absolute_cnv")
    assert resp.canonical_variation is None
    assert resp.warnings == ["absolute_cnv requires `baseline_copies`"]
