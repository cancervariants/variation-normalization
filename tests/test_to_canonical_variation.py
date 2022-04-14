"""Module for testing to_canonical_variation"""
import copy

import pytest
from ga4gh.vrsatile.pydantic.vrsatile_models import CanonicalVariation


@pytest.fixture(scope="module")
def variation1():
    """Create test fixture for NC_000013.11:20189346:GGG:GG"""
    params = {
        "_id": "ga4gh:VCC.6FxWtQdkEyVSMIOnvRj0bOEgvHgN3pRh",
        "complement": False,
        "type": "CanonicalVariation",
        "variation": {
            "_id": "ga4gh:VA.KGopzor-bEw8Ot5sAQQ5o5SVx4o7TuLN",
            "type": "Allele",
            "location": {
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
            },
            "state": {
                "type": "LiteralSequenceExpression",
                "sequence": "GG"
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
def variation3(grch38_genomic_insertion_variation):
    """Create test fixture for NC_000017.10:g.37880993_37880994insGCTTACGTGATG"""
    params = {
        "_id": "ga4gh:VCC.Wvq61Rejsn80kUCJkqSXKzCSnRwaRZNm",
        "complement": False,
        "type": "CanonicalVariation",
        "variation": grch38_genomic_insertion_variation
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
async def test_to_canonical_variation_deletion(test_query_handler, variation1):
    """Test that to_canonical_variation works correctly for deletions"""
    # https://www.ncbi.nlm.nih.gov/clinvar/variation/17014/?new_evidence=true
    q = " NC_000013.11:20189346:GGG:GG "  # 38
    resp, w = await test_query_handler.to_canonical_variation(q, fmt="spdi")
    assert resp == variation1
    assert w == []

    q = " NC_000013.10:20763485:GGG:GG "  # 37
    resp, w = await test_query_handler.to_canonical_variation(
        q, fmt="spdi", do_liftover=True)
    assert resp == variation1
    assert w == []

    q = " NC_000013.11:g.20189349del "  # 38
    resp, w = await test_query_handler.to_canonical_variation(q, fmt="hgvs")
    assert resp == variation1
    assert w == []

    q = " NC_000013.10:g.20763488del "  # 37
    resp, w = await test_query_handler.to_canonical_variation(
        q, fmt="hgvs", do_liftover=True)
    assert resp == variation1
    assert w == []


@pytest.mark.asyncio
async def test_to_canonical_variation_substitution(test_query_handler, variation2):
    """Test that to_canonical_variation works correctly for substitutions"""
    # https://www.ncbi.nlm.nih.gov/clinvar/variation/13961/
    q = "NC_000007.13:140453135:A:T"  # 37
    resp, w = await test_query_handler.to_canonical_variation(
        q, fmt="spdi", complement=True, do_liftover=True)
    assert resp == variation2
    assert w == []

    q = "NC_000007.14:140753335:A:T"  # 38
    resp, w = await test_query_handler.to_canonical_variation(
        q, fmt="spdi", complement=True)
    assert resp == variation2
    assert w == []

    resp, w = await test_query_handler.to_canonical_variation(
        q, fmt="spdi", complement=False)
    cpy_variation2 = copy.deepcopy(variation2)
    cpy_variation2.complement = False
    cpy_variation2.id = "ga4gh:VCC.W0r_NF_ecKXjgvTwcMNkyVS1pB_CXMj9"
    assert resp == cpy_variation2
    assert w == []

    q = "NC_000007.13:g.140453136A>T"  # 37
    resp, w = await test_query_handler.to_canonical_variation(
        q, fmt="hgvs", complement=True, do_liftover=True)
    assert resp == variation2
    assert w == []

    q = "NC_000007.14:g.140753336A>T"  # 38
    resp, w = await test_query_handler.to_canonical_variation(
        q, fmt="hgvs", complement=True)
    assert resp == variation2
    assert w == []

    resp, w = await test_query_handler.to_canonical_variation(
        q, fmt="hgvs", complement=False)
    cpy_variation2 = copy.deepcopy(variation2)
    cpy_variation2.complement = False
    cpy_variation2.id = "ga4gh:VCC.W0r_NF_ecKXjgvTwcMNkyVS1pB_CXMj9"
    assert resp == cpy_variation2
    assert w == []


@pytest.mark.asyncio
async def test_to_canonical_variation_duplication(test_query_handler, variation3):
    """Test that to_canonical_variation works correctly for duplications"""
    # https://www.ncbi.nlm.nih.gov/clinvar/variation/44985/
    q = "NC_000017.10:g.37880993_37880994insGCTTACGTGATG"  # 37
    resp, w = await test_query_handler.to_canonical_variation(
        q, fmt="hgvs", do_liftover=True)
    assert resp == variation3
    assert w == []

    q = "NC_000017.11:g.39724740_39724741insGCTTACGTGATG"  # 38
    resp, w = await test_query_handler.to_canonical_variation(
        q, fmt="hgvs", do_liftover=True)  # even tho it's set it wont liftover
    assert resp == variation3
    assert w == []

    q = "NC_000017.11:39724731:TACGTGATGGCT:TACGTGATGGCTTACGTGATGGCT"  # 38
    resp, w = await test_query_handler.to_canonical_variation(q, fmt="spdi")
    assert resp == variation3
    assert w == []

    q = "NC_000017.11:g.39724732_39724743dup"  # 38
    resp, w = await test_query_handler.to_canonical_variation(q, fmt="hgvs")
    assert resp == variation3
    assert w == []

    q = "NC_000017.10:g.37880985_37880996dup"  # 37
    resp, w = await test_query_handler.to_canonical_variation(
        q, fmt="hgvs", do_liftover=True)
    assert resp == variation3
    assert w == []


@pytest.mark.asyncio
async def test_to_canonical_variation_insertions(test_query_handler, variation4):
    """Test that to_canonical_variation works correctly for duplications"""
    # https://www.ncbi.nlm.nih.gov/clinvar/variation/440270/
    q = "NC_000001.10:g.2160641_2160642insCTC"  # 37
    resp, w = await test_query_handler.to_canonical_variation(
        q, fmt="hgvs", do_liftover=True)
    assert resp.dict() == variation4.dict()
    assert w == []

    q = "NC_000001.11:g.2229202_2229203insCTC"  # 37
    resp, w = await test_query_handler.to_canonical_variation(q, fmt="hgvs")
    assert resp.dict() == variation4.dict()
    assert w == []

    q = "NC_000001.10:2160640:C:CCTC"  # 37
    resp, w = await test_query_handler.to_canonical_variation(
        q, fmt="spdi", do_liftover=True)
    assert resp.dict() == variation4.dict()
    assert w == []

    q = "NC_000001.11:2229201:C:CCTC"  # 38
    resp, w = await test_query_handler.to_canonical_variation(q, fmt="spdi")
    assert resp.dict() == variation4.dict()
    assert w == []


@pytest.mark.asyncio
async def test_invalid(test_query_handler):
    """Test that invalid queries return the correct response"""
    resp, w = await test_query_handler.to_canonical_variation(
        "NC_000013.11:201845654659346:GGG:GG", fmt="spdi")
    assert resp.variation.type == "Text"
    assert w == ["start out of range (201845654659346)"]

    resp, w = await test_query_handler.to_canonical_variation(
        "NC_000013.11:2018459346:GGG:GG", fmt="spdi")
    assert resp.variation.type == "Text"
    assert w == ["Position, 2018459346, does not exist on NC_000013.11"]

    resp, w = await test_query_handler.to_canonical_variation(
        "NC_000013.1:20189346:GGG:GG", fmt="spdi")
    assert resp.variation.type == "Text"
    assert w == ["vrs-python translator raised error: seqrepo could not translate "
                 "identifier 'refseq:NC_000013.1'"]

    resp, w = await test_query_handler.to_canonical_variation(
        "NP_004324.2:p.Val600Glu", fmt="spdi")
    assert resp.variation.type == "Text"
    assert w == ["NP_004324.2:p.Val600Glu is not a valid SPDI expression"]

    resp, w = await test_query_handler.to_canonical_variation(
        "NC_000013.11:20189346:GCG:GG", fmt="spdi")
    assert resp.variation.type == "Text"
    assert w ==\
           ["Expected to find reference sequence GCG but found GGG on NC_000013.11"]

    resp, w = await test_query_handler.to_canonical_variation(
        "NC_000007.14:140753335:A:T", fmt="hgvs")
    assert resp.variation.type == "Text"
    assert w == ["Unable to tokenize 140753335"]

    resp, w = await test_query_handler.to_canonical_variation(
        "NC_000007.14:g.140753336464564654A>T", fmt="hgvs")
    assert resp.variation.type == "Text"
    assert w == ["Unable to find valid result for classifications:"
                 " {'genomic substitution'}"]

    q = "NC_000001.10:2160640:A:ACTC"
    resp, w = await test_query_handler.to_canonical_variation(
        q, fmt="spdi", do_liftover=True)
    assert resp.variation.type == "Text"
    assert w == ["Expected to find reference sequence A but found C on NC_000001.11"]
