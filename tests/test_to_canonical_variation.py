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


@pytest.mark.asyncio
async def test_canonical_spdi_to_categorical_variation(test_query_handler, variation1,
                                                       variation2):
    """Test that to_canonical_variation works correctly"""
    # https://www.ncbi.nlm.nih.gov/clinvar/variation/17014/?new_evidence=true
    q = " NC_000013.11:20189346:GGG:GG "
    resp, w = await test_query_handler.to_canonical_variation(q, fmt="spdi")
    assert resp == variation1
    assert w == []

    q = " NC_000013.11:g.20189349del "
    resp, w = await test_query_handler.to_canonical_variation(q, fmt="hgvs")
    assert resp == variation1
    assert w == []

    # https://www.ncbi.nlm.nih.gov/clinvar/variation/13961/
    q = "NC_000007.14:140753335:A:T"
    resp, w = await test_query_handler.to_canonical_variation(
        q, fmt="spdi", complement=True)
    assert resp == variation2

    resp, w = await test_query_handler.to_canonical_variation(
        q, fmt="spdi", complement=False)
    cpy_variation2 = copy.deepcopy(variation2)
    cpy_variation2.complement = False
    cpy_variation2.id = "ga4gh:VCC.W0r_NF_ecKXjgvTwcMNkyVS1pB_CXMj9"
    assert resp == cpy_variation2

    q = "NC_000007.14:g.140753336A>T"
    resp, w = await test_query_handler.to_canonical_variation(
        q, fmt="hgvs", complement=True)
    assert resp == variation2

    resp, w = await test_query_handler.to_canonical_variation(
        q, fmt="hgvs", complement=False)
    cpy_variation2 = copy.deepcopy(variation2)
    cpy_variation2.complement = False
    cpy_variation2.id = "ga4gh:VCC.W0r_NF_ecKXjgvTwcMNkyVS1pB_CXMj9"
    assert resp == cpy_variation2


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
    assert w == ["vrs-python translator raised error: Unable to parse data as "
                 "spdi variation"]

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
