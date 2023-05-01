"""Module for testing Amplification to Copy Number Change"""
import pytest
from ga4gh.vrsatile.pydantic.vrs_models import CopyNumberChange


@pytest.fixture(scope="module")
def kit_amplification():
    """Create test fixture for KIT amplification"""
    params = {
        "type": "CopyNumberChange",
        "_id": "ga4gh:CX.lMc3ItCK6GqTa1ijgmKiF7OA3i5kkw3N",
        "copy_change": "efo:0030072",
        "subject": {
            "type": "SequenceLocation",
            "_id": "ga4gh:VSL.MHuTsdf9mtx0OCn4mcv6KwJily8vIS2q",
            "sequence_id": "ga4gh:SQ.iy7Zfceb5_VGtTQzJ-v5JpPbpeifHD_V",
            "interval": {
                "start": {"type": "Number", "value": 55599320},
                "end": {"type": "Number", "value": 55599321},
                "type": "SequenceInterval"
            },
            "type": "SequenceLocation"
        }
    }
    return CopyNumberChange(**params)


def test_amplification_to_cx_cnv(test_cnv_handler, braf_amplification,
                                 prpf8_amplification, kit_amplification):
    """Test that amplification_to_cx_cnv method works correctly"""
    # Using gene normalizer
    resp = test_cnv_handler.amplification_to_cx_cnv(gene="braf")
    assert resp.copy_number_change == braf_amplification.variation
    assert resp.amplification_label == "BRAF Amplification"
    assert resp.warnings == []

    # Gene with > 1 sequence location
    resp = test_cnv_handler.amplification_to_cx_cnv(gene="PRPF8")
    assert resp.copy_number_change == prpf8_amplification.variation
    assert resp.amplification_label == "PRPF8 Amplification"
    assert resp.warnings == []

    # Gene with no location. This should NOT return a variation
    resp = test_cnv_handler.amplification_to_cx_cnv(gene="ifnr")
    assert resp.copy_number_change is None
    assert resp.amplification_label == "IFNR Amplification"
    assert resp.warnings == ["gene-normalizer could not find a priority sequence "
                             "location for gene: IFNR"]

    # Using sequence_id, start, end
    resp = test_cnv_handler.amplification_to_cx_cnv(
        gene="KIT", sequence_id="NC_000004.11", start=55599321, end=55599321)
    assert resp.copy_number_change == kit_amplification
    assert resp.amplification_label == "KIT Amplification"
    assert resp.warnings == []

    # Sequence_id not found in seqrepo
    resp = test_cnv_handler.amplification_to_cx_cnv(
        gene="BRAF", sequence_id="NC_000007", start=140453136, end=140453136)
    assert resp.copy_number_change is None
    assert resp.amplification_label == "BRAF Amplification"
    assert resp.warnings == ["SeqRepo unable to get translated identifiers for "
                             "NC_000007"]

    # pos not on valid sequence_id
    resp = test_cnv_handler.amplification_to_cx_cnv(
        gene="braf", sequence_id="NC_000007.13", start=55599321, end=9955599321)
    assert resp.copy_number_change is None
    assert resp.amplification_label == "BRAF Amplification"
    assert resp.warnings == ["End inter-residue coordinate (9955599320) is out of "
                             "index on NC_000007.13"]

    # invalid gene
    resp = test_cnv_handler.amplification_to_cx_cnv(gene="invalid")
    assert resp.copy_number_change is None
    assert resp.amplification_label is None
    assert resp.warnings == ["gene-normalizer returned no match for gene: invalid"]
