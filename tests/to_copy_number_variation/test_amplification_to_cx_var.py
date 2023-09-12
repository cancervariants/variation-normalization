"""Module for testing Amplification to Copy Number Change"""
import pytest
from ga4gh.vrs import models


@pytest.fixture(scope="module")
def kit_amplification():
    """Create test fixture for KIT amplification"""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.wQv1KnYyhMd1aKoXFrOVzT3rMNvo0OIS",
        "copyChange": "efo:0030072",
        "location": {
            "type": "SequenceLocation",
            "id": "ga4gh:SL.5UgZnBz5pAVUWzNMyC1YJBeVnAA_DGUE",
            "sequenceReference": {
                "type": "SequenceReference",
                "refgetAccession": "SQ.iy7Zfceb5_VGtTQzJ-v5JpPbpeifHD_V",
            },
            "start": 55599320,
            "end": 55599321,
        },
    }
    return models.CopyNumberChange(**params)


def test_amplification_to_cx_var(
    test_cnv_handler, braf_amplification, prpf8_amplification, kit_amplification
):
    """Test that amplification_to_cx_var method works correctly"""
    # Using gene normalizer
    resp = test_cnv_handler.amplification_to_cx_var(gene="braf")
    assert resp.copy_number_change.model_dump(
        exclude_none=True
    ) == braf_amplification.model_dump(exclude_none=True)
    assert resp.amplification_label == "BRAF Amplification"
    assert resp.warnings == []

    # Gene with > 1 sequence location
    resp = test_cnv_handler.amplification_to_cx_var(gene="PRPF8")
    assert resp.copy_number_change.model_dump(
        exclude_none=True
    ) == prpf8_amplification.model_dump(exclude_none=True)
    assert resp.amplification_label == "PRPF8 Amplification"
    assert resp.warnings == []

    # Gene with no location. This should NOT return a variation
    resp = test_cnv_handler.amplification_to_cx_var(gene="ifnr")
    assert resp.copy_number_change is None
    assert resp.amplification_label == "IFNR Amplification"
    assert resp.warnings == [
        "gene-normalizer could not find a priority sequence " "location for gene: IFNR"
    ]

    # Using sequence, start, end
    resp = test_cnv_handler.amplification_to_cx_var(
        gene="KIT", sequence="NC_000004.11", start=55599321, end=55599321
    )
    assert resp.copy_number_change.model_dump(
        exclude_none=True
    ) == kit_amplification.model_dump(exclude_none=True)
    assert resp.amplification_label == "KIT Amplification"
    assert resp.warnings == []

    # Sequence_id not found in seqrepo
    resp = test_cnv_handler.amplification_to_cx_var(
        gene="BRAF", sequence="NC_000007", start=140453136, end=140453136
    )
    assert resp.copy_number_change is None
    assert resp.amplification_label == "BRAF Amplification"
    assert resp.warnings == [
        "SeqRepo unable to get translated identifiers for " "NC_000007"
    ]

    # pos not on valid sequence
    resp = test_cnv_handler.amplification_to_cx_var(
        gene="braf", sequence="NC_000007.13", start=55599321, end=9955599321
    )
    assert resp.copy_number_change is None
    assert resp.amplification_label == "BRAF Amplification"
    assert resp.warnings == [
        "End inter-residue coordinate (9955599320) is out of " "index on NC_000007.13"
    ]

    # invalid gene
    resp = test_cnv_handler.amplification_to_cx_var(gene="invalid")
    assert resp.copy_number_change is None
    assert resp.amplification_label is None
    assert resp.warnings == ["gene-normalizer returned no match for gene: invalid"]
