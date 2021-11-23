"""Create methods used throughout tests."""
import pytest


@pytest.fixture(scope='session')
def vhl_gene_context():
    """Create a VHL gene context."""
    return {
        "id": "normalize.gene:VHL",
        "type": "GeneDescriptor",
        "label": "VHL",
        "gene_id": "hgnc:12687",
        "xrefs": [
            "ncbigene:7428",
            "ensembl:ENSG00000134086"
        ],
        "alternate_labels": [
            "HRCA1",
            "VHL1",
            "RCA1",
            "pVHL"
        ],
        "extensions": [
            {
                "type": "Extension",
                "name": "symbol_status",
                "value": "approved"
            },
            {
                "type": "Extension",
                "name": "approved_name",
                "value": "von Hippel-Lindau tumor suppressor"
            },
            {
                "type": "Extension",
                "name": "associated_with",
                "value": [
                    "ucsc:uc003bvc.4",
                    "pubmed:9671762",
                    "refseq:NM_000551",
                    "cosmic:VHL",
                    "omim:608537",
                    "vega:OTTHUMG00000128668",
                    "ccds:CCDS2598",
                    "ena.embl:L15409",
                    "orphanet:120467",
                    "ccds:CCDS2597",
                    "uniprot:P40337"
                ]
            },
            {
                "type": "Extension",
                "name": "chromosome_location",
                "value": {
                    "_id":
                        "ga4gh:VCL.S-TtMfLdsgZPVRrWEf1-jiZMyTDCt5y1",
                    "type": "ChromosomeLocation",
                    "species_id": "taxonomy:9606",
                    "chr": "3",
                    "interval": {
                        "end": "p25.3",
                        "start": "p25.3",
                        "type": "CytobandInterval"
                    }
                }
            },
            {
                "name": "previous_symbols",
                "value": [
                    "RCA1"
                ],
                "type": "Extension"
            }
        ]
    }


def assertion_checks(normalize_response, test_variation, ignore_id=False):
    """Check that normalize_response and test_variation are equal."""
    if not ignore_id:
        assert normalize_response.id == test_variation.id, "id"
    assert normalize_response.type == test_variation.type, "type"
    assert normalize_response.variation_id == \
           test_variation.variation_id, "variation_id"
    if test_variation.variation.type != "Text":
        if test_variation.variation.id:
            assert normalize_response.variation.id == \
                   test_variation.variation.id, "variation._id"
            if test_variation.variation_id:
                assert normalize_response.variation_id == \
                       normalize_response.variation.id, "variation_id == variation.id"  # noqa: E501
        assert normalize_response.variation == \
               test_variation.variation, "variation"
    else:
        if not ignore_id:
            assert normalize_response.variation.id == \
                   test_variation.variation.id
        assert normalize_response.variation.type == \
               test_variation.variation.type
        assert normalize_response.variation.definition == \
               test_variation.variation.definition
    assert normalize_response.molecule_context == \
           test_variation.molecule_context, "molecule_context"
    assert normalize_response.structural_type == \
           test_variation.structural_type, "structural_type"
    assert normalize_response.vrs_ref_allele_seq == \
           test_variation.vrs_ref_allele_seq, "vrs_ref_allele_seq"

    resp_gene_context = normalize_response.gene_context
    test_variation_context = test_variation.gene_context
    if resp_gene_context:
        assert resp_gene_context.id == \
               test_variation_context.id, "gene_context.id"
        assert resp_gene_context.label == \
               test_variation_context.label, "gene_context.label"
        assert resp_gene_context.gene_id ==\
               test_variation_context.gene_id, "gene_context.gene_id"
        assert set(resp_gene_context.xrefs) ==\
               set(test_variation_context.xrefs), "gene_context.xrefs"
        if test_variation_context.alternate_labels:
            assert set(resp_gene_context.alternate_labels) == \
                   set(test_variation_context.alternate_labels), "gene_context.alternate_labels"  # noqa: E501
        assert len(resp_gene_context.extensions) == \
               len(test_variation_context.extensions), "len gene_context.extensions"  # noqa: E501
        for resp_ext in resp_gene_context.extensions:
            for test_var in test_variation_context.extensions:
                if resp_ext.name == test_var.name:
                    if resp_ext.name == 'chromosome_location':
                        assert resp_ext.value == test_var.value, \
                            "gene_context.chromosome_location"
                    elif resp_ext.name == 'associated_with':
                        assert set(resp_ext.value) == set(test_var.value), \
                            "gene_context.associated_with"
                    else:
                        assert resp_ext.value == test_var.value,\
                            f"gene_context.{resp_ext.name}"
    else:
        assert not test_variation_context
