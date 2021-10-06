"""Create methods used throughout tests."""


def assertion_checks(normalize_response, test_variation, ignore_id=False):
    """Check that normalize_response and test_variation are equal."""
    if not ignore_id:
        assert normalize_response.id == test_variation.id
    assert normalize_response.type == test_variation.type
    assert normalize_response.variation_id == test_variation.variation_id
    assert normalize_response.variation == test_variation.variation
    assert normalize_response.molecule_context == \
           test_variation.molecule_context
    assert normalize_response.structural_type == test_variation.structural_type
    assert normalize_response.vrs_ref_allele_seq == \
           test_variation.vrs_ref_allele_seq

    resp_gene_context = normalize_response.gene_context
    test_variation_context = test_variation.gene_context
    if resp_gene_context:
        assert resp_gene_context.id == test_variation_context.id
        assert resp_gene_context.label == test_variation_context.label
        assert resp_gene_context.gene_id == test_variation_context.gene_id
        assert set(resp_gene_context.xrefs) ==\
               set(test_variation_context.xrefs)
        if test_variation_context.alternate_labels:
            assert set(resp_gene_context.alternate_labels) == \
                   set(test_variation_context.alternate_labels)
        assert len(resp_gene_context.extensions) == \
               len(test_variation_context.extensions)
        for resp_ext in resp_gene_context.extensions:
            for test_var in test_variation_context.extensions:
                if resp_ext.name == test_var.name:
                    if resp_ext.name == 'chromosome_location':
                        assert resp_ext.value == test_var.value
                    elif resp_ext.name == 'associated_with':
                        assert set(resp_ext.value) == set(test_var.value)
                    else:
                        assert resp_ext.value == test_var.value
    else:
        assert not test_variation_context
