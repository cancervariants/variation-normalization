"""Module for testing the normalize endpoint."""
import pytest
from variant.normalize import Normalize
from variant.schemas.ga4gh_vod import VariationDescriptor
from variant.to_vrs import ToVRS


@pytest.fixture(scope="module")
def test_normalize():
    """Build normalize test fixture."""
    class TestNormalize:

        def __init__(self):
            self.to_vrs = ToVRS()
            self.test_normalize = Normalize()

        def normalize(self, q):
            resp = self.test_normalize.normalize(q,
                                                 self.to_vrs.get_validations(q),  # noqa: E501
                                                 self.to_vrs.amino_acid_cache)
            return resp

    return TestNormalize()


@pytest.fixture(scope="module")
def braf_v600e():
    """Create BRAF V600E fixture."""
    params = {
        "id": "normalize.variant:BRAF%20V600E",
        "type": "VariationDescriptor",
        "value_id": "ga4gh:VA.u6sKlz0mMQvARmrlnt0Aksz6EbSkmL8z",
        "value": {
            "location": {
                "interval": {
                    "end": 600,
                    "start": 599,
                    "type": "SimpleInterval"
                },
                "sequence_id": "ga4gh:SQ.ZJwurRo2HLY018wghYjDKSfIlEH0Y8At",
                "type": "SequenceLocation"
            },
            "state": {
                "sequence": "E",
                "type": "SequenceState"
            },
            "type": "Allele"
        },
        "label": "BRAF V600E",
        "molecule_context": "protein",
        "structural_type": "SO:0001606",
        "ref_allele_seq": "V",
        "gene_context": {
            "id": "normalize.gene:BRAF",
            "type": "GeneDescriptor",
            "label": "BRAF",
            "value": {
                "gene_id": "hgnc:1097",
                "type": "Gene"
            },
            "xrefs": [
                "ncbigene:673",
                "ensembl:ENSG00000157764"
            ],
            "alternate_labels": [
                "B-Raf proto-oncogene, serine/threonine kinase",
                "BRAF1"
            ],
            "extensions": [
                {
                    "type": "Extension",
                    "name": "symbol_status",
                    "value": "approved"
                },
                {
                    "type": "Extension",
                    "name": "associated_with",
                    "value": [
                        "vega:OTTHUMG00000157457",
                        "ucsc:uc003vwc.5",
                        "ccds:CCDS5863",
                        "ccds:CCDS87555",
                        "uniprot:P15056",
                        "pubmed:2284096",
                        "pubmed:1565476",
                        "cosmic:BRAF",
                        "omim:164757",
                        "orphanet:119066",
                        "iuphar:1943",
                        "ena.embl:M95712",
                        "refseq:NM_004333"
                    ]
                },
                {
                    "type": "Extension",
                    "name": "chromosome_location",
                    "value": {
                        "_id":
                            "ga4gh:VCL.O6yCQ1cnThOrTfK9YUgMlTfM6HTqbrKw",
                        "type": "ChromosomeLocation",
                        "species_id": "taxonomy:9606",
                        "chr": "7",
                        "interval": {
                            "end": "q34",
                            "start": "q34",
                            "type": "CytobandInterval"
                        }
                    }
                }
            ]
        }
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def vhl():
    """Create VHL Tyr185Ter fixture."""
    params = {
        "id": "normalize.variant:NP_000542.1%3Ap.Tyr185Ter",
        "type": "VariationDescriptor",
        "value_id": "ga4gh:VA.5Zx8fM1_wE3T_DFPbJgEe5CD-youM0op",
        "value": {
            "location": {
                "interval": {
                    "end": 185,
                    "start": 184,
                    "type": "SimpleInterval"
                },
                "sequence_id": "ga4gh:SQ.z-Oa0pZkJ6GHJHOYM7h5mY_umc0SJzTu",
                "type": "SequenceLocation"
            },
            "state": {
                "sequence": "*",
                "type": "SequenceState"
            },
            "type": "Allele"
        },
        "label": "NP_000542.1:p.Tyr185Ter",
        "molecule_context": "protein",
        "structural_type": "SO:0001606",
        "ref_allele_seq": "Y",
        "gene_context": {
            "id": "normalize.gene:VHL",
            "type": "GeneDescriptor",
            "label": "VHL",
            "value": {
                "gene_id": "hgnc:12687",
                "type": "Gene"
            },
            "xrefs": [
                "ncbigene:7428",
                "ensembl:ENSG00000134086"
            ],
            "alternate_labels": [
                "VHL1",
                "von Hippel-Lindau tumor suppressor"
            ],
            "extensions": [
                {
                    "type": "Extension",
                    "name": "symbol_status",
                    "value": "approved"
                },
                {
                    "type": "Extension",
                    "name": "associated_with",
                    "value": [
                        "vega:OTTHUMG00000128668",
                        "ucsc:uc003bvc.4",
                        "ccds:CCDS2597",
                        "ccds:CCDS2598",
                        "uniprot:P40337",
                        "pubmed:9671762",
                        "cosmic:VHL",
                        "omim:608537",
                        "orphanet:120467",
                        "ena.embl:L15409",
                        "refseq:NM_000551"
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
                }
            ]
        }
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def kit():
    """Create NP_000213.1:p.Leu862= fixture."""
    params = {
        "id": "normalize.variant:NP_000213.1%3Ap.Leu862%3D",
        "type": "VariationDescriptor",
        "value_id": "ga4gh:VA.N9L1bGWMk2IDg9aB83D-pS-V6n-oqqxy",
        "value": {
            "location": {
                "interval": {
                    "end": 862,
                    "start": 861,
                    "type": "SimpleInterval"
                },
                "sequence_id": "ga4gh:SQ.TcMVFj5kDODDWpiy1d_1-3_gOf4BYaAB",
                "type": "SequenceLocation"
            },
            "state": {
                "sequence": "L",
                "type": "SequenceState"
            },
            "type": "Allele"
        },
        "label": "NP_000213.1:p.Leu862=",
        "molecule_context": "protein",
        "structural_type": "SO:0001606",
        "ref_allele_seq": "L",
        "gene_context": {
            "id": "normalize.gene:KIT",
            "type": "GeneDescriptor",
            "label": "KIT",
            "value": {
                "gene_id": "hgnc:6342",
                "type": "Gene"
            },
            "xrefs": [
                "ncbigene:3815",
                "ensembl:ENSG00000157404"
            ],
            "alternate_labels": [
                "C-Kit",
                "CD117",
                "SCFR",
                "PBT",
                "KIT proto-oncogene, receptor tyrosine kinase"
            ],
            "extensions": [
                {
                    "type": "Extension",
                    "name": "symbol_status",
                    "value": "approved"
                },
                {
                    "type": "Extension",
                    "name": "associated_with",
                    "value": [
                        "vega:OTTHUMG00000128713",
                        "ucsc:uc010igr.4",
                        "ccds:CCDS47058",
                        "ccds:CCDS3496",
                        "uniprot:P10721",
                        "pubmed:9027509",
                        "cosmic:KIT",
                        "omim:164920",
                        "orphanet:122862",
                        "iuphar:1805",
                        "hcdmdb:CD117",
                        "refseq:NM_000222",
                        "ena.embl:S67773"
                    ]
                },
                {
                    "type": "Extension",
                    "name": "chromosome_location",
                    "value": {
                        "_id":
                            "ga4gh:VCL.QH-9ROGxiMyAhhzVwvwcGrOfQ0kjO2yS",
                        "type": "ChromosomeLocation",
                        "species_id": "taxonomy:9606",
                        "chr": "4",
                        "interval": {
                            "end": "q12",
                            "start": "q12",
                            "type": "CytobandInterval"
                        }
                    }
                }
            ]
        }
    }
    return VariationDescriptor(**params)


def assertion_checks(normalize_response, test_variant):
    """Check that normalize_response and variant_query are equal."""
    assert normalize_response.id == test_variant.id
    assert normalize_response.type == test_variant.type
    assert normalize_response.label == test_variant.label
    assert normalize_response.value == test_variant.value
    assert normalize_response.molecule_context == test_variant.molecule_context
    assert normalize_response.structural_type == test_variant.structural_type
    assert normalize_response.ref_allele_seq == test_variant.ref_allele_seq

    resp_gene_context = normalize_response.gene_context
    test_variant_context = test_variant.gene_context
    assert resp_gene_context.id == test_variant_context.id
    assert resp_gene_context.label == test_variant_context.label
    assert resp_gene_context.value_id == test_variant_context.value_id
    assert set(resp_gene_context.xrefs) == set(test_variant_context.xrefs)
    assert set(resp_gene_context.alternate_labels) == \
           set(test_variant_context.alternate_labels)
    assert len(resp_gene_context.extensions) == \
           len(test_variant_context.extensions)
    for resp_ext in resp_gene_context.extensions:
        for test_var in test_variant_context.extensions:
            if resp_ext.name == test_var.name:
                if resp_ext.name == 'chromosome_location':
                    assert resp_ext.value == test_var.value
                elif resp_ext.name == 'associated_with':
                    assert set(resp_ext.value) == set(test_var.value)
                else:
                    assert resp_ext.value == test_var.value


def test_amino_acid_substitution(test_normalize, braf_v600e):
    """Test that amino acid substitutions normalize correctly."""
    resp = test_normalize.normalize('     BRAF      V600E    ')
    assertion_checks(resp, braf_v600e)


def test_polypeptide_truncation(test_normalize, vhl):
    """Test that polypeptide truncations normalize correctly."""
    resp = test_normalize.normalize('NP_000542.1:p.Tyr185Ter')
    assertion_checks(resp, vhl)


def test_silent_mutation(test_normalize, kit):
    """Test that silent mutations V600E normalize correctly."""
    resp = test_normalize.normalize('NP_000213.1:p.Leu862=')
    assertion_checks(resp, kit)
