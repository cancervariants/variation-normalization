"""Module for testing the normalize endpoint."""
import pytest
from variant.normalize import Normalize
from variant.schemas.ga4gh_vod import VariationDescriptor
from variant.to_vrs import ToVRS
from variant.main import normalize as normalize_get_response
from variant.main import translate as to_vrs_get_response
from datetime import datetime


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


@pytest.fixture(scope='module')
def braf_gene_context():
    """Create BRAF gene context test fixture."""
    return {
        "id": "normalize.gene:BRAF",
        "type": "GeneDescriptor",
        "label": "BRAF",
        "value": {
            "id": "hgnc:1097",
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
                    "_id": "ga4gh:VCL.O6yCQ1cnThOrTfK9YUgMlTfM6HTqbrKw",
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


@pytest.fixture(scope='module')
def vhl_gene_context():
    """Create a VHL gene context."""
    return {
        "id": "normalize.gene:VHL",
        "type": "GeneDescriptor",
        "label": "VHL",
        "value": {
            "id": "hgnc:12687",
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


@pytest.fixture(scope='module')
def erbb2_context():
    """Create test fixture for ERBB2 Gene Context."""
    return {
        "id": "normalize.gene:ERBB2",
        "type": "GeneDescriptor",
        "label": "ERBB2",
        "value": {
            "id": "hgnc:3430",
            "type": "Gene"
        },
        "xrefs": [
            "ncbigene:2064",
            "ensembl:ENSG00000141736"
        ],
        "alternate_labels": [
            "NEU",
            "HER-2",
            "HER2",
            "CD340",
            "erb-b2 receptor tyrosine kinase 2",
            "NGL"
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
                    "vega:OTTHUMG00000179300",
                    "ucsc:uc002hso.4",
                    "ccds:CCDS77016",
                    "ccds:CCDS74052",
                    "ccds:CCDS45667",
                    "ccds:CCDS32642",
                    "ccds:CCDS77017",
                    "uniprot:P04626",
                    "cosmic:ERBB2",
                    "omim:164870",
                    "iuphar:2019",
                    "hcdmdb:CD340",
                    "ena.embl:X03363",
                    "refseq:NM_004448"
                ]
            },
            {
                "type": "Extension",
                "name": "chromosome_location",
                "value": {
                    "_id": "ga4gh:VCL.pS7M3aeNymozN9LKeAwVDEB5H1nt4Kqy",
                    "type": "ChromosomeLocation",
                    "species_id": "taxonomy:9606",
                    "chr": "17",
                    "interval": {
                        "end": "q12",
                        "start": "q12",
                        "type": "CytobandInterval"
                    }
                }
            }
        ]
    }


@pytest.fixture(scope="module")
def braf_v600e(braf_gene_context):
    """Create BRAF V600E protein test fixture."""
    params = {
        "id": "normalize.variant:BRAF%20V600E",
        "type": "VariationDescriptor",
        "value_id": "ga4gh:VA.9dA0egRAIfVFDL1sdU1VP7HsBcG0-DtE",
        "value": {
            "location": {
                "interval": {
                    "end": 640,
                    "start": 639,
                    "type": "SimpleInterval"
                },
                "sequence_id": "ga4gh:SQ.WaAJ_cXXn9YpMNfhcq9lnzIvaB9ALawo",
                "type": "SequenceLocation"
            },
            "state": {
                "sequence": "E",
                "type": "SequenceState"
            },
            "type": "Allele"
        },
        "label": "NP_001361187.1:p.Val640Glu",
        "molecule_context": "protein",
        "structural_type": "SO:0001606",
        "ref_allele_seq": "V",
        "gene_context": braf_gene_context
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def vhl(vhl_gene_context):
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
        "structural_type": "SO:0001617",
        "ref_allele_seq": "Y",
        "gene_context": vhl_gene_context
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def vhl_silent(vhl_gene_context):
    """Create NP_000542.1:p.Pro61 fixture."""
    params = {
        "id": "normalize.variant:NP_000542.1%3Ap.Pro61%3D",
        "type": "VariationDescriptor",
        "label": "NP_000542.1:p.Pro61=",
        "value_id": "ga4gh:VA.LBNTm7QqFZp1alJHaFKlKuRY9cOfdHeI",
        "value": {
            "location": {
                "interval": {
                    "end": 61,
                    "start": 60,
                    "type": "SimpleInterval"
                },
                "sequence_id": "ga4gh:SQ.z-Oa0pZkJ6GHJHOYM7h5mY_umc0SJzTu",
                "type": "SequenceLocation"
            },
            "state": {
                "sequence": "P",
                "type": "SequenceState"
            },
            "type": "Allele"
        },
        "molecule_context": "protein",
        "structural_type": "SO:0001017",
        "ref_allele_seq": "P",
        "gene_context": vhl_gene_context
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope='module')
def braf_v600e_nucleotide(braf_gene_context):
    """Create a test fixture for BRAF V600E MANE select nucleotide hgvs."""
    params = {
        "id": "normalize.variant:NM_004333.4%3Ac.1799T%3EA",
        "type": "VariationDescriptor",
        "value_id": "ga4gh:VA.X_ij6wmw-fBwcoCVhHAfP7HiWUtkNfwq",
        "value": {
            "location": {
                "interval": {
                    "end": 1919,
                    "start": 1918,
                    "type": "SimpleInterval"
                },
                "sequence_id": "ga4gh:SQ.I_0feOk5bZ3VfH8ejhWQiMDe9o6o4QdR",
                "type": "SequenceLocation"
            },
            "state": {
                "sequence": "A",
                "type": "SequenceState"
            },
            "type": "Allele"
        },
        "label": "NM_001374258.1:c.1919T>A",
        "molecule_context": "transcript",
        "structural_type": "SO:0001483",
        "ref_allele_seq": "T",
        "gene_context": braf_gene_context
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope='module')
def nm_004448_coding_dna_delins(erbb2_context):
    """Create test fixture for NM_004448.4:c.2326_2327delinsCT."""
    params = {
        "id": "normalize.variant:NM_004448.4%3Ac.2326_2327delinsCT",
        "type": "VariationDescriptor",
        "value_id": "ga4gh:VA.G8pUN2zEuDtTfI8i30RNLF-gQAab4rUC",
        "value": {
            "location": {
                "interval": {
                    "end": 2327,
                    "start": 2325,
                    "type": "SimpleInterval"
                },
                "sequence_id": "ga4gh:SQ.y9b4LVMiCXpZxOg9Xt1NwRtssA03MwWM",
                "type": "SequenceLocation"
            },
            "state": {
                "sequence": "CT",
                "type": "SequenceState"
            },
            "type": "Allele"
        },
        "label": "NM_004448.4:c.2326_2327delinsCT",
        "molecule_context": "transcript",
        "structural_type": "SO:1000032",
        "ref_allele_seq": "GA",
        "gene_context": erbb2_context
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope='module')
def nc_000007_genomic_delins(braf_gene_context):
    """Create test fixture for NC_000007.13:g.140453135_140453136delinsAT."""
    params = {
        "id": "normalize.variant:NC_000007.13%3Ag.140453135_140453136delinsAT",
        "type": "VariationDescriptor",
        "value_id": "ga4gh:VA.X_ij6wmw-fBwcoCVhHAfP7HiWUtkNfwq",
        "value": {
            "location": {
                "interval": {
                    "end": 1919,
                    "start": 1918,
                    "type": "SimpleInterval"
                },
                "sequence_id": "ga4gh:SQ.I_0feOk5bZ3VfH8ejhWQiMDe9o6o4QdR",
                "type": "SequenceLocation"
            },
            "state": {
                "sequence": "A",
                "type": "SequenceState"
            },
            "type": "Allele"
        },
        "label": "NM_001374258.1:c.1919_1920delinsAT",
        "molecule_context": "genomic",
        "structural_type": "SO:1000032",
        "ref_allele_seq": "C",
        "gene_context": braf_gene_context
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope='module')
def nm_000551(vhl_gene_context):
    """Create test fixture for NM_000551.4:c.615delinsAA."""
    params = {
        "id": 'temp_id',
        "type": "VariationDescriptor",
        "value_id": "ga4gh:VA.tPF0lpCD-oDyX3SdqSDAbSXfaB_7Lo8x",
        "value": {
            "location": {
                "interval": {
                    "end": 615,
                    "start": 614,
                    "type": "SimpleInterval"
                },
                "sequence_id": "ga4gh:SQ.xBKOKptLLDr-k4hTyCetvARn16pDS_rW",
                "type": "SequenceLocation"
            },
            "state": {
                "sequence": "AA",
                "type": "SequenceState"
            },
            "type": "Allele"
        },
        "label": "NM_000551.4:c.615delinsAA",
        "molecule_context": "transcript",
        "structural_type": "SO:1000032",
        "ref_allele_seq": "G",
        "gene_context": vhl_gene_context
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope='module')
def coding_dna_silent_mutation(braf_gene_context):
    """Create test fixture for NM_004333.4:c.1799=."""
    params = {
        "id": 'normalize.variant:NM_004333.4%3Ac.1799%3D',
        "type": "VariationDescriptor",
        "value_id": "ga4gh:VA.KDlNyhYx4zWLSa_zB_db5bFOjuhmawK8",
        "value": {
            "location": {
                "interval": {
                    "end": 1799,
                    "start": 1798,
                    "type": "SimpleInterval"
                },
                "sequence_id": "ga4gh:SQ.jkiXxxRjK7uTMiW2KQFjpgvF3VQi-HhX",
                "type": "SequenceLocation"
            },
            "state": {
                "sequence": "A",
                "type": "SequenceState"
            },
            "type": "Allele"
        },
        "label": "NM_004333.4:c.1799=",
        "molecule_context": "transcript",
        "structural_type": "SO:0002073",
        "ref_allele_seq": "A",
        "gene_context": braf_gene_context
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope='module')
def nc_000007_silent_mutation():
    """Create test fixture for NC_000007.13:g.140453136=."""
    params = {
        "id": 'normalize.variant:NC_000007.13%3Ag.140453136%3D',
        "type": "VariationDescriptor",
        "value_id": "ga4gh:VA.y3OnKtG-zUM-LTJHOg7IsLVfm8x3PNLI",
        "value": {
            "location": {
                "interval": {
                    "end": 140453136,
                    "start": 140453135,
                    "type": "SimpleInterval"
                },
                "sequence_id": "ga4gh:SQ.IW78mgV5Cqf6M24hy52hPjyyo5tCCd86",
                "type": "SequenceLocation"
            },
            "state": {
                "sequence": "A",
                "type": "SequenceState"
            },
            "type": "Allele"
        },
        "label": "NC_000007.13:g.140453136=",
        "molecule_context": "genomic",
        "structural_type": "SO:0002073",
        "ref_allele_seq": "A",
        "gene_context": None
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope='module')
def amino_acid_delins():
    """Create test fixture for amino acid delins."""
    params = {
        "id": 'normalize.variant:NP_001333827.1%3Ap.Leu747_Thr751delinsPro',
        "type": "VariationDescriptor",
        "value_id": "ga4gh:VA.drLuUW5T542RCeDlVo4zbQ-_tcAiEnb6",
        "value": {
            "location": {
                "interval": {
                    "end": 751,
                    "start": 746,
                    "type": "SimpleInterval"
                },
                "sequence_id": "ga4gh:SQ.vyo55F6mA6n2LgN4cagcdRzOuh38V4mE",
                "type": "SequenceLocation"
            },
            "state": {
                "sequence": "P",
                "type": "SequenceState"
            },
            "type": "Allele"
        },
        "label": "NP_005219.2:p.Leu747_Thr751delinsPro",
        "molecule_context": "protein",
        "structural_type": "SO:1000032",
        "ref_allele_seq": "LREAT",
        "gene_context": {
            "id": "normalize.gene:EGFR",
            "type": "GeneDescriptor",
            "label": "EGFR",
            "value": {
                "id": "hgnc:3236",
                "type": "Gene"
            },
            "xrefs": [
                "ncbigene:1956",
                "ensembl:ENSG00000146648"
            ],
            "alternate_labels": [
                "ERBB1",
                "ERRP",
                "ERBB",
                "epidermal growth factor receptor"
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
                        "vega:OTTHUMG00000023661",
                        "ucsc:uc003tqk.4",
                        "ccds:CCDS87507",
                        "ccds:CCDS47587",
                        "ccds:CCDS87506",
                        "ccds:CCDS5514",
                        "ccds:CCDS5515",
                        "ccds:CCDS5516",
                        "uniprot:P00533",
                        "pubmed:1505215",
                        "cosmic:EGFR",
                        "omim:131550",
                        "orphanet:121311",
                        "iuphar:1797",
                        "refseq:NM_005228"
                    ]
                },
                {
                    "type": "Extension",
                    "name": "chromosome_location",
                    "value": {
                        "_id": "ga4gh:VCL.wgFi9e72ZIIJaOfLx5gaOeGrwP_IZoQ2",
                        "type": "ChromosomeLocation",
                        "species_id": "taxonomy:9606",
                        "chr": "7",
                        "interval": {
                            "end": "p11.2",
                            "start": "p11.2",
                            "type": "CytobandInterval"
                        }
                    }
                }
            ]
        }
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope='module')
def amino_acid_deletion_np_range(erbb2_context):
    """Create test fixture for amino acid deletion using NP accession and
    range for deletion.
    """
    params = {
        "id": 'normalize.variant:NP_004439.2%3Ap.Leu755_Thr759del',
        "type": "VariationDescriptor",
        "value_id": "ga4gh:VA.Kzk2XquE5w5Ujd_tPBLVOZcylXMP8xbW",
        "value": {
            "location": {
                "interval": {
                    "end": 759,
                    "start": 754,
                    "type": "SimpleInterval"
                },
                "sequence_id": "ga4gh:SQ.AF1UFydIo02-bMplonKSfxlWY2q6ze3m",
                "type": "SequenceLocation"
            },
            "state": {
                "sequence": "",
                "type": "SequenceState"
            },
            "type": "Allele"
        },
        "label": "NP_004439.2:p.Leu755_Thr759del",
        "molecule_context": "protein",
        "structural_type": "SO:0001604",
        "ref_allele_seq": "LRENT",
        "gene_context": erbb2_context
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope='module')
def coding_dna_deletion(erbb2_context):
    """Create test fixture for coding dna deletion range with deleted
    sequence.
    """
    params = {
        "id": 'normalize.variant:NM_004448.3%3Ac.2263_2277delTTGAGGGAAAACACA',
        "type": "VariationDescriptor",
        "value_id": "ga4gh:VA.g_DrfzkfKKB1p8LTPvYLn3pIDBgrPV0K",
        "value": {
            "location": {
                "interval": {
                    "end": 2278,
                    "start": 2263,
                    "type": "SimpleInterval"
                },
                "sequence_id": "ga4gh:SQ.y9b4LVMiCXpZxOg9Xt1NwRtssA03MwWM",
                "type": "SequenceLocation"
            },
            "state": {
                "sequence": "",
                "type": "SequenceState"
            },
            "type": "Allele"
        },
        "label": "NM_004448.4:c.2264_2278del",
        "molecule_context": "transcript",
        "structural_type": "SO:0000159",
        "ref_allele_seq": "GTGGAGCCGCTGACA",
        "gene_context": erbb2_context
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope='module')
def genomic_deletion(vhl_gene_context):
    """Create test fixture for genomic deletion range."""
    params = {
        "id": 'normalize.variant:NC_000003.11%3Ag.10188279_10188297del',
        "type": "VariationDescriptor",
        "value_id": "ga4gh:VA.s_HhEqbHB3MRwvmlgvAtBOPFFsSLUAyA",
        "value": {
            "location": {
                "interval": {
                    "end": 441,
                    "start": 421,
                    "type": "SimpleInterval"
                },
                "sequence_id": "ga4gh:SQ.xBKOKptLLDr-k4hTyCetvARn16pDS_rW",
                "type": "SequenceLocation"
            },
            "state": {
                "sequence": "C",
                "type": "SequenceState"
            },
            "type": "Allele"
        },
        "label": "NM_000551.4:c.422_440del",
        "molecule_context": "genomic",
        "structural_type": "SO:0000159",
        "ref_allele_seq": "CTCTTCAGAGATGCAGGGAC",
        "gene_context": vhl_gene_context
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
    if resp_gene_context:
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
    else:
        assert not test_variant_context


def test_amino_acid_substitution(test_normalize, braf_v600e):
    """Test that amino acid substitutions normalize correctly."""
    resp = test_normalize.normalize('     BRAF      V600E    ')
    assertion_checks(resp, braf_v600e)

    braf_id = "normalize.variant:BRAF%20V600E"

    resp = test_normalize.normalize('NP_004324.2:p.Val600Glu')
    assert resp.id == "normalize.variant:NP_004324.2%3Ap.Val600Glu"
    resp.id = braf_id
    assertion_checks(resp, braf_v600e)

    resp = test_normalize.normalize('braf v512e')
    assert resp.id == 'normalize.variant:braf%20v512e'
    resp.id = braf_id
    assertion_checks(resp, braf_v600e)

    resp = test_normalize.normalize(' NP_001365404.1:p.Val512Glu  ')
    assert resp.id == 'normalize.variant:NP_001365404.1%3Ap.Val512Glu'
    resp.id = braf_id
    assertion_checks(resp, braf_v600e)


def test_polypeptide_truncation(test_normalize, vhl):
    """Test that polypeptide truncations normalize correctly."""
    resp = test_normalize.normalize('NP_000542.1:p.Tyr185Ter')
    assertion_checks(resp, vhl)


def test_silent_mutation(test_normalize, vhl_silent):
    """Test that silent mutations normalize correctly."""
    resp = test_normalize.normalize('NP_000542.1:p.Pro61=')
    assertion_checks(resp, vhl_silent)


def test_coding_dna_and_genomic_substitution(test_normalize,
                                             braf_v600e_nucleotide):
    """Test that coding dna and genomic substitutions normalize correctly."""
    resp = test_normalize.normalize('NM_004333.4:c.1799T>A')
    assertion_checks(resp, braf_v600e_nucleotide)

    # MANE transcript
    refseq_id = 'normalize.variant:NM_004333.4%3Ac.1799T%3EA'
    refseq_label = 'NM_001374258.1:c.1919T>A'
    ensembl_label = 'ENST00000644969.2:c.1919T>A'

    # TODO: Check if this should return a different VRS object?
    resp = test_normalize.normalize('ENST00000288602.6:c.1799T>A')
    assert resp.id == 'normalize.variant:ENST00000288602.6%3Ac.1799T%3EA'
    assert resp.label == ensembl_label
    resp.id = refseq_id
    resp.label = refseq_label
    assertion_checks(resp, braf_v600e_nucleotide)

    resp = test_normalize.normalize('BRAF V600E c.1799T>A')
    assert resp.id == 'normalize.variant:BRAF%20V600E%20c.1799T%3EA'
    assert resp.label == ensembl_label
    resp.id = refseq_id
    resp.label = refseq_label
    assertion_checks(resp, braf_v600e_nucleotide)

    resp = test_normalize.normalize('BRAF V600E (c.1799T>A)')
    assert resp.id == 'normalize.variant:BRAF%20V600E%20%28c.1799T%3EA%29'
    assert resp.label == ensembl_label
    resp.id = refseq_id
    resp.label = refseq_label
    assertion_checks(resp, braf_v600e_nucleotide)

    resp = test_normalize.normalize('BRAF c.1799T>A')
    assert resp.id == 'normalize.variant:BRAF%20c.1799T%3EA'
    assert resp.label == ensembl_label
    resp.id = refseq_id
    resp.label = refseq_label
    assertion_checks(resp, braf_v600e_nucleotide)

    resp = test_normalize.normalize('NC_000007.13:g.140453136A>T')
    assert resp.id == 'normalize.variant:NC_000007.13%3Ag.140453136A%3ET'
    assert resp.ref_allele_seq == 'A'
    assert resp.molecule_context == 'genomic'
    resp.id = refseq_id
    resp.ref_allele_seq = 'T'
    resp.molecule_context = 'transcript'
    assertion_checks(resp, braf_v600e_nucleotide)

    # TODO: Issue 99
    # resp = test_normalize.normalize('BRAF V600E (g.140453136A>T)')
    # assert resp.id == 'normalize.variant:BRAF%20V600E%20%28g.140453136A%3ET%29'  # noqa: E501
    # assert resp.label == ensembl_label
    # assert resp.ref_allele_seq == 'A'
    # assert resp.molecule_context == 'genomic'
    # resp.id = refseq_id
    # resp.label = refseq_label
    # resp.ref_allele_seq = 'T'
    # resp.molecule_context = 'transcript'
    # assertion_checks(resp, braf_v600e_nucleotide)

    resp = test_normalize.normalize('BRAF g.140453136A>T')
    assert resp.id == 'normalize.variant:BRAF%20g.140453136A%3ET'
    assert resp.label == ensembl_label
    assert resp.ref_allele_seq == 'A'
    assert resp.molecule_context == 'genomic'
    resp.id = refseq_id
    resp.label = refseq_label
    resp.ref_allele_seq = 'T'
    resp.molecule_context = 'transcript'
    assertion_checks(resp, braf_v600e_nucleotide)


def assert_coding_dna_genomic_silent_mutation(resp, gene_context, start, stop):
    """Check for coding dna or genomic silent mutation."""
    assert resp
    assert resp.value_id.startswith('ga4gh:VA.')
    assert resp.value['location']['interval']['start'] == start
    assert resp.value['location']['interval']['end'] == stop
    assert resp.value['location']['sequence_id'].startswith('ga4gh:SQ.')
    assert resp.gene_context.dict(exclude_none=True) == gene_context
    assert resp.structural_type == 'SO:0002073'


def test_coding_dna_silent_mutation(test_normalize,
                                    coding_dna_silent_mutation,
                                    braf_gene_context):
    """Test that Coding DNA Silent Mutation normalizes correctly."""
    resp = test_normalize.normalize('NM_004333.4:c.1799= ')
    assertion_checks(resp, coding_dna_silent_mutation)

    # Only selects first from VRS, so this might change.
    # Update test once MANE select works on silent mutations
    resp = test_normalize.normalize('ENST00000288602.6:c.1799=')
    assert_coding_dna_genomic_silent_mutation(resp, braf_gene_context,
                                              1798, 1799)
    assert resp.id == 'normalize.variant:ENST00000288602.6%3Ac.1799%3D'
    assert resp.label == 'ENST00000288602.6:c.1799='
    assert resp.molecule_context == 'transcript'

    resp = test_normalize.normalize('BRAF    c.1799=')
    assert_coding_dna_genomic_silent_mutation(resp, braf_gene_context,
                                              1798, 1799)
    assert resp.id == 'normalize.variant:BRAF%20c.1799%3D'
    assert resp.label == 'BRAF c.1799='
    assert resp.molecule_context == 'transcript'

    resp = test_normalize.normalize('  BRAF  V600E  c.1799=  ')
    assert_coding_dna_genomic_silent_mutation(resp, braf_gene_context,
                                              1798, 1799)
    assert resp.id == 'normalize.variant:BRAF%20V600E%20c.1799%3D'
    assert resp.label == 'BRAF V600E c.1799='
    assert resp.molecule_context == 'transcript'


def test_genomic_silent_mutation(test_normalize, nc_000007_silent_mutation,
                                 braf_gene_context):
    """Test that genomic silent mutation normalizes correctly."""
    resp = test_normalize.normalize('NC_000007.13:g.140453136=')
    assertion_checks(resp, nc_000007_silent_mutation)

    resp = test_normalize.normalize('BRAF g.140453136=')
    assert_coding_dna_genomic_silent_mutation(resp, braf_gene_context,
                                              140453135, 140453136)
    assert resp.molecule_context == 'genomic'


def test_coding_dna_delins(test_normalize, nm_004448_coding_dna_delins,
                           nm_000551):
    """Test that Coding DNA DelIns normalizes correctly."""
    resp = test_normalize.normalize('    NM_004448.4:c.2326_2327delinsCT    ')
    assertion_checks(resp, nm_004448_coding_dna_delins)

    # TODO: Test ENST###.c

    resp = test_normalize.normalize('NM_000551.3:c.615delinsAA')
    nm_000551.id = 'normalize.variant:NM_000551.3%3Ac.615delinsAA'
    assertion_checks(resp, nm_000551)


def test_genomic_delins(test_normalize, nc_000007_genomic_delins,
                        nm_000551):
    """Test that Genomic DelIns normalizes correctly."""
    resp = test_normalize.normalize(
        'NC_000007.13:g.140453135_140453136delinsAT'
    )
    assertion_checks(resp, nc_000007_genomic_delins)

    resp = test_normalize.normalize('NC_000003.12:g.10149938delinsAA')
    assert resp.molecule_context == 'genomic'
    resp.molecule_context = 'transcript'
    nm_000551.id = 'normalize.variant:NC_000003.12%3Ag.10149938delinsAA'
    assertion_checks(resp, nm_000551)


def test_amino_acid_delins(test_normalize, amino_acid_delins):
    """Test that Amnio Acid DelIns normalizes correctly."""
    resp = test_normalize.normalize('NP_001333827.1:p.Leu747_Thr751delinsPro')
    assertion_checks(resp, amino_acid_delins)

    resp = test_normalize.normalize('EGFR p.Leu747_Thr751delinsPro')
    assert resp.id == 'normalize.variant:EGFR%20p.Leu747_Thr751delinsPro'
    resp.id = 'normalize.variant:NP_001333827.1%3Ap.Leu747_Thr751delinsPro'
    assertion_checks(resp, amino_acid_delins)

    resp = test_normalize.normalize('EGFR Leu747_Thr751delinsPro')
    assert resp.id == 'normalize.variant:EGFR%20Leu747_Thr751delinsPro'
    resp.id = 'normalize.variant:NP_001333827.1%3Ap.Leu747_Thr751delinsPro'
    assertion_checks(resp, amino_acid_delins)

    resp = test_normalize.normalize('EGFR L747_T751delinsP')
    assert resp.id == 'normalize.variant:EGFR%20L747_T751delinsP'
    resp.id = 'normalize.variant:NP_001333827.1%3Ap.Leu747_Thr751delinsPro'
    assertion_checks(resp, amino_acid_delins)


def test_amino_acid_deletion(test_normalize, amino_acid_deletion_np_range):
    """Test that Amino Acid Deletion normalizes correctly."""
    resp = test_normalize.normalize('NP_004439.2:p.Leu755_Thr759del')
    assertion_checks(resp, amino_acid_deletion_np_range)

    resp = test_normalize.normalize('ERBB2 p.Leu755_Thr759del')
    assert resp.id == 'normalize.variant:ERBB2%20p.Leu755_Thr759del'

    resp = test_normalize.normalize('ERBB2 Leu755_Thr759del')
    assert resp.id == 'normalize.variant:ERBB2%20Leu755_Thr759del'


def test_coding_dna_deletion(test_normalize, coding_dna_deletion):
    """Test that coding dna deletion normalizes correctly."""
    resp = \
        test_normalize.normalize('NM_004448.3:c.2263_2277delTTGAGGGAAAACACA')
    assertion_checks(resp, coding_dna_deletion)

    resp = test_normalize.normalize('ERBB2 c.2263_2277delTTGAGGGAAAACACA')
    assert resp.id == 'normalize.variant:ERBB2%20c.2263_2277delTTGAGGGAAAACACA'
    resp.id = 'normalize.variant:NM_004448.3%3Ac.2263_2277delTTGAGGGAAAACACA'
    assert resp.label == 'ENST00000269571.10:c.2264_2278del'
    resp.label = 'NM_004448.4:c.2264_2278del'
    assert resp.ref_allele_seq is None  # seqrepo can't find enst transcript
    resp.ref_allele_seq = 'GTGGAGCCGCTGACA'
    assertion_checks(resp, coding_dna_deletion)


def test_genomic_deletion(test_normalize, genomic_deletion):
    """Test that genomic deletion normalizes correctly."""
    resp = test_normalize.normalize('NC_000003.11:g.10188279_10188297del')
    assertion_checks(resp, genomic_deletion)

    resp = test_normalize.normalize('VHL g.10188279_10188297del')
    assert resp.id == 'normalize.variant:VHL%20g.10188279_10188297del'
    resp.id = 'normalize.variant:NC_000003.11%3Ag.10188279_10188297del'
    assert resp.label == 'ENST00000256474.3:c.422_440del'
    resp.label = 'NM_000551.4:c.422_440del'
    assert resp.ref_allele_seq is None  # seqrepo can't find enst transcript
    resp.ref_allele_seq = 'CTCTTCAGAGATGCAGGGAC'
    assertion_checks(resp, genomic_deletion)


def test_no_matches(test_normalize):
    """Test no matches work correctly."""
    queries = [
        "", "braf", "braf v600000932092039e", "NP_000213.1:cp.Leu862=",
        "NP_000213.1:cp.Leu862", "BRAF V600E 33", "NP_004324.2:p.Glu600Val",
        "NP_004324.2:p.Glu600Gal", "NP_004324.2839:p.Glu600Val",
        "NP_004324.2:t.Glu600Val", "this:c.54G>H", "NC_000007.13:g.4T<A",
        "test", "131", "braf z600e", "braf e600z", "Thr790Met", "p.Tyr365Ter",
        "ERBB2 G776delinsVCZ", "NP005219.2:p.Glu746_Thr751delinsValAla",
        "NP_005219.2:p.Glu746Thr751delinsValAla", "EGFR L747_L474delinsP",
        "NP_005219.2:p.Glu746_Thr751delinssValAla", "EGFR delins",
        "NM_004333.4:c.1799_1800delTGinsAT",
        "NM_173851.3(SLC30A8):c.973C>T%20(p.Arg325Trp)"
    ]
    for q in queries:
        resp = test_normalize.normalize(q)
        assert resp is None


def test_service_meta():
    """Test that service meta info populates correctly."""
    response = normalize_get_response('BRAF v600e')
    service_meta = response.service_meta_
    assert service_meta.name == "variant-normalizer"
    assert service_meta.version
    assert isinstance(service_meta.response_datetime, datetime)
    assert service_meta.url == 'https://github.com/cancervariants/variant-normalization'  # noqa: E501

    response = normalize_get_response('this-wont-normalize')
    service_meta = response.service_meta_
    assert service_meta.name == "variant-normalizer"
    assert service_meta.version
    assert isinstance(service_meta.response_datetime, datetime)
    assert service_meta.url == 'https://github.com/cancervariants/variant-normalization'  # noqa: E501

    response = to_vrs_get_response('BRAF v600e')
    service_meta = response.service_meta_
    assert service_meta.name == "variant-normalizer"
    assert service_meta.version
    assert isinstance(service_meta.response_datetime, datetime)
    assert service_meta.url == 'https://github.com/cancervariants/variant-normalization'  # noqa: E501

    response = to_vrs_get_response('this-wont-normalize')
    service_meta = response.service_meta_
    assert service_meta.name == "variant-normalizer"
    assert service_meta.version
    assert isinstance(service_meta.response_datetime, datetime)
    assert service_meta.url == 'https://github.com/cancervariants/variant-normalization'  # noqa: E501
