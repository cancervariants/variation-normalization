"""Module for testing the normalize endpoint."""
import copy
from datetime import datetime

import pytest
from ga4gh.vrsatile.pydantic.vrsatile_models import VariationDescriptor

from variation.main import normalize as normalize_get_response
from variation.main import to_vrs as to_vrs_get_response
from tests.conftest import assertion_checks


@pytest.fixture(scope="module")
def test_handler(test_query_handler):
    """Create test fixture for normalize handler"""
    return test_query_handler.normalize_handler


@pytest.fixture(scope="module")
def limk2_gene_context():
    """Create LIMK2 gene context test fixture."""
    return {
        "id": "normalize.gene:LIMK2",
        "type": "GeneDescriptor",
        "label": "LIMK2",
        "gene_id": "hgnc:6614",
        "xrefs": [
            "ncbigene:3985",
            "ensembl:ENSG00000182541"
        ],
        "extensions": [
            {
                "type": "Extension",
                "name": "symbol_status",
                "value": "approved"
            },
            {
                "name": "approved_name",
                "value": "LIM domain kinase 2",
                "type": "Extension"
            },
            {
                "type": "Extension",
                "name": "associated_with",
                "value": [
                    "refseq:NM_016733",
                    "ccds:CCDS33637",
                    "ccds:CCDS13892",
                    "ena.embl:D45906",
                    "uniprot:P53671",
                    "pubmed:10591208",
                    "vega:OTTHUMG00000151251",
                    "omim:601988",
                    "iuphar:2055",
                    "pubmed:8537403",
                    "ccds:CCDS13891",
                    "ucsc:uc003akh.4"
                ]
            },
            {
                "type": "Extension",
                "name": "hgnc_locations",
                "value": [
                    {
                        "_id": "ga4gh:VCL.IoyhTh4PxvPx8yF9P3IecXDVs_XVbDe9",
                        "type": "ChromosomeLocation",
                        "species_id": "taxonomy:9606",
                        "chr": "22",
                        "interval": {
                            "end": "q12.2",
                            "start": "q12.2",
                            "type": "CytobandInterval"
                        }
                    }
                ]
            },
            {
                "type": "Extension",
                "name": "ensembl_locations",
                "value": [
                    {
                        "_id": "ga4gh:VSL.Mat3OiuXCd7vNWAGB5lOKaSf_ESc_xK9",
                        "type": "SequenceLocation",
                        "sequence_id": "ga4gh:SQ.7B7SHsmchAR0dFcDCuSFjJAo7tX87krQ",
                        "interval": {
                            "start": {"type": "Number", "value": 31212238},
                            "end": {"type": "Number", "value": 31280080},
                            "type": "SequenceInterval"
                        }
                    }
                ]
            },
            {
                "type": "Extension",
                "name": "ncbi_locations",
                "value": [
                    {
                        "_id": "ga4gh:VCL.IoyhTh4PxvPx8yF9P3IecXDVs_XVbDe9",
                        "type": "ChromosomeLocation",
                        "species_id": "taxonomy:9606",
                        "chr": "22",
                        "interval": {
                            "end": "q12.2",
                            "start": "q12.2",
                            "type": "CytobandInterval"
                        }
                    },
                    {
                        "_id": "ga4gh:VSL.N1dI0SlDciU-8OsmdjEfSKdVDp-Q5_lf",
                        "type": "SequenceLocation",
                        "sequence_id": "ga4gh:SQ.7B7SHsmchAR0dFcDCuSFjJAo7tX87krQ",
                        "interval": {
                            "start": {"type": "Number", "value": 31212297},
                            "end": {"type": "Number", "value": 31280080},
                            "type": "SequenceInterval"
                        }
                    }
                ]
            },
            {
                "type": "Extension",
                "name": "hgnc_locus_type",
                "value": "gene with protein product"
            },
            {
                "type": "Extension",
                "name": "ncbi_gene_type",
                "value": "protein-coding"
            },
            {
                "type": "Extension",
                "name": "ensembl_biotype",
                "value": "protein_coding"
            }
        ]
    }


@pytest.fixture(scope="module")
def dis3_p63a():
    """Create DIS3 P63A test fixture."""
    params = {
        "id": "normalize.variation:DIS3%20P63A",
        "type": "VariationDescriptor",
        "variation_id": "ga4gh:VA.ueB4LJE1wba68a1w1w6Mrb-MOquy8n3C",
        "variation": {
            "_id": "ga4gh:VA.ueB4LJE1wba68a1w1w6Mrb-MOquy8n3C",
            "location": {
                "_id": "ga4gh:VSL.JcyesAEvndeQYxn4Gzi68hqC1ANZnrGN",
                "interval": {
                    "end": {"value": 63, "type": "Number"},
                    "start": {"value": 62, "type": "Number"},
                    "type": "SequenceInterval"
                },
                "sequence_id": "ga4gh:SQ.mlWsxfPKINN3o300stAI8oqN5U7P6kEu",
                "type": "SequenceLocation"
            },
            "state": {
                "sequence": "A",
                "type": "LiteralSequenceExpression"
            },
            "type": "Allele"
        },
        "molecule_context": "protein",
        "structural_type": "SO:0001606",
        "gene_context": {
            "id": "normalize.gene:DIS3",
            "type": "GeneDescriptor",
            "label": "DIS3",
            "xrefs": [
                "ensembl:ENSG00000083520",
                "ncbigene:22894"
            ],
            "alternate_labels": [
                "dis3p",
                "RRP44",
                "KIAA1008",
                "2810028N01Rik",
                "EXOSC11"
            ],
            "extensions": [
                {
                    "name": "symbol_status",
                    "value": "approved",
                    "type": "Extension"
                },
                {
                    "name": "approved_name",
                    "value": "DIS3 homolog, exosome endoribonuclease and 3'-5' exoribonuclease",  # noqa: E501
                    "type": "Extension"
                },
                {
                    "name": "hgnc_locations",
                    "value": [
                        {
                            "species_id": "taxonomy:9606",
                            "interval": {
                                "type": "CytobandInterval",
                                "start": "q21.33",
                                "end": "q21.33"
                            },
                            "_id": "ga4gh:VCL.84IPub_nKl33cWX9pNoPeGsyeVuJnyra",
                            "type": "ChromosomeLocation",
                            "chr": "13"
                        }
                    ],
                    "type": "Extension"
                },
                {
                    "type": "Extension",
                    "name": "ensembl_locations",
                    "value": [
                        {
                            "_id": "ga4gh:VSL.HZtod8n11kD7jCAbtsJLDGwLKCEhRWO1",
                            "type": "SequenceLocation",
                            "sequence_id": "ga4gh:SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT",
                            "interval": {
                                "start": {"type": "Number", "value": 72752168},
                                "end": {"type": "Number", "value": 72782096},
                                "type": "SequenceInterval"
                            }
                        }
                    ]
                },
                {
                    "type": "Extension",
                    "name": "ncbi_locations",
                    "value": [
                        {
                            "_id": "ga4gh:VCL.84IPub_nKl33cWX9pNoPeGsyeVuJnyra",
                            "type": "ChromosomeLocation",
                            "species_id": "taxonomy:9606",
                            "chr": "13",
                            "interval": {
                                "end": "q21.33",
                                "start": "q21.33",
                                "type": "CytobandInterval"
                            }
                        },
                        {
                            "_id": "ga4gh:VSL.BIRwPLT8rtyOnhd3aUXaz4xlHC4P4zA8",
                            "type": "SequenceLocation",
                            "sequence_id": "ga4gh:SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT",
                            "interval": {
                                "start": {"type": "Number", "value": 72752168},
                                "end": {"type": "Number", "value": 72781900},
                                "type": "SequenceInterval"
                            }
                        }
                    ]
                },
                {
                    "name": "associated_with",
                    "value": [
                        "vega:OTTHUMG00000017070",
                        "ccds:CCDS9447",
                        "orphanet:470196",
                        "ena.embl:AB023225",
                        "ccds:CCDS45057",
                        "omim:607533",
                        "pubmed:11935316",
                        "refseq:NM_014953",
                        "uniprot:Q9Y2L1",
                        "ccds:CCDS81772",
                        "ucsc:uc001vix.6",
                        "pubmed:9562621"
                    ],
                    "type": "Extension"
                },
                {
                    "name": "previous_symbols",
                    "value": [
                        "KIAA1008"
                    ],
                    "type": "Extension"
                },
                {
                    "type": "Extension",
                    "name": "hgnc_locus_type",
                    "value": "gene with protein product"
                },
                {
                    "type": "Extension",
                    "name": "ncbi_gene_type",
                    "value": "protein-coding"
                },
                {
                    "type": "Extension",
                    "name": "ensembl_biotype",
                    "value": "protein_coding"
                }
            ],
            "gene_id": "hgnc:20604"
        },
        "vrs_ref_allele_seq": "P"
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def vhl(vhl_gene_context):
    """Create VHL Tyr185Ter fixture."""
    params = {
        "id": "normalize.variation:NP_000542.1%3Ap.Tyr185Ter",
        "type": "VariationDescriptor",
        "variation_id": "ga4gh:VA._S0nFwX4Y2FPmv5Radf01DAsxQbxA2cc",
        "variation": {
            "_id": "ga4gh:VA._S0nFwX4Y2FPmv5Radf01DAsxQbxA2cc",
            "location": {
                "_id": "ga4gh:VSL._P3rBWI3f7OBs3a4gvZ18QJ6f6dSfqEQ",
                "interval": {
                    "end": {"value": 185, "type": "Number"},
                    "start": {"value": 184, "type": "Number"},
                    "type": "SequenceInterval"
                },
                "sequence_id": "ga4gh:SQ.z-Oa0pZkJ6GHJHOYM7h5mY_umc0SJzTu",
                "type": "SequenceLocation"
            },
            "state": {
                "sequence": "*",
                "type": "LiteralSequenceExpression"
            },
            "type": "Allele"
        },
        "molecule_context": "protein",
        "structural_type": "SO:0001617",
        "vrs_ref_allele_seq": "Y",
        "gene_context": vhl_gene_context
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def braf_v600e_nucleotide(braf_gene_context, braf_nuc_value):
    """Create a test fixture for BRAF V600E MANE select nucleotide hgvs."""
    variation = copy.deepcopy(braf_nuc_value)
    variation["_id"] = "ga4gh:VA.bsmS0aeGSffgemmymBrpTjTSW2qYeIJE"
    params = {
        "id": "normalize.variation:NM_004333.4%3Ac.1799T%3EA",
        "type": "VariationDescriptor",
        "variation_id": variation["_id"],
        "variation": variation,
        "molecule_context": "transcript",
        "structural_type": "SO:0001483",
        "vrs_ref_allele_seq": "T",
        "gene_context": braf_gene_context
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def nm_004448_coding_dna_delins(erbb2_context):
    """Create test fixture for NM_004448.4:c.2326_2327delinsCT."""
    params = {
        "id": "normalize.variation:NM_004448.4%3Ac.2326_2327delinsCT",
        "type": "VariationDescriptor",
        "variation_id": "ga4gh:VA.eMxxAEjNduAvg5U3eBZxf0nLtfcMNxqy",
        "variation": {
            "_id": "ga4gh:VA.eMxxAEjNduAvg5U3eBZxf0nLtfcMNxqy",
            "location": {
                "_id": "ga4gh:VSL.bBzTvpLChbWE2SZ7X0drm8NQj5rzNqTK",
                "interval": {
                    "end": {"value": 2502, "type": "Number"},
                    "start": {"value": 2500, "type": "Number"},
                    "type": "SequenceInterval"
                },
                "sequence_id": "ga4gh:SQ.y9b4LVMiCXpZxOg9Xt1NwRtssA03MwWM",
                "type": "SequenceLocation"
            },
            "state": {
                "sequence": "CT",
                "type": "LiteralSequenceExpression"
            },
            "type": "Allele"
        },
        "molecule_context": "transcript",
        "structural_type": "SO:1000032",
        "vrs_ref_allele_seq": "GG",
        "gene_context": erbb2_context
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def nc_000007_genomic_delins(braf_gene_context):
    """Create test fixture for NC_000007.13:g.140453135_140453136delinsAT."""
    params = {
        "id": "normalize.variation:NC_000007.13%3Ag.140453135_140453136delinsAT",  # noqa: E501
        "type": "VariationDescriptor",
        "variation_id": "ga4gh:VA.4387UZ6Yssh3XCGKjm71z_WtadpBZT3O",
        "variation": {
            "_id": "ga4gh:VA.4387UZ6Yssh3XCGKjm71z_WtadpBZT3O",
            "location": {
                "_id": "ga4gh:VSL.6PeoFwkO4ISmUjDWoYLkVsATVx8JRApd",
                "interval": {
                    "end": {"value": 2146, "type": "Number"},
                    "start": {"value": 2144, "type": "Number"},
                    "type": "SequenceInterval"
                },
                "sequence_id": "ga4gh:SQ.I_0feOk5bZ3VfH8ejhWQiMDe9o6o4QdR",
                "type": "SequenceLocation"
            },
            "state": {
                "sequence": "AT",
                "type": "LiteralSequenceExpression"
            },
            "type": "Allele"
        },
        "molecule_context": "transcript",
        "structural_type": "SO:1000032",
        "vrs_ref_allele_seq": "TG",
        "gene_context": braf_gene_context
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def nm_000551(vhl_gene_context):
    """Create test fixture for NM_000551.4:c.615delinsAA."""
    params = {
        "id": "normalize.variation:temp",
        "type": "VariationDescriptor",
        "variation_id": "ga4gh:VA._JN_AF5PO9kWKgKxB5T48cypZl7ccEsQ",
        "variation": {
            "_id": "ga4gh:VA._JN_AF5PO9kWKgKxB5T48cypZl7ccEsQ",
            "location": {
                "_id": "ga4gh:VSL.tQRFfWMPPHErWSbKvoEIh2gI2ehHIYAs",
                "interval": {
                    "end": {"value": 685, "type": "Number"},
                    "start": {"value": 684, "type": "Number"},
                    "type": "SequenceInterval"
                },
                "sequence_id": "ga4gh:SQ.xBKOKptLLDr-k4hTyCetvARn16pDS_rW",
                "type": "SequenceLocation"
            },
            "state": {
                "sequence": "AA",
                "type": "LiteralSequenceExpression"
            },
            "type": "Allele"
        },
        "molecule_context": "transcript",
        "structural_type": "SO:1000032",
        "vrs_ref_allele_seq": "C",
        "gene_context": vhl_gene_context
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def braf_nuc_value():
    """Create test fixture for BRAF V600E value on c. coordinate."""
    return {
        "location": {
            "_id": "ga4gh:VSL.WnQ4luBGVxa-tMWza5nPzFh1bIeCij3T",
            "interval": {
                "end": {"value": 2025, "type": "Number"},
                "start": {"value": 2024, "type": "Number"},
                "type": "SequenceInterval"
            },
            "sequence_id": "ga4gh:SQ.aKMPEJgmlZXt_F6gRY5cUG3THH2n-GUa",
            "type": "SequenceLocation"
        },
        "state": {
            "sequence": "A",
            "type": "LiteralSequenceExpression"
        },
        "type": "Allele"
    }


@pytest.fixture(scope="module")
def coding_dna_silent_mutation(braf_gene_context, braf_nuc_value):
    """Create test fixture for NM_004333.4:c.1799=."""
    value = copy.deepcopy(braf_nuc_value)
    value["state"]["sequence"] = "T"
    value["_id"] = "ga4gh:VA.kbgr3KgaDXI6jUR8r7A53xlNSN2cgOl8"
    params = {
        "id": "normalize.variation:NM_004333.4%3Ac.1799%3D",
        "type": "VariationDescriptor",
        "variation_id": value["_id"],
        "variation": value,
        "molecule_context": "transcript",
        "structural_type": "SO:0002073",
        "vrs_ref_allele_seq": "T",
        "gene_context": braf_gene_context
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def nc_000007_silent_mutation(braf_gene_context, braf_nuc_value):
    """Create test fixture for NC_000007.13:g.140453136=."""
    value = copy.deepcopy(braf_nuc_value)
    value["state"]["sequence"] = "T"
    value["_id"] = "ga4gh:VA.kbgr3KgaDXI6jUR8r7A53xlNSN2cgOl8"
    params = {
        "id": "normalize.variation:NC_000007.13%3Ag.140453136%3D",
        "type": "VariationDescriptor",
        "variation_id": value["_id"],
        "variation": value,
        "molecule_context": "transcript",
        "structural_type": "SO:0002073",
        "vrs_ref_allele_seq": "T",
        "gene_context": braf_gene_context
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def protein_delins(egfr_context):
    """Create test fixture for protein delins."""
    params = {
        "id": "normalize.variation:NP_001333827.1%3Ap.Leu747_Thr751delinsPro",
        "type": "VariationDescriptor",
        "variation_id": "ga4gh:VA.eDMXxJw9shlSKF3znIg5abniGoyJ3GQ4",
        "variation": {
            "_id": "ga4gh:VA.eDMXxJw9shlSKF3znIg5abniGoyJ3GQ4",
            "location": {
                "_id": "ga4gh:VSL.Mm8duqYDJyel5ZnwScnxLyGH1i9lcl3T",
                "interval": {
                    "end": {"value": 751, "type": "Number"},
                    "start": {"value": 746, "type": "Number"},
                    "type": "SequenceInterval"
                },
                "sequence_id": "ga4gh:SQ.vyo55F6mA6n2LgN4cagcdRzOuh38V4mE",
                "type": "SequenceLocation"
            },
            "state": {
                "sequence": "P",
                "type": "LiteralSequenceExpression"
            },
            "type": "Allele"
        },
        "molecule_context": "protein",
        "structural_type": "SO:1000032",
        "vrs_ref_allele_seq": "LREAT",
        "gene_context": egfr_context
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def coding_dna_deletion(erbb2_context):
    """Create test fixture for coding dna deletion range with deleted
    sequence.
    """
    params = {
        "id": "normalize.variation:NM_004448.3%3Ac.2264_2278delTGAGGGAAAACACAT",  # noqa: E501
        "type": "VariationDescriptor",
        "variation_id": "ga4gh:VA.NUCURWYivhjC4oyBtzgJZ27SaaMY08Q7",
        "variation": {
            "_id": "ga4gh:VA.NUCURWYivhjC4oyBtzgJZ27SaaMY08Q7",
            "location": {
                "_id": "ga4gh:VSL.3uPWAjsdzd8MbAqw8DV46eBLK8tQRyEs",
                "interval": {
                    "end": {"value": 2453, "type": "Number"},
                    "start": {"value": 2437, "type": "Number"},
                    "type": "SequenceInterval"
                },
                "sequence_id": "ga4gh:SQ.y9b4LVMiCXpZxOg9Xt1NwRtssA03MwWM",
                "type": "SequenceLocation"
            },
            "state": {
                "sequence": "T",
                "type": "LiteralSequenceExpression"
            },
            "type": "Allele"
        },
        "molecule_context": "transcript",
        "structural_type": "SO:0000159",
        "vrs_ref_allele_seq": "TTGAGGGAAAACACAT",
        "gene_context": erbb2_context
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def coding_dna_insertion(limk2_gene_context):
    """Create test fixture for coding DNA insertion."""
    params = {
        "id": "normalize.variation:ENST00000331728.9%3Ac.2049_2050insA",
        "type": "VariationDescriptor",
        "variation_id": "ga4gh:VA.RAcEF24tRNB_J5Hz9E6GY-qQSi4ZG932",
        "variation": {
            "_id": "ga4gh:VA.RAcEF24tRNB_J5Hz9E6GY-qQSi4ZG932",
            "location": {
                "_id": "ga4gh:VSL.WKymELRFwu4LDDmH5ci5Ip0M3XA2RObr",
                "interval": {
                    "end": {"value": 2160, "type": "Number"},
                    "start": {"value": 2160, "type": "Number"},
                    "type": "SequenceInterval"
                },
                "sequence_id": "ga4gh:SQ.7_mlQyDN-uWH0RlxTQFvFEv6ykd2D-xF",
                "type": "SequenceLocation"
            },
            "state": {
                "sequence": "A",
                "type": "LiteralSequenceExpression"
            },
            "type": "Allele"
        },
        "molecule_context": "transcript",
        "structural_type": "SO:0000667",
        "gene_context": limk2_gene_context
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_insertion(erbb2_context):
    """Create a gene insertion test fixture."""
    params = {
        "id": "normalize.variation:NC_000017.10%3Ag.37880993_37880994insGCTTACGTGATG",  # noqa: E501
        "type": "VariationDescriptor",
        "variation_id": "ga4gh:VA.nHB0_mpsq2t90S-znr81oCi2cY5CMdUe",
        "variation": {
            "_id": "ga4gh:VA.nHB0_mpsq2t90S-znr81oCi2cY5CMdUe",
            "location": {
                "_id": "ga4gh:VSL.E0o4HCXjy1EUthF1m32oj_Bc45g5YmEm",
                "interval": {
                    "end": {"value": 2500, "type": "Number"},
                    "start": {"value": 2488, "type": "Number"},
                    "type": "SequenceInterval"
                },
                "sequence_id": "ga4gh:SQ.y9b4LVMiCXpZxOg9Xt1NwRtssA03MwWM",
                "type": "SequenceLocation"
            },
            "state": {
                "sequence": "TACGTGATGGCTTACGTGATGGCT",
                "type": "LiteralSequenceExpression"
            },
            "type": "Allele"
        },
        "molecule_context": "transcript",
        "structural_type": "SO:0000667",
        "vrs_ref_allele_seq": "TACGTGATGGCT",
        "gene_context": erbb2_context
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_substitution(egfr_context):
    """Create a gene insertion test fixture."""
    params = {
        "id": "normalize.variation:NC_000007.13%3Ag.55249071C%3ET",
        "type": "VariationDescriptor",
        "variation_id": "ga4gh:VA.VkcuqgqMuSQeq8Hy0VPOGRIeyr8uSBV2",
        "variation": {
            "_id": "ga4gh:VA.VkcuqgqMuSQeq8Hy0VPOGRIeyr8uSBV2",
            "location": {
                "_id": "ga4gh:VSL.G1gIZ-om-8Exl3F0ZLxXYY8CjliwCaO1",
                "interval": {
                    "end": {"value": 2630, "type": "Number"},
                    "start": {"value": 2629, "type": "Number"},
                    "type": "SequenceInterval"
                },
                "sequence_id": "ga4gh:SQ.d_QsP29RWJi6bac7GOC9cJ9AO7s_HUMN",
                "type": "SequenceLocation"
            },
            "state": {
                "sequence": "T",
                "type": "LiteralSequenceExpression"
            },
            "type": "Allele"
        },
        "molecule_context": "transcript",
        "structural_type": "SO:0001483",
        "vrs_ref_allele_seq": "C",
        "gene_context": egfr_context
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_sub_grch38():
    """Create a genomic substitution GRCh38 test fixture."""
    params = {
        "id": "normalize.variation:NC_000007.13%3Ag.55249071C%3ET",
        "type": "VariationDescriptor",
        "variation_id": "ga4gh:VA.1ewlywoD423K7YH_K4YefZg6J_87pQTp",
        "variation": {
            "_id": "ga4gh:VA.1ewlywoD423K7YH_K4YefZg6J_87pQTp",
            "location": {
                "_id": "ga4gh:VSL.0p1nWj9-sryfUD5jvPTZZdnZeiHVHXls",
                "interval": {
                    "end": {"value": 55181378, "type": "Number"},
                    "start": {"value": 55181377, "type": "Number"},
                    "type": "SequenceInterval"
                },
                "sequence_id": "ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
                "type": "SequenceLocation"
            },
            "state": {
                "sequence": "T",
                "type": "LiteralSequenceExpression"
            },
            "type": "Allele"
        },
        "molecule_context": "genomic",
        "structural_type": "SO:0001483",
        "vrs_ref_allele_seq": "C"
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def egfr_grch38_sub(genomic_sub_grch38, egfr_context):
    """Create a genomic substitution GRCh38 test fixture."""
    params = {
        "id": "normalize.variation:NC_000007.13%3Ag.55249071C%3ET",
        "type": "VariationDescriptor",
        "variation_id": "ga4gh:VA.1ewlywoD423K7YH_K4YefZg6J_87pQTp",
        "variation": {
            "_id": "ga4gh:VA.1ewlywoD423K7YH_K4YefZg6J_87pQTp",
            "location": {
                "_id": "ga4gh:VSL.0p1nWj9-sryfUD5jvPTZZdnZeiHVHXls",
                "interval": {
                    "end": {"value": 55181378, "type": "Number"},
                    "start": {"value": 55181377, "type": "Number"},
                    "type": "SequenceInterval"
                },
                "sequence_id": "ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
                "type": "SequenceLocation"
            },
            "state": {
                "sequence": "T",
                "type": "LiteralSequenceExpression"
            },
            "type": "Allele"
        },
        "molecule_context": "genomic",
        "structural_type": "SO:0001483",
        "vrs_ref_allele_seq": "C",
        "gene_context": egfr_context
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_uncertain_del_x():
    """Create a genomic uncertain deletion on chr X test fixture."""
    params = {
        "id": "normalize.variation:NC_000023.11%3Ag.%28%3F_31120496%29_%2833339477_%3F%29del",  # noqa: E501
        "type": "VariationDescriptor",
        "variation_id": "ga4gh:VAC.oCrOJZcSJ-knvR1l6u-wpPNNMCLDCB0L",
        "variation": {
            "_id": "ga4gh:VAC.oCrOJZcSJ-knvR1l6u-wpPNNMCLDCB0L",
            "subject": {
                "location": {
                    "_id": "ga4gh:VSL.7OJ5EFgu_2C4zPFDUBgn-ziE6BZwsRcv",
                    "sequence_id": "ga4gh:SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
                    "interval": {
                        "start": {
                            "value": 31120495,
                            "comparator": "<=",
                            "type": "IndefiniteRange"
                        },
                        "end": {
                            "value": 33339477,
                            "comparator": ">=",
                            "type": "IndefiniteRange"
                        },
                        "type": "SequenceInterval"
                    },
                    "type": "SequenceLocation"
                },
                "reverse_complement": False,
                "type": "DerivedSequenceExpression"
            },
            "copies": {
                "min": 0,
                "max": 1,
                "type": "DefiniteRange"
            },
            "type": "AbsoluteCopyNumber"
        },
        "molecule_context": "genomic",
        "structural_type": "SO:0001743"
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def grch38_braf_genom_sub(braf_v600e_genomic_sub):
    """Create a genomic substitution GRCh38 test fixture for BRAF."""
    params = {
        "id": "normalize.variation:NC_000007.13%3Ag.140453136A%3ET",
        "type": "VariationDescriptor",
        "variation_id": "ga4gh:VA.fZiBjQEolbkL0AxjoTZf4SOkFy9J0ebU",
        "variation": braf_v600e_genomic_sub,
        "molecule_context": "genomic",
        "structural_type": "SO:0001483",
        "vrs_ref_allele_seq": "A"
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def grch38_braf_genom_silent_mutation():
    """Create a genomic silent mutation GRCh38 test fixture for BRAF."""
    params = {
        "id": "normalize.variation:NC_000007.13%3Ag.140453136%3D",
        "type": "VariationDescriptor",
        "variation_id": "ga4gh:VA.aMwnr5rEbtPQe5gXDDO2gZO_zSqN2RmH",
        "variation": {
            "_id": "ga4gh:VA.aMwnr5rEbtPQe5gXDDO2gZO_zSqN2RmH",
            "location": {
                "_id": "ga4gh:VSL.zga82-TpYiNmBESCfvDvAz9DyvJF98I-",
                "interval": {
                    "end": {"value": 140753336, "type": "Number"},
                    "start": {"value": 140753335, "type": "Number"},
                    "type": "SequenceInterval"
                },
                "sequence_id": "ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
                "type": "SequenceLocation"
            },
            "state": {
                "sequence": "A",
                "type": "LiteralSequenceExpression"
            },
            "type": "Allele"
        },
        "molecule_context": "genomic",
        "structural_type": "SO:0002073",
        "vrs_ref_allele_seq": "A"
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def grch38_genomic_delins1():
    """Create a test fixture for NC_000007.13:g.140453135_140453136delinsAT."""
    params = {
        "id":
            "normalize.variation:NC_000007.13%3Ag.140453135_140453136delinsAT",
        "type": "VariationDescriptor",
        "variation_id": "ga4gh:VA.mlJVnI7js6Tsb2GSLFlNRbCKE9zFRX5p",
        "variation": {
            "_id": "ga4gh:VA.mlJVnI7js6Tsb2GSLFlNRbCKE9zFRX5p",
            "location": {
                "_id": "ga4gh:VSL.b0Ldj2KcT2k0n0PZfqHCBH1YzQZYceYX",
                "interval": {
                    "end": {"value": 140753336, "type": "Number"},
                    "start": {"value": 140753334, "type": "Number"},
                    "type": "SequenceInterval"
                },
                "sequence_id": "ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
                "type": "SequenceLocation"
            },
            "state": {
                "sequence": "AT",
                "type": "LiteralSequenceExpression"
            },
            "type": "Allele"
        },
        "molecule_context": "genomic",
        "structural_type": "SO:1000032",
        "vrs_ref_allele_seq": "CA"
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def grch38_genomic_delins2():
    """Create a test fixture for NC_000003.12:g.10149938delinsAA."""
    params = {
        "id": "normalize.variation:NC_000003.12%3Ag.10149938delinsAA",
        "type": "VariationDescriptor",
        "variation_id": "ga4gh:VA.ndCPwzek-KU626kK28bQd1gWAAk2ELze",
        "variation": {
            "_id": "ga4gh:VA.ndCPwzek-KU626kK28bQd1gWAAk2ELze",
            "location": {
                "_id": "ga4gh:VSL.mVulVOKoX2frLH1XTIJfpvJa6RGdOVu_",
                "interval": {
                    "start": {"value": 10149937, "type": "Number"},
                    "end": {"value": 10149938, "type": "Number"},
                    "type": "SequenceInterval"
                },
                "sequence_id": "ga4gh:SQ.Zu7h9AggXxhTaGVsy7h_EZSChSZGcmgX",
                "type": "SequenceLocation"
            },
            "state": {
                "sequence": "AA",
                "type": "LiteralSequenceExpression"
            },
            "type": "Allele"
        },
        "molecule_context": "genomic",
        "structural_type": "SO:1000032",
        "vrs_ref_allele_seq": "C"
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def grch38_genomic_insertion(grch38_genomic_insertion_variation):
    """Create a test fixture for NC_000017.10:g.37880993_37880994insGCTTACGTGATG."""
    params = {
        "id": "normalize.variation:NC_000017.10%3Ag.37880993_37880994insGCTTACGTGATG",
        "type": "VariationDescriptor",
        "variation_id": "ga4gh:VA.tCjV190dUsV7tSjdR8qOLSQIR7Hr8VMe",
        "variation": grch38_genomic_insertion_variation,
        "molecule_context": "genomic",
        "structural_type": "SO:0000667",
        "vrs_ref_allele_seq": "TACGTGATGGCT"
    }
    return VariationDescriptor(**params)


@pytest.mark.asyncio
async def test_protein_substitution(test_handler, braf_v600e, dis3_p63a):
    """Test that protein substitutions normalize correctly."""
    resp = await test_handler.normalize("     BRAF      V600E    ")
    assertion_checks(resp.variation_descriptor, braf_v600e, "BRAF      V600E")

    braf_id = "normalize.variation:BRAF%20V600E"

    resp = await test_handler.normalize("NP_004324.2:p.Val600Glu")
    assert resp.variation_descriptor.id == \
        "normalize.variation:NP_004324.2%3Ap.Val600Glu"
    resp.variation_descriptor.id = braf_id
    assertion_checks(resp.variation_descriptor, braf_v600e, "NP_004324.2:p.Val600Glu")

    resp = await test_handler.normalize("braf v512e")
    assert resp.variation_descriptor.id == "normalize.variation:braf%20v512e"
    resp.variation_descriptor.id = braf_id
    assertion_checks(resp.variation_descriptor, braf_v600e, "braf v512e")

    resp = await test_handler.normalize(" NP_001365404.1:p.Val512Glu  ")
    assert resp.variation_descriptor.id == \
        "normalize.variation:NP_001365404.1%3Ap.Val512Glu"
    resp.variation_descriptor.id = braf_id
    assertion_checks(resp.variation_descriptor, braf_v600e,
                     "NP_001365404.1:p.Val512Glu")

    resp = await test_handler.normalize("DIS3 P63A")
    assertion_checks(resp.variation_descriptor, dis3_p63a, "DIS3 P63A")


@pytest.mark.asyncio
async def test_polypeptide_truncation(test_handler, vhl):
    """Test that polypeptide truncations normalize correctly."""
    resp = await test_handler.normalize("NP_000542.1:p.Tyr185Ter")
    assertion_checks(resp.variation_descriptor, vhl, "NP_000542.1:p.Tyr185Ter")


@pytest.mark.asyncio
async def test_silent_mutation(test_handler, vhl_silent):
    """Test that silent mutations normalize correctly."""
    resp = await test_handler.normalize("NP_000542.1:p.Pro61=")
    assertion_checks(resp.variation_descriptor, vhl_silent, "NP_000542.1:p.Pro61=")


@pytest.mark.asyncio
async def test_coding_dna_and_genomic_substitution(
        test_handler, braf_v600e_nucleotide, genomic_substitution,
        genomic_sub_grch38, egfr_grch38_sub, grch38_braf_genom_sub):
    """Test that coding dna and genomic substitutions normalize correctly."""
    resp = await test_handler.normalize("NM_004333.4:c.1799T>A")
    assertion_checks(resp.variation_descriptor, braf_v600e_nucleotide,
                     "NM_004333.4:c.1799T>A")

    # MANE transcript
    refseq_id = "normalize.variation:NM_004333.4%3Ac.1799T%3EA"

    # TODO: Check if this should return a different VRS object?
    resp = await test_handler.normalize("ENST00000288602.10:c.1799T>A")
    assert resp.variation_descriptor.id == \
        "normalize.variation:ENST00000288602.10%3Ac.1799T%3EA"
    resp.variation_descriptor.id = refseq_id
    assertion_checks(resp.variation_descriptor, braf_v600e_nucleotide,
                     "ENST00000288602.10:c.1799T>A")

    resp = await test_handler.normalize("BRAF V600E c.1799T>A")
    assert resp.variation_descriptor.id == \
        "normalize.variation:BRAF%20V600E%20c.1799T%3EA"
    resp.variation_descriptor.id = refseq_id
    assertion_checks(resp.variation_descriptor, braf_v600e_nucleotide,
                     "BRAF V600E c.1799T>A")

    resp = await test_handler.normalize("BRAF V600E (c.1799T>A)")
    assert resp.variation_descriptor.id == \
        "normalize.variation:BRAF%20V600E%20%28c.1799T%3EA%29"
    resp.variation_descriptor.id = refseq_id
    assertion_checks(resp.variation_descriptor, braf_v600e_nucleotide,
                     "BRAF V600E (c.1799T>A)")

    resp = await test_handler.normalize("BRAF c.1799T>A")
    assert resp.variation_descriptor.id == "normalize.variation:BRAF%20c.1799T%3EA"
    resp.variation_descriptor.id = refseq_id
    assertion_checks(resp.variation_descriptor, braf_v600e_nucleotide, "BRAF c.1799T>A")

    resp = await test_handler.normalize("NC_000007.13:g.140453136A>T")
    assertion_checks(resp.variation_descriptor, grch38_braf_genom_sub,
                     "NC_000007.13:g.140453136A>T")

    fixture_id = "normalize.variation:NC_000007.13%3Ag.140453136A%3ET"
    resp = await test_handler.normalize("7-140453136-A-T")  # 37
    assert resp.variation_descriptor.id == "normalize.variation:7-140453136-A-T"
    resp.variation_descriptor.id = fixture_id
    assertion_checks(resp.variation_descriptor, grch38_braf_genom_sub,
                     "7-140453136-A-T")

    resp = await test_handler.normalize("7-140753336-A-T")  # 38
    assert resp.variation_descriptor.id == "normalize.variation:7-140753336-A-T"
    resp.variation_descriptor.id = fixture_id
    assertion_checks(resp.variation_descriptor, grch38_braf_genom_sub,
                     "7-140753336-A-T")

    resp = await test_handler.normalize("BRAF V600E (g.140453136A>T)")
    assert resp.variation_descriptor.id == "normalize.variation:BRAF%20V600E%20%28g.140453136A%3ET%29"  # noqa: E501
    resp.variation_descriptor.id = refseq_id
    assertion_checks(resp.variation_descriptor, braf_v600e_nucleotide,
                     "BRAF V600E (g.140453136A>T)")

    resp = await test_handler.normalize("BRAF g.140453136A>T")
    assert resp.variation_descriptor.id == "normalize.variation:BRAF%20g.140453136A%3ET"
    resp.variation_descriptor.id = refseq_id
    assertion_checks(resp.variation_descriptor, braf_v600e_nucleotide,
                     "BRAF g.140453136A>T")

    # More than 1 gene (EGFR and EGFR-AS1)
    resp = await test_handler.normalize("NC_000007.13:g.55249071C>T")
    assertion_checks(resp.variation_descriptor, genomic_sub_grch38,
                     "NC_000007.13:g.55249071C>T")

    resp = await test_handler.normalize("EGFR g.55249071C>T")
    assert resp.variation_descriptor.id == "normalize.variation:EGFR%20g.55249071C%3ET"
    resp.variation_descriptor.id = "normalize.variation:NC_000007.13%3Ag.55249071C%3ET"
    assertion_checks(resp.variation_descriptor, genomic_substitution,
                     "EGFR g.55249071C>T")


@pytest.mark.asyncio
async def test_coding_dna_silent_mutation(test_handler,
                                          coding_dna_silent_mutation,
                                          braf_gene_context):
    """Test that Coding DNA Silent Mutation normalizes correctly."""
    resp = await test_handler.normalize("NM_004333.4:c.1799= ")
    assertion_checks(resp.variation_descriptor, coding_dna_silent_mutation,
                     "NM_004333.4:c.1799=")

    fixture_id = "normalize.variation:NM_004333.4%3Ac.1799%3D"

    resp = await test_handler.normalize("ENST00000288602.11:c.1799=")
    assert resp.variation_descriptor.id == \
        "normalize.variation:ENST00000288602.11%3Ac.1799%3D"
    resp.variation_descriptor.id = fixture_id
    assertion_checks(resp.variation_descriptor, coding_dna_silent_mutation,
                     "ENST00000288602.11:c.1799=")

    # TODO: What to do for older Ensembl transcripts that aren"t found
    #  in seqrepo or UTA
    # resp = await test_handler.normalize("ENST00000288602.6:c.1799=")
    # assert_coding_dna_genomic_silent_mutation(resp, braf_gene_context,
    #                                           1798, 1799)
    # assert resp.variation_descriptor.id == "normalize.variation:ENST00000288602.6%3Ac.1799%3D"  # noqa: E501
    # assert resp.variation_descriptor.label == "ENST00000288602.6:c.1799="
    # assert resp.variation_descriptor.molecule_context == "transcript"

    resp = await test_handler.normalize("BRAF    c.1799=")
    assert resp.variation_descriptor.id == "normalize.variation:BRAF%20c.1799%3D"
    resp.variation_descriptor.id = fixture_id
    assertion_checks(resp.variation_descriptor, coding_dna_silent_mutation,
                     "BRAF    c.1799=")

    resp = await test_handler.normalize("  BRAF  V600E  c.1799=  ")
    assert resp.variation_descriptor.id == \
        "normalize.variation:BRAF%20V600E%20c.1799%3D"
    resp.variation_descriptor.id = fixture_id
    assertion_checks(resp.variation_descriptor, coding_dna_silent_mutation,
                     "BRAF  V600E  c.1799=")


@pytest.mark.asyncio
async def test_genomic_silent_mutation(test_handler, nc_000007_silent_mutation,
                                       braf_gene_context,
                                       grch38_braf_genom_silent_mutation):
    """Test that genomic silent mutation normalizes correctly."""
    resp = await test_handler.normalize("NC_000007.13:g.140453136=")
    assertion_checks(resp.variation_descriptor, grch38_braf_genom_silent_mutation,
                     "NC_000007.13:g.140453136=")

    fixture_id = "normalize.variation:NC_000007.13%3Ag.140453136%3D"
    resp = await test_handler.normalize("7-140453136-A-A")
    assert resp.variation_descriptor.id == "normalize.variation:7-140453136-A-A"
    resp.variation_descriptor.id = fixture_id
    assertion_checks(resp.variation_descriptor, grch38_braf_genom_silent_mutation,
                     "7-140453136-A-A")

    resp = await test_handler.normalize("7-140753336-A-A")
    assert resp.variation_descriptor.id == "normalize.variation:7-140753336-A-A"
    resp.variation_descriptor.id = fixture_id
    assertion_checks(resp.variation_descriptor, grch38_braf_genom_silent_mutation,
                     "7-140753336-A-A")

    resp = await test_handler.normalize("BRAF g.140453136=")
    assert resp.variation_descriptor.id == "normalize.variation:BRAF%20g.140453136%3D"
    resp.variation_descriptor.id = fixture_id
    assertion_checks(resp.variation_descriptor, nc_000007_silent_mutation,
                     "BRAF g.140453136=")


@pytest.mark.asyncio
async def test_coding_dna_delins(test_handler, nm_004448_coding_dna_delins,
                                 nm_000551):
    """Test that Coding DNA DelIns normalizes correctly."""
    resp = await test_handler.normalize("    NM_004448.4:c.2326_2327delinsCT    ")
    assertion_checks(resp.variation_descriptor, nm_004448_coding_dna_delins,
                     "NM_004448.4:c.2326_2327delinsCT")

    # TODO: Test ENST###.c

    resp = await test_handler.normalize("NM_000551.3:c.615delinsAA")
    nm_000551.id = "normalize.variation:NM_000551.3%3Ac.615delinsAA"
    assertion_checks(resp.variation_descriptor, nm_000551, "NM_000551.3:c.615delinsAA")


@pytest.mark.asyncio
async def test_genomic_delins(test_handler, nc_000007_genomic_delins,
                              nm_000551, grch38_genomic_delins1,
                              grch38_genomic_delins2):
    """Test that Genomic DelIns normalizes correctly."""
    resp = await test_handler.normalize(
        "NC_000007.13:g.140453135_140453136delinsAT"
    )
    assertion_checks(resp.variation_descriptor, grch38_genomic_delins1,
                     "NC_000007.13:g.140453135_140453136delinsAT")

    resp = await test_handler.normalize("NC_000003.12:g.10149938delinsAA")
    assertion_checks(resp.variation_descriptor, grch38_genomic_delins2,
                     "NC_000003.12:g.10149938delinsAA")


@pytest.mark.asyncio
async def test_protein_delins(test_handler, protein_delins):
    """Test that Amnio Acid DelIns normalizes correctly."""
    resp = await test_handler.normalize("NP_001333827.1:p.Leu747_Thr751delinsPro")
    assertion_checks(resp.variation_descriptor, protein_delins,
                     "NP_001333827.1:p.Leu747_Thr751delinsPro")

    resp = await test_handler.normalize("EGFR p.Leu747_Thr751delinsPro")
    assert resp.variation_descriptor.id == \
        "normalize.variation:EGFR%20p.Leu747_Thr751delinsPro"
    resp.variation_descriptor.id = \
        "normalize.variation:NP_001333827.1%3Ap.Leu747_Thr751delinsPro"
    assertion_checks(resp.variation_descriptor, protein_delins,
                     "EGFR p.Leu747_Thr751delinsPro")

    resp = await test_handler.normalize("EGFR Leu747_Thr751delinsPro")
    assert resp.variation_descriptor.id == \
        "normalize.variation:EGFR%20Leu747_Thr751delinsPro"
    resp.variation_descriptor.id = \
        "normalize.variation:NP_001333827.1%3Ap.Leu747_Thr751delinsPro"
    assertion_checks(resp.variation_descriptor, protein_delins,
                     "EGFR Leu747_Thr751delinsPro")

    resp = await test_handler.normalize("EGFR L747_T751delinsP")
    assert resp.variation_descriptor.id == \
        "normalize.variation:EGFR%20L747_T751delinsP"
    resp.variation_descriptor.id = \
        "normalize.variation:NP_001333827.1%3Ap.Leu747_Thr751delinsPro"
    assertion_checks(resp.variation_descriptor, protein_delins, "EGFR L747_T751delinsP")


@pytest.mark.asyncio
async def test_protein_deletion(test_handler, protein_deletion_np_range):
    """Test that Protein Deletion normalizes correctly."""
    resp = await test_handler.normalize("NP_004439.2:p.Leu755_Thr759del")
    assertion_checks(resp.variation_descriptor, protein_deletion_np_range,
                     "NP_004439.2:p.Leu755_Thr759del")

    resp = await test_handler.normalize("ERBB2 p.Leu755_Thr759del")
    assert resp.variation_descriptor.id == \
        "normalize.variation:ERBB2%20p.Leu755_Thr759del"
    resp.variation_descriptor.id = \
        "normalize.variation:NP_004439.2%3Ap.Leu755_Thr759del"
    assertion_checks(resp.variation_descriptor, protein_deletion_np_range,
                     "ERBB2 p.Leu755_Thr759del")

    resp = await test_handler.normalize("ERBB2 Leu755_Thr759del")
    assert resp.variation_descriptor.id == \
        "normalize.variation:ERBB2%20Leu755_Thr759del"
    resp.variation_descriptor.id = \
        "normalize.variation:NP_004439.2%3Ap.Leu755_Thr759del"
    assertion_checks(resp.variation_descriptor, protein_deletion_np_range,
                     "ERBB2 Leu755_Thr759del")


@pytest.mark.asyncio
async def test_coding_dna_deletion(test_handler, coding_dna_deletion):
    """Test that coding dna deletion normalizes correctly."""
    # https://reg.clinicalgenome.org/redmine/projects/registry/genboree_registry/by_caid?caid=CA645372623  # noqa: E501
    q = "NM_004448.3:c.2264_2278delTGAGGGAAAACACAT"
    resp = await test_handler.normalize(q)
    assertion_checks(resp.variation_descriptor, coding_dna_deletion, q)

    q = "ERBB2 c.2264_2278delTGAGGGAAAACACAT"
    resp = await test_handler.normalize(q)
    assert resp.variation_descriptor.id == \
           "normalize.variation:ERBB2%20c.2264_2278delTGAGGGAAAACACAT"
    resp.variation_descriptor.id = \
        "normalize.variation:NM_004448.3%3Ac.2264_2278delTGAGGGAAAACACAT"
    assertion_checks(resp.variation_descriptor, coding_dna_deletion, q)


@pytest.mark.asyncio
async def test_protein_insertion(test_handler, protein_insertion):
    """Test that protein insertion normalizes correctly."""
    resp = await test_handler.normalize("NP_005219.2:p.Asp770_Asn771insGlyLeu")
    assertion_checks(resp.variation_descriptor, protein_insertion,
                     "NP_005219.2:p.Asp770_Asn771insGlyLeu")

    def change_resp(response):
        fixture_id = \
            "normalize.variation:NP_005219.2%3Ap.Asp770_Asn771insGlyLeu"
        response.id = fixture_id

    resp = await test_handler.normalize("EGFR D770_N771insGL")
    assert resp.variation_descriptor.id == "normalize.variation:EGFR%20D770_N771insGL"
    change_resp(resp.variation_descriptor)
    assertion_checks(resp.variation_descriptor, protein_insertion,
                     "EGFR D770_N771insGL")

    resp = await test_handler.normalize("EGFR p.D770_N771insGL")
    assert resp.variation_descriptor.id == "normalize.variation:EGFR%20p.D770_N771insGL"
    change_resp(resp.variation_descriptor)
    assertion_checks(resp.variation_descriptor, protein_insertion,
                     "EGFR p.D770_N771insGL")

    resp = await test_handler.normalize("EGFR Asp770_Asn771insGlyLeu")
    assert resp.variation_descriptor.id == \
        "normalize.variation:EGFR%20Asp770_Asn771insGlyLeu"
    change_resp(resp.variation_descriptor)
    assertion_checks(resp.variation_descriptor, protein_insertion,
                     "EGFR Asp770_Asn771insGlyLeu")

    resp = await test_handler.normalize("EGFR p.Asp770_Asn771insGlyLeu")
    assert resp.variation_descriptor.id == \
        "normalize.variation:EGFR%20p.Asp770_Asn771insGlyLeu"
    change_resp(resp.variation_descriptor)
    assertion_checks(resp.variation_descriptor, protein_insertion,
                     "EGFR p.Asp770_Asn771insGlyLeu")


@pytest.mark.asyncio
async def test_coding_dna_insertion(test_handler, coding_dna_insertion):
    """Test that coding dna insertion normalizes correctly."""
    resp = await test_handler.normalize("ENST00000331728.9:c.2049_2050insA")
    assertion_checks(resp.variation_descriptor, coding_dna_insertion,
                     "ENST00000331728.9:c.2049_2050insA")

    # TODO: issue-136
    # resp = await test_handler.normalize("LIMK2 c.2049_2050insA")
    # assert resp.variation_descriptor.id == "normalize.variation:LIMK2%20c.2049_2050insA"  # noqa: E501
    # resp.variation_descriptor.id = "normalize.variation:ENST00000331728.9%3Ac.2049_2050insA"  # noqa: E501
    # assertion_checks(resp.variation_descriptor, coding_dna_insertion)


@pytest.mark.asyncio
async def test_genomic_insertion(test_handler, genomic_insertion,
                                 grch38_genomic_insertion):
    """Test that genomic insertion normalizes correctly."""
    resp = await test_handler.normalize("NC_000017.10:g.37880993_37880994insGCTTACGTGATG")  # noqa: E501
    assertion_checks(resp.variation_descriptor, grch38_genomic_insertion,
                     "NC_000017.10:g.37880993_37880994insGCTTACGTGATG")

    fixture_id = \
        "normalize.variation:NC_000017.10%3Ag.37880993_37880994insGCTTACGTGATG"
    resp = await test_handler.normalize("17-37880993-G-GGCTTACGTGATG")
    assert resp.variation_descriptor.id == \
        "normalize.variation:17-37880993-G-GGCTTACGTGATG"
    resp.variation_descriptor.id = fixture_id
    assertion_checks(resp.variation_descriptor, grch38_genomic_insertion,
                     "17-37880993-G-GGCTTACGTGATG")

    resp = await test_handler.normalize(
        "ERBB2 g.37880993_37880994insGCTTACGTGATG")
    assert resp.variation_descriptor.id ==\
           "normalize.variation:ERBB2%20g.37880993_37880994insGCTTACGTGATG"
    resp.variation_descriptor.id = fixture_id
    assertion_checks(resp.variation_descriptor, genomic_insertion,
                     "ERBB2 g.37880993_37880994insGCTTACGTGATG")


@pytest.mark.asyncio
async def test_valid_queries(test_handler):
    """Test that valid queries don"t throw exceptions. Used for queries that
    revealed bugs in service.
    """
    assert await test_handler.normalize("CCND1 Y44D")

    resp = await test_handler.normalize("NC_000002.12:g.73448098_73448100delCTC")
    assert resp
    assert resp.variation_descriptor.variation.state.sequence == "CTC"
    assert resp.variation_descriptor.variation.id == \
        "ga4gh:VA.7qNgfRjiDwephtGOuwtXm04YcXDSY47x"


@pytest.mark.asyncio
async def test_no_matches(test_handler):
    """Test no matches work correctly."""
    queries = [
        "braf", "braf v600000932092039e", "NP_000213.1:cp.Leu862=",
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
        resp = await test_handler.normalize(q, untranslatable_returns_text=True)
        assert resp.variation_descriptor.type == "VariationDescriptor", q
        assert resp.variation_descriptor.variation.type == "Text", q
        assert resp.variation_descriptor.label == q.strip(), q

    resp = await test_handler.normalize("clinvar:10")
    assert resp.variation_descriptor is None

    resp = await test_handler.normalize("   ")
    assert resp.variation_descriptor is None

    resp = await test_handler.normalize("")
    assert resp.variation_descriptor is None


@pytest.mark.asyncio
async def test_service_meta():
    """Test that service meta info populates correctly."""
    response = await normalize_get_response("BRAF v600e", "default")
    service_meta = response.service_meta_
    assert service_meta.name == "variation-normalizer"
    assert service_meta.version
    assert isinstance(service_meta.response_datetime, datetime)
    assert service_meta.url == "https://github.com/cancervariants/variation-normalization"  # noqa: E501

    response = await normalize_get_response("this-wont-normalize", "default")
    service_meta = response.service_meta_
    assert service_meta.name == "variation-normalizer"
    assert service_meta.version
    assert isinstance(service_meta.response_datetime, datetime)
    assert service_meta.url == "https://github.com/cancervariants/variation-normalization"  # noqa: E501

    response = await to_vrs_get_response("BRAF v600e")
    service_meta = response.service_meta_
    assert service_meta.name == "variation-normalizer"
    assert service_meta.version
    assert isinstance(service_meta.response_datetime, datetime)
    assert service_meta.url == "https://github.com/cancervariants/variation-normalization"  # noqa: E501

    response = await to_vrs_get_response("this-wont-normalize")
    service_meta = response.service_meta_
    assert service_meta.name == "variation-normalizer"
    assert service_meta.version
    assert isinstance(service_meta.response_datetime, datetime)
    assert service_meta.url == "https://github.com/cancervariants/variation-normalization"  # noqa: E501
