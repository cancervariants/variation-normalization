"""Module for testing HGVS Dup Del mode."""
import pytest
from ga4gh.vrsatile.pydantic.vrsatile_models import VariationDescriptor

from tests.conftest import assertion_checks


@pytest.fixture(scope="module")
def test_handler(test_query_handler):
    """Create test fixture for normalize handler"""
    return test_query_handler.normalize_handler


@pytest.fixture(scope="module")
def dmd_gene_context():
    """Create test fixture for DMD gene context"""
    return {
        "id": "normalize.gene:DMD",
        "type": "GeneDescriptor",
        "label": "DMD",
        "xrefs": [
            "ensembl:ENSG00000198947",
            "ncbigene:1756"
        ],
        "alternate_labels": [
            "DXS272",
            "DXS230",
            "DXS206",
            "DXS142",
            "CMD3B",
            "DXS269",
            "BMD",
            "DXS268",
            "MRX85",
            "DXS164",
            "DXS270",
            "DXS239"
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
                "value": "dystrophin"
            },
            {
                "type": "Extension",
                "name": "hgnc_locations",
                "value": [
                    {
                        "species_id": "taxonomy:9606",
                        "interval": {
                            "type": "CytobandInterval",
                            "start": "p21.2",
                            "end": "p21.1"
                        },
                        "_id": "ga4gh:VCL.JgyIOPZJ9G6Hn6QziVAs8SQpaIWPK46H",
                        "type": "ChromosomeLocation",
                        "chr": "X"
                    }
                ]
            },
            {
                "type": "Extension",
                "name": "ensembl_locations",
                "value": [
                    {
                        "_id": "ga4gh:VSL.OrVsUUl1S5V_X_TWlFv6fZNrEwc0QRIt",
                        "type": "SequenceLocation",
                        "sequence_id": "ga4gh:SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
                        "interval": {
                            "start": {"type": "Number", "value": 31097676},
                            "end": {"type": "Number", "value": 33339609},
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
                        "_id": "ga4gh:VCL.JgyIOPZJ9G6Hn6QziVAs8SQpaIWPK46H",
                        "type": "ChromosomeLocation",
                        "species_id": "taxonomy:9606",
                        "chr": "X",
                        "interval": {
                            "end": "p21.1",
                            "start": "p21.2",
                            "type": "CytobandInterval"
                        }
                    },
                    {
                        "_id": "ga4gh:VSL.xHx4VZPL2bgxgaB8f6xhsayGPVp8GFUE",
                        "type": "SequenceLocation",
                        "sequence_id": "ga4gh:SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
                        "interval": {
                            "start": {"type": "Number", "value": 31119221},
                            "end": {"type": "Number", "value": 33339388},
                            "type": "SequenceInterval"
                        }
                    }
                ]
            },
            {
                "type": "Extension",
                "name": "associated_with",
                "value": [
                    "omim:300377",
                    "ucsc:uc004dda.2",
                    "ccds:CCDS14234",
                    "ccds:CCDS55394",
                    "pubmed:3607877",
                    "ccds:CCDS14232",
                    "ccds:CCDS55395",
                    "orphanet:121117",
                    "ccds:CCDS14233",
                    "ccds:CCDS75965",
                    "ccds:CCDS48091",
                    "vega:OTTHUMG00000021336",
                    "uniprot:P11532",
                    "ccds:CCDS14231",
                    "ena.embl:AF047505",
                    "pubmed:23900271",
                    "refseq:NM_004006",
                    "pubmed:3282674"
                ]
            },
            {
                "type": "Extension",
                "name": "previous_symbols",
                "value": [
                    "MRX85"
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
        ],
        "gene_id": "hgnc:2928"
    }


@pytest.fixture(scope="module")
def mecp2_gene_context():
    """Create test fixture for MECP2 gene context"""
    return {
        "id": "normalize.gene:MECP2",
        "type": "GeneDescriptor",
        "label": "MECP2",
        "xrefs": [
            "ensembl:ENSG00000169057",
            "ncbigene:4204"
        ],
        "alternate_labels": [
            "RTT",
            "AUTSX3",
            "RS",
            "PPMX",
            "MRX16",
            "MRXS13",
            "LOC113065",
            "RTS",
            "MRX79",
            "MRXSL"
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
                "value": "methyl-CpG binding protein 2"
            },
            {
                "type": "Extension",
                "name": "hgnc_locations",
                "value": [
                    {
                        "species_id": "taxonomy:9606",
                        "interval": {
                            "type": "CytobandInterval",
                            "start": "q28",
                            "end": "q28"
                        },
                        "_id": "ga4gh:VCL.fEBeCyej0jVKsvjw4vxyW6j1h8UVLb5S",
                        "type": "ChromosomeLocation",
                        "chr": "X"
                    }
                ]
            },
            {
                "type": "Extension",
                "name": "ensembl_locations",
                "value": [
                    {
                        "_id": "ga4gh:VSL.Z0ZN6J_WJZmzl5NPLALLbHQ0CJqYKHph",
                        "type": "SequenceLocation",
                        "sequence_id": "ga4gh:SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
                        "interval": {
                            "start": {"type": "Number", "value": 154021572},
                            "end": {"type": "Number", "value": 154137103},
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
                        "_id": "ga4gh:VCL.fEBeCyej0jVKsvjw4vxyW6j1h8UVLb5S",
                        "type": "ChromosomeLocation",
                        "species_id": "taxonomy:9606",
                        "chr": "X",
                        "interval": {
                            "end": "q28",
                            "start": "q28",
                            "type": "CytobandInterval"
                        }
                    },
                    {
                        "_id": "ga4gh:VSL.Z_sajXxiXtIMqH42SNl7mDl7gLAAfZDm",
                        "type": "SequenceLocation",
                        "sequence_id": "ga4gh:SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
                        "interval": {
                            "start": {"type": "Number", "value": 154021572},
                            "end": {"type": "Number", "value": 154097717},
                            "type": "SequenceInterval"
                        }
                    }
                ]
            },
            {
                "type": "Extension",
                "name": "associated_with",
                "value": [
                    "ccds:CCDS48193",
                    "ucsc:uc004fjv.3",
                    "ccds:CCDS14741",
                    "omim:300005",
                    "orphanet:123186",
                    "pubmed:1606614",
                    "pubmed:10508514",
                    "uniprot:P51608",
                    "refseq:NM_004992",
                    "vega:OTTHUMG00000024229",
                    "ena.embl:AF158180"
                ]
            },
            {
                "type": "Extension",
                "name": "previous_symbols",
                "value": [
                    "LOC113065",
                    "PPMX",
                    "MRX16",
                    "RTT",
                    "MRX79"
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
        ],
        "gene_id": "hgnc:6990"
    }


@pytest.fixture(scope="module")
def genomic_dup1():
    """Create test fixture containing params for genomic dup VD."""
    params = {
        "id": "normalize.variation:NC_000003.12%3Ag.49531262dup",
        "type": "VariationDescriptor",
        "variation_id": "",
        "variation": dict(),
        "molecule_context": "genomic",
        "structural_type": "SO:1000035",
        "vrs_ref_allele_seq": "GG"
    }
    return params


@pytest.fixture(scope="module")
def genomic_dup1_lse(genomic_dup1, genomic_dup1_seq_loc):
    """Create a test fixture for genomic dup LSE."""
    _id = "ga4gh:VA.aeNse-a8IJzqHiG-P5zTRYA_eVFhrJXw"
    genomic_dup1["variation_id"] = _id
    genomic_dup1["variation"] = {
        "type": "Allele",
        "_id": _id,
        "location": genomic_dup1_seq_loc,
        "state": {
            "type": "LiteralSequenceExpression",
            "sequence": "GGG"
        }
    }
    return VariationDescriptor(**genomic_dup1)


@pytest.fixture(scope="module")
def genomic_dup1_vac(genomic_dup1, genomic_dup1_38_vac):
    """Create a test fixture for genomic dup absolute CNV."""
    genomic_dup1["variation"] = genomic_dup1_38_vac
    genomic_dup1["variation_id"] = genomic_dup1["variation"]["_id"]
    return VariationDescriptor(**genomic_dup1)


@pytest.fixture(scope="module")
def genomic_dup1_vrc(genomic_dup1, genomic_dup1_seq_loc):
    """Create a test fixture for genomic dup relative CNV."""
    genomic_dup1["variation"] = {
        "type": "RelativeCopyNumber",
        "_id": "ga4gh:VRC.dhn3PLeop46_QYTlYxiIW18nrTfYqJ8N",
        "subject": genomic_dup1_seq_loc,
        "relative_copy_class": "high-level gain"
    }
    genomic_dup1["variation_id"] = genomic_dup1["variation"]["_id"]
    return VariationDescriptor(**genomic_dup1)


@pytest.fixture(scope="module")
def genomic_dup1_rse(genomic_dup1, genomic_dup1_seq_loc):
    """Create a test fixture for genomic dup RSE."""
    _id = "ga4gh:VA.lAyulP9JxvQReKWjpq0LbO50r2UTeMkl"
    genomic_dup1["variation"] = {
        "type": "Allele",
        "_id": _id,
        "location": genomic_dup1_seq_loc,
        "state": {
            "type": "RepeatedSequenceExpression",
            "seq_expr": {
                "type": "DerivedSequenceExpression",
                "location": genomic_dup1_seq_loc,
                "reverse_complement": False
            },
            "count": {
                "type": "Number",
                "value": 2
            }
        }
    }
    genomic_dup1["variation_id"] = _id
    return VariationDescriptor(**genomic_dup1)


@pytest.fixture(scope="module")
def genomic_dup1_free_text():
    """Create test fixture containing params for genomic dup VD."""
    params = {
        "id": "normalize.variation:DAG1%20g.49568695dup",
        "type": "VariationDescriptor",
        "variation_id": "",
        "variation": dict(),
        "molecule_context": "transcript",
        "structural_type": "SO:1000035",
        "vrs_ref_allele_seq": "GG",
        "gene_context": {
            "id": "normalize.gene:DAG1",
            "type": "GeneDescriptor",
            "label": "DAG1",
            "xrefs": [
                "ensembl:ENSG00000173402",
                "ncbigene:1605"
            ],
            "alternate_labels": [
                "156DAG",
                "MDDGA9",
                "AGRNR",
                "DAG",
                "LGMDR16",
                "A3a",
                "MDDGC7",
                "MDDGC9"
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
                    "value": "dystroglycan 1"
                },
                {
                    "type": "Extension",
                    "name": "hgnc_locations",
                    "value": [
                        {
                            "species_id": "taxonomy:9606",
                            "interval": {
                                "type": "CytobandInterval",
                                "start": "p21.31",
                                "end": "p21.31"
                            },
                            "_id": "ga4gh:VCL.l_F_O8hoRfwdUsaN3UScymcvqRWLeKQT",
                            "type": "ChromosomeLocation",
                            "chr": "3"
                        }
                    ]
                },
                {
                    "type": "Extension",
                    "name": "ensembl_locations",
                    "value": [
                        {
                            "_id": "ga4gh:VSL.p7nu0_gNxavS-GyUYDdGoZ00Yo72w70A",
                            "type": "SequenceLocation",
                            "sequence_id": "ga4gh:SQ.Zu7h9AggXxhTaGVsy7h_EZSChSZGcmgX",
                            "interval": {
                                "start": {"type": "Number", "value": 49468712},
                                "end": {"type": "Number", "value": 49535618},
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
                            "_id": "ga4gh:VCL.l_F_O8hoRfwdUsaN3UScymcvqRWLeKQT",
                            "type": "ChromosomeLocation",
                            "species_id": "taxonomy:9606",
                            "chr": "3",
                            "interval": {
                                "end": "p21.31",
                                "start": "p21.31",
                                "type": "CytobandInterval"
                            }
                        },
                        {
                            "_id": "ga4gh:VSL.py1JOK16I3jc9PuXWvh9Wmuj-QMIKNBg",
                            "type": "SequenceLocation",
                            "sequence_id": "ga4gh:SQ.Zu7h9AggXxhTaGVsy7h_EZSChSZGcmgX",
                            "interval": {
                                "start": {"type": "Number", "value": 49468947},
                                "end": {"type": "Number", "value": 49535615},
                                "type": "SequenceInterval"
                            }
                        }
                    ]
                },
                {
                    "type": "Extension",
                    "name": "associated_with",
                    "value": [
                        "pubmed:7774920",
                        "pubmed:1741056",
                        "merops:S72.001",
                        "uniprot:Q14118",
                        "orphanet:280347",
                        "ucsc:uc021wxz.1",
                        "vega:OTTHUMG00000156869",
                        "omim:128239",
                        "refseq:NM_001165928",
                        "ena.embl:L19711",
                        "ccds:CCDS2799"
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
            ],
            "gene_id": "hgnc:2666"
        }
    }
    return params


@pytest.fixture(scope="module")
def genomic_dup1_free_text_seq_loc():
    """Create genomic dup1 free text sequence location"""
    return {
        "_id": "ga4gh:VSL.wasOdqigAN-is7O2nEqJeDwkPlwpiOak",
        "sequence_id": "ga4gh:SQ.tpvbnWsfEGqip8gJQZnWJAF8-bWDUDKd",
        "interval": {
            "type": "SequenceInterval",
            "start": {"value": 1032, "type": "Number"},
            "end": {"value": 1034, "type": "Number"},
        },
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_dup1_free_text_lse(genomic_dup1_free_text,
                               genomic_dup1_free_text_seq_loc):
    """Create a test fixture for genomic dup LSE."""
    _id = "ga4gh:VA.eE5Kr1zJrv3PSXeBabbKTFnZxToaYxat"
    genomic_dup1_free_text["variation_id"] = _id
    genomic_dup1_free_text["variation"] = {
        "type": "Allele",
        "_id": _id,
        "location": genomic_dup1_free_text_seq_loc,
        "state": {
            "type": "LiteralSequenceExpression",
            "sequence": "GGG"
        }
    }
    return VariationDescriptor(**genomic_dup1_free_text)


@pytest.fixture(scope="module")
def genomic_dup1_free_text_vac(genomic_dup1_free_text,
                               genomic_dup1_free_text_seq_loc):
    """Create a test fixture for genomic dup absolute CNV."""
    _id = "ga4gh:VAC.qeuDGWVaGZUOf7XmF2xO1k24LvFzGVE1"
    genomic_dup1_free_text["variation"] = {
        "type": "AbsoluteCopyNumber",
        "_id": _id,
        "subject": genomic_dup1_free_text_seq_loc,
        "copies": {"type": "Number", "value": 3}
    }
    genomic_dup1_free_text["variation_id"] = _id
    return VariationDescriptor(**genomic_dup1_free_text)


@pytest.fixture(scope="module")
def genomic_dup1_free_text_rse(genomic_dup1_free_text,
                               genomic_dup1_free_text_seq_loc):
    """Create a test fixture for genomic dup RSE."""
    _id = "ga4gh:VA.VQKwP3GpeObfGc3MzvA9JNL1YwkZynKO"
    genomic_dup1_free_text["variation"] = {
        "type": "Allele",
        "_id": _id,
        "location": genomic_dup1_free_text_seq_loc,
        "state": {
            "type": "RepeatedSequenceExpression",
            "seq_expr": {
                "type": "DerivedSequenceExpression",
                "location": genomic_dup1_free_text_seq_loc,
                "reverse_complement": False
            },
            "count": {
                "type": "Number",
                "value": 2
            }
        }
    }
    genomic_dup1_free_text["variation_id"] = _id
    return VariationDescriptor(**genomic_dup1_free_text)


@pytest.fixture(scope="module")
def genomic_dup2():
    """Create test fixture containing params for genomic dup VD."""
    params = {
        "id": "normalize.variation:NC_000016.10%3Ag.2087938_2087948dup",
        "type": "VariationDescriptor",
        "variation_id": "",
        "variation": dict(),
        "molecule_context": "genomic",
        "structural_type": "SO:1000035",
        "vrs_ref_allele_seq": "AAAGGTAGGGC"
    }
    return params


@pytest.fixture(scope="module")
def genomic_dup2_lse(genomic_dup2, genomic_dup2_seq_loc):
    """Create a test fixture for genomic dup LSE."""
    _id = "ga4gh:VA.wqqxfUCrFSndedI2-4oiIuHLHHGjBFof"
    genomic_dup2["variation_id"] = _id
    genomic_dup2["variation"] = {
        "type": "Allele",
        "_id": _id,
        "location": genomic_dup2_seq_loc,
        "state": {
            "type": "LiteralSequenceExpression",
            "sequence": "AAAGGTAGGGCAAAGGTAGGGC"
        }
    }
    return VariationDescriptor(**genomic_dup2)


@pytest.fixture(scope="module")
def genomic_dup2_vac(genomic_dup2, genomic_dup2_38_vac):
    """Create a test fixture for genomic dup absolute CNV."""
    genomic_dup2["variation"] = genomic_dup2_38_vac
    genomic_dup2["variation_id"] = genomic_dup2["variation"]["_id"]
    return VariationDescriptor(**genomic_dup2)


@pytest.fixture(scope="module")
def genomic_dup2_vrc(genomic_dup2, genomic_dup2_seq_loc):
    """Create a test fixture for genomic dup relative CNV."""
    genomic_dup2["variation"] = {
        "type": "RelativeCopyNumber",
        "_id": "ga4gh:VRC.xuC326dmLBd33E4d-z9JItz0mt24JR-I",
        "subject": genomic_dup2_seq_loc,
        "relative_copy_class": "low-level gain"
    }
    genomic_dup2["variation_id"] = genomic_dup2["variation"]["_id"]
    return VariationDescriptor(**genomic_dup2)


@pytest.fixture(scope="module")
def genomic_dup2_rse(genomic_dup2, genomic_dup2_seq_loc):
    """Create a test fixture for genomic dup RSE."""
    _id = "ga4gh:VA.fXANtjCcUPJ1A4bCSgcAxSSrxoqXuL3A"
    genomic_dup2["variation"] = {
        "type": "Allele",
        "_id": _id,
        "location": genomic_dup2_seq_loc,
        "state": {
            "type": "RepeatedSequenceExpression",
            "seq_expr": {
                "type": "DerivedSequenceExpression",
                "location": genomic_dup2_seq_loc,
                "reverse_complement": False
            },
            "count": {
                "type": "Number",
                "value": 2
            }
        }
    }
    genomic_dup2["variation_id"] = _id
    return VariationDescriptor(**genomic_dup2)


@pytest.fixture(scope="module")
def genomic_dup2_free_text(dmd_gene_context):
    """Create test fixture containing params for genomic dup VD."""
    params = {
        "id": "normalize.variation:TSC2%20g.2137939_2137949dup",
        "type": "VariationDescriptor",
        "variation_id": "",
        "variation": dict(),
        "molecule_context": "transcript",
        "structural_type": "SO:1000035",
        "vrs_ref_allele_seq": "TAGA",
        "gene_context": dmd_gene_context
    }
    return params


@pytest.fixture(scope="module")
def genomic_dup2_free_text_seq_loc():
    """Create genomic dup2 free text sequence location"""
    return {
        "_id": "ga4gh:VSL.3JAa1wqyQWE510wqzNXoPptxYVXocFqj",
        "sequence_id": "ga4gh:SQ.1DeZLYHMnd-smp3GDlpRxETb9_0AokO7",
        "interval": {
            "type": "SequenceInterval",
            "start": {"value": 256, "type": "Number"},
            "end": {"value": 260, "type": "Number"},
        },
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_dup2_free_text_default(genomic_dup2_free_text,
                                   genomic_dup2_free_text_seq_loc):
    """Create a test fixture for genomic dup default and LSE."""
    _id = "ga4gh:VA.BRi89LZSxVMXaa6xVjuXIh0I_u2MyPkc"
    genomic_dup2_free_text["variation_id"] = _id
    genomic_dup2_free_text["variation"] = {
        "type": "Allele",
        "_id": _id,
        "location": genomic_dup2_free_text_seq_loc,
        "state": {
            "type": "LiteralSequenceExpression",
            "sequence": "TAGATAGA"
        }
    }
    return VariationDescriptor(**genomic_dup2_free_text)


@pytest.fixture(scope="module")
def genomic_dup2_free_text_vac(genomic_dup2_free_text,
                               genomic_dup2_free_text_seq_loc):
    """Create a test fixture for genomic dup absolute CNV."""
    _id = "ga4gh:VAC.p538XMNKHyZGEVb73xbA2DfSxSJOZG4B"
    genomic_dup2_free_text["variation"] = {
        "type": "AbsoluteCopyNumber",
        "_id": _id,
        "subject": genomic_dup2_free_text_seq_loc,
        "copies": {"type": "Number", "value": 3}
    }
    genomic_dup2_free_text["variation_id"] = _id
    return VariationDescriptor(**genomic_dup2_free_text)


@pytest.fixture(scope="module")
def genomic_dup2_free_text_rse(genomic_dup2_free_text,
                               genomic_dup2_free_text_seq_loc):
    """Create a test fixture for genomic dup RSE."""
    _id = "ga4gh:VA.Rby7K6TikhqXL9BhM8xDJHNudJlRmS3j"
    genomic_dup2_free_text["variation"] = {
        "type": "Allele",
        "_id": _id,
        "location": genomic_dup2_free_text_seq_loc,
        "state": {
            "type": "RepeatedSequenceExpression",
            "seq_expr": {
                "type": "DerivedSequenceExpression",
                "location": genomic_dup2_free_text_seq_loc,
                "reverse_complement": False
            },
            "count": {
                "type": "Number",
                "value": 2
            }
        }
    }
    genomic_dup2_free_text["variation_id"] = _id
    return VariationDescriptor(**genomic_dup2_free_text)


@pytest.fixture(scope="module")
def genomic_dup3():
    """Create test fixture containing params for genomic dup VD."""
    params = {
        "id": "normalize.variation:NC_000023.11%3Ag.%2831060227_31100351%29_%2833274278_33417151%29dup",  # noqa: E501
        "type": "VariationDescriptor",
        "variation_id": "",
        "variation": dict(),
        "molecule_context": "genomic",
        "structural_type": "SO:0001742",
        "vrs_ref_allele_seq": None
    }
    return params


@pytest.fixture(scope="module")
def genomic_dup3_vac(genomic_dup3, genomic_del3_dup3_loc):
    """Create a test fixture for genomic dup absolute cnv."""
    _id = "ga4gh:VAC.cQATJ6a1uGwXOHu-advv8lRsMgjNLKul"
    genomic_dup3["variation"] = {
        "type": "AbsoluteCopyNumber",
        "_id": _id,
        "subject": genomic_del3_dup3_loc,
        "copies": {"type": "Number", "value": 2}
    }
    genomic_dup3["variation_id"] = _id
    return VariationDescriptor(**genomic_dup3)


@pytest.fixture(scope="module")
def genomic_dup3_vrc(genomic_dup3, genomic_del3_dup3_loc):
    """Create a test fixture for genomic dup relative cnv."""
    _id = "ga4gh:VRC.APoRHUv3_v5FSbtfU5DuESzF9iuDLdYx"
    genomic_dup3["variation"] = {
        "type": "RelativeCopyNumber",
        "_id": _id,
        "subject": genomic_del3_dup3_loc,
        "relative_copy_class": "low-level gain"
    }
    genomic_dup3["variation_id"] = _id
    return VariationDescriptor(**genomic_dup3)


@pytest.fixture(scope="module")
def genomic_dup3_rse_lse(genomic_dup3):
    """Create test fixture for genomic dup rse and lse."""
    params = {
        "id": genomic_dup3["id"],
        "type": genomic_dup3["type"],
        "variation": {
            "_id": "ga4gh:VT.15sKDgSyoCPOgfrFHvSea-fHVeu7huVT",
            "type": "Text",
            "definition": "NC_000023.11:g.(31060227_31100351)_(33274278_33417151)dup"  # noqa: E501
        }
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_dup3_free_text(dmd_gene_context):
    """Create test fixture containing params for genomic dup VD."""
    params = {
        "id": "normalize.variation:DMD%20g.%2831147274_31147278%29_%2831182737_31182739%29dup",  # noqa: E501
        "type": "VariationDescriptor",
        "variation_id": "",
        "variation": dict(),
        "molecule_context": "genomic",
        "structural_type": "SO:0001742",
        "vrs_ref_allele_seq": None,
        "gene_context": dmd_gene_context
    }
    return params


@pytest.fixture(scope="module")
def genomic_dup3_free_text_subject():
    """Create test fixture for genomic dup3 free text subject"""
    return {
        "_id": "ga4gh:VSL.6JRgXRroqGleDLuwmOdHSbUK8Lm27fos",
        "sequence_id": "ga4gh:SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
        "interval": {
            "type": "SequenceInterval",
            "start": {
                "min": 31147273,
                "max": 31147277,
                "type": "DefiniteRange"
            },
            "end": {
                "min": 31182738,
                "max": 31182740,
                "type": "DefiniteRange"
            }
        },
        "type": "SequenceLocation"
    }


@pytest.fixture(scope="module")
def genomic_dup3_free_text_vrc(genomic_dup3_free_text, genomic_dup3_free_text_subject):
    """Create a test fixture for genomic dup relative cnv."""
    _id = "ga4gh:VRC.0-oCLK0W3ND8-EzZsRJUJlF5MTe3VgSD"
    genomic_dup3_free_text["variation"] = {
        "type": "RelativeCopyNumber",
        "_id": _id,
        "subject": genomic_dup3_free_text_subject,
        "relative_copy_class": "low-level gain"
    }
    genomic_dup3_free_text["variation_id"] = _id
    return VariationDescriptor(**genomic_dup3_free_text)


@pytest.fixture(scope="module")
def genomic_dup3_free_text_vac(genomic_dup3_free_text, genomic_dup3_free_text_subject):
    """Create a test fixture for genomic dup absolute cnv."""
    _id = "ga4gh:VAC.WwO4IBA5qODo__32hAMeV8cb1q5uZyqd"
    genomic_dup3_free_text["variation"] = {
        "type": "AbsoluteCopyNumber",
        "_id": _id,
        "subject": genomic_dup3_free_text_subject,
        "copies": {"type": "Number", "value": 4}
    }
    genomic_dup3_free_text["variation_id"] = _id
    return VariationDescriptor(**genomic_dup3_free_text)


@pytest.fixture(scope="module")
def genomic_dup3_free_text_rse_lse(genomic_dup3_free_text):
    """Create test fixture for genomic dup rse and lse."""
    params = {
        "id": genomic_dup3_free_text["id"],
        "type": genomic_dup3_free_text["type"],
        "variation": {
            "_id": "ga4gh:VT.F0AX-RkMN4U8KLkIE68ECGU83Y-ICWXh",
            "type": "Text",
            "definition": "DMD g.(31147274_31147278)_(31182737_31182739)dup"
        }
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_dup4():
    """Create test fixture containing params for genomic dup VD."""
    params = {
        "id": "normalize.variation:NC_000020.11%3Ag.%28%3F_30417576%29_%2831394018_%3F%29dup",  # noqa: E501
        "type": "VariationDescriptor",
        "variation_id": "",
        "variation": dict(),
        "molecule_context": "genomic",
        "structural_type": "SO:0001742",
        "vrs_ref_allele_seq": None
    }
    return params


@pytest.fixture(scope="module")
def genomic_dup4_vrc(genomic_dup4, genoimc_dup4_loc):
    """Create a test fixture for genomic dup relative cnv."""
    _id = "ga4gh:VRC.9N_ih-YzkKM91k6adSGQuGce5t4nJspV"
    genomic_dup4["variation"] = {
        "type": "RelativeCopyNumber",
        "_id": _id,
        "subject": genoimc_dup4_loc,
        "relative_copy_class": "low-level gain"
    }
    genomic_dup4["variation_id"] = _id
    return VariationDescriptor(**genomic_dup4)


@pytest.fixture(scope="module")
def genomic_dup4_vac(genomic_dup4, genoimc_dup4_loc):
    """Create a test fixture for genomic dup absolute cnv."""
    _id = "ga4gh:VAC.h594XLS8a4VA6j-ghLaghqXmof8hmF5z"
    genomic_dup4["variation"] = {
        "type": "AbsoluteCopyNumber",
        "_id": _id,
        "subject": genoimc_dup4_loc,
        "copies": {"type": "Number", "value": 3}
    }
    genomic_dup4["variation_id"] = _id
    return VariationDescriptor(**genomic_dup4)


@pytest.fixture(scope="module")
def genomic_dup4_rse_lse(genomic_dup4):
    """Create test fixture for genomic dup rse and lse."""
    params = {
        "id": genomic_dup4["id"],
        "type": genomic_dup4["type"],
        "variation": {
            "_id": "ga4gh:VT.Pga4IH82qga2iZAodjxYw9OXhB4Xa2g8",
            "type": "Text",
            "definition": "NC_000020.11:g.(?_30417576)_(31394018_?)dup"
        }
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_dup4_free_text():
    """Create test fixture containing params for genomic dup VD."""
    params = {
        "id": "normalize.variation:PRF8%20g.%28%3F_1577736%29_%281587865_%3F%29",  # noqa: E501
        "type": "VariationDescriptor",
        "variation_id": "",
        "variation": dict(),
        "molecule_context": "genomic",
        "structural_type": "SO:0001742",
        "vrs_ref_allele_seq": None,
        "gene_context": {
            "id": "normalize.gene:PRPF8",
            "type": "GeneDescriptor",
            "label": "PRPF8",
            "xrefs": [
                "ensembl:ENSG00000174231",
                "ncbigene:10594"
            ],
            "alternate_labels": [
                "PRPC8",
                "PRP8",
                "HPRP8",
                "Prp8",
                "RP13",
                "hPrp8",
                "SNRNP220"
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
                    "value": "pre-mRNA processing factor 8"
                },
                {
                    "type": "Extension",
                    "name": "hgnc_locations",
                    "value": [
                        {
                            "species_id": "taxonomy:9606",
                            "interval": {
                                "type": "CytobandInterval",
                                "start": "p13.3",
                                "end": "p13.3"
                            },
                            "_id": "ga4gh:VCL.GJ_KKaBnwZCC9_0vezbSxp_yAwM6R8c4",
                            "type": "ChromosomeLocation",
                            "chr": "17"
                        }
                    ]
                },
                {
                    "type": "Extension",
                    "name": "ensembl_locations",
                    "value": [
                        {
                            "_id": "ga4gh:VSL.REtW8dfZCgDLEvo58qhp-dkN-hHiRtDx",
                            "type": "SequenceLocation",
                            "sequence_id": "ga4gh:SQ.dLZ15tNO1Ur0IcGjwc3Sdi_0A6Yf4zm7",
                            "interval": {
                                "start": {"type": "Number", "value": 1650628},
                                "end": {"type": "Number", "value": 1684867},
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
                            "_id": "ga4gh:VCL.GJ_KKaBnwZCC9_0vezbSxp_yAwM6R8c4",
                            "type": "ChromosomeLocation",
                            "species_id": "taxonomy:9606",
                            "chr": "17",
                            "interval": {
                                "end": "p13.3",
                                "start": "p13.3",
                                "type": "CytobandInterval"
                            }
                        },
                        {
                            "_id": "ga4gh:VSL.REtW8dfZCgDLEvo58qhp-dkN-hHiRtDx",
                            "type": "SequenceLocation",
                            "sequence_id": "ga4gh:SQ.dLZ15tNO1Ur0IcGjwc3Sdi_0A6Yf4zm7",
                            "interval": {
                                "start": {"type": "Number", "value": 1650628},
                                "end": {"type": "Number", "value": 1684867},
                                "type": "SequenceInterval"
                            }
                        },
                        {
                            "_id": "ga4gh:VSL.5FvYcab11zKZuo57GyafVqW9IykgsjAh",
                            "type": "SequenceLocation",
                            "sequence_id": "ga4gh:SQ.3Wx-9rRd5d7m3WxtJ_HScX3Bz1MiQWjR",
                            "interval": {
                                "start": {"type": "Number", "value": 80656},
                                "end": {"type": "Number", "value": 114895},
                                "type": "SequenceInterval"
                            }
                        }
                    ]
                },
                {
                    "type": "Extension",
                    "name": "associated_with",
                    "value": [
                        "pubmed:10411133",
                        "ucsc:uc002fte.3",
                        "pubmed:11468273",
                        "orphanet:118066",
                        "ccds:CCDS11010",
                        "refseq:NM_006445",
                        "vega:OTTHUMG00000090553",
                        "uniprot:Q6P2Q9",
                        "ena.embl:AB007510",
                        "omim:607300"
                    ]
                },
                {
                    "type": "Extension",
                    "name": "previous_symbols",
                    "value": [
                        "RP13"
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
            ],
            "gene_id": "hgnc:17340"
        }
    }
    return params


@pytest.fixture(scope="module")
def genomic_dup4_free_text_subject():
    """Create test fixture for genomic dup4 free text subject"""
    return {
        "_id": "ga4gh:VSL.4eNCJnROfnoO-YvGnf-iGCeDHF_68g8H",
        "sequence_id": "ga4gh:SQ.dLZ15tNO1Ur0IcGjwc3Sdi_0A6Yf4zm7",
        "interval": {
            "type": "SequenceInterval",
            "start": {
                "value": 1674441,
                "comparator": "<=",
                "type": "IndefiniteRange"
            },
            "end": {
                "value": 1684571,
                "comparator": ">=",
                "type": "IndefiniteRange"
            }
        },
        "type": "SequenceLocation"
    }


@pytest.fixture(scope="module")
def genomic_dup4_free_text_vrc(genomic_dup4_free_text, genomic_dup4_free_text_subject):
    """Create a test fixture for genomic dup relative cnv."""
    _id = "ga4gh:VRC.JmhvBvSI2l_MPRhokZcjO8EPlqvU5V_g"
    genomic_dup4_free_text["variation"] = {
        "type": "RelativeCopyNumber",
        "_id": _id,
        "subject": genomic_dup4_free_text_subject,
        "relative_copy_class": "low-level gain"
    }
    genomic_dup4_free_text["variation_id"] = _id
    return VariationDescriptor(**genomic_dup4_free_text)


@pytest.fixture(scope="module")
def genomic_dup4_free_text_vac(genomic_dup4_free_text, genomic_dup4_free_text_subject):
    """Create a test fixture for genomic dup absolute cnv."""
    _id = "ga4gh:VAC.G0SbODtFctpFS9aufQGvRuGNRTR_D37E"
    genomic_dup4_free_text["variation"] = {
        "type": "AbsoluteCopyNumber",
        "_id": _id,
        "subject": genomic_dup4_free_text_subject,
        "copies": {
            "type": "Number",
            "value": 3
        }
    }
    genomic_dup4_free_text["variation_id"] = _id
    return VariationDescriptor(**genomic_dup4_free_text)


@pytest.fixture(scope="module")
def genomic_dup4_free_text_rse_lse(genomic_dup4_free_text):
    """Create test fixture for genomic dup rse and lse."""
    params = {
        "id": genomic_dup4_free_text["id"],
        "type": genomic_dup4_free_text["type"],
        "variation": {
            "_id": "ga4gh:VT.Pga4IH82qga2iZAodjxYw9OXhB4Xa2g8",
            "type": "Text",
            "definition": "NC_000020.11:g.(?_30417576)_(31394018_?)dup"
        }
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_dup5():
    """Create test fixture containing params for genomic dup VD."""
    params = {
        "id": "normalize.variation:NC_000023.11%3Ag.%28%3F_154021812%29_154092209dup",  # noqa: E501
        "type": "VariationDescriptor",
        "variation_id": "",
        "variation": dict(),
        "molecule_context": "genomic",
        "structural_type": "SO:0001742",
        "vrs_ref_allele_seq": None
    }
    return params


def genomic_dup5_abs_cnv(params, genomic_dup5_loc):
    """Create genomic dup5 aboluste cnv"""
    _id = "ga4gh:VAC.bAqgh0Ecs_NxhtrBUhhE8_GN6QCCN9Td"
    params["variation"] = {
        "type": "AbsoluteCopyNumber",
        "_id": _id,
        "subject": genomic_dup5_loc,
        "copies": {"type": "Number", "value": 3}
    }
    params["variation_id"] = _id


def genomic_dup5_rel_cnv(params, genomic_dup5_loc):
    """Create genomic dup4 relative cnv"""
    _id = "ga4gh:VRC.vWvFd5fM7Xa-unq_IO8OZUhYfHfOtSQ8"
    params["variation"] = {
        "type": "RelativeCopyNumber",
        "_id": _id,
        "subject": genomic_dup5_loc,
        "relative_copy_class": "low-level gain"
    }
    params["variation_id"] = _id


@pytest.fixture(scope="module")
def genomic_dup5_vrc(genomic_dup5, genomic_dup5_loc):
    """Create a test fixture for genomic dup5 relative cnv."""
    genomic_dup5_rel_cnv(genomic_dup5, genomic_dup5_loc)
    return VariationDescriptor(**genomic_dup5)


@pytest.fixture(scope="module")
def genomic_dup5_vac(genomic_dup5, genomic_dup5_loc):
    """Create a test fixture for genomic dup5 absolute cnv."""
    genomic_dup5_abs_cnv(genomic_dup5, genomic_dup5_loc)
    return VariationDescriptor(**genomic_dup5)


@pytest.fixture(scope="module")
def genomic_dup5_rse_lse(genomic_dup5):
    """Create test fixture for genomic dup rse and lse."""
    params = {
        "id": genomic_dup5["id"],
        "type": genomic_dup5["type"],
        "variation": {
            "_id": "ga4gh:VT.of16BEeHyU9od62SrjSCQ4LyUtbbGoKi",
            "type": "Text",
            "definition": "NC_000023.11:g.(?_154021812)_154092209dup"
        }
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_dup5_free_text(mecp2_gene_context):
    """Create test fixture containing params for genomic dup VD."""
    params = {
        "id": "normalize.variation:MECP2%20g.%28%3F_154021812%29_154092209dup",
        "type": "VariationDescriptor",
        "variation_id": "",
        "variation": dict(),
        "molecule_context": "genomic",
        "structural_type": "SO:0001742",
        "vrs_ref_allele_seq": None,
        "gene_context": mecp2_gene_context
    }
    return params


@pytest.fixture(scope="module")
def genomic_dup5_free_text_vrc(genomic_dup5_free_text, genomic_dup5_loc):
    """Create a test fixture for genomic dup relative cnv."""
    genomic_dup5_rel_cnv(genomic_dup5_free_text, genomic_dup5_loc)
    return VariationDescriptor(**genomic_dup5_free_text)


@pytest.fixture(scope="module")
def genomic_dup5_free_text_vac(genomic_dup5_free_text, genomic_dup5_loc):
    """Create a test fixture for genomic dup absolute cnv."""
    genomic_dup5_abs_cnv(genomic_dup5_free_text, genomic_dup5_loc)
    return VariationDescriptor(**genomic_dup5_free_text)


@pytest.fixture(scope="module")
def genomic_dup5_free_text_rse_lse(genomic_dup5_free_text):
    """Create test fixture for genomic dup rse and lse."""
    params = {
        "id": genomic_dup5_free_text["id"],
        "type": genomic_dup5_free_text["type"],
        "variation": {
            "_id": "ga4gh:VT.Kw18bSFpQp9xdKg88fqW7zUx4_VXFIiW",
            "type": "Text",
            "definition": "MECP2 g.(?_154021812)_154092209dup"
        }
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_dup6():
    """Create test fixture containing params for genomic dup VD."""
    params = {
        "id": "normalize.variation:NC_000023.11%3Ag.154021812_%28154092209_%3F%29dup",  # noqa: E501
        "type": "VariationDescriptor",
        "variation_id": "",
        "variation": dict(),
        "molecule_context": "genomic",
        "structural_type": "SO:0001742",
        "vrs_ref_allele_seq": None
    }
    return params


def genomic_dup6_rel_cnv(params, genoimc_dup6_loc):
    """Create genomic dup6 relative cnv"""
    _id = "ga4gh:VRC.YPZLS_Cld8CEkEveTafSmWfH_QyB3okT"
    params["variation"] = {
        "type": "RelativeCopyNumber",
        "_id": _id,
        "subject": genoimc_dup6_loc,
        "relative_copy_class": "low-level gain"
    }
    params["variation_id"] = _id


def genomic_dup6_abs_cnv(params, genoimc_dup6_loc):
    """Create genomic dup6 absolute cnv"""
    _id = "ga4gh:VAC.ZkgR6TD7VypzVrLAYFnSb-D7DXp62Yfn"
    params["variation"] = {
        "type": "AbsoluteCopyNumber",
        "_id": _id,
        "subject": genoimc_dup6_loc,
        "copies": {"type": "Number", "value": 2}
    }
    params["variation_id"] = _id


@pytest.fixture(scope="module")
def genomic_dup6_vrc(genomic_dup6, genoimc_dup6_loc):
    """Create a test fixture for genomic dup relative cnv."""
    genomic_dup6_rel_cnv(genomic_dup6, genoimc_dup6_loc)
    return VariationDescriptor(**genomic_dup6)


@pytest.fixture(scope="module")
def genomic_dup6_vac(genomic_dup6, genoimc_dup6_loc):
    """Create a test fixture for genomic dup absolute cnv."""
    genomic_dup6_abs_cnv(genomic_dup6, genoimc_dup6_loc)
    return VariationDescriptor(**genomic_dup6)


@pytest.fixture(scope="module")
def genomic_dup6_rse_lse(genomic_dup6):
    """Create test fixture for genomic dup rse and lse."""
    params = {
        "id": genomic_dup6["id"],
        "type": genomic_dup6["type"],
        "variation": {
            "_id": "ga4gh:VT.2k5AWTbGJxvLVT6bUW0pUMq6XGAcEjXW",
            "type": "Text",
            "definition": "NC_000023.11:g.154021812_(154092209_?)dup"
        }
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_dup6_free_text(mecp2_gene_context):
    """Create test fixture containing params for genomic dup VD."""
    params = {
        "id": "normalize.variation:MECP2%20g.154021812_%28154092209_%3F%29dup",
        "type": "VariationDescriptor",
        "variation_id": "",
        "variation": dict(),
        "molecule_context": "genomic",
        "structural_type": "SO:0001742",
        "vrs_ref_allele_seq": None,
        "gene_context": mecp2_gene_context
    }
    return params


@pytest.fixture(scope="module")
def genomic_dup6_free_text_vrc(genomic_dup6_free_text, genoimc_dup6_loc):
    """Create a test fixture for genomic dup relative cnv."""
    genomic_dup6_rel_cnv(genomic_dup6_free_text, genoimc_dup6_loc)
    return VariationDescriptor(**genomic_dup6_free_text)


@pytest.fixture(scope="module")
def genomic_dup6_free_text_vac(genomic_dup6_free_text, genoimc_dup6_loc):
    """Create a test fixture for genomic dup absolute cnv."""
    genomic_dup6_abs_cnv(genomic_dup6_free_text, genoimc_dup6_loc)
    return VariationDescriptor(**genomic_dup6_free_text)


@pytest.fixture(scope="module")
def genomic_dup6_free_text_rse_lse(genomic_dup6_free_text):
    """Create test fixture for genomic dup rse and lse."""
    params = {
        "id": genomic_dup6_free_text["id"],
        "type": genomic_dup6_free_text["type"],
        "variation": {
            "_id": "ga4gh:VT.LbAqiLmJs1t9-FgEKD0-KDKwzvM3AAlz",
            "type": "Text",
            "definition": "MECP2 g.154021812_(154092209_?)dup"
        }
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_del1():
    """Create test fixture containing params for genomic del VD."""
    params = {
        "id": "normalize.variation:NC_000003.12%3Ag.10149811del",
        "type": "VariationDescriptor",
        "variation_id": "",
        "variation": dict(),
        "molecule_context": "genomic",
        "structural_type": "SO:0000159",
        "vrs_ref_allele_seq": "T"
    }
    return params


@pytest.fixture(scope="module")
def genomic_del1_lse(genomic_del1, genomic_del1_seq_loc):
    """Create a test fixture for genomic del LSE."""
    _id = "ga4gh:VA.jUeT1n4AuBzwtt5TT-Iaac1KasATWjKE"
    genomic_del1["variation_id"] = _id
    genomic_del1["variation"] = {
        "type": "Allele",
        "_id": _id,
        "location": genomic_del1_seq_loc,
        "state": {
            "type": "LiteralSequenceExpression",
            "sequence": ""
        }
    }
    return VariationDescriptor(**genomic_del1)


@pytest.fixture(scope="module")
def genomic_del1_vac(genomic_del1, genomic_del1_38_vac):
    """Create a test fixture for genomic del absolute CNV."""
    genomic_del1["variation"] = genomic_del1_38_vac
    genomic_del1["variation_id"] = genomic_del1["variation"]["_id"]
    return VariationDescriptor(**genomic_del1)


@pytest.fixture(scope="module")
def genomic_del1_vrc(genomic_del1, genomic_del1_seq_loc):
    """Create a test fixture for genomic del relative CNV."""
    genomic_del1["variation"] = {
        "type": "RelativeCopyNumber",
        "_id": "ga4gh:VRC.FHImCphKfSBeobf9HO6qpu_Bm5U9VfHz",
        "subject": genomic_del1_seq_loc,
        "relative_copy_class": "copy neutral"
    }
    genomic_del1["variation_id"] = genomic_del1["variation"]["_id"]
    return VariationDescriptor(**genomic_del1)


@pytest.fixture(scope="module")
def genomic_del1_rse(genomic_del1, genomic_del1_seq_loc):
    """Create a test fixture for genomic del RSE."""
    _id = "ga4gh:VA.6fIEZ3R2W4wIaltUX1jyw9ap5YN6oGDT"
    genomic_del1["variation"] = {
        "type": "Allele",
        "_id": _id,
        "location": genomic_del1_seq_loc,
        "state": {
            "type": "RepeatedSequenceExpression",
            "seq_expr": {
                "type": "DerivedSequenceExpression",
                "location": genomic_del1_seq_loc,
                "reverse_complement": False
            },
            "count": {
                "type": "Number",
                "value": 0
            }
        }
    }
    genomic_del1["variation_id"] = _id
    return VariationDescriptor(**genomic_del1)


@pytest.fixture(scope="module")
def genomic_del1_free_text(vhl_gene_context):
    """Create test fixture containing params for genomic del VD."""
    params = {
        "id": "normalize.variation:VHL%20g.10191495del",
        "type": "VariationDescriptor",
        "variation_id": "",
        "variation": dict(),
        "molecule_context": "transcript",
        "structural_type": "SO:0000159",
        "vrs_ref_allele_seq": "T",
        "gene_context": vhl_gene_context
    }
    return params


@pytest.fixture(scope="module")
def genomic_del1_free_text_seq_loc():
    """Create genomic del1 free text sequence location"""
    return {
        "_id": "ga4gh:VSL.90XXYrpPCTvaFcyb7L4W4EcE9OexpmNv",
        "sequence_id": "ga4gh:SQ.xBKOKptLLDr-k4hTyCetvARn16pDS_rW",
        "interval": {
            "type": "SequenceInterval",
            "start": {"value": 557, "type": "Number"},
            "end": {"value": 558, "type": "Number"},
        },
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_del1_free_text_lse(genomic_del1_free_text,
                               genomic_del1_free_text_seq_loc):
    """Create a test fixture for genomic del LSE."""
    _id = "ga4gh:VA.DdtLZ_d22R0O0VU020WcCLvNhXNZtU2j"
    genomic_del1_free_text["variation_id"] = _id
    genomic_del1_free_text["variation"] = {
        "type": "Allele",
        "_id": _id,
        "location": genomic_del1_free_text_seq_loc,
        "state": {
            "type": "LiteralSequenceExpression",
            "sequence": ""
        }
    }
    return VariationDescriptor(**genomic_del1_free_text)


@pytest.fixture(scope="module")
def genomic_del1_free_text_vac(genomic_del1_free_text,
                               genomic_del1_free_text_seq_loc):
    """Create a test fixture for genomic del absolute CNV."""
    _id = "ga4gh:VAC.Wgvw8a4LXRY4d1jopC5tZjUlaEKci5Ai"
    genomic_del1_free_text["variation"] = {
        "type": "AbsoluteCopyNumber",
        "_id": _id,
        "subject": genomic_del1_free_text_seq_loc,
        "copies": {"type": "Number", "value": 1}
    }
    genomic_del1_free_text["variation_id"] = _id
    return VariationDescriptor(**genomic_del1_free_text)


@pytest.fixture(scope="module")
def genomic_del1_free_text_rse(genomic_del1_free_text,
                               genomic_del1_free_text_seq_loc):
    """Create a test fixture for genomic del RSE."""
    _id = "ga4gh:VA.o8kDqsCKM-cakyb_Pa5HWXLFxKqHtZA4"
    genomic_del1_free_text["variation"] = {
        "type": "Allele",
        "_id": _id,
        "location": genomic_del1_free_text_seq_loc,
        "state": {
            "type": "RepeatedSequenceExpression",
            "seq_expr": {
                "type": "DerivedSequenceExpression",
                "location": genomic_del1_free_text_seq_loc,
                "reverse_complement": False
            },
            "count": {
                "type": "Number",
                "value": 0
            }
        }
    }
    genomic_del1_free_text["variation_id"] = _id
    return VariationDescriptor(**genomic_del1_free_text)


@pytest.fixture(scope="module")
def genomic_del2():
    """Create test fixture containing params for genomic del VD."""
    params = {
        "id": "normalize.variation:NC_000003.12%3Ag.10146595_10146613del",
        "type": "VariationDescriptor",
        "variation_id": "",
        "variation": dict(),
        "molecule_context": "genomic",
        "structural_type": "SO:0000159",
        "vrs_ref_allele_seq": "ATGTTGACGGACAGCCTAT"
    }
    return params


@pytest.fixture(scope="module")
def genomic_del2_lse(genomic_del2, genomic_del2_seq_loc):
    """Create a test fixture for genomic del LSE."""
    _id = "ga4gh:VA.CSWNhR5w_geMmJTxkbO3UCLCvT0S2Ypx"
    genomic_del2["variation_id"] = _id
    genomic_del2["variation"] = {
        "type": "Allele",
        "_id": _id,
        "location": genomic_del2_seq_loc,
        "state": {
            "type": "LiteralSequenceExpression",
            "sequence": ""
        }
    }
    return VariationDescriptor(**genomic_del2)


@pytest.fixture(scope="module")
def genomic_del2_vac(genomic_del2, genomic_del2_38_vac):
    """Create a test fixture for genomic del absolute CNV."""
    genomic_del2["variation"] = genomic_del2_38_vac
    genomic_del2["variation_id"] = genomic_del2["variation"]["_id"]
    return VariationDescriptor(**genomic_del2)


@pytest.fixture(scope="module")
def genomic_del2_vrc(genomic_del2, genomic_del2_seq_loc):
    """Create a test fixture for genomic del relative CNV."""
    genomic_del2["variation"] = {
        "type": "RelativeCopyNumber",
        "_id": "ga4gh:VRC.WtEs7KGUIjxWaV9Wx0lCCgDyyWVl-ykM",
        "subject": genomic_del2_seq_loc,
        "relative_copy_class": "complete loss"
    }
    genomic_del2["variation_id"] = genomic_del2["variation"]["_id"]
    return VariationDescriptor(**genomic_del2)


@pytest.fixture(scope="module")
def genomic_del2_rse(genomic_del2, genomic_del2_seq_loc):
    """Create a test fixture for genomic del RSE."""
    _id = "ga4gh:VA.aQeEhbisBWYrzVbf3-VPOZtGJu1vKmfx"
    genomic_del2["variation"] = {
        "type": "Allele",
        "_id": _id,
        "location": genomic_del2_seq_loc,
        "state": {
            "type": "RepeatedSequenceExpression",
            "seq_expr": {
                "type": "DerivedSequenceExpression",
                "location": genomic_del2_seq_loc,
                "reverse_complement": False
            },
            "count": {
                "type": "Number",
                "value": 0
            }
        }
    }
    genomic_del2["variation_id"] = _id
    return VariationDescriptor(**genomic_del2)


@pytest.fixture(scope="module")
def genomic_del2_free_text(vhl_gene_context):
    """Create test fixture containing params for genomic del VD."""
    params = {
        "id": "normalize.variation:VHL%20g.10188279_10188297del",
        "type": "VariationDescriptor",
        "variation_id": "",
        "variation": dict(),
        "molecule_context": "transcript",
        "structural_type": "SO:0000159",
        "vrs_ref_allele_seq": "ATGTTGACGGACAGCCTAT",
        "gene_context": vhl_gene_context
    }
    return params


@pytest.fixture(scope="module")
def genomic_del2_free_text_seq_loc():
    """Create genomic del2 free text sequence location"""
    return {
        "_id": "ga4gh:VSL.9fIfzZxIhfm4AlUhBlU9PswkG8ei57lR",
        "sequence_id": "ga4gh:SQ.xBKOKptLLDr-k4hTyCetvARn16pDS_rW",
        "interval": {
            "type": "SequenceInterval",
            "start": {"value": 491, "type": "Number"},
            "end": {"value": 510, "type": "Number"},
        },
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_del2_free_text_default(genomic_del2_free_text,
                                   genomic_del2_free_text_seq_loc):
    """Create a test fixture for genomic del default and LSE."""
    _id = "ga4gh:VA.V0TeIIZTlBnFTIc64hqxzjbhAH3I4VZI"
    genomic_del2_free_text["variation_id"] = _id
    genomic_del2_free_text["variation"] = {
        "type": "Allele",
        "_id": _id,
        "location": genomic_del2_free_text_seq_loc,
        "state": {
            "type": "LiteralSequenceExpression",
            "sequence": ""
        }
    }
    return VariationDescriptor(**genomic_del2_free_text)


@pytest.fixture(scope="module")
def genomic_del2_free_text_cnv(genomic_del2_free_text,
                               genomic_del2_free_text_seq_loc):
    """Create a test fixture for genomic del CNV."""
    _id = "ga4gh:VAC.UKlDhLtC7JSqSCJA0cuusgWhVcFd0hSS"
    genomic_del2_free_text["variation"] = {
        "type": "AbsoluteCopyNumber",
        "_id": _id,
        "subject": genomic_del2_free_text_seq_loc,
        "copies": {"type": "Number", "value": 1}
    }
    genomic_del2_free_text["variation_id"] = _id
    return VariationDescriptor(**genomic_del2_free_text)


@pytest.fixture(scope="module")
def genomic_del2_free_text_rse(genomic_del2_free_text,
                               genomic_del2_free_text_seq_loc):
    """Create a test fixture for genomic del RSE."""
    _id = "ga4gh:VA.uED5jM7zwbFLiXfCufVuwIs2ufkPF2KJ"
    genomic_del2_free_text["variation"] = {
        "type": "Allele",
        "_id": _id,
        "location": genomic_del2_free_text_seq_loc,
        "state": {
            "type": "RepeatedSequenceExpression",
            "seq_expr": {
                "type": "DerivedSequenceExpression",
                "location": genomic_del2_free_text_seq_loc,
                "reverse_complement": False
            },
            "count": {
                "type": "Number",
                "value": 0
            }
        }
    }
    genomic_del2_free_text["variation_id"] = _id
    return VariationDescriptor(**genomic_del2_free_text)


@pytest.fixture(scope="module")
def genomic_del3():
    """Create test fixture containing params for genomic del VD."""
    params = {
        "id": "normalize.variation:NC_000023.11%3Ag.%2831060227_31100351%29_%2833274278_33417151%29del",  # noqa: E501
        "type": "VariationDescriptor",
        "variation_id": "",
        "variation": dict(),
        "molecule_context": "genomic",
        "structural_type": "SO:0001743",
        "vrs_ref_allele_seq": None
    }
    return params


@pytest.fixture(scope="module")
def genomic_del3_vrc(genomic_del3, genomic_del3_dup3_loc):
    """Create a test fixture for genomic del relative cnv."""
    _id = "ga4gh:VRC.7uGeBoQVduNHyz3dDmnTDaVsDhCMAaZe"
    genomic_del3["variation"] = {
        "type": "RelativeCopyNumber",
        "_id": _id,
        "subject": genomic_del3_dup3_loc,
        "relative_copy_class": "partial loss"
    }
    genomic_del3["variation_id"] = _id
    return VariationDescriptor(**genomic_del3)


@pytest.fixture(scope="module")
def genomic_del3_vac(genomic_del3, genomic_del3_dup3_loc):
    """Create a test fixture for genomic del absolute cnv."""
    _id = "ga4gh:VAC.cQATJ6a1uGwXOHu-advv8lRsMgjNLKul"
    genomic_del3["variation"] = {
        "type": "AbsoluteCopyNumber",
        "_id": _id,
        "subject": genomic_del3_dup3_loc,
        "copies": {"type": "Number", "value": 2}
    }
    genomic_del3["variation_id"] = _id
    return VariationDescriptor(**genomic_del3)


@pytest.fixture(scope="module")
def genomic_del3_rse_lse(genomic_del3):
    """Create test fixture for genomic del rse and lse."""
    params = {
        "id": genomic_del3["id"],
        "type": genomic_del3["type"],
        "variation": {
            "_id": "ga4gh:VT.tmA3mpMy9HKUweaB8aYsq6uuejEx9iK7",
            "type": "Text",
            "definition": "NC_000023.11:g.(31060227_31100351)_(33274278_33417151)del"  # noqa: E501
        }
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_del3_free_text():
    """Create test fixture containing params for genomic del VD."""
    params = {
        "id": "normalize.variation:EFNB1%20g.%2868839265_68839268%29_%2868841120_68841125%29del",  # noqa: E501
        "type": "VariationDescriptor",
        "variation_id": "",
        "variation": dict(),
        "molecule_context": "genomic",
        "structural_type": "SO:0001743",
        "vrs_ref_allele_seq": None,
        "gene_context": {
            "id": "normalize.gene:EFNB1",
            "type": "GeneDescriptor",
            "label": "EFNB1",
            "xrefs": [
                "ensembl:ENSG00000090776",
                "ncbigene:1947"
            ],
            "alternate_labels": [
                "EPLG2",
                "Elk-L",
                "CFND",
                "CFNS",
                "EFB1",
                "LERK2",
                "EFL3"
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
                    "value": "ephrin B1"
                },
                {
                    "type": "Extension",
                    "name": "hgnc_locations",
                    "value": [
                        {
                            "species_id": "taxonomy:9606",
                            "interval": {
                                "type": "CytobandInterval",
                                "start": "q13.1",
                                "end": "q13.1"
                            },
                            "_id": "ga4gh:VCL.2INIrDKtMs_uh9lw8BWws2AMpzlbMaBB",
                            "type": "ChromosomeLocation",
                            "chr": "X"
                        }
                    ]
                },
                {
                    "type": "Extension",
                    "name": "ensembl_locations",
                    "value": [
                        {
                            "_id": "ga4gh:VSL.BemPsrGqpo2gqkjdItXADblyNGSNkUwB",
                            "type": "SequenceLocation",
                            "sequence_id": "ga4gh:SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
                            "interval": {
                                "start": {"type": "Number", "value": 68829020},
                                "end": {"type": "Number", "value": 68842160},
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
                            "_id": "ga4gh:VCL.2INIrDKtMs_uh9lw8BWws2AMpzlbMaBB",
                            "type": "ChromosomeLocation",
                            "species_id": "taxonomy:9606",
                            "chr": "X",
                            "interval": {
                                "end": "q13.1",
                                "start": "q13.1",
                                "type": "CytobandInterval"
                            }
                        },
                        {
                            "_id": "ga4gh:VSL.BemPsrGqpo2gqkjdItXADblyNGSNkUwB",
                            "type": "SequenceLocation",
                            "sequence_id": "ga4gh:SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
                            "interval": {
                                "start": {"type": "Number", "value": 68829020},
                                "end": {"type": "Number", "value": 68842160},
                                "type": "SequenceInterval"
                            }
                        }
                    ]
                },
                {
                    "type": "Extension",
                    "name": "associated_with",
                    "value": [
                        "uniprot:P98172",
                        "ena.embl:U09303",
                        "ucsc:uc004dxd.5",
                        "omim:300035",
                        "pubmed:16526919",
                        "refseq:NM_004429",
                        "ccds:CCDS14391",
                        "orphanet:121305",
                        "vega:OTTHUMG00000021751",
                        "iuphar:4913",
                        "pubmed:7774950"
                    ]
                },
                {
                    "type": "Extension",
                    "name": "previous_symbols",
                    "value": [
                        "CFNS",
                        "EPLG2"
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
            ],
            "gene_id": "hgnc:3226"
        }
    }
    return params


@pytest.fixture(scope="module")
def genomic_del3_free_text_subject():
    """Create test fixture for genomic del3 free text subject"""
    return {
        "_id": "ga4gh:VSL.gqWO-oN2bMIXm_YuZR4_beT57QN-kRGJ",
        "sequence_id": "ga4gh:SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
        "interval": {
            "type": "SequenceInterval",
            "start": {
                "min": 68839264,
                "max": 68839267,
                "type": "DefiniteRange"
            },
            "end": {
                "min": 68841121,
                "max": 68841126,
                "type": "DefiniteRange"
            }
        },
        "type": "SequenceLocation"
    }


@pytest.fixture(scope="module")
def genomic_del3_free_text_vrc(genomic_del3_free_text, genomic_del3_free_text_subject):
    """Create a test fixture for genomic del relative cnv."""
    _id = "ga4gh:VRC.hHXnhBGUcR870shryKPt8qibKpaRjDE4"
    genomic_del3_free_text["variation"] = {
        "type": "RelativeCopyNumber",
        "_id": _id,
        "subject": genomic_del3_free_text_subject,
        "relative_copy_class": "partial loss"
    }
    genomic_del3_free_text["variation_id"] = _id
    return VariationDescriptor(**genomic_del3_free_text)


@pytest.fixture(scope="module")
def genomic_del3_free_text_vac(genomic_del3_free_text, genomic_del3_free_text_subject):
    """Create a test fixture for genomic del absolute cnv."""
    _id = "ga4gh:VAC.14iSDwh4eU37jPuiBVZ_tdtRdQAqbjqG"
    genomic_del3_free_text["variation"] = {
        "type": "AbsoluteCopyNumber",
        "_id": _id,
        "subject": genomic_del3_free_text_subject,
        "copies": {"type": "Number", "value": 2}
    }
    genomic_del3_free_text["variation_id"] = _id
    return VariationDescriptor(**genomic_del3_free_text)


@pytest.fixture(scope="module")
def genomic_del3_free_text_rse_lse(genomic_del3_free_text):
    """Create test fixture for genomic del rse and lse."""
    params = {
        "id": genomic_del3_free_text["id"],
        "type": genomic_del3_free_text["type"],
        "variation": {
            "_id": "ga4gh:VT.9mGg0U_Z7NZCFV3jrLdGxSQU03g7z3Z1",
            "type": "Text",
            "definition": "EFNB1 g.(68839265_68839268)_(68841120_68841125)del"
        }
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_del4():
    """Create test fixture containing params for genomic del VD."""
    params = {
        "id": "normalize.variation:NC_000023.11%3Ag.%28%3F_31120496%29_%2833339477_%3F%29del",  # noqa: E501
        "type": "VariationDescriptor",
        "variation_id": "",
        "variation": dict(),
        "molecule_context": "genomic",
        "structural_type": "SO:0001743",
        "vrs_ref_allele_seq": None
    }
    return params


@pytest.fixture(scope="module")
def genomic_del4_vrc(genomic_del4, genomic_del4_seq_loc):
    """Create a test fixture for genomic del relative cnv."""
    _id = "ga4gh:VRC.dtwRjvChZ6LuyDXyWTnVGQJidyfJsQfe"
    genomic_del4["variation"] = {
        "type": "RelativeCopyNumber",
        "_id": _id,
        "subject": genomic_del4_seq_loc,
        "relative_copy_class": "partial loss"
    }
    genomic_del4["variation_id"] = _id
    return VariationDescriptor(**genomic_del4)


@pytest.fixture(scope="module")
def genomic_del4_vac(genomic_del4, genomic_del4_seq_loc):
    """Create a test fixture for genomic del absolute cnv."""
    _id = "ga4gh:VAC.9AmYiW27hEz5V2fCisXqbFJpYeyGaZBv"
    genomic_del4["variation"] = {
        "type": "AbsoluteCopyNumber",
        "_id": _id,
        "subject": genomic_del4_seq_loc,
        "copies": {"type": "Number", "value": 1}
    }
    genomic_del4["variation_id"] = _id
    return VariationDescriptor(**genomic_del4)


@pytest.fixture(scope="module")
def genomic_del4_rse_lse(genomic_del4):
    """Create test fixture for genomic dup rse and lse."""
    params = {
        "id": genomic_del4["id"],
        "type": genomic_del4["type"],
        "variation": {
            "_id": "ga4gh:VT.whBY5P24WVxF1wneDcI8x8btqorJUWXQ",
            "type": "Text",
            "definition": "NC_000023.11:g.(?_31120496)_(33339477_?)del"
        }
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_del4_free_text():
    """Create test fixture containing params for genomic del VD."""
    params = {
        "id": "normalize.variation:COL4A4%20g.%28%3F_227022028%29_%28227025830_%3F%29del",  # noqa: E501
        "type": "VariationDescriptor",
        "variation_id": "",
        "variation": dict(),
        "molecule_context": "genomic",
        "structural_type": "SO:0001743",
        "vrs_ref_allele_seq": None,
        "gene_context": {
            "id": "normalize.gene:COL4A4",
            "type": "GeneDescriptor",
            "label": "COL4A4",
            "xrefs": [
                "ensembl:ENSG00000081052",
                "ncbigene:1286"
            ],
            "alternate_labels": [
                "BFH",
                "ATS2",
                "CA44"
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
                    "value": "collagen type IV alpha 4 chain"
                },
                {
                    "type": "Extension",
                    "name": "hgnc_locations",
                    "value": [
                        {
                            "species_id": "taxonomy:9606",
                            "interval": {
                                "type": "CytobandInterval",
                                "start": "q36.3",
                                "end": "q36.3"
                            },
                            "_id": "ga4gh:VCL.1raDfW4j_diAb62KX4wnjRGD3A6va_BB",
                            "type": "ChromosomeLocation",
                            "chr": "2"
                        }
                    ]
                },
                {
                    "type": "Extension",
                    "name": "ensembl_locations",
                    "value": [
                        {
                            "_id": "ga4gh:VSL.0O5Ozp39UHxxgO5XK4gRcrkX_VwlFkmM",
                            "type": "SequenceLocation",
                            "sequence_id": "ga4gh:SQ.pnAqCRBrTsUoBghSD1yp_jXWSmlbdh4g",
                            "interval": {
                                "start": {"type": "Number", "value": 227002713},
                                "end": {"type": "Number", "value": 227164453},
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
                            "_id": "ga4gh:VCL.1raDfW4j_diAb62KX4wnjRGD3A6va_BB",
                            "type": "ChromosomeLocation",
                            "species_id": "taxonomy:9606",
                            "chr": "2",
                            "interval": {
                                "end": "q36.3",
                                "start": "q36.3",
                                "type": "CytobandInterval"
                            }
                        },
                        {
                            "_id": "ga4gh:VSL.yQ4kA-x2w3oTE6CbiKVJlldbzdrOi0yU",
                            "type": "SequenceLocation",
                            "sequence_id": "ga4gh:SQ.pnAqCRBrTsUoBghSD1yp_jXWSmlbdh4g",
                            "interval": {
                                "start": {"type": "Number", "value": 226967359},
                                "end": {"type": "Number", "value": 227164488},
                                "type": "SequenceInterval"
                            }
                        }
                    ]
                },
                {
                    "type": "Extension",
                    "name": "associated_with",
                    "value": [
                        "omim:120131",
                        "ucsc:uc061teu.1",
                        "ccds:CCDS42828",
                        "orphanet:120720",
                        "pubmed:1639407",
                        "vega:OTTHUMG00000149892",
                        "refseq:NM_000092",
                        "uniprot:P53420"
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
            ],
            "gene_id": "hgnc:2206"
        }
    }
    return params


@pytest.fixture(scope="module")
def genomic_del4_free_text_subject():
    """Create test fixture for genomic del4 free text subject"""
    return {
        "_id": "ga4gh:VSL.s4_6D986zFS0HIBuEDFl5aq2-VCl45h1",
        "sequence_id": "ga4gh:SQ.pnAqCRBrTsUoBghSD1yp_jXWSmlbdh4g",
        "interval": {
            "type": "SequenceInterval",
            "start": {
                "value": 227022027,
                "comparator": "<=",
                "type": "IndefiniteRange"
            },
            "end": {
                "value": 227025830,
                "comparator": ">=",
                "type": "IndefiniteRange"
            }
        },
        "type": "SequenceLocation"
    }


@pytest.fixture(scope="module")
def genomic_del4_free_text_vrc(genomic_del4_free_text, genomic_del4_free_text_subject):
    """Create a test fixture for genomic del relative cnv."""
    _id = "ga4gh:VRC.xNHxMmWtqOkKQS7hpJlCN4jjGgBsdk7b"
    genomic_del4_free_text["variation"] = {
        "type": "RelativeCopyNumber",
        "_id": _id,
        "subject": genomic_del4_free_text_subject,
        "relative_copy_class": "partial loss"
    }
    genomic_del4_free_text["variation_id"] = _id
    return VariationDescriptor(**genomic_del4_free_text)


@pytest.fixture(scope="module")
def genomic_del4_free_text_vac(genomic_del4_free_text, genomic_del4_free_text_subject):
    """Create a test fixture for genomic del absolute cnv."""
    _id = "ga4gh:VAC.DrUaDS4OCFqKSvvXhJse2PCduwnZMwkU"
    genomic_del4_free_text["variation"] = {
        "type": "AbsoluteCopyNumber",
        "_id": _id,
        "subject": genomic_del4_free_text_subject,
        "copies": {"type": "Number", "value": 1}
    }
    genomic_del4_free_text["variation_id"] = _id
    return VariationDescriptor(**genomic_del4_free_text)


@pytest.fixture(scope="module")
def genomic_del4_free_text_rse_lse(genomic_del4_free_text):
    """Create test fixture for genomic dup rse and lse."""
    params = {
        "id": genomic_del4_free_text["id"],
        "type": genomic_del4_free_text["type"],
        "variation": {
            "_id": "ga4gh:VT.lT0rFYhOGFLA9MYA8ypnCf5q-CkV8dJv",
            "type": "Text",
            "definition": "COL4A4 g.(?_227022028)_(227025830_?)del"
        }
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_uncertain_del_2():
    """Create a genomic uncertain deletion on chr 2 test fixture."""
    params = {
        "id": "normalize.variation:NC_000002.12%3Ag.%28%3F_110104900%29_%28110207160_%3F%29del",  # noqa: E501
        "type": "VariationDescriptor",
        "variation_id": "ga4gh:VRC.GYVDxHLsPXu54Ri8x2kFqRWhGBg1RgFc",
        "variation": {
            "_id": "ga4gh:VRC.GYVDxHLsPXu54Ri8x2kFqRWhGBg1RgFc",
            "subject": {
                "_id": "ga4gh:VSL.75GQmJvq7dyP9-wom8Jffjk0Q9Le7Q9O",
                "sequence_id": "ga4gh:SQ.pnAqCRBrTsUoBghSD1yp_jXWSmlbdh4g",
                "interval": {
                    "start": {
                        "value": 110104899,
                        "comparator": "<=",
                        "type": "IndefiniteRange"
                    },
                    "end": {
                        "value": 110207160,
                        "comparator": ">=",
                        "type": "IndefiniteRange"
                    },
                    "type": "SequenceInterval"
                },
                "type": "SequenceLocation"
            },
            "relative_copy_class": "partial loss",
            "type": "RelativeCopyNumber"
        },
        "molecule_context": "genomic",
        "structural_type": "SO:0001743"
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_uncertain_del_y():
    """Create a genomic uncertain deletion on chr Y test fixture."""
    params = {
        "id": "normalize.variation:NC_000024.10%3Ag.%28%3F_14076802%29_%2857165209_%3F%29del",  # noqa: E501
        "type": "VariationDescriptor",
        "variation_id": "ga4gh:VRC.L0F1kbaN2aaGbLwum29WEWNwji5tit2E",
        "variation": {
            "_id": "ga4gh:VRC.L0F1kbaN2aaGbLwum29WEWNwji5tit2E",
            "subject": {
                "_id": "ga4gh:VSL.1xIN_RumlXTIsdTWvyJznzuzxJlwUfiD",
                "sequence_id": "ga4gh:SQ.8_liLu1aycC0tPQPFmUaGXJLDs5SbPZ5",
                "interval": {
                    "start": {
                        "value": 14076801,
                        "comparator": "<=",
                        "type": "IndefiniteRange"
                    },
                    "end": {
                        "value": 57165209,
                        "comparator": ">=",
                        "type": "IndefiniteRange"
                    },
                    "type": "SequenceInterval"
                },
                "type": "SequenceLocation"
            },
            "relative_copy_class": "partial loss",
            "type": "RelativeCopyNumber"
        },
        "molecule_context": "genomic",
        "structural_type": "SO:0001743"
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_del5():
    """Create test fixture containing params for genomic del VD."""
    params = {
        "id": "normalize.variation:NC_000023.11%3Ag.%28%3F_18575354%29_18653629del",  # noqa: E501
        "type": "VariationDescriptor",
        "variation_id": "",
        "variation": dict(),
        "molecule_context": "genomic",
        "structural_type": "SO:0001743",
        "vrs_ref_allele_seq": None
    }
    return params


def genomic_del5_abs_cnv(params, genomic_del5_seq_loc):
    """Create genomic del5 absolute cnv"""
    _id = "ga4gh:VAC.foWmgooK2J3amtQE6xSGSAnlGzzBsIjy"
    params["variation"] = {
        "type": "AbsoluteCopyNumber",
        "_id": _id,
        "subject": genomic_del5_seq_loc,
        "copies": {"type": "Number", "value": 3}
    }
    params["variation_id"] = _id


def genomic_del5_rel_cnv(params, genomic_del5_seq_loc):
    """Create genomic del5 relative cnv"""
    _id = "ga4gh:VRC.fxZ8m34viXiB5fjxuw2xXYItk5Xi5cve"
    params["variation"] = {
        "type": "RelativeCopyNumber",
        "_id": _id,
        "subject": genomic_del5_seq_loc,
        "relative_copy_class": "partial loss"
    }
    params["variation_id"] = _id


@pytest.fixture(scope="module")
def genomic_del5_vrc(genomic_del5, genomic_del5_seq_loc):
    """Create a test fixture for genomic del relative cnv."""
    genomic_del5_rel_cnv(genomic_del5, genomic_del5_seq_loc)
    return VariationDescriptor(**genomic_del5)


@pytest.fixture(scope="module")
def genomic_del5_vac(genomic_del5, genomic_del5_seq_loc):
    """Create a test fixture for genomic del absolute cnv."""
    genomic_del5_abs_cnv(genomic_del5, genomic_del5_seq_loc)
    return VariationDescriptor(**genomic_del5)


@pytest.fixture(scope="module")
def genomic_del5_rse_lse(genomic_del5):
    """Create test fixture for genomic del rse and lse."""
    params = {
        "id": genomic_del5["id"],
        "type": genomic_del5["type"],
        "variation": {
            "_id": "ga4gh:VT.xCLHh3GpCebrP6KDMsWZRdIiW7Sti27H",
            "type": "Text",
            "definition": "NC_000023.11:g.(?_18575354)_18653629del"
        }
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_del5_free_text():
    """Create test fixture containing params for genomic del VD."""
    params = {
        "id": "normalize.variation:CDKL5%20g.%28%3F_18575354%29_18653629del",
        "type": "VariationDescriptor",
        "variation_id": "",
        "variation": dict(),
        "molecule_context": "genomic",
        "structural_type": "SO:0001743",
        "vrs_ref_allele_seq": None,
        "gene_context": {
            "id": "normalize.gene:CDKL5",
            "type": "GeneDescriptor",
            "label": "CDKL5",
            "xrefs": [
                "ensembl:ENSG00000008086",
                "ncbigene:6792"
            ],
            "alternate_labels": [
                "STK9",
                "ISSX",
                "EIEE2",
                "CFAP247",
                "DEE2"
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
                    "value": "cyclin dependent kinase like 5"
                },
                {
                    "type": "Extension",
                    "name": "hgnc_locations",
                    "value": [
                        {
                            "species_id": "taxonomy:9606",
                            "interval": {
                                "type": "CytobandInterval",
                                "start": "p22.13",
                                "end": "p22.13"
                            },
                            "_id": "ga4gh:VCL.BzhQOPmaVZVLol6JOVltNZrsv0XRekWR",
                            "type": "ChromosomeLocation",
                            "chr": "X"
                        }
                    ]
                },
                {
                    "type": "Extension",
                    "name": "ensembl_locations",
                    "value": [
                        {
                            "_id": "ga4gh:VSL.JQbEJe-xBEW39qI60yfBUC3siCPrf5NK",
                            "type": "SequenceLocation",
                            "sequence_id": "ga4gh:SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
                            "interval": {
                                "start": {"type": "Number", "value": 18425582},
                                "end": {"type": "Number", "value": 18653629},
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
                            "_id": "ga4gh:VCL.BzhQOPmaVZVLol6JOVltNZrsv0XRekWR",
                            "type": "ChromosomeLocation",
                            "species_id": "taxonomy:9606",
                            "chr": "X",
                            "interval": {
                                "end": "p22.13",
                                "start": "p22.13",
                                "type": "CytobandInterval"
                            }
                        },
                        {
                            "_id": "ga4gh:VSL._lU2syu-OtdipRW6XFXBBpZnsjdoKKN0",
                            "type": "SequenceLocation",
                            "sequence_id": "ga4gh:SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
                            "interval": {
                                "start": {"type": "Number", "value": 18425607},
                                "end": {"type": "Number", "value": 18653629},
                                "type": "SequenceInterval"
                            }
                        }
                    ]
                },
                {
                    "type": "Extension",
                    "name": "associated_with",
                    "value": [
                        "orphanet:119297",
                        "ena.embl:Y15057",
                        "pubmed:16935860",
                        "ccds:CCDS83458",
                        "omim:300203",
                        "vega:OTTHUMG00000021214",
                        "uniprot:O76039",
                        "pubmed:9721213",
                        "refseq:NM_003159",
                        "ccds:CCDS14186",
                        "iuphar:1986",
                        "ucsc:uc004cyn.4"
                    ]
                },
                {
                    "type": "Extension",
                    "name": "previous_symbols",
                    "value": [
                        "STK9"
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
            ],
            "gene_id": "hgnc:11411"
        }
    }
    return params


@pytest.fixture(scope="module")
def genomic_del5_free_text_vrc(genomic_del5_free_text, genomic_del5_seq_loc):
    """Create a test fixture for genomic del relative cnv."""
    genomic_del5_rel_cnv(genomic_del5_free_text, genomic_del5_seq_loc)
    return VariationDescriptor(**genomic_del5_free_text)


@pytest.fixture(scope="module")
def genomic_del5_free_text_vac(genomic_del5_free_text, genomic_del5_seq_loc):
    """Create a test fixture for genomic del absolute cnv."""
    genomic_del5_abs_cnv(genomic_del5_free_text, genomic_del5_seq_loc)
    return VariationDescriptor(**genomic_del5_free_text)


@pytest.fixture(scope="module")
def genomic_del5_free_text_rse_lse(genomic_del5_free_text):
    """Create test fixture for genomic del rse and lse."""
    params = {
        "id": genomic_del5_free_text["id"],
        "type": genomic_del5_free_text["type"],
        "variation": {
            "_id": "ga4gh:VT.xCLHh3GpCebrP6KDMsWZRdIiW7Sti27H",
            "type": "Text",
            "definition": "NC_000023.11:g.(?_18575354)_18653629del"
        }
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_del6():
    """Create test fixture containing params for genomic del VD."""
    params = {
        "id": "normalize.variation:NC_000006.12%3Ag.133462764_%28133464858_%3F%29del",  # noqa: E501
        "type": "VariationDescriptor",
        "variation_id": "",
        "variation": dict(),
        "molecule_context": "genomic",
        "structural_type": "SO:0001743",
        "vrs_ref_allele_seq": None
    }
    return params


def genomic_del6_rel_cnv(params, genomic_del6_seq_loc):
    """Create genomic del6 relative cnv"""
    _id = "ga4gh:VRC.mO5A42JUswyNIbJPLuE-2fhquQEZpXxn"
    params["variation"] = {
        "type": "RelativeCopyNumber",
        "_id": _id,
        "subject": genomic_del6_seq_loc,
        "relative_copy_class": "partial loss"
    }
    params["variation_id"] = _id


def genomic_del6_abs_cnv(params, genomic_del6_seq_loc):
    """Create genomic del6 absolute cnv"""
    _id = "ga4gh:VAC.6RkHgDOiRMZKMKgI6rmG9C3T6WuMhcex"
    params["variation"] = {
        "type": "AbsoluteCopyNumber",
        "_id": _id,
        "subject": genomic_del6_seq_loc,
        "copies": {"type": "Number", "value": 1}
    }
    params["variation_id"] = _id


@pytest.fixture(scope="module")
def genomic_del6_vrc(genomic_del6, genomic_del6_seq_loc):
    """Create a test fixture for genomic del relative cnv."""
    genomic_del6_rel_cnv(genomic_del6, genomic_del6_seq_loc)
    return VariationDescriptor(**genomic_del6)


@pytest.fixture(scope="module")
def genomic_del6_vac(genomic_del6, genomic_del6_seq_loc):
    """Create a test fixture for genomic del absolute cnv."""
    genomic_del6_abs_cnv(genomic_del6, genomic_del6_seq_loc)
    return VariationDescriptor(**genomic_del6)


@pytest.fixture(scope="module")
def genomic_del6_rse_lse(genomic_del6):
    """Create test fixture for genomic del rse and lse."""
    params = {
        "id": genomic_del6["id"],
        "type": genomic_del6["type"],
        "variation": {
            "_id": "ga4gh:VT.Df49jbB-kZ2LSm180uA9wn4TT_p215yX",
            "type": "Text",
            "definition": "NC_000006.12:g.133462764_(133464858_?)del"
        }
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_del6_free_text():
    """Create test fixture containing params for genomic del VD."""
    params = {
        "id": "normalize.variation:EYA4%20g.133462764_%28133464858_%3F%29del",
        "type": "VariationDescriptor",
        "variation_id": "",
        "variation": dict(),
        "molecule_context": "genomic",
        "structural_type": "SO:0001743",
        "vrs_ref_allele_seq": None,
        "gene_context": {
            "id": "normalize.gene:EYA4",
            "type": "GeneDescriptor",
            "label": "EYA4",
            "xrefs": [
                "ensembl:ENSG00000112319",
                "ncbigene:2070"
            ],
            "alternate_labels": [
                "CMD1J",
                "DFNA10"
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
                    "value": "EYA transcriptional coactivator "
                             "and phosphatase 4"
                },
                {
                    "type": "Extension",
                    "name": "hgnc_locations",
                    "value": [
                        {
                            "species_id": "taxonomy:9606",
                            "interval": {
                                "type": "CytobandInterval",
                                "start": "q23.2",
                                "end": "q23.2"
                            },
                            "_id": "ga4gh:VCL.h3RPutxMQk5_6dIPcN5GL8KRahoTi9fm",
                            "type": "ChromosomeLocation",
                            "chr": "6"
                        }
                    ]
                },
                {
                    "type": "Extension",
                    "name": "ensembl_locations",
                    "value": [
                        {
                            "_id": "ga4gh:VSL.hswgyI183l3KGlAaOOGotbOU95oFNFxo",
                            "type": "SequenceLocation",
                            "sequence_id": "ga4gh:SQ.0iKlIQk2oZLoeOG9P1riRU6hvL5Ux8TV",
                            "interval": {
                                "start": {"type": "Number", "value": 133240513},
                                "end": {"type": "Number", "value": 133532128},
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
                            "_id": "ga4gh:VCL.h3RPutxMQk5_6dIPcN5GL8KRahoTi9fm",
                            "type": "ChromosomeLocation",
                            "species_id": "taxonomy:9606",
                            "chr": "6",
                            "interval": {
                                "end": "q23.2",
                                "start": "q23.2",
                                "type": "CytobandInterval"
                            }
                        },
                        {
                            "_id": "ga4gh:VSL.D6_ZpKkLD2mGUBoCwA8pEWWzUCpdXZg_",
                            "type": "SequenceLocation",
                            "sequence_id": "ga4gh:SQ.0iKlIQk2oZLoeOG9P1riRU6hvL5Ux8TV",
                            "interval": {
                                "start": {"type": "Number", "value": 133240592},
                                "end": {"type": "Number", "value": 133532128},
                                "type": "SequenceInterval"
                            }
                        }
                    ]
                },
                {
                    "type": "Extension",
                    "name": "associated_with",
                    "value": [
                        "pubmed:9887327",
                        "ccds:CCDS43506",
                        "ccds:CCDS75521",
                        "pubmed:11159937",
                        "ucsc:uc011ecs.3",
                        "orphanet:121654",
                        "ccds:CCDS75523",
                        "omim:603550",
                        "refseq:NM_004100",
                        "ccds:CCDS5165",
                        "ena.embl:Y17114",
                        "vega:OTTHUMG00000015602",
                        "uniprot:O95677"
                    ]
                },
                {
                    "type": "Extension",
                    "name": "previous_symbols",
                    "value": [
                        "CMD1J",
                        "DFNA10"
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
            ],
            "gene_id": "hgnc:3522"
        }
    }
    return params


@pytest.fixture(scope="module")
def genomic_del6_free_text_vrc(genomic_del6_free_text, genomic_del6_seq_loc):
    """Create a test fixture for genomic del relative cnv."""
    genomic_del6_rel_cnv(genomic_del6_free_text, genomic_del6_seq_loc)
    return VariationDescriptor(**genomic_del6_free_text)


@pytest.fixture(scope="module")
def genomic_del6_free_text_vac(genomic_del6_free_text, genomic_del6_seq_loc):
    """Create a test fixture for genomic del absolute cnv."""
    genomic_del6_abs_cnv(genomic_del6_free_text, genomic_del6_seq_loc)
    return VariationDescriptor(**genomic_del6_free_text)


@pytest.fixture(scope="module")
def genomic_del6_free_text_rse_lse(genomic_del6_free_text):
    """Create test fixture for genomic del rse and lse."""
    params = {
        "id": genomic_del6_free_text["id"],
        "type": genomic_del6_free_text["type"],
        "variation": {
            "_id": "ga4gh:VT.a3kXhodtO3tgsdPlEL39Ql4jOuCpOc0s",
            "type": "Text",
            "definition": "EYA4 g.133462764_(133464858_?)del"
        }
    }
    return VariationDescriptor(**params)


@pytest.mark.asyncio
async def assert_text_variation(query_list, test_handler):
    """Make assertion checks for invalid queries"""
    for q in query_list:
        resp = await test_handler.normalize(q, "default",
                                            untranslatable_returns_text=True)
        assert resp.variation_descriptor.label == q.strip()
        assert (resp.variation_descriptor.variation.type == "Text"), q


@pytest.mark.asyncio
async def test_genomic_dup1(test_handler, genomic_dup1_lse,
                            genomic_dup1_vac, genomic_dup1_vrc, genomic_dup1_rse,
                            genomic_dup1_free_text_lse,
                            genomic_dup1_free_text_vac, genomic_dup1_free_text_rse):
    """Test that genomic duplication works correctly."""
    q = "NC_000003.12:g.49531262dup"  # 38
    resp = await test_handler.normalize(q, "default")
    assertion_checks(resp.variation_descriptor, genomic_dup1_lse, q)

    resp = await test_handler.normalize(q, "default")
    assertion_checks(resp.variation_descriptor, genomic_dup1_lse, q)

    resp = await test_handler.normalize(q, "absolute_cnv", baseline_copies=2)
    assertion_checks(resp.variation_descriptor, genomic_dup1_vac, q)

    resp = await test_handler.normalize(q, "relative_cnv",
                                        relative_copy_class="high-level gain")
    assertion_checks(resp.variation_descriptor, genomic_dup1_vrc, q)

    resp = await test_handler.normalize(q, "repeated_seq_expr")
    assertion_checks(resp.variation_descriptor, genomic_dup1_rse, q)

    resp = await test_handler.normalize(q, "literal_seq_expr")
    assertion_checks(resp.variation_descriptor, genomic_dup1_lse, q)

    q = "NC_000003.11:g.49568695dup"  # 37
    resp = await test_handler.normalize(q, "default")
    assertion_checks(resp.variation_descriptor, genomic_dup1_lse, q, ignore_id=True)

    resp = await test_handler.normalize(q, "default")
    assertion_checks(resp.variation_descriptor, genomic_dup1_lse, q, ignore_id=True)

    resp = await test_handler.normalize(q, "absolute_cnv", baseline_copies=2)
    assertion_checks(resp.variation_descriptor, genomic_dup1_vac, q, ignore_id=True)

    resp = await test_handler.normalize(q, "relative_cnv",
                                        relative_copy_class="high-level gain")
    assertion_checks(resp.variation_descriptor, genomic_dup1_vrc, q, ignore_id=True)

    resp = await test_handler.normalize(q, "repeated_seq_expr")
    assertion_checks(resp.variation_descriptor, genomic_dup1_rse, q, ignore_id=True)

    resp = await test_handler.normalize(q, "literal_seq_expr")
    assertion_checks(resp.variation_descriptor, genomic_dup1_lse, q, ignore_id=True)

    # Free Text
    for q in [
        "DAG1 g.49568695dup",  # 37
        "DAG1 g.49531262dup"  # 38
    ]:
        resp = await test_handler.normalize(q, "default")
        assertion_checks(resp.variation_descriptor, genomic_dup1_free_text_lse, q,
                         ignore_id=True)

        resp = await test_handler.normalize(q, "default")
        assertion_checks(resp.variation_descriptor, genomic_dup1_free_text_lse, q,
                         ignore_id=True)

        resp = await test_handler.normalize(q, "absolute_cnv", baseline_copies=2)
        assertion_checks(resp.variation_descriptor, genomic_dup1_free_text_vac, q,
                         ignore_id=True)

        resp = await test_handler.normalize(q, "repeated_seq_expr")
        assertion_checks(resp.variation_descriptor, genomic_dup1_free_text_rse, q,
                         ignore_id=True)

        resp = await test_handler.normalize(q, "literal_seq_expr")
        assertion_checks(resp.variation_descriptor, genomic_dup1_free_text_lse, q,
                         ignore_id=True)

    # Invalid
    invalid_queries = [
        "NC_000007.13:g.159138670dup",
        "NC_000007.14:g.159345976dup",
        "BRAF g.140219337dup", "BRAF g.141024929dup"
    ]
    await assert_text_variation(invalid_queries, test_handler)


@pytest.mark.asyncio
async def test_genomic_dup2(test_handler, genomic_dup2_lse, genomic_dup2_vac,
                            genomic_dup2_vrc, genomic_dup2_rse,
                            genomic_dup2_free_text_default, genomic_dup2_free_text_vac,
                            genomic_dup2_free_text_rse):
    """Test that genomic duplication works correctly."""
    q = "NC_000016.10:g.2087938_2087948dup"  # 38
    resp = await test_handler.normalize(q, "default")
    assertion_checks(resp.variation_descriptor, genomic_dup2_lse, q)

    resp = await test_handler.normalize(q, "absolute_cnv", baseline_copies=2)
    assertion_checks(resp.variation_descriptor, genomic_dup2_vac, q)

    resp = await test_handler.normalize(q, "relative_cnv")
    assertion_checks(resp.variation_descriptor, genomic_dup2_vrc, q)

    resp = await test_handler.normalize(q, "repeated_seq_expr")
    assertion_checks(resp.variation_descriptor, genomic_dup2_rse, q)

    resp = await test_handler.normalize(q, "literal_seq_expr")
    assertion_checks(resp.variation_descriptor, genomic_dup2_lse, q)

    q = "NC_000016.9:g.2137939_2137949dup"  # 37
    resp = await test_handler.normalize(q, "default")
    assertion_checks(resp.variation_descriptor, genomic_dup2_lse, q, ignore_id=True)

    resp = await test_handler.normalize(q, "absolute_cnv", baseline_copies=2)
    assertion_checks(resp.variation_descriptor, genomic_dup2_vac, q, ignore_id=True)

    resp = await test_handler.normalize(q, "relative_cnv")
    assertion_checks(resp.variation_descriptor, genomic_dup2_vrc, q, ignore_id=True)

    resp = await test_handler.normalize(q, "repeated_seq_expr")
    assertion_checks(resp.variation_descriptor, genomic_dup2_rse, q, ignore_id=True)

    resp = await test_handler.normalize(q, "literal_seq_expr")
    assertion_checks(resp.variation_descriptor, genomic_dup2_lse, q, ignore_id=True)

    # Free text
    for q in [
        "DMD g.33229407_33229410dup",  # 37
        "DMD g.33211290_33211293dup"  # 38
    ]:
        resp = await test_handler.normalize(q, "default")
        assertion_checks(resp.variation_descriptor, genomic_dup2_free_text_default, q,
                         ignore_id=True)

        resp = await test_handler.normalize(q, "absolute_cnv", baseline_copies=2)
        assertion_checks(resp.variation_descriptor, genomic_dup2_free_text_vac, q,
                         ignore_id=True)

        resp = await test_handler.normalize(q, "repeated_seq_expr")
        assertion_checks(resp.variation_descriptor, genomic_dup2_free_text_rse, q,
                         ignore_id=True)

        resp = await test_handler.normalize(q, "literal_seq_expr")
        assertion_checks(resp.variation_descriptor, genomic_dup2_free_text_default, q,
                         ignore_id=True)

    # Invalid
    invalid_queries = [
        "NC_000007.13:g.140413127_159138670dup",
        "NC_000007.14:g.140413127_159345976dup",
        "BRAF g.140219337_140924929dup", "BRAF g.140719326_141024929dup"
    ]
    await assert_text_variation(invalid_queries, test_handler)


@pytest.mark.asyncio
async def test_genomic_dup3(test_handler, genomic_dup3_vrc, genomic_dup3_vac,
                            genomic_dup3_rse_lse, genomic_dup3_free_text_vac,
                            genomic_dup3_free_text_vrc,
                            genomic_dup3_free_text_rse_lse):
    """Test that genomic duplication works correctly."""
    q = "NC_000023.11:g.(31060227_31100351)_(33274278_33417151)dup"  # 38
    resp = await test_handler.normalize(q, "default")
    assertion_checks(resp.variation_descriptor, genomic_dup3_vrc, q)

    resp = await test_handler.normalize(q, "absolute_cnv", baseline_copies=1)
    assertion_checks(resp.variation_descriptor, genomic_dup3_vac, q)

    resp = await test_handler.normalize(q, "relative_cnv",
                                        relative_copy_class="low-level gain")
    assertion_checks(resp.variation_descriptor, genomic_dup3_vrc, q)

    resp = await test_handler.normalize(q, "repeated_seq_expr",
                                        untranslatable_returns_text=True)
    assertion_checks(resp.variation_descriptor, genomic_dup3_rse_lse, q)

    resp = await test_handler.normalize(q, "literal_seq_expr",
                                        untranslatable_returns_text=True)
    assertion_checks(resp.variation_descriptor, genomic_dup3_rse_lse, q)

    q = "NC_000023.10:g.(31078344_31118468)_(33292395_33435268)dup"  # 37
    resp = await test_handler.normalize(q, "default")
    assertion_checks(resp.variation_descriptor, genomic_dup3_vrc, q, ignore_id=True)

    resp = await test_handler.normalize(q, "absolute_cnv", baseline_copies=1)
    assertion_checks(resp.variation_descriptor, genomic_dup3_vac, q, ignore_id=True)

    resp = await test_handler.normalize(q, "relative_cnv")
    assertion_checks(resp.variation_descriptor, genomic_dup3_vrc, q, ignore_id=True)

    genomic_dup3_rse_lse.variation.definition = q
    resp = await test_handler.normalize(q, "repeated_seq_expr",
                                        untranslatable_returns_text=True)
    assertion_checks(resp.variation_descriptor, genomic_dup3_rse_lse, q, ignore_id=True)

    resp = await test_handler.normalize(q, "literal_seq_expr",
                                        untranslatable_returns_text=True)
    assertion_checks(resp.variation_descriptor, genomic_dup3_rse_lse, q, ignore_id=True)

    # Free Text
    for q in [
        # TODO:  issue-176
        # "DMD g.(31165391_31165395)_(31200854_31200856)dup",
        "DMD g.(31147274_31147278)_(31182737_31182739)dup"  # 38
    ]:
        resp = await test_handler.normalize(q, "default")
        assertion_checks(resp.variation_descriptor, genomic_dup3_free_text_vrc, q,
                         ignore_id=True)

        resp = await test_handler.normalize(q, "absolute_cnv", baseline_copies=3)
        assertion_checks(resp.variation_descriptor, genomic_dup3_free_text_vac, q,
                         ignore_id=True)

        genomic_dup3_rse_lse.variation.definition = q
        resp = await test_handler.normalize(q, "repeated_seq_expr",
                                            untranslatable_returns_text=True)
        assertion_checks(resp.variation_descriptor, genomic_dup3_free_text_rse_lse, q,
                         ignore_id=True)

        resp = await test_handler.normalize(q, "literal_seq_expr",
                                            untranslatable_returns_text=True)
        assertion_checks(resp.variation_descriptor, genomic_dup3_free_text_rse_lse, q,
                         ignore_id=True)

    # Invalid
    invalid_queries = [
        "NC_000023.10:g.(31119221_31119227)_(31119300_155270562)dup",
        "NC_000023.11:g.(31119221_31119227)_(31119300_156040899)dup",
        "DMD g.(31060227_31100351)_(33274278_33417151)dup"
    ]
    await assert_text_variation(invalid_queries, test_handler)


@pytest.mark.asyncio
async def test_genomic_dup4(test_handler, genomic_dup4_vac, genomic_dup4_vrc,
                            genomic_dup4_rse_lse, genomic_dup4_free_text_vac,
                            genomic_dup4_free_text_vrc,
                            genomic_dup4_free_text_rse_lse):
    """Test that genomic duplication works correctly."""
    q = "NC_000020.11:g.(?_30417576)_(31394018_?)dup"  # 38
    resp = await test_handler.normalize(q, "default")
    assertion_checks(resp.variation_descriptor, genomic_dup4_vrc, q)

    resp = await test_handler.normalize(q, "absolute_cnv", baseline_copies=2)
    assertion_checks(resp.variation_descriptor, genomic_dup4_vac, q)

    resp = await test_handler.normalize(q, "relative_cnv")
    assertion_checks(resp.variation_descriptor, genomic_dup4_vrc, q)

    resp = await test_handler.normalize(q, "repeated_seq_expr",
                                        untranslatable_returns_text=True)
    assertion_checks(resp.variation_descriptor, genomic_dup4_rse_lse, q)

    resp = await test_handler.normalize(q, "literal_seq_expr",
                                        untranslatable_returns_text=True)
    assertion_checks(resp.variation_descriptor, genomic_dup4_rse_lse, q)

    q = "NC_000020.10:g.(?_29652252)_(29981821_?)dup"  # 37
    resp = await test_handler.normalize(q, "default")
    assertion_checks(resp.variation_descriptor, genomic_dup4_vrc, q, ignore_id=True)

    resp = await test_handler.normalize(q, "absolute_cnv", baseline_copies=2)
    assertion_checks(resp.variation_descriptor, genomic_dup4_vac, q, ignore_id=True)

    resp = await test_handler.normalize(q, "relative_cnv")
    assertion_checks(resp.variation_descriptor, genomic_dup4_vrc, q, ignore_id=True)

    genomic_dup4_rse_lse.variation.definition = q
    resp = await test_handler.normalize(q, "repeated_seq_expr",
                                        untranslatable_returns_text=True)
    assertion_checks(resp.variation_descriptor, genomic_dup4_rse_lse, q, ignore_id=True)

    resp = await test_handler.normalize(q, "literal_seq_expr",
                                        untranslatable_returns_text=True)
    assertion_checks(resp.variation_descriptor, genomic_dup4_rse_lse, q, ignore_id=True)

    # Free Text
    for q in [
        "PRPF8 g.(?_1577736)_(1587865_?)dup",  # 37
        "PRPF8 g.(?_1674442)_(1684571_?)dup"  # 38
    ]:
        resp = await test_handler.normalize(q, "default")
        assertion_checks(resp.variation_descriptor, genomic_dup4_free_text_vrc, q,
                         ignore_id=True)

        resp = await test_handler.normalize(q, "absolute_cnv", baseline_copies=2)
        assertion_checks(resp.variation_descriptor, genomic_dup4_free_text_vac, q,
                         ignore_id=True)

        genomic_dup4_rse_lse.variation.definition = q
        resp = await test_handler.normalize(q, "repeated_seq_expr",
                                            untranslatable_returns_text=True)
        genomic_dup4_free_text_rse_lse.variation.definition = q
        assertion_checks(resp.variation_descriptor, genomic_dup4_free_text_rse_lse, q,
                         ignore_id=True)

        resp = await test_handler.normalize(q, "literal_seq_expr",
                                            untranslatable_returns_text=True)
        assertion_checks(resp.variation_descriptor, genomic_dup4_free_text_rse_lse, q,
                         ignore_id=True)

    # Invalid
    invalid_queries = [
        "NC_000020.10:g.(?_29652252)_(63025530_?)dup",
        "NC_000020.11:g.(?_29652252)_(64444169_?)dup",
        "PRPF8 g.(?_1650628)_(1684571_?)dup"
    ]
    await assert_text_variation(invalid_queries, test_handler)


@pytest.mark.asyncio
async def test_genomic_dup5(test_handler, genomic_dup5_vac, genomic_dup5_vrc,
                            genomic_dup5_rse_lse, genomic_dup5_free_text_vac,
                            genomic_dup5_free_text_vrc,
                            genomic_dup5_free_text_rse_lse):
    """Test that genomic duplication works correctly."""
    q = "NC_000023.11:g.(?_154021812)_154092209dup"  # 38
    resp = await test_handler.normalize(q, "default")
    assertion_checks(resp.variation_descriptor, genomic_dup5_vrc, q)

    resp = await test_handler.normalize(q, "absolute_cnv", baseline_copies=2)
    assertion_checks(resp.variation_descriptor, genomic_dup5_vac, q)

    resp = await test_handler.normalize(q, "relative_cnv")
    assertion_checks(resp.variation_descriptor, genomic_dup5_vrc, q)

    resp = await test_handler.normalize(q, "repeated_seq_expr",
                                        untranslatable_returns_text=True)
    assertion_checks(resp.variation_descriptor, genomic_dup5_rse_lse, q)

    resp = await test_handler.normalize(q, "literal_seq_expr",
                                        untranslatable_returns_text=True)
    assertion_checks(resp.variation_descriptor, genomic_dup5_rse_lse, q)

    q = "NC_000023.10:g.(?_153287263)_153357667dup"  # 37
    resp = await test_handler.normalize(q, "default")
    assertion_checks(resp.variation_descriptor, genomic_dup5_vrc, q, ignore_id=True)

    resp = await test_handler.normalize(q, "absolute_cnv", baseline_copies=2)
    assertion_checks(resp.variation_descriptor, genomic_dup5_vac, q, ignore_id=True)

    resp = await test_handler.normalize(q, "relative_cnv")
    assertion_checks(resp.variation_descriptor, genomic_dup5_vrc, q, ignore_id=True)

    genomic_dup5_rse_lse.variation.definition = q
    resp = await test_handler.normalize(q, "repeated_seq_expr",
                                        untranslatable_returns_text=True)
    assertion_checks(resp.variation_descriptor, genomic_dup5_rse_lse, q, ignore_id=True)

    resp = await test_handler.normalize(q, "literal_seq_expr",
                                        untranslatable_returns_text=True)
    assertion_checks(resp.variation_descriptor, genomic_dup5_rse_lse, q, ignore_id=True)

    # Free Text
    for q in [
        "MECP2 g.(?_153287263)_153357667dup",  # 37
        "MECP2 g.(?_154021812)_154092209dup"  # 38
    ]:
        resp = await test_handler.normalize(q, "default")
        assertion_checks(resp.variation_descriptor, genomic_dup5_free_text_vrc, q,
                         ignore_id=True)

        resp = await test_handler.normalize(q, "absolute_cnv", baseline_copies=2)
        assertion_checks(resp.variation_descriptor, genomic_dup5_free_text_vac, q,
                         ignore_id=True)

        genomic_dup5_free_text_rse_lse.variation.definition = q
        resp = await test_handler.normalize(q, "repeated_seq_expr",
                                            untranslatable_returns_text=True)
        assertion_checks(resp.variation_descriptor, genomic_dup5_free_text_rse_lse, q,
                         ignore_id=True)

        resp = await test_handler.normalize(q, "literal_seq_expr",
                                            untranslatable_returns_text=True)
        assertion_checks(resp.variation_descriptor, genomic_dup5_free_text_rse_lse, q,
                         ignore_id=True)

    # Invalid
    for q in [
        "NC_000023.10:g.(?_153287263)_155270561dup",
        "NC_000023.11:g.(?_154021812)_156040896dup",
        "MECP2 g.(?_154021812)_154097733dup"  # 37
        "MECP2 g.(?_154021572)_154092209dup",  # 38
    ]:
        resp = await test_handler.normalize(q, "default",
                                            untranslatable_returns_text=True)
        assert resp.variation_descriptor.variation.type == "Text"


@pytest.mark.asyncio
async def test_genomic_dup6(test_handler, genomic_dup6_vac, genomic_dup6_vrc,
                            genomic_dup6_rse_lse, genomic_dup6_free_text_vac,
                            genomic_dup6_free_text_vrc,
                            genomic_dup6_free_text_rse_lse):
    """Test that genomic duplication works correctly."""
    q = "NC_000023.11:g.154021812_(154092209_?)dup"  # 38
    resp = await test_handler.normalize(q, "default")
    assertion_checks(resp.variation_descriptor, genomic_dup6_vrc, q)

    resp = await test_handler.normalize(q, "absolute_cnv", baseline_copies=1)
    assertion_checks(resp.variation_descriptor, genomic_dup6_vac, q)

    resp = await test_handler.normalize(q, "relative_cnv")
    assertion_checks(resp.variation_descriptor, genomic_dup6_vrc, q)

    resp = await test_handler.normalize(q, "repeated_seq_expr",
                                        untranslatable_returns_text=True)
    assertion_checks(resp.variation_descriptor, genomic_dup6_rse_lse, q)

    resp = await test_handler.normalize(q, "literal_seq_expr",
                                        untranslatable_returns_text=True)
    assertion_checks(resp.variation_descriptor, genomic_dup6_rse_lse, q)

    q = "NC_000023.10:g.153287263_(153357667_?)dup"  # 37
    resp = await test_handler.normalize(q, "default")
    assertion_checks(resp.variation_descriptor, genomic_dup6_vrc, q, ignore_id=True)

    resp = await test_handler.normalize(q, "absolute_cnv", baseline_copies=1)
    assertion_checks(resp.variation_descriptor, genomic_dup6_vac, q, ignore_id=True)

    resp = await test_handler.normalize(q, "relative_cnv")
    assertion_checks(resp.variation_descriptor, genomic_dup6_vrc, q, ignore_id=True)

    genomic_dup6_rse_lse.variation.definition = q
    resp = await test_handler.normalize(q, "repeated_seq_expr",
                                        untranslatable_returns_text=True)
    assertion_checks(resp.variation_descriptor, genomic_dup6_rse_lse, q, ignore_id=True)

    resp = await test_handler.normalize(q, "literal_seq_expr",
                                        untranslatable_returns_text=True)
    assertion_checks(resp.variation_descriptor, genomic_dup6_rse_lse, q, ignore_id=True)

    # Free Text
    for q in [
        "MECP2 g.153287263_(153357667_?)dup",  # 37
        "MECP2 g.154021812_(154092209_?)dup"  # 38
    ]:
        resp = await test_handler.normalize(q, "default")
        assertion_checks(resp.variation_descriptor, genomic_dup6_free_text_vrc, q,
                         ignore_id=True)

        resp = await test_handler.normalize(q, "absolute_cnv", baseline_copies=1)
        assertion_checks(resp.variation_descriptor, genomic_dup6_free_text_vac, q,
                         ignore_id=True)

        genomic_dup6_free_text_rse_lse.variation.definition = q
        resp = await test_handler.normalize(q, "repeated_seq_expr",
                                            untranslatable_returns_text=True)
        assertion_checks(resp.variation_descriptor, genomic_dup6_free_text_rse_lse, q,
                         ignore_id=True)

        resp = await test_handler.normalize(q, "literal_seq_expr",
                                            untranslatable_returns_text=True)
        assertion_checks(resp.variation_descriptor, genomic_dup6_free_text_rse_lse, q,
                         ignore_id=True)

    # Invalid
    for q in [
        "NC_000023.10:g.153287263_(155270561_?)dup",
        "NC_000023.11:g.154021812_(156040896_?)dup",
        "MECP2 g.154021812_(154097733_?)dup"  # 37
        "MECP2 g.154021572_(154092209_?)dup",  # 38
    ]:
        resp = await test_handler.normalize(q, "default",
                                            untranslatable_returns_text=True)
        assert resp.variation_descriptor.variation.type == "Text"


@pytest.mark.asyncio
async def test_genomic_del1(test_handler, genomic_del1_lse, genomic_del1_vac,
                            genomic_del1_vrc, genomic_del1_rse,
                            genomic_del1_free_text_lse, genomic_del1_free_text_vac,
                            genomic_del1_free_text_rse):
    """Test that genomic deletion works correctly."""
    q = "NC_000003.12:g.10149811del"  # 38
    resp = await test_handler.normalize(q, "default")
    assertion_checks(resp.variation_descriptor, genomic_del1_lse, q)

    resp = await test_handler.normalize(q, "absolute_cnv", baseline_copies=2)
    assertion_checks(resp.variation_descriptor, genomic_del1_vac, q)

    resp = await test_handler.normalize(q, "relative_cnv",
                                        relative_copy_class="copy neutral")
    assertion_checks(resp.variation_descriptor, genomic_del1_vrc, q)

    resp = await test_handler.normalize(q, "repeated_seq_expr")
    assertion_checks(resp.variation_descriptor, genomic_del1_rse, q)

    resp = await test_handler.normalize(q, "literal_seq_expr")
    assertion_checks(resp.variation_descriptor, genomic_del1_lse, q)

    q = "NC_000003.11:g.10191495del"  # 37
    resp = await test_handler.normalize(q, "default")
    assertion_checks(resp.variation_descriptor, genomic_del1_lse, q, ignore_id=True)

    resp = await test_handler.normalize(q, "absolute_cnv", baseline_copies=2)
    assertion_checks(resp.variation_descriptor, genomic_del1_vac, q, ignore_id=True)

    resp = await test_handler.normalize(q, "relative_cnv",
                                        relative_copy_class="copy neutral")
    assertion_checks(resp.variation_descriptor, genomic_del1_vrc, q, ignore_id=True)

    resp = await test_handler.normalize(q, "repeated_seq_expr")
    assertion_checks(resp.variation_descriptor, genomic_del1_rse, q, ignore_id=True)

    resp = await test_handler.normalize(q, "literal_seq_expr")
    assertion_checks(resp.variation_descriptor, genomic_del1_lse, q, ignore_id=True)

    # Free text
    for q in [
        "VHL g.10191495del",  # 37
        "VHL g.10149811del"  # 38
    ]:
        resp = await test_handler.normalize(q, "default")
        assertion_checks(resp.variation_descriptor, genomic_del1_free_text_lse, q,
                         ignore_id=True)

        resp = await test_handler.normalize(q, "absolute_cnv", baseline_copies=2)
        assertion_checks(resp.variation_descriptor, genomic_del1_free_text_vac, q,
                         ignore_id=True)

        resp = await test_handler.normalize(q, "repeated_seq_expr")
        assertion_checks(resp.variation_descriptor, genomic_del1_free_text_rse, q,
                         ignore_id=True)

        resp = await test_handler.normalize(q, "literal_seq_expr")
        assertion_checks(resp.variation_descriptor, genomic_del1_free_text_lse, q,
                         ignore_id=True)

    # gnomad vcf
    q = "3-10149810-CT-C"  # 38
    resp = await test_handler.normalize(q)
    assertion_checks(resp.variation_descriptor, genomic_del1_lse, q, ignore_id=True)

    resp = await test_handler.normalize(q, "absolute_cnv")
    assertion_checks(resp.variation_descriptor, genomic_del1_lse, q, ignore_id=True)

    q = "3-10191494-CT-C"  # 37
    resp = await test_handler.normalize(q)
    assertion_checks(resp.variation_descriptor, genomic_del1_lse, q, ignore_id=True)

    resp = await test_handler.normalize(q, "repeated_seq_expr")
    assertion_checks(resp.variation_descriptor, genomic_del1_lse, q, ignore_id=True)

    # Invalid
    invalid_queries = [
        "NC_000003.11:g.198022431del",
        "NC_000003.12:g.198295567del",
        "BRAF g.140413127del", "BRAF g.141024929del"
    ]
    await assert_text_variation(invalid_queries, test_handler)


@pytest.mark.asyncio
async def test_genomic_del2(test_handler, genomic_del2_lse, genomic_del2_vac,
                            genomic_del2_vrc, genomic_del2_rse,
                            genomic_del2_free_text_default, genomic_del2_free_text_cnv,
                            genomic_del2_free_text_rse):
    """Test that genomic deletion works correctly."""
    q = "NC_000003.12:g.10146595_10146613del"  # 38
    resp = await test_handler.normalize(q, "default")
    assertion_checks(resp.variation_descriptor, genomic_del2_lse, q)

    resp = await test_handler.normalize(q, "absolute_cnv", baseline_copies=2)
    assertion_checks(resp.variation_descriptor, genomic_del2_vac, q)

    resp = await test_handler.normalize(q, "relative_cnv",
                                        relative_copy_class="complete loss")
    assertion_checks(resp.variation_descriptor, genomic_del2_vrc, q)

    resp = await test_handler.normalize(q, "repeated_seq_expr")
    assertion_checks(resp.variation_descriptor, genomic_del2_rse, q)

    resp = await test_handler.normalize(q, "literal_seq_expr")
    assertion_checks(resp.variation_descriptor, genomic_del2_lse, q)

    q = "NC_000003.11:g.10188279_10188297del"  # 37
    resp = await test_handler.normalize(q, "default")
    assertion_checks(resp.variation_descriptor, genomic_del2_lse, q, ignore_id=True)

    resp = await test_handler.normalize(q, "absolute_cnv", baseline_copies=2)
    assertion_checks(resp.variation_descriptor, genomic_del2_vac, q, ignore_id=True)

    resp = await test_handler.normalize(q, "relative_cnv",
                                        relative_copy_class="complete loss")
    assertion_checks(resp.variation_descriptor, genomic_del2_vrc, q, ignore_id=True)

    resp = await test_handler.normalize(q, "repeated_seq_expr")
    assertion_checks(resp.variation_descriptor, genomic_del2_rse, q, ignore_id=True)

    resp = await test_handler.normalize(q, "literal_seq_expr")
    assertion_checks(resp.variation_descriptor, genomic_del2_lse, q, ignore_id=True)

    # Free text
    for q in [
        "VHL g.10188279_10188297del",  # 37
        "VHL g.10146595_10146613del"  # 38
    ]:
        resp = await test_handler.normalize(q, "default")
        assertion_checks(resp.variation_descriptor, genomic_del2_free_text_default, q,
                         ignore_id=True)

        resp = await test_handler.normalize(q, "absolute_cnv", baseline_copies=2)
        assertion_checks(resp.variation_descriptor, genomic_del2_free_text_cnv, q,
                         ignore_id=True)

        resp = await test_handler.normalize(q, "repeated_seq_expr")
        assertion_checks(resp.variation_descriptor, genomic_del2_free_text_rse, q,
                         ignore_id=True)

        resp = await test_handler.normalize(q, "literal_seq_expr")
        assertion_checks(resp.variation_descriptor, genomic_del2_free_text_default, q,
                         ignore_id=True)

    # gnomad vcf
    q = "3-10146594-AATGTTGACGGACAGCCTAT-A"
    resp = await test_handler.normalize(q, "default")
    assertion_checks(resp.variation_descriptor, genomic_del2_lse, q, ignore_id=True)

    q = "3-10188278-AATGTTGACGGACAGCCTAT-A"
    resp = await test_handler.normalize(q)
    assertion_checks(resp.variation_descriptor, genomic_del2_lse, q, ignore_id=True)

    # Invalid
    invalid_queries = [
        "NC_000003.12:g.10146595_198295580del",
        "NC_000003.11:g.198022435_198022437del",
        "BRAF g.140413127_140419136del", "BRAF g.140719326_141024929del"
    ]
    await assert_text_variation(invalid_queries, test_handler)


@pytest.mark.asyncio
async def test_genomic_del3(test_handler, genomic_del3_vac, genomic_del3_vrc,
                            genomic_del3_rse_lse, genomic_del3_free_text_vac,
                            genomic_del3_free_text_vrc,
                            genomic_del3_free_text_rse_lse):
    """Test that genomic deletion works correctly."""
    q = "NC_000023.11:g.(31060227_31100351)_(33274278_33417151)del"  # 38
    resp = await test_handler.normalize(q, "default")
    assertion_checks(resp.variation_descriptor, genomic_del3_vrc, q)

    resp = await test_handler.normalize(q, "absolute_cnv", baseline_copies=3)
    assertion_checks(resp.variation_descriptor, genomic_del3_vac, q)

    resp = await test_handler.normalize(q, "relative_cnv")
    assertion_checks(resp.variation_descriptor, genomic_del3_vrc, q)

    resp = await test_handler.normalize(q, "repeated_seq_expr",
                                        untranslatable_returns_text=True)
    assertion_checks(resp.variation_descriptor, genomic_del3_rse_lse, q)

    resp = await test_handler.normalize(q, "literal_seq_expr",
                                        untranslatable_returns_text=True)
    assertion_checks(resp.variation_descriptor, genomic_del3_rse_lse, q)

    q = "NC_000023.10:g.(31078344_31118468)_(33292395_33435268)del"  # 37
    resp = await test_handler.normalize(q, "default")
    assertion_checks(resp.variation_descriptor, genomic_del3_vrc, q, ignore_id=True)

    resp = await test_handler.normalize(q, "absolute_cnv", baseline_copies=3)
    assertion_checks(resp.variation_descriptor, genomic_del3_vac, q, ignore_id=True)

    resp = await test_handler.normalize(q, "relative_cnv")
    assertion_checks(resp.variation_descriptor, genomic_del3_vrc, q, ignore_id=True)

    genomic_del3_rse_lse.variation.definition = q
    resp = await test_handler.normalize(q, "repeated_seq_expr",
                                        untranslatable_returns_text=True)
    assertion_checks(resp.variation_descriptor, genomic_del3_rse_lse, q, ignore_id=True)

    resp = await test_handler.normalize(q, "literal_seq_expr",
                                        untranslatable_returns_text=True)
    assertion_checks(resp.variation_descriptor, genomic_del3_rse_lse, q, ignore_id=True)

    # Free Text
    for q in [
        "EFNB1 g.(68059108_68059111)_(68060963_68060968)del",  # 37
        "EFNB1 g.(68839265_68839268)_(68841120_68841125)del"  # 38
    ]:
        resp = await test_handler.normalize(q, "default")
        assertion_checks(resp.variation_descriptor, genomic_del3_free_text_vrc, q,
                         ignore_id=True)

        resp = await test_handler.normalize(q, "absolute_cnv", baseline_copies=3)
        assertion_checks(resp.variation_descriptor, genomic_del3_free_text_vac, q,
                         ignore_id=True)

        genomic_del3_free_text_rse_lse.variation.definition = q
        resp = await test_handler.normalize(q, "repeated_seq_expr",
                                            untranslatable_returns_text=True)
        assertion_checks(resp.variation_descriptor, genomic_del3_free_text_rse_lse, q,
                         ignore_id=True)

        resp = await test_handler.normalize(q, "literal_seq_expr",
                                            untranslatable_returns_text=True)
        assertion_checks(resp.variation_descriptor, genomic_del3_free_text_rse_lse, q,
                         ignore_id=True)

    # Invalid
    invalid_queries = [
        "NC_000023.11:g.(156040880_156040883)_(156040896_156040899)del",
        "NC_000023.10:g.(155270550_155270555)_(155270560_155270562)del",
        "EFNB1 g.(68048863_68048870)_(68842150_68842152)del",  # 37
        "EFNB1 g.(68829022_68829030)_(68842150_68842161)del"  # 38
    ]
    await assert_text_variation(invalid_queries, test_handler)


@pytest.mark.asyncio
async def test_genomic_del4(test_handler, genomic_del4_vac, genomic_del4_vrc,
                            genomic_del4_rse_lse, genomic_uncertain_del_2,
                            genomic_uncertain_del_y, genomic_del4_free_text_vac,
                            genomic_del4_free_text_rse_lse, genomic_del4_free_text_vrc):
    """Test that genomic deletion works correctly."""
    q = "NC_000023.11:g.(?_31120496)_(33339477_?)del"  # 38
    resp = await test_handler.normalize(q, "default")
    assertion_checks(resp.variation_descriptor, genomic_del4_vrc, q)

    resp = await test_handler.normalize(q, "absolute_cnv", baseline_copies=2)
    assertion_checks(resp.variation_descriptor, genomic_del4_vac, q)

    resp = await test_handler.normalize(q, "relative_cnv")
    assertion_checks(resp.variation_descriptor, genomic_del4_vrc, q)

    resp = await test_handler.normalize(q, "repeated_seq_expr",
                                        untranslatable_returns_text=True)
    assertion_checks(resp.variation_descriptor, genomic_del4_rse_lse, q)

    resp = await test_handler.normalize(q, "literal_seq_expr",
                                        untranslatable_returns_text=True)
    assertion_checks(resp.variation_descriptor, genomic_del4_rse_lse, q)

    q = "NC_000023.10:g.(?_31138613)_(33357594_?)del"  # 37
    resp = await test_handler.normalize(q, "default")
    assertion_checks(resp.variation_descriptor, genomic_del4_vrc, q, ignore_id=True)

    resp = await test_handler.normalize(q, "absolute_cnv", baseline_copies=2)
    assertion_checks(resp.variation_descriptor, genomic_del4_vac, q, ignore_id=True)

    resp = await test_handler.normalize(q, "relative_cnv")
    assertion_checks(resp.variation_descriptor, genomic_del4_vrc, q, ignore_id=True)

    genomic_del4_rse_lse.variation.definition = q
    resp = await test_handler.normalize(q, "repeated_seq_expr",
                                        untranslatable_returns_text=True)
    assertion_checks(resp.variation_descriptor, genomic_del4_rse_lse, q, ignore_id=True)

    resp = await test_handler.normalize(q, "literal_seq_expr",
                                        untranslatable_returns_text=True)
    assertion_checks(resp.variation_descriptor, genomic_del4_rse_lse, q, ignore_id=True)

    q = "NC_000002.12:g.(?_110104900)_(110207160_?)del"
    resp = await test_handler.normalize(q)
    assertion_checks(resp.variation_descriptor, genomic_uncertain_del_2, q)

    q = "NC_000024.10:g.(?_14076802)_(57165209_?)del"
    resp = await test_handler.normalize(q)
    assertion_checks(resp.variation_descriptor, genomic_uncertain_del_y, q)

    # Free Text
    for q in [
        # TODO:  issue-176
        # "COL4A4 g.(?_227886744)_(227890546_?)del",  # 37
        "COL4A4 g.(?_227022028)_(227025830_?)del"  # 38
    ]:
        resp = await test_handler.normalize(q, "default")
        assertion_checks(resp.variation_descriptor, genomic_del4_free_text_vrc, q,
                         ignore_id=True)

        resp = await test_handler.normalize(q, "absolute_cnv", baseline_copies=2)
        assertion_checks(resp.variation_descriptor, genomic_del4_free_text_vac, q,
                         ignore_id=True)

        resp = await test_handler.normalize(q, "repeated_seq_expr",
                                            untranslatable_returns_text=True)
        assertion_checks(resp.variation_descriptor, genomic_del4_free_text_rse_lse, q,
                         ignore_id=True)

        resp = await test_handler.normalize(q, "literal_seq_expr",
                                            untranslatable_returns_text=True)
        assertion_checks(resp.variation_descriptor, genomic_del4_free_text_rse_lse, q,
                         ignore_id=True)

    # Invalid
    invalid_queries = [
        "NC_000023.11:g.(?_156040899)_(156040900_?)del",
        "NC_000024.10:g.(?_155270565)_(155270568_?)del",
        "COL4A4 g.(?_227002710)_(227003710_?)del",
        "COL4A4 g.(?_227867430)_(228029276_?)del",
    ]
    await assert_text_variation(invalid_queries, test_handler)


@pytest.mark.asyncio
async def test_genomic_del5(test_handler, genomic_del5_vac, genomic_del5_vrc,
                            genomic_del5_rse_lse, genomic_del5_free_text_vac,
                            genomic_del5_free_text_vrc,
                            genomic_del5_free_text_rse_lse):
    """Test that genomic deletion works correctly."""
    q = "NC_000023.11:g.(?_18575354)_18653629del"  # 38
    resp = await test_handler.normalize(q, "default")
    assertion_checks(resp.variation_descriptor, genomic_del5_vrc, q)

    resp = await test_handler.normalize(q, "absolute_cnv", baseline_copies=4)
    assertion_checks(resp.variation_descriptor, genomic_del5_vac, q)

    resp = await test_handler.normalize(q, "relative_cnv")
    assertion_checks(resp.variation_descriptor, genomic_del5_vrc, q)

    resp = await test_handler.normalize(q, "repeated_seq_expr",
                                        untranslatable_returns_text=True)
    assertion_checks(resp.variation_descriptor, genomic_del5_rse_lse, q)

    resp = await test_handler.normalize(q, "literal_seq_expr",
                                        untranslatable_returns_text=True)
    assertion_checks(resp.variation_descriptor, genomic_del5_rse_lse, q)

    q = "NC_000023.10:g.(?_18593474)_18671749del"  # 37
    resp = await test_handler.normalize(q, "default")
    assertion_checks(resp.variation_descriptor, genomic_del5_vrc, q, ignore_id=True)

    resp = await test_handler.normalize(q, "absolute_cnv", baseline_copies=4)
    assertion_checks(resp.variation_descriptor, genomic_del5_vac, q, ignore_id=True)

    resp = await test_handler.normalize(q, "relative_cnv")
    assertion_checks(resp.variation_descriptor, genomic_del5_vrc, q, ignore_id=True)

    genomic_del5_rse_lse.variation.definition = q
    resp = await test_handler.normalize(q, "repeated_seq_expr",
                                        untranslatable_returns_text=True)
    assertion_checks(resp.variation_descriptor, genomic_del5_rse_lse, q, ignore_id=True)

    resp = await test_handler.normalize(q, "literal_seq_expr",
                                        untranslatable_returns_text=True)
    assertion_checks(resp.variation_descriptor, genomic_del5_rse_lse, q, ignore_id=True)

    # Free text
    for q in [
        # TODO:  issue-176
        # "CDKL5 g.(?_18593474)_18671749del",
        "CDKL5 g.(?_18575354)_18653629del"
    ]:
        resp = await test_handler.normalize(q, "default")
        assertion_checks(resp.variation_descriptor, genomic_del5_free_text_vrc, q,
                         ignore_id=True)

        resp = await test_handler.normalize(q, "absolute_cnv", baseline_copies=4)
        assertion_checks(resp.variation_descriptor, genomic_del5_free_text_vac, q,
                         ignore_id=True)

        genomic_del5_free_text_rse_lse.variation.definition = q
        resp = await test_handler.normalize(q, "repeated_seq_expr",
                                            untranslatable_returns_text=True)
        assertion_checks(resp.variation_descriptor, genomic_del5_free_text_rse_lse, q,
                         ignore_id=True)

        resp = await test_handler.normalize(q, "literal_seq_expr",
                                            untranslatable_returns_text=True)
        assertion_checks(resp.variation_descriptor, genomic_del5_free_text_rse_lse, q,
                         ignore_id=True)

    # Invalid
    invalid_queries = [
        "NC_000023.10:g.(?_155270550)_155270570del",
        "NC_000023.11:g.(?_18593474)_18671749del"
        "CDKL5  g.(?_18443702)_18671700del",  # 37
        "CDKL5  g.(?_18425585)_18653631del",  # 38
        "CDKL5  g.(?_18425582)_18653500del"  # 38
    ]
    await assert_text_variation(invalid_queries, test_handler)


@pytest.mark.asyncio
async def test_genomic_del6(test_handler, genomic_del6_vac, genomic_del6_vrc,
                            genomic_del6_rse_lse, genomic_del6_free_text_vac,
                            genomic_del6_free_text_vrc,
                            genomic_del6_free_text_rse_lse):
    """Test that genomic deletion works correctly."""
    q = "NC_000006.12:g.133462764_(133464858_?)del"  # 38
    resp = await test_handler.normalize(q, "default")
    assertion_checks(resp.variation_descriptor, genomic_del6_vrc, q)

    resp = await test_handler.normalize(q, "absolute_cnv", baseline_copies=2)
    assertion_checks(resp.variation_descriptor, genomic_del6_vac, q)

    resp = await test_handler.normalize(q, "relative_cnv")
    assertion_checks(resp.variation_descriptor, genomic_del6_vrc, q)

    resp = await test_handler.normalize(q, "repeated_seq_expr",
                                        untranslatable_returns_text=True)
    assertion_checks(resp.variation_descriptor, genomic_del6_rse_lse, q)

    resp = await test_handler.normalize(q, "literal_seq_expr",
                                        untranslatable_returns_text=True)
    assertion_checks(resp.variation_descriptor, genomic_del6_rse_lse, q)

    q = "NC_000006.11:g.133783902_(133785996_?)del"  # 37
    resp = await test_handler.normalize(q, "default")
    assertion_checks(resp.variation_descriptor, genomic_del6_vrc, q, ignore_id=True)

    resp = await test_handler.normalize(q, "absolute_cnv", baseline_copies=2)
    assertion_checks(resp.variation_descriptor, genomic_del6_vac, q, ignore_id=True)

    resp = await test_handler.normalize(q, "relative_cnv")
    assertion_checks(resp.variation_descriptor, genomic_del6_vrc, q, ignore_id=True)

    genomic_del6_rse_lse.variation.definition = q
    resp = await test_handler.normalize(q, "repeated_seq_expr",
                                        untranslatable_returns_text=True)
    assertion_checks(resp.variation_descriptor, genomic_del6_rse_lse, q, ignore_id=True)

    resp = await test_handler.normalize(q, "literal_seq_expr",
                                        untranslatable_returns_text=True)
    assertion_checks(resp.variation_descriptor, genomic_del6_rse_lse, q, ignore_id=True)

    # Free text
    for q in [
        # TODO:  issue-176
        # "EYA4 g.133783902_(133785996_?)del",  # 37
        "EYA4 g.133462764_(133464858_?)del"  # 38
    ]:
        resp = await test_handler.normalize(q, "default")
        assertion_checks(resp.variation_descriptor, genomic_del6_free_text_vrc, q,
                         ignore_id=True)

        resp = await test_handler.normalize(q, "absolute_cnv", baseline_copies=2)
        assertion_checks(resp.variation_descriptor, genomic_del6_free_text_vac, q,
                         ignore_id=True)

        genomic_del6_rse_lse.variation.definition = q
        resp = await test_handler.normalize(q, "repeated_seq_expr",
                                            untranslatable_returns_text=True)
        assertion_checks(resp.variation_descriptor, genomic_del6_free_text_rse_lse, q,
                         ignore_id=True)

        resp = await test_handler.normalize(q, "literal_seq_expr",
                                            untranslatable_returns_text=True)
        assertion_checks(resp.variation_descriptor, genomic_del6_free_text_rse_lse, q,
                         ignore_id=True)

    # Invalid
    invalid_queries = [
        "NC_000006.11:g.171115069_(171115080_?)del",
        "NC_000006.12:g.170805981_(170805989_?)del"
        "EYA4 g.133561700_(133853270_?)del",  # 37
        "EYA4 g.133561651_(133561708_?)del",  # 37
        "EYA4 g.133240513_(133240600_?)del",  # 38
        "EYA4 g.133240515_(133532130_?)del"  # 38
    ]
    await assert_text_variation(invalid_queries, test_handler)


@pytest.mark.asyncio
async def test_parameters(test_handler):
    """Check that valid and invalid parameters work as intended."""
    resp = await test_handler.normalize("7-140453136-A-T")
    assert resp.variation_descriptor
    assert resp.warnings == []

    q = "NC_000003.12:g.49531262dup"
    resp = await test_handler.normalize(q, "")
    assert resp.variation_descriptor
    assert resp.warnings == []

    resp = await test_handler.normalize(q, None)
    assert resp.variation_descriptor
    assert resp.warnings == []

    resp = await test_handler.normalize(q, " absolute_CnV ", baseline_copies=2)
    assert resp.variation_descriptor
    assert resp.warnings == []

    resp = await test_handler.normalize(q, " absolute_CnV ")
    assert resp.variation_descriptor is None
    assert resp.warnings == ["absolute_cnv mode requires `baseline_copies`"]

    resp = await test_handler.normalize(q, " relative_CnV ")
    assert resp.variation_descriptor
    assert resp.warnings == []
