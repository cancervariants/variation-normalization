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
                "name": "chromosome_location",
                "value": {
                    "species_id": "taxonomy:9606",
                    "start": "p21.2",
                    "end": "p21.1",
                    "id": "ga4gh:CL.7D3ErxYNxYRlDyaNuBIV9EKeo--IT2iS",
                    "type": "ChromosomeLocation",
                    "chr": "X"
                }
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
        "gene": "hgnc:2928"
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
                "name": "chromosome_location",
                "value": {
                    "species_id": "taxonomy:9606",
                    "start": "q28",
                    "end": "q28",
                    "id": "ga4gh:CL.p5Va-YpCTrSTYWyJrpR-rvnxO1YWPIDY",
                    "type": "ChromosomeLocation",
                    "chr": "X"
                }
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
        "gene": "hgnc:6990"
    }


@pytest.fixture(scope="module")
def genomic_dup1():
    """Create test fixture containing params for genomic dup VD."""
    params = {
        "id": "normalize.variation:NC_000003.12%3Ag.49531262dup",
        "type": "VariationDescriptor",
        "variation": dict(),
        "molecule_context": "genomic",
        "structural_type": "SO:1000035",
        "vrs_ref_allele_seq": "GG"
    }
    return params


@pytest.fixture(scope="module")
def genomic_dup1_lse(genomic_dup1, genomic_dup1_seq_loc):
    """Create a test fixture for genomic dup LSE."""
    _id = "ga4gh:VA.J4tgSfdlqvfFtIFW2QY_ux7RKzFco2pd"
    genomic_dup1["variation"] = {
        "type": "Allele",
        "id": _id,
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
    return VariationDescriptor(**genomic_dup1)


@pytest.fixture(scope="module")
def genomic_dup1_vrc(genomic_dup1, genomic_dup1_seq_loc):
    """Create a test fixture for genomic dup relative CNV."""
    genomic_dup1["variation"] = {
        "type": "RelativeCopyNumber",
        "id": "ga4gh:RCN.qSnPCsf9ylOaecpSBjTkxWFqxmmyYWtL",
        "location": genomic_dup1_seq_loc,
        "relative_copy_class": "high-level gain"
    }
    return VariationDescriptor(**genomic_dup1)


@pytest.fixture(scope="module")
def genomic_dup1_rse(genomic_dup1, genomic_dup1_seq_loc):
    """Create a test fixture for genomic dup RSE."""
    _id = "ga4gh:VA.C32Z4YjYHzsEZwzbW85O-jz_CBWl1Blu"
    genomic_dup1["variation"] = {
        "type": "Allele",
        "id": _id,
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
    return VariationDescriptor(**genomic_dup1)


@pytest.fixture(scope="module")
def genomic_dup1_free_text():
    """Create test fixture containing params for genomic dup VD."""
    params = {
        "id": "normalize.variation:DAG1%20g.49568695dup",
        "type": "VariationDescriptor",
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
                    "name": "chromosome_location",
                    "value": {
                        "species_id": "taxonomy:9606",
                        "start": "p21.31",
                        "end": "p21.31",
                        "id": "ga4gh:CL.soO7dXM3WndwwT134GGncAwGwcwBqpnS",
                        "type": "ChromosomeLocation",
                        "chr": "3"
                    }
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
            "gene": "hgnc:2666"
        }
    }
    return params


@pytest.fixture(scope="module")
def genomic_dup1_free_text_seq_loc():
    """Create genomic dup1 free text sequence location"""
    return {
        "id": "ga4gh:SL.CnJQs8qYCvlsG8NM2XFEU5IjvjQaCVx1",
        "sequence_id": "ga4gh:SQ.tpvbnWsfEGqip8gJQZnWJAF8-bWDUDKd",
        "start": {"value": 1032, "type": "Number"},
        "end": {"value": 1034, "type": "Number"},
        "type": "SequenceLocation"
    }


@pytest.fixture(scope="module")
def genomic_dup1_free_text_lse(genomic_dup1_free_text,
                               genomic_dup1_free_text_seq_loc):
    """Create a test fixture for genomic dup LSE."""
    _id = "ga4gh:VA.qkuCWe2_OQeWHXjHUa8NQpDitjhYpNhm"
    genomic_dup1_free_text["variation"] = {
        "type": "Allele",
        "id": _id,
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
    _id = "ga4gh:ACN.TpprCQKgSaGlbBdbAAv4RoOS5yvE4jf2"
    genomic_dup1_free_text["variation"] = {
        "type": "AbsoluteCopyNumber",
        "id": _id,
        "location": genomic_dup1_free_text_seq_loc,
        "copies": {"type": "Number", "value": 3}
    }
    return VariationDescriptor(**genomic_dup1_free_text)


@pytest.fixture(scope="module")
def genomic_dup1_free_text_rse(genomic_dup1_free_text,
                               genomic_dup1_free_text_seq_loc):
    """Create a test fixture for genomic dup RSE."""
    _id = "ga4gh:VA.sY_rNabMagcE7AQW2gAhOII3kvArGjAx"
    genomic_dup1_free_text["variation"] = {
        "type": "Allele",
        "id": _id,
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
    return VariationDescriptor(**genomic_dup1_free_text)


@pytest.fixture(scope="module")
def genomic_dup2():
    """Create test fixture containing params for genomic dup VD."""
    params = {
        "id": "normalize.variation:NC_000016.10%3Ag.2087938_2087948dup",
        "type": "VariationDescriptor",
        "variation": dict(),
        "molecule_context": "genomic",
        "structural_type": "SO:1000035",
        "vrs_ref_allele_seq": "AAAGGTAGGGC"
    }
    return params


@pytest.fixture(scope="module")
def genomic_dup2_lse(genomic_dup2, genomic_dup2_seq_loc):
    """Create a test fixture for genomic dup LSE."""
    _id = "ga4gh:VA.9veMbNiqeGwZKcYUlMrJsQ4vknTDmNrG"
    genomic_dup2["variation"] = {
        "type": "Allele",
        "id": _id,
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
    return VariationDescriptor(**genomic_dup2)


@pytest.fixture(scope="module")
def genomic_dup2_vrc(genomic_dup2, genomic_dup2_seq_loc):
    """Create a test fixture for genomic dup relative CNV."""
    genomic_dup2["variation"] = {
        "type": "RelativeCopyNumber",
        "id": "ga4gh:RCN.zizBAY460E_wbra9Y6oCxojd8YCbIfnx",
        "location": genomic_dup2_seq_loc,
        "relative_copy_class": "low-level gain"
    }
    return VariationDescriptor(**genomic_dup2)


@pytest.fixture(scope="module")
def genomic_dup2_rse(genomic_dup2, genomic_dup2_seq_loc):
    """Create a test fixture for genomic dup RSE."""
    _id = "ga4gh:VA.OSONaDOtYL2WSJp6FuXzL7mXji_mJMZL"
    genomic_dup2["variation"] = {
        "type": "Allele",
        "id": _id,
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
    return VariationDescriptor(**genomic_dup2)


@pytest.fixture(scope="module")
def genomic_dup2_free_text(dmd_gene_context):
    """Create test fixture containing params for genomic dup VD."""
    params = {
        "id": "normalize.variation:TSC2%20g.2137939_2137949dup",
        "type": "VariationDescriptor",
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
        "id": "ga4gh:SL.Ij3eMZP3euBYYAAcxBgayh0miq6wGvcN",
        "sequence_id": "ga4gh:SQ.1DeZLYHMnd-smp3GDlpRxETb9_0AokO7",
        "start": {"value": 256, "type": "Number"},
        "end": {"value": 260, "type": "Number"},
        "type": "SequenceLocation"
    }


@pytest.fixture(scope="module")
def genomic_dup2_free_text_default(genomic_dup2_free_text,
                                   genomic_dup2_free_text_seq_loc):
    """Create a test fixture for genomic dup default and LSE."""
    _id = "ga4gh:VA.ofKjBdAtvrfIAng0zLYiScDbfzIMfZSm"
    genomic_dup2_free_text["variation"] = {
        "type": "Allele",
        "id": _id,
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
    _id = "ga4gh:ACN.McPtQ3LSldeGzC6-qMEGJXIjPLGNEyky"
    genomic_dup2_free_text["variation"] = {
        "type": "AbsoluteCopyNumber",
        "id": _id,
        "location": genomic_dup2_free_text_seq_loc,
        "copies": {"type": "Number", "value": 3}
    }
    return VariationDescriptor(**genomic_dup2_free_text)


@pytest.fixture(scope="module")
def genomic_dup2_free_text_rse(genomic_dup2_free_text,
                               genomic_dup2_free_text_seq_loc):
    """Create a test fixture for genomic dup RSE."""
    _id = "ga4gh:VA.Ys6Td_en3Nc72EMtr6fUEGp2WxDI-CNs"
    genomic_dup2_free_text["variation"] = {
        "type": "Allele",
        "id": _id,
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
    return VariationDescriptor(**genomic_dup2_free_text)


@pytest.fixture(scope="module")
def genomic_dup3():
    """Create test fixture containing params for genomic dup VD."""
    params = {
        "id": "normalize.variation:NC_000023.11%3Ag.%2831060227_31100351%29_%2833274278_33417151%29dup",  # noqa: E501
        "type": "VariationDescriptor",
        "variation": dict(),
        "molecule_context": "genomic",
        "structural_type": "SO:0001742",
        "vrs_ref_allele_seq": None
    }
    return params


@pytest.fixture(scope="module")
def genomic_dup3_vac(genomic_dup3, genomic_del3_dup3_loc):
    """Create a test fixture for genomic dup absolute cnv."""
    _id = "ga4gh:ACN.2ZUQcccwvtoGZ5LZZRUoDZp6218Y6sQK"
    genomic_dup3["variation"] = {
        "type": "AbsoluteCopyNumber",
        "id": _id,
        "location": genomic_del3_dup3_loc,
        "copies": {"type": "Number", "value": 2}
    }
    return VariationDescriptor(**genomic_dup3)


@pytest.fixture(scope="module")
def genomic_dup3_vrc(genomic_dup3, genomic_del3_dup3_loc):
    """Create a test fixture for genomic dup relative cnv."""
    _id = "ga4gh:RCN.4m1TD2538i7v4NQ_OeJ-pIEhlRpYZ3_y"
    genomic_dup3["variation"] = {
        "type": "RelativeCopyNumber",
        "id": _id,
        "location": genomic_del3_dup3_loc,
        "relative_copy_class": "low-level gain"
    }
    return VariationDescriptor(**genomic_dup3)


@pytest.fixture(scope="module")
def genomic_dup3_rse_lse(genomic_dup3):
    """Create test fixture for genomic dup rse and lse."""
    params = {
        "id": genomic_dup3["id"],
        "type": genomic_dup3["type"],
        "variation": {
            "id": "ga4gh:VT.15sKDgSyoCPOgfrFHvSea-fHVeu7huVT",
            "type": "Text",
            "definition": "NC_000023.11:g.(31060227_31100351)_(33274278_33417151)dup"
        }
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_dup3_free_text(dmd_gene_context):
    """Create test fixture containing params for genomic dup VD."""
    params = {
        "id": "normalize.variation:DMD%20g.%2831147274_31147278%29_%2831182737_31182739%29dup",  # noqa: E501
        "type": "VariationDescriptor",
        "variation": dict(),
        "molecule_context": "genomic",
        "structural_type": "SO:0001742",
        "vrs_ref_allele_seq": None,
        "gene_context": dmd_gene_context
    }
    return params


@pytest.fixture(scope="module")
def genomic_dup3_free_text_subject():
    """Create test fixture for genomic dup3 free text location"""
    return {
        "id": "ga4gh:SL.tdjgxR-PHAmEUyQJHQSPo1CI0Vgn19Gg",
        "sequence_id": "ga4gh:SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
        "start": {"min": 31147273, "max": 31147277, "type": "DefiniteRange"},
        "end": {"min": 31182738, "max": 31182740, "type": "DefiniteRange"},
        "type": "SequenceLocation"
    }


@pytest.fixture(scope="module")
def genomic_dup3_free_text_vrc(genomic_dup3_free_text, genomic_dup3_free_text_subject):
    """Create a test fixture for genomic dup relative cnv."""
    _id = "ga4gh:RCN.ReWfNwAnchjIPJQEXM038T9M3OsOO7yK"
    genomic_dup3_free_text["variation"] = {
        "type": "RelativeCopyNumber",
        "id": _id,
        "location": genomic_dup3_free_text_subject,
        "relative_copy_class": "low-level gain"
    }
    return VariationDescriptor(**genomic_dup3_free_text)


@pytest.fixture(scope="module")
def genomic_dup3_free_text_vac(genomic_dup3_free_text, genomic_dup3_free_text_subject):
    """Create a test fixture for genomic dup absolute cnv."""
    _id = "ga4gh:ACN.5B1zY9uirjKRQ5gT_6C0Z0FhdJHwG_yS"
    genomic_dup3_free_text["variation"] = {
        "type": "AbsoluteCopyNumber",
        "id": _id,
        "location": genomic_dup3_free_text_subject,
        "copies": {"type": "Number", "value": 4}
    }
    return VariationDescriptor(**genomic_dup3_free_text)


@pytest.fixture(scope="module")
def genomic_dup3_free_text_rse_lse(genomic_dup3_free_text):
    """Create test fixture for genomic dup rse and lse."""
    params = {
        "id": genomic_dup3_free_text["id"],
        "type": genomic_dup3_free_text["type"],
        "variation": {
            "id": "ga4gh:VT.F0AX-RkMN4U8KLkIE68ECGU83Y-ICWXh",
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
        "variation": dict(),
        "molecule_context": "genomic",
        "structural_type": "SO:0001742",
        "vrs_ref_allele_seq": None
    }
    return params


@pytest.fixture(scope="module")
def genomic_dup4_vrc(genomic_dup4, genoimc_dup4_loc):
    """Create a test fixture for genomic dup relative cnv."""
    _id = "ga4gh:RCN.uSGvLYzpzivqDzhuKR44DHc5imZJSmoV"
    genomic_dup4["variation"] = {
        "type": "RelativeCopyNumber",
        "id": _id,
        "location": genoimc_dup4_loc,
        "relative_copy_class": "low-level gain"
    }
    return VariationDescriptor(**genomic_dup4)


@pytest.fixture(scope="module")
def genomic_dup4_vac(genomic_dup4, genoimc_dup4_loc):
    """Create a test fixture for genomic dup absolute cnv."""
    _id = "ga4gh:ACN.OBSQ6SY9waWJ7YOyYA_qK0te_mTfOT4A"
    genomic_dup4["variation"] = {
        "type": "AbsoluteCopyNumber",
        "id": _id,
        "location": genoimc_dup4_loc,
        "copies": {"type": "Number", "value": 3}
    }
    return VariationDescriptor(**genomic_dup4)


@pytest.fixture(scope="module")
def genomic_dup4_rse_lse(genomic_dup4):
    """Create test fixture for genomic dup rse and lse."""
    params = {
        "id": genomic_dup4["id"],
        "type": genomic_dup4["type"],
        "variation": {
            "id": "ga4gh:VT.Pga4IH82qga2iZAodjxYw9OXhB4Xa2g8",
            "type": "Text",
            "definition": "NC_000020.11:g.(?_30417576)_(31394018_?)dup"
        }
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_dup4_free_text():
    """Create test fixture containing params for genomic dup VD."""
    params = {
        "id": "normalize.variation:PRF8%20g.%28%3F_1577736%29_%281587865_%3F%29",
        "type": "VariationDescriptor",
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
                    "name": "chromosome_location",
                    "value": {
                        "species_id": "taxonomy:9606",
                        "start": "p13.3",
                        "end": "p13.3",
                        "id": "ga4gh:CL.cgXboJas91hqg8zGKLMZYN0RdB_RCj6k",
                        "type": "ChromosomeLocation",
                        "chr": "17"
                    }
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
            "gene": "hgnc:17340"
        }
    }
    return params


@pytest.fixture(scope="module")
def genomic_dup4_free_text_subject():
    """Create test fixture for genomic dup4 free text location"""
    return {
        "id": "ga4gh:SL.C177kUIMvRZlQySAXkDN9K4pg5HWqbo0",
        "sequence_id": "ga4gh:SQ.dLZ15tNO1Ur0IcGjwc3Sdi_0A6Yf4zm7",
        "start": {"value": 1674441, "comparator": "<=", "type": "IndefiniteRange"},
        "end": {"value": 1684571, "comparator": ">=", "type": "IndefiniteRange"},
        "type": "SequenceLocation"
    }


@pytest.fixture(scope="module")
def genomic_dup4_free_text_vrc(genomic_dup4_free_text, genomic_dup4_free_text_subject):
    """Create a test fixture for genomic dup relative cnv."""
    _id = "ga4gh:RCN.Cf9r8JDpUC18VkkEv44jf8b8WMqUIRFu"
    genomic_dup4_free_text["variation"] = {
        "type": "RelativeCopyNumber",
        "id": _id,
        "location": genomic_dup4_free_text_subject,
        "relative_copy_class": "low-level gain"
    }
    return VariationDescriptor(**genomic_dup4_free_text)


@pytest.fixture(scope="module")
def genomic_dup4_free_text_vac(genomic_dup4_free_text, genomic_dup4_free_text_subject):
    """Create a test fixture for genomic dup absolute cnv."""
    _id = "ga4gh:ACN.g4-wQni3sHHQNArM4DaIYjh2YuZFKkiv"
    genomic_dup4_free_text["variation"] = {
        "type": "AbsoluteCopyNumber",
        "id": _id,
        "location": genomic_dup4_free_text_subject,
        "copies": {"type": "Number", "value": 3}
    }
    return VariationDescriptor(**genomic_dup4_free_text)


@pytest.fixture(scope="module")
def genomic_dup4_free_text_rse_lse(genomic_dup4_free_text):
    """Create test fixture for genomic dup rse and lse."""
    params = {
        "id": genomic_dup4_free_text["id"],
        "type": genomic_dup4_free_text["type"],
        "variation": {
            "id": "ga4gh:VT.Pga4IH82qga2iZAodjxYw9OXhB4Xa2g8",
            "type": "Text",
            "definition": "NC_000020.11:g.(?_30417576)_(31394018_?)dup"
        }
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_dup5():
    """Create test fixture containing params for genomic dup VD."""
    params = {
        "id": "normalize.variation:NC_000023.11%3Ag.%28%3F_154021812%29_154092209dup",
        "type": "VariationDescriptor",
        "variation": dict(),
        "molecule_context": "genomic",
        "structural_type": "SO:0001742",
        "vrs_ref_allele_seq": None
    }
    return params


def genomic_dup5_abs_cnv(params, genomic_dup5_loc):
    """Create genomic dup5 aboluste cnv"""
    _id = "ga4gh:ACN.BoEis2HKjiA0yZKs2yvvX-xOSkLxKHlC"
    params["variation"] = {
        "type": "AbsoluteCopyNumber",
        "id": _id,
        "location": genomic_dup5_loc,
        "copies": {"type": "Number", "value": 3}
    }


def genomic_dup5_rel_cnv(params, genomic_dup5_loc):
    """Create genomic dup4 relative cnv"""
    _id = "ga4gh:RCN.tr_brFSOfykm3I3ufLQTS9pV8KKjqLhK"
    params["variation"] = {
        "type": "RelativeCopyNumber",
        "id": _id,
        "location": genomic_dup5_loc,
        "relative_copy_class": "low-level gain"
    }


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
            "id": "ga4gh:VT.of16BEeHyU9od62SrjSCQ4LyUtbbGoKi",
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
            "id": "ga4gh:VT.Kw18bSFpQp9xdKg88fqW7zUx4_VXFIiW",
            "type": "Text",
            "definition": "MECP2 g.(?_154021812)_154092209dup"
        }
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_dup6():
    """Create test fixture containing params for genomic dup VD."""
    params = {
        "id": "normalize.variation:NC_000023.11%3Ag.154021812_%28154092209_%3F%29dup",
        "type": "VariationDescriptor",
        "variation": dict(),
        "molecule_context": "genomic",
        "structural_type": "SO:0001742",
        "vrs_ref_allele_seq": None
    }
    return params


def genomic_dup6_rel_cnv(params, genoimc_dup6_loc):
    """Create genomic dup6 relative cnv"""
    _id = "ga4gh:RCN.c5Uq3TFDNpQSyDlbXb0BRonw9AYHHk0H"
    params["variation"] = {
        "type": "RelativeCopyNumber",
        "id": _id,
        "location": genoimc_dup6_loc,
        "relative_copy_class": "low-level gain"
    }


def genomic_dup6_abs_cnv(params, genoimc_dup6_loc):
    """Create genomic dup6 absolute cnv"""
    _id = "ga4gh:ACN.6q8D5d_ie9MO0HdNEJmRJmGMg8C5LdAM"
    params["variation"] = {
        "type": "AbsoluteCopyNumber",
        "id": _id,
        "location": genoimc_dup6_loc,
        "copies": {"type": "Number", "value": 2}
    }


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
            "id": "ga4gh:VT.2k5AWTbGJxvLVT6bUW0pUMq6XGAcEjXW",
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
            "id": "ga4gh:VT.LbAqiLmJs1t9-FgEKD0-KDKwzvM3AAlz",
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
        "variation": dict(),
        "molecule_context": "genomic",
        "structural_type": "SO:0000159",
        "vrs_ref_allele_seq": "T"
    }
    return params


@pytest.fixture(scope="module")
def genomic_del1_lse(genomic_del1, genomic_del1_seq_loc):
    """Create a test fixture for genomic del LSE."""
    _id = "ga4gh:VA.FVRzUwKTV-A-8zvxPUyREBR9mCunjIPO"
    genomic_del1["variation"] = {
        "type": "Allele",
        "id": _id,
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
    return VariationDescriptor(**genomic_del1)


@pytest.fixture(scope="module")
def genomic_del1_vrc(genomic_del1, genomic_del1_seq_loc):
    """Create a test fixture for genomic del relative CNV."""
    genomic_del1["variation"] = {
        "type": "RelativeCopyNumber",
        "id": "ga4gh:RCN.z7MU8QUSR_aeWG7MP161H4jwPGoyo1No",
        "location": genomic_del1_seq_loc,
        "relative_copy_class": "copy neutral"
    }
    return VariationDescriptor(**genomic_del1)


@pytest.fixture(scope="module")
def genomic_del1_rse(genomic_del1, genomic_del1_seq_loc):
    """Create a test fixture for genomic del RSE."""
    _id = "ga4gh:VA.Mubi52bBVfOHkfemgLXVg1vtl6WLfyxe"
    genomic_del1["variation"] = {
        "type": "Allele",
        "id": _id,
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
    return VariationDescriptor(**genomic_del1)


@pytest.fixture(scope="module")
def genomic_del1_free_text(vhl_gene_context):
    """Create test fixture containing params for genomic del VD."""
    params = {
        "id": "normalize.variation:VHL%20g.10191495del",
        "type": "VariationDescriptor",
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
        "id": "ga4gh:SL.C8wsPU7c4uq-YG88CXZzEldP888a1FMm",
        "sequence_id": "ga4gh:SQ.xBKOKptLLDr-k4hTyCetvARn16pDS_rW",
        "start": {"value": 557, "type": "Number"},
        "end": {"value": 558, "type": "Number"},
        "type": "SequenceLocation"
    }


@pytest.fixture(scope="module")
def genomic_del1_free_text_lse(genomic_del1_free_text,
                               genomic_del1_free_text_seq_loc):
    """Create a test fixture for genomic del LSE."""
    _id = "ga4gh:VA.VbssoYULETSOoFI0FwjGypt2YinbTEOc"
    genomic_del1_free_text["variation"] = {
        "type": "Allele",
        "id": _id,
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
    _id = "ga4gh:ACN.3dlkcT_WnurdradNYgs8gHW_pyP-FLYA"
    genomic_del1_free_text["variation"] = {
        "type": "AbsoluteCopyNumber",
        "id": _id,
        "location": genomic_del1_free_text_seq_loc,
        "copies": {"type": "Number", "value": 1}
    }
    return VariationDescriptor(**genomic_del1_free_text)


@pytest.fixture(scope="module")
def genomic_del1_free_text_rse(genomic_del1_free_text,
                               genomic_del1_free_text_seq_loc):
    """Create a test fixture for genomic del RSE."""
    _id = "ga4gh:VA.FMSY98r0s6L1GF5YIELVx0bmmpFmJ4n-"
    genomic_del1_free_text["variation"] = {
        "type": "Allele",
        "id": _id,
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
    return VariationDescriptor(**genomic_del1_free_text)


@pytest.fixture(scope="module")
def genomic_del2():
    """Create test fixture containing params for genomic del VD."""
    params = {
        "id": "normalize.variation:NC_000003.12%3Ag.10146595_10146613del",
        "type": "VariationDescriptor",
        "variation": dict(),
        "molecule_context": "genomic",
        "structural_type": "SO:0000159",
        "vrs_ref_allele_seq": "ATGTTGACGGACAGCCTAT"
    }
    return params


@pytest.fixture(scope="module")
def genomic_del2_lse(genomic_del2, genomic_del2_seq_loc):
    """Create a test fixture for genomic del LSE."""
    _id = "ga4gh:VA.UgJSDSWAaJFwhRm56LA0Gez47_PYqv0k"
    genomic_del2["variation"] = {
        "type": "Allele",
        "id": _id,
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
    return VariationDescriptor(**genomic_del2)


@pytest.fixture(scope="module")
def genomic_del2_vrc(genomic_del2, genomic_del2_seq_loc):
    """Create a test fixture for genomic del relative CNV."""
    genomic_del2["variation"] = {
        "type": "RelativeCopyNumber",
        "id": "ga4gh:RCN._19twDngCtP3U-8ED2Kly2HM53I_7CV7",
        "location": genomic_del2_seq_loc,
        "relative_copy_class": "complete loss"
    }
    return VariationDescriptor(**genomic_del2)


@pytest.fixture(scope="module")
def genomic_del2_rse(genomic_del2, genomic_del2_seq_loc):
    """Create a test fixture for genomic del RSE."""
    _id = "ga4gh:VA._QLDzH5kHqog6-RaKOLg36EmEtWP_qE7"
    genomic_del2["variation"] = {
        "type": "Allele",
        "id": _id,
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
    return VariationDescriptor(**genomic_del2)


@pytest.fixture(scope="module")
def genomic_del2_free_text(vhl_gene_context):
    """Create test fixture containing params for genomic del VD."""
    params = {
        "id": "normalize.variation:VHL%20g.10188279_10188297del",
        "type": "VariationDescriptor",
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
        "id": "ga4gh:SL.0ietJcxUtWPRKr_9XtqQoH3cN4XfCrDM",
        "sequence_id": "ga4gh:SQ.xBKOKptLLDr-k4hTyCetvARn16pDS_rW",
        "start": {"value": 491, "type": "Number"},
        "end": {"value": 510, "type": "Number"},
        "type": "SequenceLocation"
    }


@pytest.fixture(scope="module")
def genomic_del2_free_text_default(genomic_del2_free_text,
                                   genomic_del2_free_text_seq_loc):
    """Create a test fixture for genomic del default and LSE."""
    _id = "ga4gh:VA.6wTYBh0btGq6SlXDu4V7iEK9UrehXS-6"
    genomic_del2_free_text["variation"] = {
        "type": "Allele",
        "id": _id,
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
    _id = "ga4gh:ACN.frpOM82FmXnicV0mpJ7R1oGLzJsGjj8F"
    genomic_del2_free_text["variation"] = {
        "type": "AbsoluteCopyNumber",
        "id": _id,
        "location": genomic_del2_free_text_seq_loc,
        "copies": {"type": "Number", "value": 1}
    }
    return VariationDescriptor(**genomic_del2_free_text)


@pytest.fixture(scope="module")
def genomic_del2_free_text_rse(genomic_del2_free_text,
                               genomic_del2_free_text_seq_loc):
    """Create a test fixture for genomic del RSE."""
    _id = "ga4gh:VA.bXvkvrOgc7R9KG2MM--H7dRJEDv-CEVa"
    genomic_del2_free_text["variation"] = {
        "type": "Allele",
        "id": _id,
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
    return VariationDescriptor(**genomic_del2_free_text)


@pytest.fixture(scope="module")
def genomic_del3():
    """Create test fixture containing params for genomic del VD."""
    params = {
        "id": "normalize.variation:NC_000023.11%3Ag.%2831060227_31100351%29_%2833274278_33417151%29del",  # noqa: E501
        "type": "VariationDescriptor",
        "variation": dict(),
        "molecule_context": "genomic",
        "structural_type": "SO:0001743",
        "vrs_ref_allele_seq": None
    }
    return params


@pytest.fixture(scope="module")
def genomic_del3_vrc(genomic_del3, genomic_del3_dup3_loc):
    """Create a test fixture for genomic del relative cnv."""
    _id = "ga4gh:RCN.dYerC8FSiqcexo8X1n3XUKpAckoAsfOK"
    genomic_del3["variation"] = {
        "type": "RelativeCopyNumber",
        "id": _id,
        "location": genomic_del3_dup3_loc,
        "relative_copy_class": "partial loss"
    }
    return VariationDescriptor(**genomic_del3)


@pytest.fixture(scope="module")
def genomic_del3_vac(genomic_del3, genomic_del3_dup3_loc):
    """Create a test fixture for genomic del absolute cnv."""
    _id = "ga4gh:ACN.2ZUQcccwvtoGZ5LZZRUoDZp6218Y6sQK"
    genomic_del3["variation"] = {
        "type": "AbsoluteCopyNumber",
        "id": _id,
        "location": genomic_del3_dup3_loc,
        "copies": {"type": "Number", "value": 2}
    }
    return VariationDescriptor(**genomic_del3)


@pytest.fixture(scope="module")
def genomic_del3_rse_lse(genomic_del3):
    """Create test fixture for genomic del rse and lse."""
    params = {
        "id": genomic_del3["id"],
        "type": genomic_del3["type"],
        "variation": {
            "id": "ga4gh:VT.tmA3mpMy9HKUweaB8aYsq6uuejEx9iK7",
            "type": "Text",
            "definition": "NC_000023.11:g.(31060227_31100351)_(33274278_33417151)del"
        }
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_del3_free_text():
    """Create test fixture containing params for genomic del VD."""
    params = {
        "id": "normalize.variation:EFNB1%20g.%2868839265_68839268%29_%2868841120_68841125%29del",  # noqa: E501
        "type": "VariationDescriptor",
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
                    "name": "chromosome_location",
                    "value": {
                        "species_id": "taxonomy:9606",
                        "start": "q13.1",
                        "end": "q13.1",
                        "id": "ga4gh:CL.gZnektnDFlImvvJ_gmIzLx8Efyqv2l29",
                        "type": "ChromosomeLocation",
                        "chr": "X"
                    }
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
            "gene": "hgnc:3226"
        }
    }
    return params


@pytest.fixture(scope="module")
def genomic_del3_free_text_subject():
    """Create test fixture for genomic del3 free text location"""
    return {
        "id": "ga4gh:SL.rqRKvcOMBF9hB5DMhRZ50hu7pnDlCMpi",
        "sequence_id": "ga4gh:SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
        "start": {"min": 68839264, "max": 68839267, "type": "DefiniteRange"},
        "end": {"min": 68841121, "max": 68841126, "type": "DefiniteRange"},
        "type": "SequenceLocation"
    }


@pytest.fixture(scope="module")
def genomic_del3_free_text_vrc(genomic_del3_free_text, genomic_del3_free_text_subject):
    """Create a test fixture for genomic del relative cnv."""
    _id = "ga4gh:RCN.6c0tRlHyFYGSeDEmSyn0nrZJLQDkwVsG"
    genomic_del3_free_text["variation"] = {
        "type": "RelativeCopyNumber",
        "id": _id,
        "location": genomic_del3_free_text_subject,
        "relative_copy_class": "partial loss"
    }
    return VariationDescriptor(**genomic_del3_free_text)


@pytest.fixture(scope="module")
def genomic_del3_free_text_vac(genomic_del3_free_text, genomic_del3_free_text_subject):
    """Create a test fixture for genomic del absolute cnv."""
    _id = "ga4gh:ACN.ZIaBXgnKLD-AkQhs3W3SMuhMmk80DSlQ"
    genomic_del3_free_text["variation"] = {
        "type": "AbsoluteCopyNumber",
        "id": _id,
        "location": genomic_del3_free_text_subject,
        "copies": {"type": "Number", "value": 2}
    }
    return VariationDescriptor(**genomic_del3_free_text)


@pytest.fixture(scope="module")
def genomic_del3_free_text_rse_lse(genomic_del3_free_text):
    """Create test fixture for genomic del rse and lse."""
    params = {
        "id": genomic_del3_free_text["id"],
        "type": genomic_del3_free_text["type"],
        "variation": {
            "id": "ga4gh:VT.9mGg0U_Z7NZCFV3jrLdGxSQU03g7z3Z1",
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
        "variation": dict(),
        "molecule_context": "genomic",
        "structural_type": "SO:0001743",
        "vrs_ref_allele_seq": None
    }
    return params


@pytest.fixture(scope="module")
def genomic_del4_vrc(genomic_del4, genomic_del4_seq_loc):
    """Create a test fixture for genomic del relative cnv."""
    _id = "ga4gh:RCN.BwZOFAfo5u8TcwbR3DMi8qbIImv96VQU"
    genomic_del4["variation"] = {
        "type": "RelativeCopyNumber",
        "id": _id,
        "location": genomic_del4_seq_loc,
        "relative_copy_class": "partial loss"
    }
    return VariationDescriptor(**genomic_del4)


@pytest.fixture(scope="module")
def genomic_del4_vac(genomic_del4, genomic_del4_seq_loc):
    """Create a test fixture for genomic del absolute cnv."""
    _id = "ga4gh:ACN.C9lafWMOwiuIc8SbJ3Ie32zYFiJwS9j7"
    genomic_del4["variation"] = {
        "type": "AbsoluteCopyNumber",
        "id": _id,
        "location": genomic_del4_seq_loc,
        "copies": {"type": "Number", "value": 1}
    }
    return VariationDescriptor(**genomic_del4)


@pytest.fixture(scope="module")
def genomic_del4_rse_lse(genomic_del4):
    """Create test fixture for genomic dup rse and lse."""
    params = {
        "id": genomic_del4["id"],
        "type": genomic_del4["type"],
        "variation": {
            "id": "ga4gh:VT.whBY5P24WVxF1wneDcI8x8btqorJUWXQ",
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
                    "name": "chromosome_location",
                    "value": {
                        "species_id": "taxonomy:9606",
                        "start": "q36.3",
                        "end": "q36.3",
                        "id": "ga4gh:CL.9SOUKUdu11HvIXHPui35LFZpUp9TECF6",
                        "type": "ChromosomeLocation",
                        "chr": "2"
                    }
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
            "gene": "hgnc:2206"
        }
    }
    return params


@pytest.fixture(scope="module")
def genomic_del4_free_text_subject():
    """Create test fixture for genomic del4 free text location"""
    return {
        "id": "ga4gh:SL.TxgMX9W0YT_v1ix8uLviaLGdcqMd6WRg",
        "sequence_id": "ga4gh:SQ.pnAqCRBrTsUoBghSD1yp_jXWSmlbdh4g",
        "start": {"value": 227022027, "comparator": "<=", "type": "IndefiniteRange"},
        "end": {"value": 227025830, "comparator": ">=", "type": "IndefiniteRange"},
        "type": "SequenceLocation"
    }


@pytest.fixture(scope="module")
def genomic_del4_free_text_vrc(genomic_del4_free_text, genomic_del4_free_text_subject):
    """Create a test fixture for genomic del relative cnv."""
    _id = "ga4gh:RCN.XqfrZ9k9mwDO0cM9duK7ooOih0iR1H2Q"
    genomic_del4_free_text["variation"] = {
        "type": "RelativeCopyNumber",
        "id": _id,
        "location": genomic_del4_free_text_subject,
        "relative_copy_class": "partial loss"
    }
    return VariationDescriptor(**genomic_del4_free_text)


@pytest.fixture(scope="module")
def genomic_del4_free_text_vac(genomic_del4_free_text, genomic_del4_free_text_subject):
    """Create a test fixture for genomic del absolute cnv."""
    _id = "ga4gh:ACN.xE7GFDRDsq5DHaRhcJ8xgdC6FYawHoxP"
    genomic_del4_free_text["variation"] = {
        "type": "AbsoluteCopyNumber",
        "id": _id,
        "location": genomic_del4_free_text_subject,
        "copies": {"type": "Number", "value": 1}
    }
    return VariationDescriptor(**genomic_del4_free_text)


@pytest.fixture(scope="module")
def genomic_del4_free_text_rse_lse(genomic_del4_free_text):
    """Create test fixture for genomic dup rse and lse."""
    params = {
        "id": genomic_del4_free_text["id"],
        "type": genomic_del4_free_text["type"],
        "variation": {
            "id": "ga4gh:VT.lT0rFYhOGFLA9MYA8ypnCf5q-CkV8dJv",
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
        "variation": {
            "id": "ga4gh:RCN.7EM-Wsg_7mmAE1LW8cmRI3QwKhJCA24a",
            "location": {
                "id": "ga4gh:SL.gUeB872FGVaphqoSAfI0gz4KXJvpZKL_",
                "sequence_id": "ga4gh:SQ.pnAqCRBrTsUoBghSD1yp_jXWSmlbdh4g",
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
        "variation": {
            "id": "ga4gh:RCN.2q7DKevv8nUh87Sl00Z7l50h047Ti2at",
            "location": {
                "id": "ga4gh:SL.ykRzA8IFueiCG7oznnN4teL2nXXBshHV",
                "sequence_id": "ga4gh:SQ.8_liLu1aycC0tPQPFmUaGXJLDs5SbPZ5",
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
        "id": "normalize.variation:NC_000023.11%3Ag.%28%3F_18575354%29_18653629del",
        "type": "VariationDescriptor",
        "variation": dict(),
        "molecule_context": "genomic",
        "structural_type": "SO:0001743",
        "vrs_ref_allele_seq": None
    }
    return params


def genomic_del5_abs_cnv(params, genomic_del5_seq_loc):
    """Create genomic del5 absolute cnv"""
    _id = "ga4gh:ACN.UFPa9YtomnJ00-rlfsQPzc2G05EXgFM3"
    params["variation"] = {
        "type": "AbsoluteCopyNumber",
        "id": _id,
        "location": genomic_del5_seq_loc,
        "copies": {"type": "Number", "value": 3}
    }


def genomic_del5_rel_cnv(params, genomic_del5_seq_loc):
    """Create genomic del5 relative cnv"""
    _id = "ga4gh:RCN.9rG3a5u3JODwQGVrv1IgAjG6SZgdvraH"
    params["variation"] = {
        "type": "RelativeCopyNumber",
        "id": _id,
        "location": genomic_del5_seq_loc,
        "relative_copy_class": "partial loss"
    }


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
            "id": "ga4gh:VT.xCLHh3GpCebrP6KDMsWZRdIiW7Sti27H",
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
                    "name": "chromosome_location",
                    "value": {
                        "species_id": "taxonomy:9606",
                        "start": "p22.13",
                        "end": "p22.13",
                        "id": "ga4gh:CL.so_qfQ0i7DCprrv5baguJALdalff7oAW",
                        "type": "ChromosomeLocation",
                        "chr": "X"
                    }
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
            "gene": "hgnc:11411"
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
            "id": "ga4gh:VT.xCLHh3GpCebrP6KDMsWZRdIiW7Sti27H",
            "type": "Text",
            "definition": "NC_000023.11:g.(?_18575354)_18653629del"
        }
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_del6():
    """Create test fixture containing params for genomic del VD."""
    params = {
        "id": "normalize.variation:NC_000006.12%3Ag.133462764_%28133464858_%3F%29del",
        "type": "VariationDescriptor",
        "variation": dict(),
        "molecule_context": "genomic",
        "structural_type": "SO:0001743",
        "vrs_ref_allele_seq": None
    }
    return params


def genomic_del6_rel_cnv(params, genomic_del6_seq_loc):
    """Create genomic del6 relative cnv"""
    _id = "ga4gh:RCN.zsagK87b_RdK4_QZGMnbNl39LCuJUjAr"
    params["variation"] = {
        "type": "RelativeCopyNumber",
        "id": _id,
        "location": genomic_del6_seq_loc,
        "relative_copy_class": "partial loss"
    }


def genomic_del6_abs_cnv(params, genomic_del6_seq_loc):
    """Create genomic del6 absolute cnv"""
    _id = "ga4gh:ACN.ZnnJNutwCrHNzFQamAWXMbLC7PfILmqA"
    params["variation"] = {
        "type": "AbsoluteCopyNumber",
        "id": _id,
        "location": genomic_del6_seq_loc,
        "copies": {"type": "Number", "value": 1}
    }


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
            "id": "ga4gh:VT.Df49jbB-kZ2LSm180uA9wn4TT_p215yX",
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
                    "name": "chromosome_location",
                    "value": {
                        "species_id": "taxonomy:9606",
                        "start": "q23.2",
                        "end": "q23.2",
                        "id": "ga4gh:CL.mVXPq-PjEnHwc-7o9vpkOzOhSAghUu4d",
                        "type": "ChromosomeLocation",
                        "chr": "6"
                    }
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
            "gene": "hgnc:3522"
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
            "id": "ga4gh:VT.a3kXhodtO3tgsdPlEL39Ql4jOuCpOc0s",
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
