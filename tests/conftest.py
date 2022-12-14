"""Create methods used throughout tests."""
import asyncio

import pytest
from ga4gh.vrsatile.pydantic.vrsatile_models import VariationDescriptor

from variation.query import QueryHandler


@pytest.fixture(scope="session")
def event_loop(request):
    """Create an instance of the default event loop for each test case."""
    loop = asyncio.get_event_loop_policy().new_event_loop()
    yield loop
    loop.close()


@pytest.fixture(scope="session")
def test_query_handler():
    """Build normalize test fixture."""
    return QueryHandler()


@pytest.fixture(scope="session")
def test_cnv_handler(test_query_handler):
    """Create test fixture for copy number variation handler"""
    return test_query_handler.to_copy_number_handler


@pytest.fixture(scope="session")
def vhl_gene_context():
    """Create a VHL gene context."""
    return {
        "id": "normalize.gene:VHL",
        "type": "GeneDescriptor",
        "label": "VHL",
        "gene": "hgnc:12687",
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
                    "iuphar:3204",
                    "orphanet:120467",
                    "ccds:CCDS2597",
                    "uniprot:P40337"
                ]
            },
            {
                "type": "Extension",
                "name": "hgnc_locations",
                "value": [
                    {
                        "start": "p25.3",
                        "species_id": "taxonomy:9606",
                        "end": "p25.3",
                        "id": "ga4gh:CL.idJ3P9Ld6SC4XmgZwV5zvHCqzzFyaXBA",
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
                        "sequence_id": "ga4gh:SQ.Zu7h9AggXxhTaGVsy7h_EZSChSZGcmgX",
                        "start": {"type": "Number", "value": 10141777},
                        "end": {"type": "Number", "value": 10153667},
                        "id": "ga4gh:SL.HqLrCh14JJXTqj2RXkGNEpni9eVnG_P9",
                        "type": "SequenceLocation"
                    }
                ]
            },
            {
                "type": "Extension",
                "name": "ncbi_locations",
                "value": [
                    {
                        "start": "p25.3",
                        "species_id": "taxonomy:9606",
                        "end": "p25.3",
                        "id": "ga4gh:CL.idJ3P9Ld6SC4XmgZwV5zvHCqzzFyaXBA",
                        "type": "ChromosomeLocation",
                        "chr": "3"
                    },
                    {
                        "sequence_id": "ga4gh:SQ.Zu7h9AggXxhTaGVsy7h_EZSChSZGcmgX",
                        "start": {"type": "Number", "value": 10141777},
                        "end": {"type": "Number", "value": 10153667},
                        "id": "ga4gh:SL.HqLrCh14JJXTqj2RXkGNEpni9eVnG_P9",
                        "type": "SequenceLocation"
                    }
                ]
            },
            {
                "name": "previous_symbols",
                "value": [
                    "RCA1"
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
        ]
    }


@pytest.fixture(scope="session")
def braf_ncbi_seq_loc():
    """Create test fixture for BRAF ncbi priority sequence location"""
    return {
        "sequence_id": "ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
        "start": {"type": "Number", "value": 140713327},
        "end": {"type": "Number", "value": 140924929},
        "id": "ga4gh:SL.po-AExwyqkstDx3JWYn6plIlxn5eojv4",
        "type": "SequenceLocation"
    }


@pytest.fixture(scope="session")
def braf_gene_context(braf_ncbi_seq_loc):
    """Create BRAF gene context test fixture."""
    return {
        "id": "normalize.gene:BRAF",
        "type": "GeneDescriptor",
        "label": "BRAF",
        "gene": "hgnc:1097",
        "xrefs": [
            "ncbigene:673",
            "ensembl:ENSG00000157764"
        ],
        "alternate_labels": [
            "BRAF1",
            "RAFB1",
            "B-raf",
            "B-RAF1",
            "NS7",
            "BRAF-1"
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
                "value": "B-Raf proto-oncogene, serine/threonine kinase"
            },
            {
                "type": "Extension",
                "name": "associated_with",
                "value": [
                    "ucsc:uc003vwc.5",
                    "pubmed:1565476",
                    "omim:164757",
                    "vega:OTTHUMG00000157457",
                    "ccds:CCDS5863",
                    "iuphar:1943",
                    "ccds:CCDS87555",
                    "orphanet:119066",
                    "refseq:NM_004333",
                    "ena.embl:M95712",
                    "pubmed:2284096",
                    "uniprot:P15056",
                    "cosmic:BRAF"
                ]
            },
            {
                "type": "Extension",
                "name": "hgnc_locations",
                "value": [
                    {
                        "start": "q34",
                        "species_id": "taxonomy:9606",
                        "end": "q34",
                        "id": "ga4gh:CL.ZZZYpOwuW1BLLJXc_Dm4eVZ5E0smVYCc",
                        "type": "ChromosomeLocation",
                        "chr": "7"
                    }
                ]
            },
            {
                "type": "Extension",
                "name": "ensembl_locations",
                "value": [
                    {
                        "sequence_id": "ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
                        "start": {"type": "Number", "value": 140719326},
                        "end": {"type": "Number", "value": 140924929},
                        "id": "ga4gh:SL.MPDI-H-JK1AMOnLdvVO6zQjbAXmSnzQ4",
                        "type": "SequenceLocation"
                    }
                ]
            },
            {
                "type": "Extension",
                "name": "ncbi_locations",
                "value": [
                    {
                        "start": "q34",
                        "species_id": "taxonomy:9606",
                        "end": "q34",
                        "id": "ga4gh:CL.ZZZYpOwuW1BLLJXc_Dm4eVZ5E0smVYCc",
                        "type": "ChromosomeLocation",
                        "chr": "7"
                    },
                    braf_ncbi_seq_loc
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


@pytest.fixture(scope="session")
def egfr_context():
    """Create EGFR gene context test fixture"""
    return {
        "id": "normalize.gene:EGFR",
        "type": "GeneDescriptor",
        "label": "EGFR",
        "gene": "hgnc:3236",
        "xrefs": [
            "ncbigene:1956",
            "ensembl:ENSG00000146648"
        ],
        "alternate_labels": [
            "HER1",
            "NISBD2",
            "ERBB",
            "PIG61",
            "mENA",
            "ERBB1",
            "ERRP"
        ],
        "extensions": [
            {
                "type": "Extension",
                "name": "symbol_status",
                "value": "approved"
            },
            {
                "name": "approved_name",
                "value": "epidermal growth factor receptor",
                "type": "Extension"
            },
            {
                "type": "Extension",
                "name": "associated_with",
                "value": [
                    "ccds:CCDS5516",
                    "ccds:CCDS5514",
                    "iuphar:1797",
                    "uniprot:P00533",
                    "vega:OTTHUMG00000023661",
                    "ucsc:uc003tqk.4",
                    "ccds:CCDS5515",
                    "refseq:NM_005228",
                    "ccds:CCDS87506",
                    "ccds:CCDS47587",
                    "pubmed:1505215",
                    "cosmic:EGFR",
                    "ccds:CCDS87507",
                    "omim:131550",
                    "orphanet:121311"
                ]
            },
            {
                "type": "Extension",
                "name": "hgnc_locations",
                "value": [
                    {
                        "start": "p11.2",
                        "species_id": "taxonomy:9606",
                        "end": "p11.2",
                        "id": "ga4gh:CL.UTlCnUItQ_jX44tlpLrVx277hUtRvDNG",
                        "type": "ChromosomeLocation",
                        "chr": "7"
                    }
                ]
            },
            {
                "type": "Extension",
                "name": "ensembl_locations",
                "value": [
                    {
                        "sequence_id": "ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
                        "start": {"type": "Number", "value": 55019016},
                        "end": {"type": "Number", "value": 55211628},
                        "id": "ga4gh:SL.OBv25qadiUzZU3eJ2ImcA2XQjI6wprss",
                        "type": "SequenceLocation"
                    }
                ]
            },
            {
                "type": "Extension",
                "name": "ncbi_locations",
                "value": [
                    {
                        "start": "p11.2",
                        "species_id": "taxonomy:9606",
                        "end": "p11.2",
                        "id": "ga4gh:CL.UTlCnUItQ_jX44tlpLrVx277hUtRvDNG",
                        "type": "ChromosomeLocation",
                        "chr": "7"
                    },
                    {
                        "sequence_id": "ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
                        "start": {"type": "Number", "value": 55019016},
                        "end": {"type": "Number", "value": 55211628},
                        "id": "ga4gh:SL.OBv25qadiUzZU3eJ2ImcA2XQjI6wprss",
                        "type": "SequenceLocation"
                    }
                ]
            },
            {
                "name": "previous_symbols",
                "value": [
                    "ERBB"
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
        ]
    }


@pytest.fixture(scope="session")
def erbb2_context():
    """Create test fixture for ERBB2 Gene Context."""
    return {
        "id": "normalize.gene:ERBB2",
        "type": "GeneDescriptor",
        "label": "ERBB2",
        "gene": "hgnc:3430",
        "xrefs": [
            "ncbigene:2064",
            "ensembl:ENSG00000141736"
        ],
        "alternate_labels": [
            "NGL",
            "CD340",
            "HER2",
            "NEU",
            "TKR1",
            "HER-2",
            "HER-2/neu",
            "VSCN2",
            "MLN 19",
            "MLN-19",
            "c-ERB-2",
            "c-ERB2",
            "p185(erbB2)"
        ],
        "extensions": [
            {
                "type": "Extension",
                "name": "symbol_status",
                "value": "approved"
            },
            {
                "name": "approved_name",
                "value": "erb-b2 receptor tyrosine kinase 2",
                "type": "Extension"
            },
            {
                "type": "Extension",
                "name": "associated_with",
                "value": [
                    "ucsc:uc002hso.4",
                    "ena.embl:X03363",
                    "ccds:CCDS77017",
                    "vega:OTTHUMG00000179300",
                    "ccds:CCDS77016",
                    "uniprot:P04626",
                    "refseq:NM_004448",
                    "ccds:CCDS74052",
                    "hcdmdb:CD340",
                    "omim:164870",
                    "ccds:CCDS32642",
                    "ccds:CCDS45667",
                    "cosmic:ERBB2",
                    "iuphar:2019",
                    "pubmed:1675005",
                    "pubmed:2885835",
                    "pubmed:2903500"
                ]
            },
            {
                "type": "Extension",
                "name": "hgnc_locations",
                "value": [
                    {
                        "start": "q12",
                        "species_id": "taxonomy:9606",
                        "end": "q12",
                        "id": "ga4gh:CL.cD8l0i9fES2sth5nVc5qUGTD-QW-xDyT",
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
                        "sequence_id": "ga4gh:SQ.dLZ15tNO1Ur0IcGjwc3Sdi_0A6Yf4zm7",
                        "start": {"type": "Number", "value": 39687913},
                        "end": {"type": "Number", "value": 39730426},
                        "id": "ga4gh:SL.A-GcJQazJLUu9Ad20lpwAu5vbRjxeWSn",
                        "type": "SequenceLocation"
                    }
                ]
            },
            {
                "type": "Extension",
                "name": "ncbi_locations",
                "value": [
                    {
                        "start": "q12",
                        "species_id": "taxonomy:9606",
                        "end": "q12",
                        "id": "ga4gh:CL.cD8l0i9fES2sth5nVc5qUGTD-QW-xDyT",
                        "type": "ChromosomeLocation",
                        "chr": "17"
                    },
                    {
                        "sequence_id": "ga4gh:SQ.dLZ15tNO1Ur0IcGjwc3Sdi_0A6Yf4zm7",
                        "start": {"type": "Number", "value": 39688093},
                        "end": {"type": "Number", "value": 39728658},
                        "id": "ga4gh:SL.IAvrcyN2dV0U9wV7zY9lbD-4cV0xejdP",
                        "type": "SequenceLocation"
                    }
                ]
            },
            {
                "name": "previous_symbols",
                "value": [
                    "NGL"
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
        ]
    }


@pytest.fixture(scope="session")
def prpf8_ncbi_seq_loc():
    """Create test fixture for PRPF8 ncbi priority sequence location"""
    return {
        "sequence_id": "ga4gh:SQ.dLZ15tNO1Ur0IcGjwc3Sdi_0A6Yf4zm7",
        "start": {"type": "Number", "value": 1650628},
        "end": {"type": "Number", "value": 1684867},
        "id": "ga4gh:SL.tiIyPIF9TqW6xrAkH1T0Z_oKIGQf-G3A",
        "type": "SequenceLocation"
    }


@pytest.fixture(scope="session")
def prpf8_gene_context(prpf8_ncbi_seq_loc):
    """Create test fixture for PRPF8 gene context"""
    return {
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
                        "start": "p13.3",
                        "species_id": "taxonomy:9606",
                        "end": "p13.3",
                        "id": "ga4gh:CL.cgXboJas91hqg8zGKLMZYN0RdB_RCj6k",
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
                        "sequence_id": "ga4gh:SQ.dLZ15tNO1Ur0IcGjwc3Sdi_0A6Yf4zm7",
                        "start": {"type": "Number", "value": 1650628},
                        "end": {"type": "Number", "value": 1684867},
                        "id": "ga4gh:SL.tiIyPIF9TqW6xrAkH1T0Z_oKIGQf-G3A",
                        "type": "SequenceLocation"
                    }
                ]
            },
            {
                "type": "Extension",
                "name": "ncbi_locations",
                "value": [
                    {
                        "start": "p13.3",
                        "species_id": "taxonomy:9606",
                        "end": "p13.3",
                        "id": "ga4gh:CL.cgXboJas91hqg8zGKLMZYN0RdB_RCj6k",
                        "type": "ChromosomeLocation",
                        "chr": "17"
                    },
                    prpf8_ncbi_seq_loc,
                    {
                        "sequence_id": "ga4gh:SQ.3Wx-9rRd5d7m3WxtJ_HScX3Bz1MiQWjR",
                        "start": {"type": "Number", "value": 80656},
                        "end": {"type": "Number", "value": 114895},
                        "id": "ga4gh:SL.Z1Vjxb1xSfKYknyz4MmOARdEqpCMuIFR",
                        "type": "SequenceLocation"
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
        "gene": "hgnc:17340"
    }


@pytest.fixture(scope="session")
def braf_600loc():
    """Create test fixture for BRAF 600 location"""
    return {
        "id": "ga4gh:SL.xfBTztcmMstx8jrrdgPiE_BUoLHLFMMS",
        "end": {"value": 600, "type": "Number"},
        "start": {"value": 599, "type": "Number"},
        "sequence_id": "ga4gh:SQ.cQvw4UsHHRRlogxbWCB8W-mKD4AraM9y",
        "type": "SequenceLocation"
    }


@pytest.fixture(scope="session")
def braf_v600e(braf_gene_context, braf_600loc):
    """Create BRAF V600E protein test fixture."""
    params = {
        "id": "normalize.variation:BRAF%20V600E",
        "type": "VariationDescriptor",
        "variation": {
            "id": "ga4gh:VA.h313H4CQh6pogbbSJ3H5pI1cPoh9YMm_",
            "location": braf_600loc,
            "state": {
                "sequence": "E",
                "type": "LiteralSequenceExpression"
            },
            "type": "Allele"
        },
        "molecule_context": "protein",
        "structural_type": "SO:0001606",
        "vrs_ref_allele_seq": "V",
        "gene_context": braf_gene_context
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="session")
def vhl_silent(vhl_gene_context):
    """Create NP_000542.1:p.Pro61 fixture."""
    params = {
        "id": "normalize.variation:NP_000542.1%3Ap.Pro61%3D",
        "type": "VariationDescriptor",
        "variation": {
            "id": "ga4gh:VA.hFvpLOEfU4Qtxd1_pdSx-_XnIJng9Xnb",
            "location": {
                "id": "ga4gh:SL.Kx99ER-oHCNP3RwwatOYn9IN5LRRxiy-",
                "end": {"value": 61, "type": "Number"},
                "start": {"value": 60, "type": "Number"},
                "sequence_id": "ga4gh:SQ.z-Oa0pZkJ6GHJHOYM7h5mY_umc0SJzTu",
                "type": "SequenceLocation"
            },
            "state": {
                "sequence": "P",
                "type": "LiteralSequenceExpression"
            },
            "type": "Allele"
        },
        "molecule_context": "protein",
        "structural_type": "SO:0001017",
        "vrs_ref_allele_seq": "P",
        "gene_context": vhl_gene_context
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="session")
def protein_insertion(egfr_context):
    """Create test fixture for NP protein insertion."""
    params = {
        "id": "normalize.variation:NP_005219.2%3Ap.Asp770_Asn771insGlyLeu",
        "type": "VariationDescriptor",
        "variation": {
            "id": "ga4gh:VA.Daydg17safvIKs_ENTrvKpfoSooKgImP",
            "location": {
                "id": "ga4gh:SL.ozw2OUd_hkRUcAo4zUM_jH40Wlbd_lb0",
                "end": {"value": 770, "type": "Number"},
                "start": {"value": 770, "type": "Number"},
                "sequence_id": "ga4gh:SQ.vyo55F6mA6n2LgN4cagcdRzOuh38V4mE",
                "type": "SequenceLocation"
            },
            "state": {
                "sequence": "GL",
                "type": "LiteralSequenceExpression"
            },
            "type": "Allele"
        },
        "molecule_context": "protein",
        "structural_type": "SO:0001605",
        "gene_context": egfr_context
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="session")
def protein_deletion_np_range(erbb2_context):
    """Create test fixture for protein deletion using NP accession and
    range for deletion.
    """
    params = {
        "id": "normalize.variation:NP_004439.2%3Ap.Leu755_Thr759del",
        "type": "VariationDescriptor",
        "variation": {
            "id": "ga4gh:VA.gb8REjJatcpIgb1B3LKiIDI4DhJd70Bk",
            "location": {
                "id": "ga4gh:SL.Ca3e2urICwddGClFXppmmMeGr3Zqg-i8",
                "end": {"value": 759, "type": "Number"},
                "start": {"value": 754, "type": "Number"},
                "sequence_id": "ga4gh:SQ.AF1UFydIo02-bMplonKSfxlWY2q6ze3m",
                "type": "SequenceLocation"
            },
            "state": {
                "sequence": "",
                "type": "LiteralSequenceExpression"
            },
            "type": "Allele"
        },
        "molecule_context": "protein",
        "structural_type": "SO:0001604",
        "vrs_ref_allele_seq": "LRENT",
        "gene_context": erbb2_context
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="session")
def braf_v600e_genomic_sub():
    """Create test fixture for NC_000007.14:g.140753336A>T"""
    return {
        "id": "ga4gh:VA.3xFHF399HGbG1JUf5uwcj3oWVKZJ70oX",
        "location": {
            "id": "ga4gh:SL.WBLxdkoypnRME6b8tJtlOWqZKU1ruqY1",
            "end": {"value": 140753336, "type": "Number"},
            "start": {"value": 140753335, "type": "Number"},
            "sequence_id": "ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
            "type": "SequenceLocation"
        },
        "state": {
            "sequence": "T",
            "type": "LiteralSequenceExpression"
        },
        "type": "Allele"
    }


@pytest.fixture(scope="session")
def genomic_dup1_seq_loc():
    """Create test fixture containing genomic dup1 sequence location"""
    return {
        "id": "ga4gh:SL.KefUQwlqEBGtzoNO-MzOozx7_H1uP-fD",
        "sequence_id": "ga4gh:SQ.Zu7h9AggXxhTaGVsy7h_EZSChSZGcmgX",
        "start": {"value": 49531260, "type": "Number"},
        "end": {"value": 49531262, "type": "Number"},
        "type": "SequenceLocation"
    }


@pytest.fixture(scope="session")
def genomic_dup1_38_vac(genomic_dup1_seq_loc):
    """Create test fixture for absolute copy number dup1 on GRCh38"""
    return {
        "type": "AbsoluteCopyNumber",
        "id": "ga4gh:ACN._C6yp4SRGVfuRmMiJShIKYCK3dSX0vNF",
        "location": genomic_dup1_seq_loc,
        "copies": {"type": "Number", "value": 3}
    }


@pytest.fixture(scope="session")
def genomic_dup2_seq_loc():
    """Create genomic dup2 sequence location"""
    return {
        "id": "ga4gh:SL.Efhss2fsPGP9kmSO2eoXKUiMS-GZGBgh",
        "sequence_id": "ga4gh:SQ.yC_0RBj3fgBlvgyAuycbzdubtLxq-rE0",
        "start": {"value": 2087937, "type": "Number"},
        "end": {"value": 2087948, "type": "Number"},
        "type": "SequenceLocation"
    }


@pytest.fixture(scope="session")
def genomic_dup2_38_vac(genomic_dup2_seq_loc):
    """Create test fixture for absolute copy number dup2 on GRCh38"""
    return {
        "type": "AbsoluteCopyNumber",
        "id": "ga4gh:ACN.pUVPSxSZ5n9__GwA2FOQZvNtF_XtIzh1",
        "location": genomic_dup2_seq_loc,
        "copies": {"type": "Number", "value": 3}
    }


@pytest.fixture(scope="session")
def genomic_del3_dup3_loc():
    """Create genomic del3 dup3 sequence location"""
    return {
        "id": "ga4gh:SL.RANaZSqxDM1hoJeXZAfSXErh6XqPrijo",
        "sequence_id": "ga4gh:SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
        "start": {"min": 31060226, "max": 31100350, "type": "DefiniteRange"},
        "end": {"min": 33274279, "max": 33417152, "type": "DefiniteRange"},
        "type": "SequenceLocation"
    }


@pytest.fixture(scope="session")
def genoimc_dup4_loc():
    """Create genoimc dup4 sequence location"""
    return {
        "id": "ga4gh:SL.uD8efGXIXdiMNyHX4MogVF0jA28jIWb4",
        "sequence_id": "ga4gh:SQ.-A1QmD_MatoqxvgVxBLZTONHz9-c7nQo",
        "start": {"value": 30417575, "comparator": "<=", "type": "IndefiniteRange"},
        "end": {"value": 31394018, "comparator": ">=", "type": "IndefiniteRange"},
        "type": "SequenceLocation"
    }


@pytest.fixture(scope="session")
def genomic_dup5_loc():
    """Create genoimc dup5 sequence location"""
    return {
        "id": "ga4gh:SL.GzmuP1MBA9qILR8fVFhp4BdUEcaLwKaR",
        "sequence_id": "ga4gh:SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
        "start": {"value": 154021811, "comparator": "<=", "type": "IndefiniteRange"},
        "end": {"value": 154092209, "type": "Number"},
        "type": "SequenceLocation"
    }


@pytest.fixture(scope="session")
def genoimc_dup6_loc():
    """Create genoimc dup6 sequence location"""
    return {
        "id": "ga4gh:SL.8j0dwTvx7zKHVk2JCDt1eqoEt-o993hg",
        "sequence_id": "ga4gh:SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
        "start": {"value": 154021811, "type": "Number"},
        "end": {"value": 154092209, "comparator": ">=", "type": "IndefiniteRange"},
        "type": "SequenceLocation"
    }


@pytest.fixture(scope="session")
def genomic_del1_seq_loc():
    """Create genomic del1 sequence location"""
    return {
        "id": "ga4gh:SL.UI2JMyUJaYVNqYd1MI-sAxPsxnolalcw",
        "sequence_id": "ga4gh:SQ.Zu7h9AggXxhTaGVsy7h_EZSChSZGcmgX",
        "start": {"value": 10149810, "type": "Number"},
        "end": {"value": 10149811, "type": "Number"},
        "type": "SequenceLocation"
    }


@pytest.fixture(scope="session")
def genomic_del1_38_vac(genomic_del1_seq_loc):
    """Create test fixture for absolute copy number del1 on GRCh38"""
    return {
        "type": "AbsoluteCopyNumber",
        "id": "ga4gh:ACN.-B_iqSMiz6RCkwNC4PeL0DGUPHz0ZuAO",
        "location": genomic_del1_seq_loc,
        "copies": {"type": "Number", "value": 1}
    }


@pytest.fixture(scope="session")
def genomic_del2_seq_loc():
    """Create genomic del2 sequence location"""
    return {
        "id": "ga4gh:SL.rTvq-fnraKFAvs9aMgduEsi4Z7RTaxO5",
        "sequence_id": "ga4gh:SQ.Zu7h9AggXxhTaGVsy7h_EZSChSZGcmgX",
        "start": {"value": 10146594, "type": "Number"},
        "end": {"value": 10146613, "type": "Number"},
        "type": "SequenceLocation"
    }


@pytest.fixture(scope="session")
def genomic_del2_38_vac(genomic_del2_seq_loc):
    """Create test fixture for absolute copy number del1 on GRCh38"""
    return {
        "type": "AbsoluteCopyNumber",
        "id": "ga4gh:ACN.GXkCtgRsF6xmzYowbGt6RFTYck07nCWW",
        "location": genomic_del2_seq_loc,
        "copies": {"type": "Number", "value": 1}
    }


@pytest.fixture(scope="session")
def genomic_del4_seq_loc():
    """Create genomic del4 sequence location"""
    return {
        "id": "ga4gh:SL.dRc1d9ymsXhbb439OQE830RBELZ4aMXi",
        "sequence_id": "ga4gh:SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
        "start": {"value": 31120495, "comparator": "<=", "type": "IndefiniteRange"},
        "end": {"value": 33339477, "comparator": ">=", "type": "IndefiniteRange"},
        "type": "SequenceLocation"
    }


@pytest.fixture(scope="session")
def genomic_del5_seq_loc():
    """Create genomic del5 sequence location"""
    return {
        "id": "ga4gh:SL.t3PI8S74FRjr39sAp8l3SiJrbcRGRaMx",
        "sequence_id": "ga4gh:SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
        "start": {"value": 18575353, "comparator": "<=", "type": "IndefiniteRange"},
        "end": {"value": 18653629, "type": "Number"},
        "type": "SequenceLocation"
    }


@pytest.fixture(scope="session")
def genomic_del6_seq_loc():
    """Create genomic del6 sequence location"""
    return {
        "id": "ga4gh:SL.sawJf1SKy79nxruXFK4lb4pt2g6AU1Fy",
        "sequence_id": "ga4gh:SQ.0iKlIQk2oZLoeOG9P1riRU6hvL5Ux8TV",
        "start": {"value": 133462763, "type": "Number"},
        "end": {"value": 133464858, "comparator": ">=", "type": "IndefiniteRange"},
        "type": "SequenceLocation"
    }


@pytest.fixture(scope="session")
def grch38_genomic_insertion_seq_loc():
    """Create test fixture for GRCh38 genomic insertion seq location"""
    return {
        "id": "ga4gh:SL.IQbKwy5bXuJrZxrcbeazxmYMCJmEjgsW",
        "end": {"value": 39724743, "type": "Number"},
        "start": {"value": 39724731, "type": "Number"},
        "sequence_id": "ga4gh:SQ.dLZ15tNO1Ur0IcGjwc3Sdi_0A6Yf4zm7",
        "type": "SequenceLocation"
    }


@pytest.fixture(scope="session")
def grch38_genomic_insertion_variation(grch38_genomic_insertion_seq_loc):
    """Create a test fixture for NC_000017.10:g.37880993_37880994insGCTTACGTGATG"""
    return {
        "id": "ga4gh:VA.Zx2gPVfZkaplqGNrv1rzPvj3rjx0soSd",
        "location": grch38_genomic_insertion_seq_loc,
        "state": {
            "sequence": "TACGTGATGGCTTACGTGATGGCT",
            "type": "LiteralSequenceExpression"
        },
        "type": "Allele"
    }


@pytest.fixture(scope="session")
def braf_amplification(braf_ncbi_seq_loc, braf_gene_context):
    """Create test fixture for BRAF Amplification"""
    params = {
        "id": "normalize.variation:BRAF%20Amplification",
        "type": "VariationDescriptor",
        "variation": {
            "id": "ga4gh:RCN.tXX8oMzsJx3r9ZlqQlzk_K8Luz-bswdT",
            "location": braf_ncbi_seq_loc,
            "relative_copy_class": "EFO:0030072",
            "type": "RelativeCopyNumber"
        },
        "molecule_context": "genomic",
        "structural_type": "SO:0001880",
        "gene_context": braf_gene_context
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="session")
def prpf8_amplification(prpf8_ncbi_seq_loc, prpf8_gene_context):
    """Create test fixture for PRPF8 Amplification"""
    params = {
        "id": "normalize.variation:PRPF8%20AMPLIFICATION",
        "type": "VariationDescriptor",
        "variation": {
            "id": "ga4gh:RCN.DW0vRfIA0aI4AR0epEh_k-qrB2pdpZVw",
            "location": prpf8_ncbi_seq_loc,
            "relative_copy_class": "EFO:0030072",
            "type": "RelativeCopyNumber"
        },
        "molecule_context": "genomic",
        "structural_type": "SO:0001880",
        "gene_context": prpf8_gene_context
    }
    return VariationDescriptor(**params)


def assertion_checks(normalize_response, test_variation, label, ignore_id=False):
    """Check that normalize_response and test_variation are equal."""
    if not ignore_id:
        assert normalize_response.id == test_variation.id, "id"
    assert normalize_response.label == label
    assert normalize_response.type == test_variation.type, "type"
    if test_variation.variation.type != "Text":
        if test_variation.variation.id:
            assert normalize_response.variation.id == \
                   test_variation.variation.id, "variation.id"
        assert normalize_response.variation == test_variation.variation, "variation"
    else:
        if not ignore_id:
            assert normalize_response.variation.id == test_variation.variation.id
        assert normalize_response.variation.type == test_variation.variation.type
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
        if not ignore_id:
            assert resp_gene_context.id == test_variation_context.id, "gene_context.id"
        assert resp_gene_context.label == \
               test_variation_context.label, "gene_context.label"
        assert resp_gene_context.gene ==\
               test_variation_context.gene, "gene_context.gene"
        assert set(resp_gene_context.xrefs) ==\
               set(test_variation_context.xrefs), "gene_context.xrefs"
        if test_variation_context.alternate_labels:
            assert set(resp_gene_context.alternate_labels) == \
                   set(test_variation_context.alternate_labels), "gene_context.alternate_labels"  # noqa: E501
        assert len(resp_gene_context.extensions) == \
               len(test_variation_context.extensions), "len gene_context.extensions"
        for resp_ext in resp_gene_context.extensions:
            for test_var in test_variation_context.extensions:
                if resp_ext.name == test_var.name:
                    if resp_ext.name == "chromosome_location":
                        assert resp_ext.value == test_var.value, \
                            "gene_context.chromosome_location"
                    elif resp_ext.name == "associated_with":
                        assert set(resp_ext.value) == set(test_var.value), \
                            "gene_context.associated_with"
                    else:
                        if isinstance(test_var.value, list) and isinstance(test_var.value[0], str):  # noqa: E501
                            assert set(resp_ext.value) == set(test_var.value), \
                                f"gene_context.{resp_ext.name}"
                        else:
                            assert resp_ext.value == test_var.value,\
                                f"gene_context.{resp_ext.name}"
    else:
        assert not test_variation_context
