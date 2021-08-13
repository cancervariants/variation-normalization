"""Module for testing the normalize endpoint."""
import pytest
from ga4gh.vrs.dataproxy import SeqRepoDataProxy
from ga4gh.vrs.extras.translator import Translator
from variation.normalize import Normalize
from variation.schemas.ga4gh_vrsatile import VariationDescriptor
from variation.to_vrs import ToVRS
from variation.main import normalize as normalize_get_response
from variation.main import translate as to_vrs_get_response
from variation.classifiers import Classify
from variation.tokenizers import Tokenize
from variation.validators import Validate
from variation.translators import Translate
from variation.data_sources import SeqRepoAccess, TranscriptMappings, \
    UTA, MANETranscriptMappings
from variation.mane_transcript import MANETranscript
from variation.tokenizers import GeneSymbol
from variation.tokenizers.caches import GeneSymbolCache, AminoAcidCache
from datetime import datetime
import copy


@pytest.fixture(scope="module")
def test_normalize():
    """Build normalize test fixture."""
    class TestNormalize:

        def __init__(self):
            tokenizer = Tokenize()
            classifier = Classify()
            seqrepo_access = SeqRepoAccess()
            transcript_mappings = TranscriptMappings()
            gene_symbol = GeneSymbol(GeneSymbolCache())
            amino_acid_cache = AminoAcidCache()
            uta = UTA()
            mane_transcript_mappings = MANETranscriptMappings()
            dp = SeqRepoDataProxy(seqrepo_access.seq_repo_client)
            tlr = Translator(data_proxy=dp)
            mane_transcript = MANETranscript(seqrepo_access,
                                             transcript_mappings,
                                             mane_transcript_mappings, uta)
            validator = Validate(seqrepo_access, transcript_mappings,
                                 gene_symbol, mane_transcript, uta,
                                 dp, tlr, amino_acid_cache)
            translator = Translate()

            self.to_vrs = ToVRS(tokenizer, classifier, seqrepo_access,
                                transcript_mappings, gene_symbol,
                                amino_acid_cache, uta,
                                mane_transcript_mappings, mane_transcript,
                                validator, translator)
            self.test_normalize = Normalize(seqrepo_access, uta)

        def normalize(self, q):
            validations, warnings = self.to_vrs.get_validations(
                q, normalize_endpoint=True
            )
            resp = \
                self.test_normalize.normalize(q,
                                              validations,
                                              warnings)
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
            "BRAF1",
            "RAFB1",
            "B-raf",
            "B-RAF1",
            "NS7"
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
            "NGL",
            "CD340",
            "HER2",
            "NEU",
            "TKR1",
            "HER-2",
            "HER-2/neu",
            "VSCN2",
            "MLN 19"
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
                    "iuphar:2019"
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
            },
            {
                "name": "previous_symbols",
                "value": [
                    "NGL"
                ],
                "type": "Extension"
            }
        ]
    }


@pytest.fixture(scope='module')
def egfr_context():
    """Create EGFR gene context test fixture"""
    return {
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
            },
            {
                "name": "previous_symbols",
                "value": [
                    "ERBB"
                ],
                "type": "Extension"
            }
        ]
    }


@pytest.fixture(scope='module')
def limk2_gene_context():
    """Create LIMK2 gene context test fixture."""
    return {
        "id": "normalize.gene:LIMK2",
        "type": "GeneDescriptor",
        "label": "LIMK2",
        "value": {
            "id": "hgnc:6614",
            "type": "Gene"
        },
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
                "name": "chromosome_location",
                "value": {
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
            }
        ]
    }


@pytest.fixture(scope="module")
def braf_v600e(braf_gene_context):
    """Create BRAF V600E protein test fixture."""
    params = {
        "id": "normalize.variation:BRAF%20V600E",
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
        "molecule_context": "protein",
        "structural_type": "SO:0001606",
        "vrs_ref_allele_seq": "V",
        "gene_context": braf_gene_context
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def vhl(vhl_gene_context):
    """Create VHL Tyr185Ter fixture."""
    params = {
        "id": "normalize.variation:NP_000542.1%3Ap.Tyr185Ter",
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
        "molecule_context": "protein",
        "structural_type": "SO:0001617",
        "vrs_ref_allele_seq": "Y",
        "gene_context": vhl_gene_context
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def vhl_silent(vhl_gene_context):
    """Create NP_000542.1:p.Pro61 fixture."""
    params = {
        "id": "normalize.variation:NP_000542.1%3Ap.Pro61%3D",
        "type": "VariationDescriptor",
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
        "vrs_ref_allele_seq": "P",
        "gene_context": vhl_gene_context
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope='module')
def braf_v600e_nucleotide(braf_gene_context, braf_nuc_value):
    """Create a test fixture for BRAF V600E MANE select nucleotide hgvs."""
    params = {
        "id": "normalize.variation:NM_004333.4%3Ac.1799T%3EA",
        "type": "VariationDescriptor",
        "value_id": "ga4gh:VA.19rEOp0IBkrDkUA4gwwM-4Gde08-kBb1",
        "value": braf_nuc_value,
        "molecule_context": "transcript",
        "structural_type": "SO:0001483",
        "vrs_ref_allele_seq": "T",
        "gene_context": braf_gene_context
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope='module')
def nm_004448_coding_dna_delins(erbb2_context):
    """Create test fixture for NM_004448.4:c.2326_2327delinsCT."""
    params = {
        "id": "normalize.variation:NM_004448.4%3Ac.2326_2327delinsCT",
        "type": "VariationDescriptor",
        "value_id": "ga4gh:VA.sSFX2CO2DPTvE4MqnJ5VifnaQOGS0CVb",
        "value": {
            "location": {
                "interval": {
                    "end": 2502,
                    "start": 2500,
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
        "molecule_context": "transcript",
        "structural_type": "SO:1000032",
        "vrs_ref_allele_seq": "GG",
        "gene_context": erbb2_context
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope='module')
def nc_000007_genomic_delins(braf_gene_context):
    """Create test fixture for NC_000007.13:g.140453135_140453136delinsAT."""
    params = {
        "id": "normalize.variation:NC_000007.13%3Ag.140453135_140453136delinsAT",  # noqa: E501
        "type": "VariationDescriptor",
        "value_id": "ga4gh:VA._GzAG8_K8YwcYQk6bEvINNGM_hEViytU",
        "value": {
            "location": {
                "interval": {
                    "end": 2146,
                    "start": 2144,
                    "type": "SimpleInterval"
                },
                "sequence_id": "ga4gh:SQ.I_0feOk5bZ3VfH8ejhWQiMDe9o6o4QdR",
                "type": "SequenceLocation"
            },
            "state": {
                "sequence": "AT",
                "type": "SequenceState"
            },
            "type": "Allele"
        },
        "molecule_context": "transcript",
        "structural_type": "SO:1000032",
        "vrs_ref_allele_seq": "TG",
        "gene_context": braf_gene_context
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope='module')
def nm_000551(vhl_gene_context):
    """Create test fixture for NM_000551.4:c.615delinsAA."""
    params = {
        "id": 'normalize.variation:temp',
        "type": "VariationDescriptor",
        "value_id": "ga4gh:VA.SjJnUcJL1EyRFUb6f8PSJA4u3fyin2Wj",
        "value": {
            "location": {
                "interval": {
                    "end": 685,
                    "start": 684,
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
        "molecule_context": "transcript",
        "structural_type": "SO:1000032",
        "vrs_ref_allele_seq": "C",
        "gene_context": vhl_gene_context
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope='module')
def braf_nuc_value():
    """Create test fixture for BRAF V600E value on c. coordinate."""
    return {
        "location": {
            "interval": {
                "end": 2145,
                "start": 2144,
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
    }


@pytest.fixture(scope='module')
def coding_dna_silent_mutation(braf_gene_context, braf_nuc_value):
    """Create test fixture for NM_004333.4:c.1799=."""
    value = copy.deepcopy(braf_nuc_value)
    value['state']['sequence'] = 'T'
    params = {
        "id": 'normalize.variation:NM_004333.4%3Ac.1799%3D',
        "type": "VariationDescriptor",
        "value_id": "ga4gh:VA.bVNMOANetNE2z4PZ1j0DmwUL1rULmqkN",
        "value": value,
        "molecule_context": "transcript",
        "structural_type": "SO:0002073",
        "vrs_ref_allele_seq": "T",
        "gene_context": braf_gene_context
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope='module')
def nc_000007_silent_mutation(braf_gene_context, braf_nuc_value):
    """Create test fixture for NC_000007.13:g.140453136=."""
    value = copy.deepcopy(braf_nuc_value)
    value['state']['sequence'] = 'T'
    params = {
        "id": 'normalize.variation:NC_000007.13%3Ag.140453136%3D',
        "type": "VariationDescriptor",
        "value_id": "ga4gh:VA.bVNMOANetNE2z4PZ1j0DmwUL1rULmqkN",
        "value": value,
        "molecule_context": "transcript",
        "structural_type": "SO:0002073",
        "vrs_ref_allele_seq": "T",
        "gene_context": braf_gene_context
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope='module')
def amino_acid_delins(egfr_context):
    """Create test fixture for amino acid delins."""
    params = {
        "id": 'normalize.variation:NP_001333827.1%3Ap.Leu747_Thr751delinsPro',
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
        "molecule_context": "protein",
        "structural_type": "SO:1000032",
        "vrs_ref_allele_seq": "LREAT",
        "gene_context": egfr_context
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope='module')
def amino_acid_deletion_np_range(erbb2_context):
    """Create test fixture for amino acid deletion using NP accession and
    range for deletion.
    """
    params = {
        "id": 'normalize.variation:NP_004439.2%3Ap.Leu755_Thr759del',
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
        "molecule_context": "protein",
        "structural_type": "SO:0001604",
        "vrs_ref_allele_seq": "LRENT",
        "gene_context": erbb2_context
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope='module')
def coding_dna_deletion(erbb2_context):
    """Create test fixture for coding dna deletion range with deleted
    sequence.
    """
    params = {
        "id": 'normalize.variation:NM_004448.3%3Ac.2263_2277delTTGAGGGAAAACACA',  # noqa: E501
        "type": "VariationDescriptor",
        "value_id": "ga4gh:VA.FqrZbBlsAwpXWOXiBq2glfhCvLqp4xLC",
        "value": {
            "location": {
                "interval": {
                    "end": 2453,
                    "start": 2437,
                    "type": "SimpleInterval"
                },
                "sequence_id": "ga4gh:SQ.y9b4LVMiCXpZxOg9Xt1NwRtssA03MwWM",
                "type": "SequenceLocation"
            },
            "state": {
                "sequence": "T",
                "type": "SequenceState"
            },
            "type": "Allele"
        },
        "molecule_context": "transcript",
        "structural_type": "SO:0000159",
        "vrs_ref_allele_seq": "TTGAGGGAAAACACAT",
        "gene_context": erbb2_context
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope='module')
def genomic_deletion(vhl_gene_context):
    """Create test fixture for genomic deletion range."""
    params = {
        "id": 'normalize.variation:NC_000003.11%3Ag.10188279_10188297del',
        "type": "VariationDescriptor",
        "value_id": "ga4gh:VA.uagNswLQY5rgN2c30_J3-45UMpIySM4C",
        "value": {
            "location": {
                "interval": {
                    "end": 510,
                    "start": 491,
                    "type": "SimpleInterval"
                },
                "sequence_id": "ga4gh:SQ.xBKOKptLLDr-k4hTyCetvARn16pDS_rW",
                "type": "SequenceLocation"
            },
            "state": {
                "sequence": "",
                "type": "SequenceState"
            },
            "type": "Allele"
        },
        "molecule_context": "transcript",
        "structural_type": "SO:0000159",
        "vrs_ref_allele_seq": "ATGTTGACGGACAGCCTAT",
        "gene_context": vhl_gene_context
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope='module')
def amino_acid_insertion(egfr_context):
    """Create test fixture for NP amino acid insertion."""
    params = {
        "id": 'normalize.variation:NP_005219.2%3Ap.Asp770_Asn771insGlyLeu',
        "type": "VariationDescriptor",
        "value_id": "ga4gh:VA.d3dLtsYaLYE2Yh_GENvPUtTVZWlwLnJw",
        "value": {
            "location": {
                "interval": {
                    "end": 770,
                    "start": 770,
                    "type": "SimpleInterval"
                },
                "sequence_id": "ga4gh:SQ.vyo55F6mA6n2LgN4cagcdRzOuh38V4mE",
                "type": "SequenceLocation"
            },
            "state": {
                "sequence": "GL",
                "type": "SequenceState"
            },
            "type": "Allele"
        },
        "molecule_context": "protein",
        "structural_type": "SO:0001605",
        "gene_context": egfr_context
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope='module')
def coding_dna_insertion(limk2_gene_context):
    """Create test fixture for coding DNA insertion."""
    params = {
        "id": 'normalize.variation:ENST00000331728.9%3Ac.2049_2050insA',
        "type": "VariationDescriptor",
        "value_id": "ga4gh:VA.195Sg1AkyM4uQOhxLhBhANe2BUbnbEcR",
        "value": {
            "location": {
                "interval": {
                    "end": 2160,
                    "start": 2160,
                    "type": "SimpleInterval"
                },
                "sequence_id": "ga4gh:SQ.7_mlQyDN-uWH0RlxTQFvFEv6ykd2D-xF",
                "type": "SequenceLocation"
            },
            "state": {
                "sequence": "A",
                "type": "SequenceState"
            },
            "type": "Allele"
        },
        "molecule_context": "transcript",
        "structural_type": "SO:0000667",
        "gene_context": limk2_gene_context
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope='module')
def genomic_insertion(erbb2_context):
    """Create a gene insertion test fixture."""
    params = {
        "id": 'normalize.variation:NC_000017.10%3Ag.37880993_37880994insGCTTACGTGATG',  # noqa: E501
        "type": "VariationDescriptor",
        "value_id": "ga4gh:VA.qk5UNMFwxqQQWjO6CGMk3tryHBN3Sm_P",
        "value": {
            "location": {
                "interval": {
                    "end": 2500,
                    "start": 2488,
                    "type": "SimpleInterval"
                },
                "sequence_id": "ga4gh:SQ.y9b4LVMiCXpZxOg9Xt1NwRtssA03MwWM",
                "type": "SequenceLocation"
            },
            "state": {
                "sequence": "TACGTGATGGCTTACGTGATGGCT",
                "type": "SequenceState"
            },
            "type": "Allele"
        },
        "molecule_context": "transcript",
        "structural_type": "SO:0000667",
        "vrs_ref_allele_seq": "TACGTGATGGCT",
        "gene_context": erbb2_context
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope='module')
def genomic_substitution(egfr_context):
    """Create a gene insertion test fixture."""
    params = {
        "id": 'normalize.variation:NC_000007.13%3Ag.55249071C%3ET',
        "type": "VariationDescriptor",
        "value_id": "ga4gh:VA.c8ePmPEstWMCAJmtg3FuPb10XDr1G_8E",
        "value": {
            "location": {
                "interval": {
                    "end": 2630,
                    "start": 2629,
                    "type": "SimpleInterval"
                },
                "sequence_id": "ga4gh:SQ.d_QsP29RWJi6bac7GOC9cJ9AO7s_HUMN",
                "type": "SequenceLocation"
            },
            "state": {
                "sequence": "T",
                "type": "SequenceState"
            },
            "type": "Allele"
        },
        "molecule_context": "transcript",
        "structural_type": "SO:0001483",
        "vrs_ref_allele_seq": "C",
        "gene_context": egfr_context
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope='module')
def genomic_sub_grch38():
    """Create a genomic substitution GRCh38 test fixture."""
    params = {
        "id": 'normalize.variation:NC_000007.13%3Ag.55249071C%3ET',
        "type": "VariationDescriptor",
        "value_id": "ga4gh:VA.vWT6m5QcrdIJ37MfeQRsEO0avQiufIEx",
        "value": {
            "location": {
                "interval": {
                    "end": 55181378,
                    "start": 55181377,
                    "type": "SimpleInterval"
                },
                "sequence_id": "ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
                "type": "SequenceLocation"
            },
            "state": {
                "sequence": "T",
                "type": "SequenceState"
            },
            "type": "Allele"
        },
        "molecule_context": "genomic",
        "structural_type": "SO:0001483",
        "vrs_ref_allele_seq": "C"
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope='module')
def egfr_grch38_sub(genomic_sub_grch38, egfr_context):
    """Create a genomic substitution GRCh38 test fixture."""
    params = {
        "id": 'normalize.variation:NC_000007.13%3Ag.55249071C%3ET',
        "type": "VariationDescriptor",
        "value_id": "ga4gh:VA.vWT6m5QcrdIJ37MfeQRsEO0avQiufIEx",
        "value": {
            "location": {
                "interval": {
                    "end": 55181378,
                    "start": 55181377,
                    "type": "SimpleInterval"
                },
                "sequence_id": "ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
                "type": "SequenceLocation"
            },
            "state": {
                "sequence": "T",
                "type": "SequenceState"
            },
            "type": "Allele"
        },
        "molecule_context": "genomic",
        "structural_type": "SO:0001483",
        "vrs_ref_allele_seq": "C",
        "gene_context": egfr_context
    }
    return VariationDescriptor(**params)


def assertion_checks(normalize_response, test_variation):
    """Check that normalize_response and test_variation are equal."""
    assert normalize_response.id == test_variation.id
    assert normalize_response.type == test_variation.type
    assert normalize_response.value == test_variation.value
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
        assert resp_gene_context.value_id == test_variation_context.value_id
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


def test_amino_acid_substitution(test_normalize, braf_v600e):
    """Test that amino acid substitutions normalize correctly."""
    resp = test_normalize.normalize('     BRAF      V600E    ')
    assertion_checks(resp, braf_v600e)

    braf_id = "normalize.variation:BRAF%20V600E"

    resp = test_normalize.normalize('NP_004324.2:p.Val600Glu')
    assert resp.id == "normalize.variation:NP_004324.2%3Ap.Val600Glu"
    resp.id = braf_id
    assertion_checks(resp, braf_v600e)

    resp = test_normalize.normalize('braf v512e')
    assert resp.id == 'normalize.variation:braf%20v512e'
    resp.id = braf_id
    assertion_checks(resp, braf_v600e)

    resp = test_normalize.normalize(' NP_001365404.1:p.Val512Glu  ')
    assert resp.id == 'normalize.variation:NP_001365404.1%3Ap.Val512Glu'
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


def test_coding_dna_and_genomic_substitution(
        test_normalize, braf_v600e_nucleotide, genomic_substitution,
        genomic_sub_grch38, egfr_grch38_sub):
    """Test that coding dna and genomic substitutions normalize correctly."""
    resp = test_normalize.normalize('NM_004333.4:c.1799T>A')
    assertion_checks(resp, braf_v600e_nucleotide)

    # MANE transcript
    refseq_id = 'normalize.variation:NM_004333.4%3Ac.1799T%3EA'

    # TODO: Check if this should return a different VRS object?
    resp = test_normalize.normalize('ENST00000288602.10:c.1799T>A')
    assert resp.id == 'normalize.variation:ENST00000288602.10%3Ac.1799T%3EA'
    resp.id = refseq_id
    assertion_checks(resp, braf_v600e_nucleotide)

    resp = test_normalize.normalize('BRAF V600E c.1799T>A')
    assert resp.id == 'normalize.variation:BRAF%20V600E%20c.1799T%3EA'
    resp.id = refseq_id
    assertion_checks(resp, braf_v600e_nucleotide)

    resp = test_normalize.normalize('BRAF V600E (c.1799T>A)')
    assert resp.id == 'normalize.variation:BRAF%20V600E%20%28c.1799T%3EA%29'
    resp.id = refseq_id
    assertion_checks(resp, braf_v600e_nucleotide)

    resp = test_normalize.normalize('BRAF c.1799T>A')
    assert resp.id == 'normalize.variation:BRAF%20c.1799T%3EA'
    resp.id = refseq_id
    assertion_checks(resp, braf_v600e_nucleotide)

    resp = test_normalize.normalize('NC_000007.13:g.140453136A>T')
    assert resp.id == 'normalize.variation:NC_000007.13%3Ag.140453136A%3ET'
    resp.id = refseq_id
    assertion_checks(resp, braf_v600e_nucleotide)

    # TODO: Issue 99
    resp = test_normalize.normalize('BRAF V600E (g.140453136A>T)')
    assert resp.id == 'normalize.variation:BRAF%20V600E%20%28g.140453136A%3ET%29'  # noqa: E501
    resp.id = refseq_id
    assertion_checks(resp, braf_v600e_nucleotide)

    resp = test_normalize.normalize('BRAF g.140453136A>T')
    assert resp.id == 'normalize.variation:BRAF%20g.140453136A%3ET'
    resp.id = refseq_id
    assertion_checks(resp, braf_v600e_nucleotide)

    # More than 1 gene (EGFR and EGFR-AS1)
    resp = test_normalize.normalize('NC_000007.13:g.55249071C>T')
    assertion_checks(resp, genomic_sub_grch38)

    resp = test_normalize.normalize('EGFR g.55249071C>T')
    assert resp.id == 'normalize.variation:EGFR%20g.55249071C%3ET'
    resp.id = 'normalize.variation:NC_000007.13%3Ag.55249071C%3ET'
    assertion_checks(resp, genomic_substitution)


def test_coding_dna_silent_mutation(test_normalize,
                                    coding_dna_silent_mutation,
                                    braf_gene_context):
    """Test that Coding DNA Silent Mutation normalizes correctly."""
    resp = test_normalize.normalize('NM_004333.4:c.1799= ')
    assertion_checks(resp, coding_dna_silent_mutation)

    fixture_id = 'normalize.variation:NM_004333.4%3Ac.1799%3D'

    resp = test_normalize.normalize('ENST00000288602.11:c.1799=')
    assert resp.id == 'normalize.variation:ENST00000288602.11%3Ac.1799%3D'
    resp.id = fixture_id
    assertion_checks(resp, coding_dna_silent_mutation)

    # TODO: What to do for older Ensembl transcripts that aren't found
    #  in seqrepo or UTA
    # resp = test_normalize.normalize('ENST00000288602.6:c.1799=')
    # assert_coding_dna_genomic_silent_mutation(resp, braf_gene_context,
    #                                           1798, 1799)
    # assert resp.id == 'normalize.variation:ENST00000288602.6%3Ac.1799%3D'
    # assert resp.label == 'ENST00000288602.6:c.1799='
    # assert resp.molecule_context == 'transcript'

    resp = test_normalize.normalize('BRAF    c.1799=')
    assert resp.id == 'normalize.variation:BRAF%20c.1799%3D'
    resp.id = fixture_id
    assertion_checks(resp, coding_dna_silent_mutation)

    resp = test_normalize.normalize('  BRAF  V600E  c.1799=  ')
    assert resp.id == 'normalize.variation:BRAF%20V600E%20c.1799%3D'
    resp.id = fixture_id
    assertion_checks(resp, coding_dna_silent_mutation)


def test_genomic_silent_mutation(test_normalize, nc_000007_silent_mutation,
                                 braf_gene_context):
    """Test that genomic silent mutation normalizes correctly."""
    resp = test_normalize.normalize('NC_000007.13:g.140453136=')
    assertion_checks(resp, nc_000007_silent_mutation)

    resp = test_normalize.normalize('BRAF g.140453136=')
    assert resp.id == 'normalize.variation:BRAF%20g.140453136%3D'
    resp.id = 'normalize.variation:NC_000007.13%3Ag.140453136%3D'
    assertion_checks(resp, nc_000007_silent_mutation)


def test_coding_dna_delins(test_normalize, nm_004448_coding_dna_delins,
                           nm_000551):
    """Test that Coding DNA DelIns normalizes correctly."""
    resp = test_normalize.normalize('    NM_004448.4:c.2326_2327delinsCT    ')
    assertion_checks(resp, nm_004448_coding_dna_delins)

    # TODO: Test ENST###.c

    resp = test_normalize.normalize('NM_000551.3:c.615delinsAA')
    nm_000551.id = 'normalize.variation:NM_000551.3%3Ac.615delinsAA'
    assertion_checks(resp, nm_000551)


def test_genomic_delins(test_normalize, nc_000007_genomic_delins,
                        nm_000551):
    """Test that Genomic DelIns normalizes correctly."""
    resp = test_normalize.normalize(
        'NC_000007.13:g.140453135_140453136delinsAT'
    )
    assertion_checks(resp, nc_000007_genomic_delins)

    resp = test_normalize.normalize('NC_000003.12:g.10149938delinsAA')
    nm_000551.id = 'normalize.variation:NC_000003.12%3Ag.10149938delinsAA'
    assertion_checks(resp, nm_000551)


def test_amino_acid_delins(test_normalize, amino_acid_delins):
    """Test that Amnio Acid DelIns normalizes correctly."""
    resp = test_normalize.normalize('NP_001333827.1:p.Leu747_Thr751delinsPro')
    assertion_checks(resp, amino_acid_delins)

    resp = test_normalize.normalize('EGFR p.Leu747_Thr751delinsPro')
    assert resp.id == 'normalize.variation:EGFR%20p.Leu747_Thr751delinsPro'
    resp.id = 'normalize.variation:NP_001333827.1%3Ap.Leu747_Thr751delinsPro'
    assertion_checks(resp, amino_acid_delins)

    resp = test_normalize.normalize('EGFR Leu747_Thr751delinsPro')
    assert resp.id == 'normalize.variation:EGFR%20Leu747_Thr751delinsPro'
    resp.id = 'normalize.variation:NP_001333827.1%3Ap.Leu747_Thr751delinsPro'
    assertion_checks(resp, amino_acid_delins)

    resp = test_normalize.normalize('EGFR L747_T751delinsP')
    assert resp.id == 'normalize.variation:EGFR%20L747_T751delinsP'
    resp.id = 'normalize.variation:NP_001333827.1%3Ap.Leu747_Thr751delinsPro'
    assertion_checks(resp, amino_acid_delins)


def test_amino_acid_deletion(test_normalize, amino_acid_deletion_np_range):
    """Test that Amino Acid Deletion normalizes correctly."""
    resp = test_normalize.normalize('NP_004439.2:p.Leu755_Thr759del')
    assertion_checks(resp, amino_acid_deletion_np_range)

    resp = test_normalize.normalize('ERBB2 p.Leu755_Thr759del')
    assert resp.id == 'normalize.variation:ERBB2%20p.Leu755_Thr759del'

    resp = test_normalize.normalize('ERBB2 Leu755_Thr759del')
    assert resp.id == 'normalize.variation:ERBB2%20Leu755_Thr759del'


def test_coding_dna_deletion(test_normalize, coding_dna_deletion):
    """Test that coding dna deletion normalizes correctly."""
    resp = \
        test_normalize.normalize('NM_004448.3:c.2263_2277delTTGAGGGAAAACACA')
    assertion_checks(resp, coding_dna_deletion)

    resp = test_normalize.normalize('ERBB2 c.2263_2277delTTGAGGGAAAACACA')
    assert resp.id == \
           'normalize.variation:ERBB2%20c.2263_2277delTTGAGGGAAAACACA'
    resp.id = 'normalize.variation:NM_004448.3%3Ac.2263_2277delTTGAGGGAAAACACA'
    assertion_checks(resp, coding_dna_deletion)


def test_genomic_deletion(test_normalize, genomic_deletion):
    """Test that genomic deletion normalizes correctly."""
    resp = test_normalize.normalize('NC_000003.11:g.10188279_10188297del')
    assertion_checks(resp, genomic_deletion)

    resp = test_normalize.normalize('VHL g.10188279_10188297del')
    assert resp.id == 'normalize.variation:VHL%20g.10188279_10188297del'
    resp.id = 'normalize.variation:NC_000003.11%3Ag.10188279_10188297del'
    assertion_checks(resp, genomic_deletion)


def test_amino_acid_insertion(test_normalize, amino_acid_insertion):
    """Test that amino acid insertion normalizes correctly."""
    resp = test_normalize.normalize('NP_005219.2:p.Asp770_Asn771insGlyLeu')
    assertion_checks(resp, amino_acid_insertion)

    def change_resp(response):
        fixture_id = \
            'normalize.variation:NP_005219.2%3Ap.Asp770_Asn771insGlyLeu'
        response.id = fixture_id

    resp = test_normalize.normalize('EGFR D770_N771insGL')
    assert resp.id == 'normalize.variation:EGFR%20D770_N771insGL'
    change_resp(resp)
    assertion_checks(resp, amino_acid_insertion)

    resp = test_normalize.normalize('EGFR p.D770_N771insGL')
    assert resp.id == 'normalize.variation:EGFR%20p.D770_N771insGL'
    change_resp(resp)
    assertion_checks(resp, amino_acid_insertion)

    resp = test_normalize.normalize('EGFR Asp770_Asn771insGlyLeu')
    assert resp.id == 'normalize.variation:EGFR%20Asp770_Asn771insGlyLeu'
    change_resp(resp)
    assertion_checks(resp, amino_acid_insertion)

    resp = test_normalize.normalize('EGFR p.Asp770_Asn771insGlyLeu')
    assert resp.id == 'normalize.variation:EGFR%20p.Asp770_Asn771insGlyLeu'
    change_resp(resp)
    assertion_checks(resp, amino_acid_insertion)


def test_coding_dna_insertion(test_normalize, coding_dna_insertion):
    """Test that coding dna insertion normalizes correctly."""
    resp = test_normalize.normalize('ENST00000331728.9:c.2049_2050insA')
    assertion_checks(resp, coding_dna_insertion)

    resp = test_normalize.normalize('LIMK2 c.2049_2050insA')
    assert resp.id == 'normalize.variation:LIMK2%20c.2049_2050insA'
    resp.id = 'normalize.variation:ENST00000331728.9%3Ac.2049_2050insA'
    assertion_checks(resp, coding_dna_insertion)


def test_genomic_insertion(test_normalize, genomic_insertion):
    """Test that genomic insertion normalizes correctly."""
    resp = test_normalize.normalize('NC_000017.10:g.37880993_37880994insGCTTACGTGATG')  # noqa: E501
    assertion_checks(resp, genomic_insertion)

    resp = test_normalize.normalize('ERBB2 g.37880993_37880994insGCTTACGTGATG')
    assert resp.id ==\
           'normalize.variation:ERBB2%20g.37880993_37880994insGCTTACGTGATG'
    resp.id = \
        'normalize.variation:NC_000017.10%3Ag.37880993_37880994insGCTTACGTGATG'
    assertion_checks(resp, genomic_insertion)


def test_no_matches(test_normalize):
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
        resp = test_normalize.normalize(q)
        assert resp.type == 'VariationDescriptor'
        assert resp.value.type == 'Text'

    resp = test_normalize.normalize('clinvar:10')
    assert resp.type == 'VariationDescriptor'
    assert resp.value.definition == 'clinvar:10'

    resp = test_normalize.normalize('   ')
    assert resp is None

    resp = test_normalize.normalize('')
    assert resp is None

    resp = test_normalize.normalize(None)
    assert resp is None


def test_service_meta():
    """Test that service meta info populates correctly."""
    response = normalize_get_response('BRAF v600e')
    service_meta = response.service_meta_
    assert service_meta.name == "variation-normalizer"
    assert service_meta.version
    assert isinstance(service_meta.response_datetime, datetime)
    assert service_meta.url == 'https://github.com/cancervariants/variation-normalization'  # noqa: E501

    response = normalize_get_response('this-wont-normalize')
    service_meta = response.service_meta_
    assert service_meta.name == "variation-normalizer"
    assert service_meta.version
    assert isinstance(service_meta.response_datetime, datetime)
    assert service_meta.url == 'https://github.com/cancervariants/variation-normalization'  # noqa: E501

    response = to_vrs_get_response('BRAF v600e')
    service_meta = response.service_meta_
    assert service_meta.name == "variation-normalizer"
    assert service_meta.version
    assert isinstance(service_meta.response_datetime, datetime)
    assert service_meta.url == 'https://github.com/cancervariants/variation-normalization'  # noqa: E501

    response = to_vrs_get_response('this-wont-normalize')
    service_meta = response.service_meta_
    assert service_meta.name == "variation-normalizer"
    assert service_meta.version
    assert isinstance(service_meta.response_datetime, datetime)
    assert service_meta.url == 'https://github.com/cancervariants/variation-normalization'  # noqa: E501
