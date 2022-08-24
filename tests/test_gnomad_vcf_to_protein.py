"""Module for testing gnomad_vcf_to_protein works correctly"""
import pytest
from ga4gh.vrsatile.pydantic.vrsatile_models import VariationDescriptor

from tests.conftest import assertion_checks


@pytest.fixture(scope="module")
def test_handler(test_query_handler):
    """Create test fixture for gnomad vcf to protein handler"""
    return test_query_handler.gnomad_vcf_to_protein_handler


@pytest.fixture(scope="module")
def mmel1_gene_context():
    """Create MMEL1 gene context test fixture"""
    return {
        "id": "normalize.gene:HGNC%3A14668",
        "type": "GeneDescriptor",
        "label": "MMEL1",
        "xrefs": [
            "ensembl:ENSG00000142606",
            "ncbigene:79258"
        ],
        "alternate_labels": [
            "SEP",
            "NL1",
            "MMEL2",
            "MMEL1",
            "MELL1",
            "NL2",
            "NEPII",
            "NEP2"
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
                "value": "membrane metalloendopeptidase like 1"
            },
            {
                "type": "Extension",
                "name": "chromosome_location",
                "value": {
                    "species_id": "taxonomy:9606",
                    "start": "p36.32",
                    "end": "p36.32",
                    "id": "ga4gh:CL.GvnxCN3xVRhtAGEabC3Tl1p3wsHM7-J7",
                    "type": "ChromosomeLocation",
                    "chr": "1"
                }
            },
            {
                "type": "Extension",
                "name": "associated_with",
                "value": [
                    "uniprot:Q495T6",
                    "ccds:CCDS30569",
                    "ucsc:uc001ajy.3",
                    "merops:M13.008",
                    "pubmed:11560781",
                    "pubmed:18539150",
                    "vega:OTTHUMG00000000846",
                    "orphanet:239881",
                    "refseq:NM_033467",
                    "ena.embl:AF336981",
                    "omim:618104"
                ]
            },
            {
                "type": "Extension",
                "name": "previous_symbols",
                "value": [
                    "MMEL2",
                    "MELL1",
                    "MMEL1"
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
        "gene_id": "hgnc:14668"
    }


@pytest.fixture(scope="module")
def cdk11a_gene_context():
    """Create test fixture for CDK11A"""
    return {
        "id": "normalize.gene:HGNC%3A1730",
        "type": "GeneDescriptor",
        "label": "CDK11A",
        "xrefs": [
            "ensembl:ENSG00000008128",
            "ncbigene:728642"
        ],
        "alternate_labels": [
            "LOC100294398",
            "LOC100134433",
            "CDC2L3",
            "CDK11-p46",
            "CDC2L2",
            "CDK11-p58",
            "CDK11-p110",
            "p58GTA",
            "PITSLRE"
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
                "value": "cyclin dependent kinase 11A"
            },
            {
                "type": "Extension",
                "name": "chromosome_location",
                "value": {
                    "species_id": "taxonomy:9606",
                    "start": "p36.33",
                    "end": "p36.33",
                    "id": "ga4gh:CL.08UTaK9ZRrooFYQRV-WD5k2Z4iOyTmMG",
                    "type": "ChromosomeLocation",
                    "chr": "1"
                }
            },
            {
                "type": "Extension",
                "name": "associated_with",
                "value": [
                    "pubmed:9750192",
                    "pubmed:19884882",
                    "vega:OTTHUMG00000000703",
                    "ccds:CCDS44043",
                    "ccds:CCDS81254",
                    "ccds:CCDS44042",
                    "iuphar:1963",
                    "refseq:NM_024011",
                    "ucsc:uc009vks.4",
                    "ena.embl:AF067522",
                    "omim:116951",
                    "pubmed:7920654",
                    "uniprot:Q9UQ88",
                    "ccds:CCDS81253"
                ]
            },
            {
                "type": "Extension",
                "name": "previous_symbols",
                "value": [
                    "LOC100294398",
                    "LOC100134433",
                    "CDC2L2",
                    "CDC2L3"
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
        "gene_id": "hgnc:1730"
    }


@pytest.fixture(scope="module")
def lrp8_gene_context():
    """Create LRP8 gene context test fixture"""
    return {
        "id": "normalize.gene:HGNC%3A6700",
        "type": "GeneDescriptor",
        "label": "LRP8",
        "xrefs": [
            "ensembl:ENSG00000157193",
            "ncbigene:7804"
        ],
        "alternate_labels": [
            "APOER2",
            "LRP-8",
            "HSZ75190",
            "LRP8",
            "MIPS",
            "MCI1"
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
                "value": "LDL receptor related protein 8"
            },
            {
                "type": "Extension",
                "name": "chromosome_location",
                "value": {
                    "species_id": "taxonomy:9606",
                    "start": "p32.3",
                    "end": "p32.3",
                    "id": "ga4gh:CL.u-8S9gYnhTxY-5dMXCjo_ZMuDMVbbobs",
                    "type": "ChromosomeLocation",
                    "chr": "1"
                }
            },
            {
                "type": "Extension",
                "name": "associated_with",
                "value": [
                    "ucsc:uc001cvi.4",
                    "ccds:CCDS30720",
                    "ccds:CCDS579",
                    "uniprot:Q14114",
                    "ccds:CCDS580",
                    "vega:OTTHUMG00000008924",
                    "pubmed:9079678",
                    "pubmed:8626535",
                    "ccds:CCDS578",
                    "refseq:NM_004631",
                    "ena.embl:D50678",
                    "omim:602600"
                ]
            },
            {
                "type": "Extension",
                "name": "previous_symbols",
                "value": [
                    "HSZ75190",
                    "MIPS",
                    "LRP8"
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
        "gene_id": "hgnc:6700"
    }


@pytest.fixture(scope="module")
def atad3a_gene_context():
    """Create ATAD3A gene context test fixture"""
    return {
        "id": "normalize.gene:HGNC%3A25567",
        "type": "GeneDescriptor",
        "label": "ATAD3A",
        "xrefs": [
            "ensembl:ENSG00000197785",
            "ncbigene:55210"
        ],
        "alternate_labels": [
            "HAYOS",
            "FLJ10709",
            "PHRINL"
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
                "value": "ATPase family AAA domain containing 3A"
            },
            {
                "type": "Extension",
                "name": "chromosome_location",
                "value": {
                    "species_id": "taxonomy:9606",
                    "start": "p36.33",
                    "end": "p36.33",
                    "id": "ga4gh:CL.08UTaK9ZRrooFYQRV-WD5k2Z4iOyTmMG",
                    "type": "ChromosomeLocation",
                    "chr": "1"
                }
            },
            {
                "type": "Extension",
                "name": "associated_with",
                "value": [
                    "ccds:CCDS53259",
                    "ucsc:uc001aga.3",
                    "ena.embl:AK025865",
                    "vega:OTTHUMG00000000575",
                    "ccds:CCDS31",
                    "refseq:NM_018188",
                    "omim:612316",
                    "ccds:CCDS53260",
                    "orphanet:469968",
                    "pubmed:28158749",
                    "uniprot:Q9NVI7"
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
        "gene_id": "hgnc:25567"
    }


@pytest.fixture(scope="module")
def mmel1_l30m(mmel1_gene_context):
    """Create test fixture for MMEL1 L30M"""
    params = {
        "id": "normalize.variation:1-2629397-G-T",
        "type": "VariationDescriptor",
        "variation_id": "ga4gh:VA.TNZKL_HZqeU1anmA0vGWnvV_RjsNb1mB",
        "variation": {
            "id": "ga4gh:VA.TNZKL_HZqeU1anmA0vGWnvV_RjsNb1mB",
            "location": {
                "id": "ga4gh:SL.faQG3PE4Z3J2sCLd2mvOi5GtbcNeTs3K",
                "end": {"value": 30, "type": "Number"},
                "start": {"value": 29, "type": "Number"},
                "sequence_id": "ga4gh:SQ.iQ8F_pnsiQOLohiV2qh3OWRZiftUt8jZ",
                "type": "SequenceLocation"
            },
            "state": {
                "sequence": "M",
                "type": "LiteralSequenceExpression"
            },
            "type": "Allele"
        },
        "molecule_context": "protein",
        "structural_type": "SO:0001606",
        "vrs_ref_allele_seq": "L",
        "gene_context": mmel1_gene_context
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def cdk11a_e314del(cdk11a_gene_context):
    """Create test fixture for CDK11A Glu314del"""
    params = {
        "id": "normalize.variation:1-1708855-TTCC-T",
        "type": "VariationDescriptor",
        "variation_id": "ga4gh:VA.onGUSj8yycu9JLUsygAcvY5T2eRZGeN-",
        "variation": {
            "id": "ga4gh:VA.onGUSj8yycu9JLUsygAcvY5T2eRZGeN-",
            "location": {
                "id": "ga4gh:SL.PfUIK14UgNsxxUrI_6FFv68SC_tvItrP",
                "end": {"value": 321, "type": "Number"},
                "start": {"value": 308, "type": "Number"},
                "sequence_id": "ga4gh:SQ.N728VSRRMHJ1SrhJgKqJOCaa3l5Z4sqm",
                "type": "SequenceLocation"
            },
            "state": {
                "sequence": "EEEEEEEEEEEE",
                "type": "LiteralSequenceExpression"
            },
            "type": "Allele"
        },
        "molecule_context": "protein",
        "structural_type": "SO:0001604",
        "vrs_ref_allele_seq": "EEEEEEEEEEEEE",
        "gene_context": cdk11a_gene_context
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def protein_insertion2(lrp8_gene_context):
    """Create test fixture for LRP8 p.Gln25_Leu26insArg"""
    params = {
        "id": "normalize.variation:1-53327836-A-ACGC",
        "type": "VariationDescriptor",
        "variation_id": "ga4gh:VA.pt-iZhrxYtI6FjWk8bxxUXWKsH_rFj2t",
        "variation": {
            "id": "ga4gh:VA.pt-iZhrxYtI6FjWk8bxxUXWKsH_rFj2t",
            "location": {
                "id": "ga4gh:SL.lfzCaOsWzXwPOpCa9boheo2N-OGLlsMV",
                "end": {"value": 25, "type": "Number"},
                "start": {"value": 25, "type": "Number"},
                "sequence_id": "ga4gh:SQ.qgIh8--4F6IpxRwX_lVtD2BhepH5B5Ef",
                "type": "SequenceLocation"
            },
            "state": {
                "sequence": "R",
                "type": "LiteralSequenceExpression"
            },
            "type": "Allele"
        },
        "molecule_context": "protein",
        "structural_type": "SO:0001605",
        "gene_context": lrp8_gene_context
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def atad3a_loc():
    """Create test fixture for ATAD3A location"""
    return {
        "id": "ga4gh:SL.gPvwSGKYW_3l4spkahhv6nA4VO3Oe9ye",
        "end": {"value": 7, "type": "Number"},
        "start": {"value": 6, "type": "Number"},
        "sequence_id": "ga4gh:SQ.MHPOY_7fv8V9SktyvaTxulVFSK6XCxM8",
        "type": "SequenceLocation"
    }


@pytest.fixture(scope="module")
def atad3a_i7v(atad3a_gene_context, atad3a_loc):
    """Create test fixture for ATAD3A Ile7Val"""
    params = {
        "id": "normalize.variation:1-1512287-A-G",
        "type": "VariationDescriptor",
        "variation_id": "ga4gh:VA.0n_eZ3A0nCDSf3D8tBPsmg5-sOU2eYS1",
        "variation": {
            "id": "ga4gh:VA.0n_eZ3A0nCDSf3D8tBPsmg5-sOU2eYS1",
            "location": atad3a_loc,
            "state": {
                "sequence": "V",
                "type": "LiteralSequenceExpression"
            },
            "type": "Allele"
        },
        "molecule_context": "protein",
        "structural_type": "SO:0001606",
        "vrs_ref_allele_seq": "I",
        "gene_context": atad3a_gene_context
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def atad3a_i7t(atad3a_gene_context, atad3a_loc):
    """Create test fixture for ATAD3A Ile7Thr"""
    params = {
        "id": "normalize.variation:1-1512288-T-C",
        "type": "VariationDescriptor",
        "variation_id": "ga4gh:VA.xvlGODCInBKado5IkXueBbkHY248GQeG",
        "variation": {
            "id": "ga4gh:VA.xvlGODCInBKado5IkXueBbkHY248GQeG",
            "location": atad3a_loc,
            "state": {
                "sequence": "T",
                "type": "LiteralSequenceExpression"
            },
            "type": "Allele"
        },
        "molecule_context": "protein",
        "structural_type": "SO:0001606",
        "vrs_ref_allele_seq": "I",
        "gene_context": atad3a_gene_context
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def atad3a_i7m(atad3a_gene_context, atad3a_loc):
    """Create test fixture for ATAD3A Ile7Met"""
    params = {
        "id": "normalize.variation:1-1512289-T-G",
        "type": "VariationDescriptor",
        "variation_id": "ga4gh:VA.SB4A2_i6ldjbr02igJtvMBOMxVxw_371",
        "variation": {
            "id": "ga4gh:VA.SB4A2_i6ldjbr02igJtvMBOMxVxw_371",
            "location": atad3a_loc,
            "state": {
                "sequence": "M",
                "type": "LiteralSequenceExpression"
            },
            "type": "Allele"
        },
        "molecule_context": "protein",
        "structural_type": "SO:0001606",
        "vrs_ref_allele_seq": "I",
        "gene_context": atad3a_gene_context
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="session")
def braf_v600l(braf_gene_context, braf_600loc):
    """Create test fixture for BRAF Val600Leu."""
    params = {
        "id": "normalize.variation:7-140753337-C-A",
        "type": "VariationDescriptor",
        "variation_id": "ga4gh:VA.HVHeMnSZ0nIsCxqI9xE35RaiBhcSo9IS",
        "variation": {
            "id": "ga4gh:VA.HVHeMnSZ0nIsCxqI9xE35RaiBhcSo9IS",
            "location": braf_600loc,
            "state": {
                "sequence": "L",
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
def braf_600_silent_mutation(braf_gene_context, braf_600loc):
    """Create test fixture for BRAF Val600=."""
    params = {
        "id": "normalize.variation:7-140753335-C-A",
        "type": "VariationDescriptor",
        "variation_id": "ga4gh:VA.dbEtYnjAkO0rz_D2y_xaLusONuEW27Qd",
        "variation": {
            "id": "ga4gh:VA.dbEtYnjAkO0rz_D2y_xaLusONuEW27Qd",
            "location": braf_600loc,
            "state": {
                "sequence": "V",
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


@pytest.fixture(scope="module")
def kras_g12d():
    """Fixture for KRAS G12C"""
    return {
        "id": "ga4gh:VA.jSiuIG2uZEbml8MfCtZi3J9jPBCBiaPw",
        "type": "Allele",
        "location": {
            "id": "ga4gh:SL.NZ2azyPR9iXnUHx3Kv2Eam1DeprkcoQR",
            "type": "SequenceLocation",
            "sequence_id": "ga4gh:SQ.fytWhQSNGnA-86vDiQCxTSzgkk_WfQRS",
            "start": {"type": "Number", "value": 11},
            "end": {"type": "Number", "value": 12}
        },
        "state": {
            "type": "LiteralSequenceExpression",
            "sequence": "D"
        }
    }


@pytest.mark.asyncio
async def test_substitution(test_handler, braf_v600e, braf_v600l,
                            braf_600_silent_mutation, mmel1_l30m, atad3a_i7v,
                            atad3a_i7t, atad3a_i7m, kras_g12d):
    """Test that substitution queries return correct response"""
    # Reading Frame 1, Negative Strand
    resp = await test_handler.gnomad_vcf_to_protein("7-140753337-C-A")
    assertion_checks(resp.variation_descriptor, braf_v600l, "7-140753337-C-A",
                     ignore_id=True)
    assert resp.warnings == []

    # Reading Frame 2, Negative Strand
    resp = await test_handler.gnomad_vcf_to_protein("7-140753336-A-T")
    assertion_checks(resp.variation_descriptor, braf_v600e, "7-140753336-A-T",
                     ignore_id=True)
    assert resp.warnings == []

    # Reading Frame 3, Negative Strand
    resp = await test_handler.gnomad_vcf_to_protein("7-140753335-C-A")
    assertion_checks(resp.variation_descriptor, braf_600_silent_mutation,
                     "7-140753335-C-A", ignore_id=True)
    assert resp.warnings == []

    # Reading Frame 3, Negative Strand
    resp = await test_handler.gnomad_vcf_to_protein("1-2629397-G-T")
    assertion_checks(resp.variation_descriptor, mmel1_l30m, "1-2629397-G-T")
    assert resp.warnings == []

    # Reading Frame 1, Positive Strand
    resp = await test_handler.gnomad_vcf_to_protein("1-1512287-A-G")
    assertion_checks(resp.variation_descriptor, atad3a_i7v, "1-1512287-A-G")
    assert resp.warnings == []

    # Reading Frame 2, Positive Strand
    resp = await test_handler.gnomad_vcf_to_protein("1-1512288-T-C")
    assertion_checks(resp.variation_descriptor, atad3a_i7t, "1-1512288-T-C")
    assert resp.warnings == []

    # Reading Frame 3, Positive Strand
    resp = await test_handler.gnomad_vcf_to_protein("1-1512289-T-G")
    assertion_checks(resp.variation_descriptor, atad3a_i7m, "1-1512289-T-G",
                     ignore_id=True)
    assert resp.warnings == []

    resp = await test_handler.gnomad_vcf_to_protein("12-25245350-C-T")
    assert resp
    variation = resp.variation_descriptor.dict()["variation"]
    assert variation == kras_g12d


@pytest.mark.asyncio
async def test_silent_mutation(test_handler, vhl_silent):
    """Test that silent queries return correct response"""
    resp = await test_handler.gnomad_vcf_to_protein("3-10183714-C-C")
    assertion_checks(resp.variation_descriptor, vhl_silent, "3-10183714-C-C",
                     ignore_id=True)
    assert resp.warnings == []


@pytest.mark.asyncio
async def test_insertion(test_handler, protein_insertion,
                         protein_insertion2):
    """Test that insertion queries return correct response"""
    resp = await test_handler.gnomad_vcf_to_protein("7-55181319-C-CGGGTTG")
    assertion_checks(resp.variation_descriptor, protein_insertion,
                     "7-55181319-C-CGGGTTG", ignore_id=True)
    assert resp.warnings == []

    resp = await test_handler.gnomad_vcf_to_protein("1-53327836-A-ACGC")
    assertion_checks(resp.variation_descriptor, protein_insertion2, "1-53327836-A-ACGC")
    assert resp.warnings == []


@pytest.mark.asyncio
async def test_deletion(test_handler, protein_deletion_np_range,
                        cdk11a_e314del):
    """Test that deletion queries return correct response"""
    resp = await test_handler.gnomad_vcf_to_protein(
        "17-39723966-TTGAGGGAAAACACAT-T")
    assertion_checks(resp.variation_descriptor, protein_deletion_np_range,
                     "17-39723966-TTGAGGGAAAACACAT-T", ignore_id=True)
    assert resp.warnings == []

    resp = await test_handler.gnomad_vcf_to_protein("1-1708855-TTCC-T")
    assertion_checks(resp.variation_descriptor, cdk11a_e314del, "1-1708855-TTCC-T")
    assert resp.warnings == []


@pytest.mark.asyncio
async def test_invalid(test_handler):
    """Test that invalid queries return correct response"""
    resp = await test_handler.gnomad_vcf_to_protein("dummy")
    assert resp.variation_descriptor is None

    resp = await test_handler.gnomad_vcf_to_protein("BRAF V600E",
                                                    untranslatable_returns_text=True)
    assert resp.variation_descriptor.variation.type == "Text"
    assert resp.variation_descriptor.label == "BRAF V600E"
    assert resp.warnings == ["BRAF V600E is not a supported gnomad vcf query"]

    resp = await test_handler.gnomad_vcf_to_protein("7-140753336-T-G",
                                                    untranslatable_returns_text=True)
    assert resp.variation_descriptor.variation.type == "Text"
    assert resp.variation_descriptor.label == "7-140753336-T-G"
    assert set(resp.warnings) == {
        "Unable to get transcripts given NC_000007.13 and 140753336",
        "Expected T but found A on NC_000007.14 at position 140753336",
    }

    resp = await test_handler.gnomad_vcf_to_protein("20-2-TC-TG",
                                                    untranslatable_returns_text=True)
    assert resp.variation_descriptor.variation.type == "Text"
    assert resp.variation_descriptor.label == "20-2-TC-TG"
    assert resp.warnings == ["Unable to get protein variation for 20-2-TC-TG"]
