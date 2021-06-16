"""Module for testing MANE Transcript Mapping class."""
import pytest
from variation.data_sources import MANETranscriptMappings


@pytest.fixture(scope='module')
def test_mane_transcript_mappings():
    """Build MANE transcript mappings test fixture."""
    class TestMANETranscriptMappings:

        def __init__(self):
            self.test_mane_transcript_mappings = MANETranscriptMappings()
            self.test_df = \
                self.test_mane_transcript_mappings._load_mane_transcript_data()

        def get_row_from_gene(self, gene_symbol):
            return self.test_mane_transcript_mappings.get_gene_mane_data(
                gene_symbol)

    return TestMANETranscriptMappings()


@pytest.fixture(scope='module')
def braf():
    """Create test fixture for BRAF MANE Transcript data."""
    return {
        "#NCBI_GeneID": "GeneID:673",
        "Ensembl_Gene": "ENSG00000157764.14",
        "HGNC_ID": "HGNC:1097",
        "symbol": "BRAF",
        "name": "B-Raf proto-oncogene, serine/threonine kinase",
        "RefSeq_nuc": "NM_001374258.1",
        "RefSeq_prot": "NP_001361187.1",
        "Ensembl_nuc": "ENST00000644969.2",
        "Ensembl_prot": "ENSP00000496776.1",
        "MANE_status": "MANE Select",
        "GRCh38_chr": "7",
        "chr_start": 140719337,
        "chr_end": 140924929,
        "chr_strand": "-"
    }


@pytest.fixture(scope='module')
def ercc6_plus_clinical():
    """Create test fixture for ERCC6 MANE Plus Clinical Transcript data."""
    return {
        "#NCBI_GeneID": "GeneID:2074",
        "Ensembl_Gene": "ENSG00000225830.16",
        "HGNC_ID": "HGNC:3438",
        "symbol": "ERCC6",
        "name": "ERCC excision repair 6, chromatin remodeling factor",
        "RefSeq_nuc": "NM_001277058.2",
        "RefSeq_prot": "NP_001263987.1",
        "Ensembl_nuc": "ENST00000447839.7",
        "Ensembl_prot": "ENSP00000387966.2",
        "MANE_status": "MANE Plus Clinical",
        "GRCh38_chr": "10",
        "chr_start": 49515198,
        "chr_end": 49539121,
        "chr_strand": "-"
    }


@pytest.fixture(scope='module')
def ercc6_select():
    """Create test fixture for ERCC6 MANE Select Transcript data."""
    return {
        "#NCBI_GeneID": "GeneID:2074",
        "Ensembl_Gene": "ENSG00000225830.16",
        "HGNC_ID": "HGNC:3438",
        "symbol": "ERCC6",
        "name": "ERCC excision repair 6, chromatin remodeling factor",
        "RefSeq_nuc": "NM_000124.4",
        "RefSeq_prot": "NP_000115.1",
        "Ensembl_nuc": "ENST00000355832.10",
        "Ensembl_prot": "ENSP00000348089.5",
        "MANE_status": "MANE Select",
        "GRCh38_chr": "10",
        "chr_start": 49454470,
        "chr_end": 49539121,
        "chr_strand": "-"
    }


def test_get_gene_mane_data(test_mane_transcript_mappings, braf, ercc6_select,
                            ercc6_plus_clinical):
    """Test that get_gene_mane_data method works correctly."""
    # MANE Select
    actual = test_mane_transcript_mappings.get_row_from_gene('BRAF')
    assert len(actual) == 1
    actual = actual[0]
    assert actual == braf

    actual = test_mane_transcript_mappings.get_row_from_gene('braf')
    assert len(actual) == 1
    actual = actual[0]
    assert actual == braf

    # MANE Select and MANE Plus Clinical
    actual = test_mane_transcript_mappings.get_row_from_gene('ERCC6')
    assert len(actual) == 2
    assert actual[0] == ercc6_plus_clinical
    assert actual[1] == ercc6_select

    actual = test_mane_transcript_mappings.get_row_from_gene('ercc6')
    assert actual[0] == ercc6_plus_clinical
    assert actual[1] == ercc6_select

    # No Matches
    actual = test_mane_transcript_mappings.get_row_from_gene('BRAFF')
    assert actual is None

    actual = test_mane_transcript_mappings.get_row_from_gene('')
    assert actual is None
