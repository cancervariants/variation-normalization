"""Module for testing MANE Transcript class."""
import pytest
from variation.mane_transcript import MANETranscript
from variation.data_sources import TranscriptMappings, \
    MANETranscriptMappings, UTA, SeqRepoAccess
import hgvs.parser
import copy


@pytest.fixture(scope='module')
def test_mane_transcript():
    """Build mane transcript test fixture."""
    class TestMANETranscript:

        def __init__(self):
            self.test_mane_transcript = MANETranscript(
                SeqRepoAccess(), TranscriptMappings(),
                MANETranscriptMappings(), UTA()
            )
            self.test_hgvs_parser = hgvs.parser.Parser()

        def _get_reading_frame(self, pos):
            return self.test_mane_transcript._get_reading_frame(pos)

        def _p_to_c_pos(self, start, end):
            return self.test_mane_transcript._p_to_c_pos(start, end)

        def _p_to_c(self, ac, start_pos, end_pos):
            return self.test_mane_transcript._p_to_c(ac, start_pos, end_pos)

        def _c_to_g(self, ac, pos):
            return self.test_mane_transcript._c_to_g(ac, pos)

        def _g_to_mane_c(self, alt_tx_data, mane_data):
            return self.test_mane_transcript._g_to_mane_c(
                alt_tx_data, mane_data
            )

        def _get_mane_p(self, mane_data, mane_c_pos_range):
            return self.test_mane_transcript._get_mane_p(
                mane_data, mane_c_pos_range
            )

        def get_mane_transcript(self, ac, start_pos, end_pos,
                                start_annotation_layer):
            return self.test_mane_transcript.get_mane_transcript(
                ac, start_pos, end_pos, start_annotation_layer
            )

        def g_to_mane_c(self, ac, start_pos, end_pos):
            return self.test_mane_transcript.g_to_mane_c(ac, start_pos,
                                                         end_pos)

    return TestMANETranscript()


@pytest.fixture(scope='module')
def braf_mane_data():
    """Create test fixture for BRAF MANE data."""
    return {
        '#NCBI_GeneID': 'GeneID:673',
        'Ensembl_Gene': 'ENSG00000157764.14',
        'HGNC_ID': 'HGNC:1097',
        'symbol': 'BRAF',
        'name': 'B-Raf proto-oncogene, serine/threonine kinase',
        'RefSeq_nuc': 'NM_001374258.1',
        'RefSeq_prot': 'NP_001361187.1',
        'Ensembl_nuc': 'ENST00000644969.2',
        'Ensembl_prot': 'ENSP00000496776.1',
        'MANE_status': 'MANE Select',
        'GRCh38_chr': '7',
        'chr_start': 140719337,
        'chr_end': 140924929,
        'chr_strand': '-'
    }


@pytest.fixture(scope='module')
def nm_004333v6_g():
    """Create test fixture for NM_004333.6 genomic data."""
    return {
        'gene': 'BRAF',
        'tx_ac': 'NM_004333.6',
        'tx_pos_range': (1967, 2086),
        'alt_ac': 'NC_000007.14',
        'alt_pos_change_range': (140753331, 140753333),
        'alt_pos_range': (140753274, 140753393),
        'pos_change': (57, 60),
        'strand': '-',
        'alt_aln_method': 'splign',
        'tx_exon_id': 7649345,
        'alt_exon_id': 9507338,
        'coding_start_site': 226
    }


@pytest.fixture(scope='module')
def braf_v600e_mane_p():
    """Create test fixture for BRAF V600E MANE Transcript on p coordinate."""
    return {
        'refseq': 'NP_001361187.1',
        'ensembl': 'ENSP00000496776.1',
        'pos': (640, 640),
        'status': 'MANE Select',
        'strand': '-',
        'gene': 'BRAF'
    }


@pytest.fixture(scope='module')
def egfr_l858r_mane_p():
    """Create test fixture for EGFR L858R MANE Transcript on p coordinate."""
    return {
        'refseq': 'NP_005219.2',
        'ensembl': 'ENSP00000275493.2',
        'pos': (858, 858),
        'status': 'MANE Select',
        'strand': '+',
        'gene': 'EGFR'
    }


@pytest.fixture(scope='module')
def braf_v600e_mane_c():
    """Create test fixture for BRAF V600E MANE Transcript on c coordinate."""
    return {
        'refseq': 'NM_001374258.1',
        'ensembl': 'ENST00000644969.2',
        'pos': (1919, 1919),
        'status': 'MANE Select',
        'strand': '-',
        'coding_start_site': 226,
        'coding_end_site': 2650,
        'gene': 'BRAF'
    }


@pytest.fixture(scope='module')
def egfr_l858r_mane_c():
    """Create test fixture for EGFR L858R MANE Transcript on c coordinate."""
    return {
        'refseq': 'NM_005228.5',
        'ensembl': 'ENST00000275493.7',
        'pos': (2573, 2573),
        'status': 'MANE Select',
        'strand': '+',
        'coding_start_site': 261,
        'coding_end_site': 3894,
        'gene': 'EGFR'
    }


def test_get_reading_frame(test_mane_transcript):
    """Test that _get_reading_frame works correctly."""
    rf = test_mane_transcript._get_reading_frame(1797)
    assert rf == 3

    rf = test_mane_transcript._get_reading_frame(1798)
    assert rf == 1

    rf = test_mane_transcript._get_reading_frame(1799)
    assert rf == 2

    rf = test_mane_transcript._get_reading_frame(1800)
    assert rf == 3

    rf = test_mane_transcript._get_reading_frame(2573)
    assert rf == 2


def test_p_to_c_pos(test_mane_transcript):
    """Test that _p_to_c_pos method works correctly."""
    # https://civicdb.org/events/genes/5/summary/variants/12/summary#variant
    c_pos = test_mane_transcript._p_to_c_pos(600, 600)
    assert c_pos == (1798, 1800)

    # https://civicdb.org/events/genes/19/summary/variants/33/summary#variant
    c_pos = test_mane_transcript._p_to_c_pos(858, 858)
    assert c_pos == (2572, 2574)

    # https://civicdb.org/events/genes/29/summary/variants/72/summary#variant
    c_pos = test_mane_transcript._p_to_c_pos(576, 576)
    assert c_pos == (1726, 1728)

    # https://civicdb.org/events/genes/38/summary/variants/99/summary#variant
    c_pos = test_mane_transcript._p_to_c_pos(842, 842)
    assert c_pos == (2524, 2526)

    # https://civicdb.org/events/genes/30/summary/variants/78/summary#variant
    c_pos = test_mane_transcript._p_to_c_pos(12, 12)
    assert c_pos == (34, 36)

    c_pos = test_mane_transcript._p_to_c_pos(747, 751)
    assert c_pos == (2239, 2253)


def test_p_to_c(test_mane_transcript):
    """Test that _p_to_c method works correctly."""
    # Amino Acid Substitution
    expected_pos = 1798, 1800
    ac, pos = test_mane_transcript._p_to_c('NP_004324.2', 600, 600)
    assert ac == 'NM_004333.6'
    assert pos == expected_pos

    ac, pos = test_mane_transcript._p_to_c('ENSP00000288602.7', 600, 600)
    assert ac == 'ENST00000288602.11'
    assert pos == expected_pos

    expected_pos = 2572, 2574
    ac, pos = test_mane_transcript._p_to_c('NP_005219.2', 858, 858)
    assert ac == 'NM_005228.5'
    assert pos == expected_pos

    ac, pos = test_mane_transcript._p_to_c('ENSP00000275493.2', 858, 858)
    assert ac == 'ENST00000275493.7'
    assert pos == expected_pos

    # Polypeptide Truncation
    expected_pos = 553, 555
    ac, pos = test_mane_transcript._p_to_c('NP_000542.1', 185, 185)
    assert ac == 'NM_000551.4'
    assert pos == expected_pos

    ac, pos = test_mane_transcript._p_to_c('ENSP00000256474.3', 185, 185)
    assert ac == 'ENST00000256474.3'
    assert pos == expected_pos

    # Silent Mutation
    expected_pos = 181, 183
    ac, pos = test_mane_transcript._p_to_c('NP_000542.1', 61, 61)
    assert ac == 'NM_000551.4'
    assert pos == expected_pos


def test_c_to_g(test_mane_transcript, nm_004333v6_g):
    """Test that _c_to_g method works correctly."""
    tx_ac = 'NM_004333.6'
    g = test_mane_transcript._c_to_g(tx_ac, (1798, 1800))
    assert g == nm_004333v6_g


def test__g_to_mane_c(test_mane_transcript, braf_mane_data, nm_004333v6_g,
                      braf_v600e_mane_c):
    """Test that _g_to_mane_c method works correctly."""
    mane_c = test_mane_transcript._g_to_mane_c(
        nm_004333v6_g, braf_mane_data
    )
    expected = copy.deepcopy(braf_v600e_mane_c)
    expected['pos'] = (1918, 1920)
    assert mane_c == expected


def test_get_mane_p(test_mane_transcript, braf_mane_data, braf_v600e_mane_p):
    """Test that _get_mane_p method works correctly."""
    mane_p = test_mane_transcript._get_mane_p(braf_mane_data, (1918, 1920))
    assert mane_p == braf_v600e_mane_p


def test_p_to_mane_p(test_mane_transcript, braf_v600e_mane_p,
                     egfr_l858r_mane_p):
    """Test that p_to_mane_p method works correctly."""
    # BRAF V600E RefSeq Accessions
    mane_p = test_mane_transcript.get_mane_transcript('NP_004324.2', 600,
                                                      None, 'p')
    assert mane_p == braf_v600e_mane_p

    mane_p = test_mane_transcript.get_mane_transcript('NP_004324.2', 600,
                                                      600, 'p')
    assert mane_p == braf_v600e_mane_p

    # BRAF V600E Ensembl Accessions
    mane_p = test_mane_transcript.get_mane_transcript('ENSP00000288602.7',
                                                      600, None, 'p')
    assert mane_p == braf_v600e_mane_p

    mane_p = test_mane_transcript.get_mane_transcript('ENSP00000288602.7',
                                                      600, 600, 'p')
    assert mane_p == braf_v600e_mane_p

    # EGFR L858R RefSeq Accessions
    mane_p = test_mane_transcript.get_mane_transcript('NP_005219.2', 858,
                                                      None, 'p')
    assert mane_p == egfr_l858r_mane_p

    mane_p = test_mane_transcript.get_mane_transcript('NP_005219.2', 858,
                                                      858, 'p')
    assert mane_p == egfr_l858r_mane_p

    # EGFR L858R Ensembl Accessions
    mane_p = test_mane_transcript.get_mane_transcript('ENSP00000275493.2',
                                                      858, None, 'p')
    assert mane_p == egfr_l858r_mane_p

    mane_p = test_mane_transcript.get_mane_transcript('ENSP00000275493.2',
                                                      858, 858, 'p')
    assert mane_p == egfr_l858r_mane_p

    assert test_mane_transcript.get_mane_transcript('NP_004439.2', 755,
                                                    759, 'p')


def test_c_to_mane_c(test_mane_transcript, braf_v600e_mane_c,
                     egfr_l858r_mane_c):
    """Test that c_to_mane_p method works correctly."""
    # BRAF V600E RefSeq Accessions
    mane_c = test_mane_transcript.get_mane_transcript('NM_004333.4', 1799,
                                                      None, 'c')
    assert mane_c == braf_v600e_mane_c

    mane_c = test_mane_transcript.get_mane_transcript('NM_004333.5', 1799,
                                                      None, 'C')
    assert mane_c == braf_v600e_mane_c

    mane_c = test_mane_transcript.get_mane_transcript('NM_004333.6', 1799,
                                                      None, 'c')
    assert mane_c == braf_v600e_mane_c

    # BRAF V600E Ensembl Accessions
    mane_c = test_mane_transcript.get_mane_transcript('ENST00000288602.10',
                                                      1799, None, 'c')
    assert mane_c == braf_v600e_mane_c

    mane_c = test_mane_transcript.get_mane_transcript('ENST00000288602.11',
                                                      1799, None, 'c')
    assert mane_c == braf_v600e_mane_c

    # EGFR L858R RefSeq Accessions
    mane_c = test_mane_transcript.get_mane_transcript('NM_005228.3',
                                                      2573, None, 'c')
    assert mane_c == egfr_l858r_mane_c

    mane_c = test_mane_transcript.get_mane_transcript('NM_005228.4',
                                                      2573, None, 'c')
    assert mane_c == egfr_l858r_mane_c

    mane_c = test_mane_transcript.get_mane_transcript('NM_005228.5', 2573,
                                                      2573, 'c')
    assert mane_c == egfr_l858r_mane_c

    # EGFR L858R Ensembl Accessions
    mane_c = test_mane_transcript.get_mane_transcript('ENST00000275493.7',
                                                      2573, None, 'c')
    assert mane_c == egfr_l858r_mane_c

    mane_c = test_mane_transcript.get_mane_transcript('ENST00000275493.6',
                                                      2573, None, 'c')
    assert mane_c == egfr_l858r_mane_c


def test_g_to_mane_c(test_mane_transcript, egfr_l858r_mane_c,
                     braf_v600e_mane_c):
    """Test that g_to_mane_c method works correctly."""
    mane_c = test_mane_transcript.g_to_mane_c('NC_000007.13', 55259515, None)
    assert mane_c == egfr_l858r_mane_c

    mane_c = test_mane_transcript.g_to_mane_c('NC_000007.13', 55259515,
                                              55259515)
    assert mane_c == egfr_l858r_mane_c

    mane_c = test_mane_transcript.g_to_mane_c('NC_000007.13', 140453136, None)
    assert mane_c == braf_v600e_mane_c


def test_no_matches(test_mane_transcript):
    """Test that invalid queries return None."""
    # Invalid ENST version
    mane_c = test_mane_transcript.get_mane_transcript(
        'ENST00000275493.15645', 2573, None, 'c'
    )
    assert mane_c is None
