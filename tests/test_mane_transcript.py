"""Module for testing MANE Transcript class."""
import pytest
from variant.mane_transcript import MANETranscript
from variant.data_sources import TranscriptMappings
from variant.tokenizers.caches import AminoAcidCache
import hgvs.parser


@pytest.fixture(scope='module')
def test_mane_transcript():
    """Build mane transcript test fixture."""
    class TestMANETranscript:

        def __init__(self):
            self.test_mane_transcript = MANETranscript(TranscriptMappings(),
                                                       AminoAcidCache())
            self.test_hgvs_parser = hgvs.parser.Parser()

        def p_to_c(self, transcript, token):
            return self.test_mane_transcript.p_to_c(transcript, token)
    return TestMANETranscript()


def test_p_to_c(test_mane_transcript):
    """Test that p_to_c method works correctly."""
    # Amino Acid Substitution
    expected_pos = 1798, 1800
    ac, pos = test_mane_transcript.p_to_c('NP_004324.2', 600)
    assert ac == 'NM_004333.4'
    assert pos == expected_pos

    ac, pos = test_mane_transcript.p_to_c('ENSP00000288602.7', 600)
    assert ac == 'ENST00000288602.11'
    assert pos == expected_pos

    expected_pos = 2572, 2574
    ac, pos = test_mane_transcript.p_to_c('NP_005219.2', 858)
    assert ac == 'NM_005228.3'
    assert pos == expected_pos

    ac, pos = test_mane_transcript.p_to_c('ENSP00000275493.2', 858)
    assert ac == 'ENST00000275493.7'
    assert pos == expected_pos

    # Polypeptide Truncation
    expected_pos = 554, 556
    ac, pos = test_mane_transcript.p_to_c('NP_000542.1', 185)
    assert ac == 'NM_000551.3'
    assert pos == expected_pos

    ac, pos = test_mane_transcript.p_to_c('ENSP00000256474.3', 185)
    assert ac == 'ENST00000256474.3'
    assert pos == expected_pos

    # Silent Mutation
    expected_pos = 182, 184
    ac, pos = test_mane_transcript.p_to_c('NP_000542.1', 61)
    assert ac == 'NM_000551.3'
    assert pos == expected_pos
