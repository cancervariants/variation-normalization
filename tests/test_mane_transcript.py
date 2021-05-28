"""Module for testing MANE Transcript class."""
import pytest
from variant.mane_transcript import MANETranscript
from variant.data_sources import TranscriptMappings
from variant.tokenizers.caches import AminoAcidCache
from variant.schemas.token_response_schema import AminoAcidSubstitutionToken,\
    PolypeptideTruncationToken, SilentMutationToken
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


@pytest.fixture(scope='module')
def p_braf_v600e():
    """Create test fixture for BRAF V600E Amino Acid Substitution Token."""
    return AminoAcidSubstitutionToken(
        token='Val600E',
        match_type=5,
        input_string='BRAF V600E',
        ref_protein='Val',
        alt_protein='Glu',
        position=600
    )


@pytest.fixture(scope='module')
def p_egfr_l858r():
    """Create test fixture for EGFR L858R Amino Acid Substitution Token."""
    return AminoAcidSubstitutionToken(
        token='Leu858Arg',
        match_type=5,
        input_string='EGFR L858R',
        ref_protein='Leu',
        alt_protein='Arg',
        position=858
    )


@pytest.fixture(scope='module')
def p_vhl_tyr185ter():
    """Create test fixture for VHL Y185* Polypeptide Truncation Token."""
    return PolypeptideTruncationToken(
        token='Tyr185Ter',
        match_type=5,
        input_string='VHL Y185*',
        ref_protein='Tyr',
        alt_protein='Ter',
        position=185
    )


@pytest.fixture(scope='module')
def p_vhl_p61p():
    """Create test fixture for VHL Pro61= Silent Mutation Token."""
    return SilentMutationToken(
        token='Pro61=',
        match_type=5,
        input_string='VHL Pro61=',
        ref_protein='Pro',
        position=61
    )


def test_p_to_c(test_mane_transcript, p_braf_v600e, p_egfr_l858r,
                p_vhl_tyr185ter, p_vhl_p61p):
    """Test that p_to_c method works correctly."""
    # Amino Acid Substitution
    ac, pos = test_mane_transcript.p_to_c('NP_004324.2:p.Val600Glu',
                                          p_braf_v600e)
    assert ac == 'NM_004333.4'
    assert pos == 1799

    ac, pos = test_mane_transcript.p_to_c('ENSP00000288602.7:p.Val600Glu',
                                          p_braf_v600e)
    assert ac == 'ENST00000288602.11'
    assert pos == 1799

    ac, pos = test_mane_transcript.p_to_c('NP_005219.2:p.Leu858Arg',
                                          p_egfr_l858r)
    assert ac == 'NM_005228.3'
    assert pos == 2573

    ac, pos = test_mane_transcript.p_to_c('ENSP00000275493.2:p.Leu858Arg',
                                          p_egfr_l858r)
    assert ac == 'ENST00000275493.7'
    assert pos == 2573

    # Polypeptide Truncation
    ac, pos = test_mane_transcript.p_to_c('NP_000542.1:p.Tyr185Ter',
                                          p_vhl_tyr185ter)
    assert ac == 'NM_000551.3'
    assert pos == 555

    ac, pos = test_mane_transcript.p_to_c('ENSP00000256474.3:p.Tyr185Ter',
                                          p_vhl_tyr185ter)
    assert ac == 'ENST00000256474.3'
    assert pos == 555

    # Silent Mutation
    ac, pos = test_mane_transcript.p_to_c('NP_000542.1:p.Pro61=', p_vhl_p61p)
    assert ac == 'NM_000551.3'
    assert pos == 183
