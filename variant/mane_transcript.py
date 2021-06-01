"""Module for retrieving MANE transcript."""
from typing import Optional, Tuple
from variant.schemas.token_response_schema import ReferenceSequence
from variant.data_sources import CodonTable
import hgvs.parser
import logging


logger = logging.getLogger('variant')
logger.setLevel(logging.DEBUG)


class MANETranscript:
    """Class for retrieving MANE transcripts."""

    def __init__(self, transcript_mappings, amino_acid_cache) -> None:
        """Initialize the MANETranscript class.

        :param TranscriptMappings transcript_mappings: Access to transcript
            accession mappings
        :param AminoAcidCache amino_acid_cache: Access to amino acid codes
            and conversions
        """
        self.hgvs_parser = hgvs.parser.Parser()
        self.transcript_mappings = transcript_mappings
        self.codon_table = CodonTable(amino_acid_cache)

    def p_to_c(self, ac, token) -> Optional[Tuple[str, int, int]]:
        """Convert protein (p.) annotation to cDNA (c.) annotation.

        :param str ac: Transcript accession
        :param Token token: Classification token
        :return: [cDNA transcript accession, cDNA position, reading frame]
        """
        if token.reference_sequence != ReferenceSequence.PROTEIN:
            logger.debug(f"{token} does not have a "
                         f"protein reference sequence.")

        # TODO: Check version mappings 1 to 1 relationship
        if ac.startswith('NP_'):
            ac = self.transcript_mappings.np_to_nm[ac.split(':')[0]]
        elif ac.startswith('ENSP'):
            ac = \
                self.transcript_mappings.ensp_to_enst[ac.split(':')[0]]
        else:
            return None

        pos = self._p_to_c_pos(token.position)
        rf = self._get_reading_frame(token.position)
        return ac, pos, rf

    def _get_reading_frame(self, pos) -> int:
        """Return reading frame number.

        :param int pos: Position
        :return: Reading frame
        """
        pos_mod_3 = (pos - 1) % 3
        return pos_mod_3 + 1

    def _p_to_c_pos(self, p_pos) -> int:
        """Return cDNA position given a protein position.

        :param int p_pos: Protein position
        :return: cDNA position
        """
        pos_mod_3 = p_pos % 3
        pos = p_pos * 3
        if pos_mod_3 == 0:
            pos -= 1
        return pos
