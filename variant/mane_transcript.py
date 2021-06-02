"""Module for retrieving MANE transcript."""
from typing import Optional, Tuple
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

    def p_to_c(self, ac, pos) -> Optional[Tuple[str, Tuple[int, int]]]:
        """Convert protein (p.) annotation to cDNA (c.) annotation.

        :param str ac: Transcript accession
        :param int pos: Protein position where change occurred
        :return: [cDNA transcript accession, cDNA pos start, cDNA pos end]
        """
        # TODO: Check version mappings 1 to 1 relationship
        if ac.startswith('NP_'):
            ac = self.transcript_mappings.np_to_nm[ac]
        elif ac.startswith('ENSP'):
            ac = \
                self.transcript_mappings.ensp_to_enst[ac]
        else:
            return None

        pos = self._p_to_c_pos(pos)
        return ac, pos

    def _get_reading_frame(self, pos) -> int:
        """Return reading frame number.

        :param int pos: Position
        :return: Reading frame
        """
        pos_mod_3 = (pos - 1) % 3
        return pos_mod_3

    def _p_to_c_pos(self, p_pos) -> Tuple[int, int]:
        """Return cDNA position given a protein position.

        :param int p_pos: Protein position
        :return: cDNA position
        """
        pos_mod_3 = p_pos % 3
        pos = p_pos * 3
        if pos_mod_3 == 0:
            pos -= 1
        return pos - 1, pos + 1
