"""Module for retrieving MANE transcript."""
from typing import Optional, Tuple, Dict
import hgvs.parser
import logging


logger = logging.getLogger('variation')
logger.setLevel(logging.DEBUG)


class MANETranscript:
    """Class for retrieving MANE transcripts."""

    def __init__(self, transcript_mappings, amino_acid_cache,
                 mane_transcript_mappings, uta) -> None:
        """Initialize the MANETranscript class.

        :param TranscriptMappings transcript_mappings: Access to transcript
            accession mappings
        :param AminoAcidCache amino_acid_cache: Access to amino acid codes
            and conversions
        """
        self.hgvs_parser = hgvs.parser.Parser()
        self.transcript_mappings = transcript_mappings
        self.mane_transcript_mappings = mane_transcript_mappings
        self.uta = uta

    def _get_preferred_annotation(self, gene_symbol) -> Dict:
        """Get preferred annotation for a gene symbol.

        :param str gene_symbol: Gene symbol
        :return: MANE transcript data
        """
        # TODO: If MANE Transcript not found, select longest transcript
        return self.mane_transcript_mappings.get_gene_mane_data(gene_symbol)

    def _p_to_g(self, ac, pos) -> Optional[Tuple[str, str, Tuple[int, int]]]:
        """Convert protein annotation to genomic annotation.

        :param str ac: Protein accession
        :param int pos: Protein change position
        :return:
            [Gene Symbol, NC accession, [Genomic start pos, Genomic end pos]]
        """
        c = self._p_to_c(ac, pos)
        if not c:
            return None
        c_ac, pos = c

        return self._c_to_g(c_ac, pos)

    def _p_to_c(self, ac, pos) -> Optional[Tuple[str, Tuple[int, int]]]:
        """Convert protein (p.) annotation to cDNA (c.) annotation.

        :param str ac: Transcript accession
        :param int pos: Protein position where change occurred
        :return: [cDNA transcript accession, cDNA pos start, cDNA pos end]
        """
        # TODO: Check version mappings 1 to 1 relationship
        ac = self.uta.p_to_c_ac(ac)
        if ac:
            ac = ac[-1][1]
        else:
            if ac.startswith('NP_'):
                ac = self.transcript_mappings.np_to_nm[ac]
            elif ac.startswith('ENSP'):
                ac = \
                    self.transcript_mappings.ensp_to_enst[ac]
            else:
                logger.warning(f"Unable to find accession: {ac}")
                return None
        pos = self._p_to_c_pos(pos)
        return ac, pos

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

    def _c_to_g(self, ac, pos):
        """Get g. annotation from c. annotation.

        :param str ac: cDNA accession
        :param tuple pos: [cDNA pos start, cDNA pos end]
        """
        alt_tx_data = self.uta.get_alt_tx_data(ac, pos)
        if not alt_tx_data:
            return None

        gene_symbol, nc_ac, alt_pos, strand = alt_tx_data
        alt_pos, strand = self.uta.liftover_to_38(nc_ac, alt_pos,
                                                  strand=strand)
        return gene_symbol, nc_ac, alt_pos, strand

    def _get_reading_frame(self, pos) -> int:
        """Return reading frame number.

        :param int pos: Position
        :return: Reading frame
        """
        pos_mod_3 = (pos - 1) % 3
        return pos_mod_3
