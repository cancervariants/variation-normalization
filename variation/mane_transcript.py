"""Module for retrieving MANE transcript."""
from typing import Optional, Tuple, Dict
import hgvs.parser
import logging


logger = logging.getLogger('variation')
logger.setLevel(logging.DEBUG)


class MANETranscript:
    """Class for retrieving MANE transcripts."""

    def __init__(self, transcript_mappings,
                 mane_transcript_mappings, uta) -> None:
        """Initialize the MANETranscript class.

        :param TranscriptMappings transcript_mappings: Access to transcript
            accession mappings and conversions
        :param MANETranscriptMappings mane_transcript_mappings: Access to
            MANE Transcript accession mapping data
        :param UTA uta: UTA instance to give access to query methods for
            transcript alignments
        """
        self.hgvs_parser = hgvs.parser.Parser()
        self.transcript_mappings = transcript_mappings
        self.mane_transcript_mappings = mane_transcript_mappings
        self.uta = uta

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
        :return: cDNA position start, cDNA position end
        """
        pos_mod_3 = p_pos % 3
        pos = p_pos * 3
        if pos_mod_3 == 0:
            pos -= 1
        return pos - 1, pos + 1

    def _p_to_c(self, ac, pos) -> Optional[Tuple[str, Tuple[int, int]]]:
        """Convert protein (p.) annotation to cDNA (c.) annotation.

        :param str ac: Transcript accession
        :param int pos: Protein position where change occurred
        :return: [cDNA transcript accession, [cDNA pos start, cDNA pos end]]
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

    def _c_to_g(self, ac, pos) -> Optional[Dict]:
        """Get g. annotation from c. annotation.

        :param str ac: cDNA accession
        :param tuple pos: [cDNA pos start, cDNA pos end]
        :return: Gene, Transcript accession and position change,
            Altered transcript accession and position change, Strand
        """
        alt_tx_data = self.uta.get_alt_tx_data(ac, pos)
        if not alt_tx_data:
            return None

        self.uta.liftover_to_38(alt_tx_data)
        return alt_tx_data

    def _g_to_mane_c(self, nc_ac, mane_c_ac, genomic_change_range,
                     alt_tx_data, mane_data) -> Optional[Dict]:
        """Get MANE Transcript c. annotation from g. annotation.

        :param str nc_ac: NC accession
        :param str mane_c_ac: MANE Transcript c. accession
        :param tuple[int, int] genomic_change_range: [genomic change start
            position, genomic change end position]
        :param dict alt_tx_data: Transcript data
        :param dict mane_data: MANE Transcript data (Transcript accessions,
            gene, and location information)
        :return: MANE Transcripts accessions for RefSeq and Ensembl c.
            coordinates, and position where change occurred on these accessions
        """
        result = self.uta.get_mane_tx_c_data(mane_c_ac, nc_ac,
                                             genomic_change_range)
        if not result:
            logger.warning(f"Unable to find MANE Transcript {mane_c_ac} "
                           f"position change.")
            return None
        else:
            result = result[-1]

        mane_tx_pos_range = result[6], result[7]
        mane_c_pos_change = (
            mane_tx_pos_range[0] + alt_tx_data['pos_change'][0],
            mane_tx_pos_range[1] - alt_tx_data['pos_change'][1]
        )

        return dict(
            refseq=mane_data['RefSeq_nuc'],
            ensembl=mane_data['Ensembl_nuc'],
            pos=mane_c_pos_change
        )

    def _mane_c_to_mane_p(self, mane_data, mane_c_pos_range) -> Dict:
        """Translate MANE Transcript c. annotation to p. annotation

        :param dict mane_data: MANE Transcript data
        :param tuple[int, int] mane_c_pos_range: Position change range
            on MANE Transcript c. coordinate
        :return: MANE transcripts accessions and position change on
            p. coordinate
        """
        return dict(
            refseq=mane_data['RefSeq_prot'],
            ensembl=mane_data['Ensembl_prot'],
            pos=int(mane_c_pos_range[1] / 3)
        )

    def p_to_mane_p(self, ac, pos) -> Dict:
        """Return MANE Transcript on the p. coordinate.
        p->c->g->GRCh38->MANE c.->MANE p.

        :param str ac: Protein accession
        :param int pos: Protein change position
        :return: MANE transcripts with position change on p. coordinate
        """
        c = self._p_to_c(ac, pos)
        if not c:
            return None
        c_ac, pos = c
        g = self._c_to_g(c_ac, pos)
        mane_data = self.mane_transcript_mappings.get_gene_mane_data(g['gene'])
        mane_c = self._g_to_mane_c(g['alt_ac'], mane_data['RefSeq_nuc'],
                                   g['alt_pos_range'], g, mane_data)
        return self._mane_c_to_mane_p(mane_data, mane_c['pos'])
