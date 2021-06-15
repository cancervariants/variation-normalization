"""Module for retrieving MANE transcript."""
from typing import Optional, Tuple, Dict
import hgvs.parser
import logging
import math
from pydantic.types import StrictBool

logger = logging.getLogger('variation')
logger.setLevel(logging.DEBUG)


# TODO:
#  ENST queries
#  Validation:
#     Check correct reading frame
#     Checks references match (should this be in validators?)
#     Exon Structure
#  g -> MANE c


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

    def _p_to_c(self, ac, start_pos, end_pos)\
            -> Optional[Tuple[str, Tuple[int, int]]]:
        """Convert protein (p.) annotation to cDNA (c.) annotation.

        :param str ac: Transcript accession
        :param int start_pos: Protein start position
        :param int end_pos: Protein end position
        :return: [cDNA transcript accession, [cDNA pos start, cDNA pos end]]
        """
        # TODO: Check version mappings 1 to 1 relationship
        temp_ac = self.uta.p_to_c_ac(ac)
        if temp_ac:
            ac = temp_ac[-1][1]
        else:
            try:
                if ac.startswith('NP_'):
                    ac = self.transcript_mappings.np_to_nm[ac]
                elif ac.startswith('ENSP'):
                    ac = \
                        self.transcript_mappings.ensp_to_enst[ac]
                else:
                    logger.warning(f"Unable to find accession: {ac}")
                    return None
            except KeyError:
                logger.warning(f"{ac} not found in transcript_mappings")
                return None

        pos = self._p_to_c_pos(start_pos)
        if end_pos is not None:
            end_pos = self._p_to_c_pos(end_pos)
            pos = pos[0], end_pos[1]
        return ac, pos

    def _c_to_g(self, ac, pos) -> Optional[Dict]:
        """Get g. annotation from c. annotation.

        :param str ac: cDNA accession
        :param tuple pos: [cDNA pos start, cDNA pos end]
        :return: Gene, Transcript accession and position change,
            Altered transcript accession and position change, Strand
        """
        # UTA does not store ENST versions
        if ac.startswith('ENST'):
            ac = ac.split('.')[0]

        coding_start_site = self.uta.get_coding_start_site(ac)
        pos = pos[0] + coding_start_site, pos[1] + coding_start_site

        alt_tx_data = self.uta.get_alt_tx_data(ac, pos)
        if not alt_tx_data:
            return None

        self.uta.liftover_to_38(alt_tx_data)
        return alt_tx_data

    def _g_to_mane_c(self, g, mane_data) -> Optional[Dict]:
        """Get MANE Transcript c. annotation from g. annotation.

        :param dict g: Genomic data
        :param dict mane_data: MANE Transcript data (Transcript accessions,
            gene, and location information)
        :return: MANE Transcripts accessions for RefSeq and Ensembl c.
            coordinates, and position where change occurred on these accessions
        """
        mane_c_ac = mane_data['RefSeq_nuc']
        result = self.uta.get_mane_tx_c_data(mane_c_ac, g['alt_ac'],
                                             g['alt_pos_range'])
        if not result:
            logger.warning(f"Unable to find MANE Transcript {mane_c_ac} "
                           f"position change.")
            return None
        else:
            result = result[-1]

        coding_start_site = \
            self.uta.get_coding_start_site(mane_data['RefSeq_nuc'])

        mane_tx_pos_range = result[6], result[7]
        mane_c_pos_change = (
            mane_tx_pos_range[0] + g['pos_change'][0] - coding_start_site,
            mane_tx_pos_range[1] - g['pos_change'][1] - coding_start_site
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
            pos=(math.ceil(mane_c_pos_range[0] / 3),
                 math.floor(mane_c_pos_range[1] / 3))  # TODO: Check
        )

    def _validate_reading_frames(self, ac, start_pos, end_pos,
                                 mane_transcript) -> StrictBool:
        """Return whether reading frames are the same after translation.

        :param str ac: Query accession
        :param int start_pos: Original start position change
        :param int end_pos: Original end position change
        :param dict mane_transcript: Ensembl and RefSeq transcripts with
            corresponding position change
        """
        start_og_rf = self._get_reading_frame(start_pos)
        start_mane_rf = self._get_reading_frame(mane_transcript['pos'][0])
        if start_og_rf != start_mane_rf:
            logger.warning(f"{ac} original start reading frame "
                           f"({start_og_rf}) does not match"
                           f" {mane_transcript['refseq']}, "
                           f"{mane_transcript['ensembl']} MANE start "
                           f"reading frame ({start_mane_rf})")
            return False

        if end_pos is not None:
            end_og_rf = self._get_reading_frame(end_pos)
            end_mane_rf = self._get_reading_frame(mane_transcript['pos'][1])
            if end_og_rf != end_mane_rf:
                logger.warning(f"{ac} original end reading frame "
                               f"({start_og_rf}) does not match"
                               f" {mane_transcript['refseq']}, "
                               f"{mane_transcript['ensembl']} MANE end "
                               f"reading frame ({start_mane_rf})")
                return False
        return True

    def p_to_mane_p(self, ac, start_pos, end_pos) -> Optional[Dict]:
        """Return MANE Transcript on the p. coordinate.
        p->c->g->GRCh38->MANE c.->MANE p.

        :param str ac: Protein accession
        :param int start_pos: Protein start position
        :param int end_pos: Protein end position
        :return: MANE transcripts with position change on p. coordinate
        """
        c = self._p_to_c(ac, start_pos, end_pos)
        if not c:
            return None
        c_ac, pos = c
        g = self._c_to_g(c_ac, pos)
        mane_data = self.mane_transcript_mappings.get_gene_mane_data(g['gene'])
        mane_c = self._g_to_mane_c(g, mane_data)
        mane_p = self._mane_c_to_mane_p(mane_data, mane_c['pos'])

        # Check reading frames
        valid_reading_frames = self._validate_reading_frames(ac, start_pos,
                                                             end_pos, mane_p)
        if not valid_reading_frames:
            return None
        return mane_p

    def c_to_mane_c(self, ac, pos) -> Optional[Dict]:
        """Return MANE Transcript on the c. coordinate.
        c->g->GRCh38->MANE c.

        :param str ac: Transcript accession on c. coordinate
        :param int pos: cDNA change position
        :return: MANE Transcripts with cDNA change on c. coordinate
        """
        pos = pos - 1, pos + 1
        g = self._c_to_g(ac, pos)
        mane_data = self.mane_transcript_mappings.get_gene_mane_data(g['gene'])
        mane_c = self._g_to_mane_c(g, mane_data)

        # Check reading frames
        valid_reading_frames = self._validate_reading_frames(ac, pos[0],
                                                             pos[1], mane_c)
        if not valid_reading_frames:
            return None

        return mane_c
