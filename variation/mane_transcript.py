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
#     Exon Structure
#  g -> MANE c


class MANETranscript:
    """Class for retrieving MANE transcripts."""

    def __init__(self, seqrepo_access, transcript_mappings,
                 mane_transcript_mappings, uta) -> None:
        """Initialize the MANETranscript class.

        :param TranscriptMappings transcript_mappings: Access to transcript
            accession mappings and conversions
        :param MANETranscriptMappings mane_transcript_mappings: Access to
            MANE Transcript accession mapping data
        :param UTA uta: UTA instance to give access to query methods for
            transcript alignments
        """
        self.seqrepo_access = seqrepo_access
        self.hgvs_parser = hgvs.parser.Parser()
        self.transcript_mappings = transcript_mappings
        self.mane_transcript_mappings = mane_transcript_mappings
        self.uta = uta

    def _get_reading_frame(self, pos) -> int:
        """Return reading frame number.

        :param int pos: Position
        :return: Reading frame
        """
        pos_mod_3 = pos % 3
        if pos_mod_3 == 0:
            pos_mod_3 = 3
        return pos_mod_3

    def _p_to_c_pos(self, p_pos) -> Tuple[int, int]:
        """Return cDNA position given a protein position.

        :param int p_pos: Protein position
        :return: cDNA position start, cDNA position end
        """
        pos = p_pos * 3 - 1
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
            if not self.transcript_mappings.ensembl_transcript_version_to_gene_symbol.get(ac):  # noqa: E501
                try:
                    self.seqrepo_access.seq_repo_client.fetch(ac)
                except KeyError:
                    logger.warning(f"Ensembl transcript {ac} not found")
                    return None

            temp_ac = ac.split('.')[0]
        else:
            temp_ac = ac

        coding_start_site = self.uta.get_coding_start_site(temp_ac)
        if coding_start_site is None:
            logger.warning(f"Accession {temp_ac} not found in UTA")
            return None

        pos = pos[0] + coding_start_site, pos[1] + coding_start_site

        genomic_tx_data = self.uta.get_genomic_tx_data(ac, pos)
        if not genomic_tx_data:
            return None

        og_alt_exon_id = genomic_tx_data['alt_exon_id']
        self.uta.liftover_to_38(genomic_tx_data)
        liftover_alt_exon_id = genomic_tx_data['alt_exon_id']

        # Validation check: Exon structure
        if og_alt_exon_id != liftover_alt_exon_id:
            logger.warning(f"Original alt_exon_id {og_alt_exon_id} "
                           f"does not match liftover alt_exon_id "
                           f"{liftover_alt_exon_id}")
            return None

        return genomic_tx_data

    def _get_mane_c(self, mane_data, mane_c_pos_change):
        """Return MANE Transcript data on c. coordinate.

        :param dict mane_data: MANE Transcript data (transcript accessions,
            gene, and location information)
        :param tuple[int, int] mane_c_pos_change: Start and end positions
            for change on c. coordinate
        """
        return dict(
            refseq=mane_data['RefSeq_nuc'],
            ensembl=mane_data['Ensembl_nuc'],
            pos=mane_c_pos_change,
            strand=mane_data['chr_strand'],
            mane_status=mane_data['MANE_status']
        )

    def _get_mane_p(self, mane_data, mane_c_pos_range) -> Dict:
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
                 math.floor(mane_c_pos_range[1] / 3)),  # TODO: Check
            strand=mane_data['chr_strand'],
            mane_status=mane_data['MANE_status']
        )

    def _g_to_mane_c(self, g, mane_data) -> Optional[Dict]:
        """Get MANE Transcript c. annotation from g. annotation.

        :param dict g: Genomic data
        :param dict mane_data: MANE Transcript data (Transcript accessions,
            gene, and location information)
        :return: MANE Transcripts accessions for RefSeq and Ensembl c.
            coordinates, and position where change occurred on these accessions
        """
        mane_c_ac = mane_data['RefSeq_nuc']
        result = self.uta.get_tx_exon_aln_v_data(
            mane_c_ac, g['alt_pos_change_range'][0],
            g['alt_pos_change_range'][1], alt_ac=g['alt_ac'], use_tx_pos=False
        )

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

        return self._get_mane_c(mane_data, mane_c_pos_change)

    def _validate_reading_frames(self, ac, start_pos, end_pos,
                                 mane_transcript) -> StrictBool:
        """Return whether reading frames are the same after translation.

        :param str ac: Query accession
        :param int start_pos: Original start position change
        :param int end_pos: Original end position change
        :param dict mane_transcript: Ensembl and RefSeq transcripts with
            corresponding position change
        """
        for pos, mane_pos_index in [(start_pos, 0), (end_pos, 1)]:
            if pos is not None:
                og_rf = self._get_reading_frame(pos)
                mane_rf = self._get_reading_frame(
                    mane_transcript['pos'][mane_pos_index]
                )

                if og_rf != mane_rf:
                    logger.warning(f"{ac} original reading frame ({og_rf}) "
                                   f"does not match MANE "
                                   f"{mane_transcript['ensembl']}, "
                                   f"{mane_transcript['refseq']} reading "
                                   f"frame ({mane_rf})")
                    return False
            else:
                if mane_pos_index == 0:
                    logger.warning(f"{ac} must having start position")
                    return False
        return True

    def _validate_references(self, ac, start_pos, end_pos,
                             mane_transcript) -> StrictBool:
        """Return whether or not reference changes are the same.

        :param str ac: Query accession
        :param int start_pos: Original start position change
        :param int end_pos: Origin end position change
        :param dict mane_transcript: Ensembl and RefSeq transcripts with
            corresponding position change
        """
        if start_pos == end_pos:
            ref = self.seqrepo_access.sequence_at_position(ac, start_pos)
        else:
            ref = self.seqrepo_access.get_sequence(ac, start_pos, end_pos)
        if not ref:
            return False

        mane_start_pos = mane_transcript['pos'][0]
        mane_end_pos = mane_transcript['pos'][1]
        if mane_start_pos == mane_end_pos:
            mane_ref = self.seqrepo_access.sequence_at_position(
                mane_transcript['refseq'],
                mane_start_pos
            )
        else:
            mane_ref = self.seqrepo_access.get_sequence(
                mane_transcript['refseq'],
                mane_transcript['pos'][0],
                mane_transcript['pos'][1]
            )
        if not mane_ref:
            return False

        if ref != mane_ref:
            logger.warning(f"Original accession, {ac}, ref {ref} does not "
                           f"match MANE accession, {mane_transcript['refseq']}"
                           f", ref {mane_ref}")
            return False

        return True

    def get_mane_transcript(self, ac, start_pos, end_pos,
                            start_annotation_layer) -> Optional[Dict]:
        """Return mane transcript.

        :param str ac: Accession
        :param int start_pos: Start position change
        :param int end_pos: End position change
        :param start_annotation_layer: Annotation layer we are starting from.
            Must be either `p`, `c`, or `g`.
        :return: MANE transcript
        """
        anno = start_annotation_layer.lower()
        if end_pos is None:
            end_pos = start_pos
        if anno in ['p', 'c']:
            # Get accession and position on c. coordinate
            if anno == 'p':
                c = self._p_to_c(ac, start_pos, end_pos)
                if not c:
                    return None
                c_ac, c_pos = c
            else:
                c_ac = ac
                c_pos = start_pos, end_pos

            # Go from c -> g annotation (liftover as well)
            g = self._c_to_g(c_ac, c_pos)
            if g is None:
                return None

            # Go from g -> mane transcript
            mane_data = \
                self.mane_transcript_mappings.get_gene_mane_data(g['gene'])
            if not mane_data:
                return None
            mane_data_len = len(mane_data)

            # Transcript Priority (Must pass validation checks):
            #  1. MANE Select
            #  2. MANE Plus Clinical
            #  3. Longest Compatible Remaining
            for i in range(mane_data_len):
                index = mane_data_len - i - 1
                current_mane_data = mane_data[index]

                mane = self._g_to_mane_c(g, current_mane_data)

                valid_reading_frame = self._validate_reading_frames(
                    c_ac, c_pos[0], c_pos[1], mane
                )
                if not valid_reading_frame:
                    continue

                if anno == 'p':
                    mane = self._get_mane_p(current_mane_data, mane['pos'])

                # TODO: Fix
                valid_references = self._validate_references(
                    ac, start_pos, end_pos, mane
                )
                if not valid_references:
                    continue

                return mane
            return None
        elif anno == 'g':
            # TODO: Uncomment below once working
            # self.g_to_mane_c(ac, start_pos, end_pos)
            pass
        else:
            logger.warning(f"Annotation layer not supported: {anno}")

    def g_to_mane_c(self, ac, start_pos, end_pos):
        """Return MANE Transcript on the c. coordinate.
        g->GRCh38->MANE c.

        :param str ac: Transcript accession on g. coordinate
        :param int start_pos: genomic change start position
        :param int end_pos: genomic change end position
        :return: MANE Transcripts with cDNA change on c. coordinate
        """
        if end_pos is None:
            end_pos = start_pos

        gene_symbol = self.uta.get_gene_from_ac(ac, start_pos, end_pos)
        if not gene_symbol:
            return None

        mane_data =\
            self.mane_transcript_mappings.get_gene_mane_data(gene_symbol)
        if not mane_data:
            return None
        mane_data_len = len(mane_data)
        for i in range(mane_data_len):
            index = mane_data_len - i - 1
            current_mane_data = mane_data[index]

            mane_c_ac = current_mane_data['RefSeq_nuc']
            mane_tx_genomic_data = self.uta.get_mane_c_genomic_data(
                mane_c_ac, ac, start_pos, end_pos
            )

            self.uta.liftover_to_38(mane_tx_genomic_data)

            mane_tx_genomic_data = self.uta.get_mane_c_genomic_data(
                mane_c_ac, mane_tx_genomic_data['alt_ac'],
                mane_tx_genomic_data['alt_pos_change_range'][0],
                mane_tx_genomic_data['alt_pos_change_range'][1]
            )

            coding_start_site = mane_tx_genomic_data['coding_start_site']
            alt_pos_change = (
                mane_tx_genomic_data['alt_pos_change_range'][0] - mane_tx_genomic_data['alt_pos_range'][0],  # noqa: E501
                mane_tx_genomic_data['alt_pos_range'][1] - mane_tx_genomic_data['alt_pos_change_range'][1]  # noqa: E501
            )

            mane_c_pos_change = (
                mane_tx_genomic_data['tx_pos_range'][0] + alt_pos_change[0] - coding_start_site,  # noqa: E501
                mane_tx_genomic_data['tx_pos_range'][1] - alt_pos_change[1] - coding_start_site  # noqa: E501
            )

            return self._get_mane_c(current_mane_data, mane_c_pos_change)
