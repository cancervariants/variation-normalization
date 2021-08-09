"""Module for retrieving MANE Transcript from variation on p/c/g coordinate.

Steps:
1. Map annotation layer to genome
2. Liftover to preferred genome
    We want to liftover to GRCh38. We do not support getting MANE transcripts
    for GRCh36 and earlier assemblies.
3. Select preferred compatible annotation
4. Map back to correct annotation layer
"""
from typing import Optional, Tuple, Dict
import hgvs.parser
import logging
import math
from pydantic.types import StrictBool
from pyliftover import LiftOver
from variation.data_sources.uta import GRCH_TO_HG

logger = logging.getLogger('variation')
logger.setLevel(logging.DEBUG)


class MANETranscript:
    """Class for retrieving MANE transcripts."""

    def __init__(self, seqrepo_access, transcript_mappings,
                 mane_transcript_mappings, uta) -> None:
        """Initialize the MANETranscript class.

        :param SeqRepoAccess seqrepo_access: Access to seqrepo queries
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
        Only used on c. coordinate

        :param int pos: Position
        :return: Reading frame
        """
        pos_mod_3 = pos % 3
        if pos_mod_3 == 0:
            pos_mod_3 = 3
        return pos_mod_3

    def _p_to_c_pos(self, start, end) -> Tuple[int, int]:
        """Return cDNA position given a protein position.

        :param int start: Start protein position
        :param int end: End protein position
        :return: cDNA position start, cDNA position end
        """
        start_pos = start * 3 - 1
        if end != start:
            end_pos = end * 3 - 1
        else:
            end_pos = start_pos

        return start_pos - 1, end_pos + 1

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
            ac = temp_ac[-1][0]
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

        pos = self._p_to_c_pos(start_pos, end_pos)
        return ac, pos

    def _c_to_g(self, ac, pos) -> Optional[Dict]:
        """Get g. annotation from c. annotation.

        :param str ac: cDNA accession
        :param tuple pos: [cDNA pos start, cDNA pos end]
        :return: Gene, Transcript accession and position change,
            Altered transcript accession and position change, Strand
        """
        # UTA does not store ENST versions
        # So we want to make sure version is valid
        if ac.startswith('ENST'):
            if not self.transcript_mappings.ensembl_transcript_version_to_gene_symbol.get(ac):  # noqa: E501
                try:
                    self.seqrepo_access.seq_repo_client.fetch(ac)
                except KeyError:
                    logger.warning(f"Ensembl transcript not found: {ac}")
                    return None

            temp_ac = ac.split('.')[0]
        else:
            temp_ac = ac

        # c. coordinate does not contain cds start, so we need to add it
        cds_start_end = self.uta.get_cds_start_end(temp_ac)
        if not cds_start_end:
            logger.warning(f"Accession {temp_ac} not found in UTA")
            return None
        coding_start_site = cds_start_end[0]
        pos = pos[0] + coding_start_site, pos[1] + coding_start_site

        genomic_tx_data = self.uta.get_genomic_tx_data(ac, pos)
        if not genomic_tx_data:
            return None
        genomic_tx_data['coding_start_site'] = coding_start_site

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

    def _get_mane_c(self, mane_data, mane_c_pos_change, cds_start_end):
        """Return MANE Transcript data on c. coordinate.

        :param dict mane_data: MANE Transcript data (transcript accessions,
            gene, and location information)
        :param tuple[int, int] mane_c_pos_change: Start and end positions
            for change on c. coordinate
        :param DictRow[int, int]: Coding start and end site for MANE transcript
        """
        cds_start = cds_start_end[0]
        cds_end = cds_start_end[1]
        lt_cds_start = (mane_c_pos_change[0] < cds_start and mane_c_pos_change[1] < cds_start)  # noqa: E501
        gt_cds_end = (mane_c_pos_change[1] > cds_end and mane_c_pos_change[1] > cds_end)  # noqa: E501

        if lt_cds_start or gt_cds_end:
            logger.info(f"{mane_data['RefSeq_nuc']} with position"
                        f" {mane_c_pos_change} is not within CDS start/end")

        return dict(
            gene=mane_data['symbol'],
            refseq=mane_data['RefSeq_nuc'],
            ensembl=mane_data['Ensembl_nuc'],
            coding_start_site=cds_start,
            coding_end_site=cds_end,
            pos=mane_c_pos_change,
            strand=mane_data['chr_strand'],
            status=mane_data['MANE_status']
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
            gene=mane_data['symbol'],
            refseq=mane_data['RefSeq_prot'],
            ensembl=mane_data['Ensembl_prot'],
            pos=(math.ceil(mane_c_pos_range[0] / 3),
                 math.floor(mane_c_pos_range[1] / 3)),
            strand=mane_data['chr_strand'],
            status=mane_data['MANE_status']
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

        cds_start_end = self.uta.get_cds_start_end(mane_data['RefSeq_nuc'])
        if not cds_start_end:
            return None
        coding_start_site = cds_start_end[0]

        g_pos = g['alt_pos_change_range'][0], g['alt_pos_change_range'][1]
        g_pos_change = g_pos[0] - result[8], result[9] - g_pos[1]

        mane_tx_pos_range = result[6], result[7]
        mane_c_pos_change = (
            mane_tx_pos_range[0] + g_pos_change[0] - coding_start_site,
            mane_tx_pos_range[1] - g_pos_change[1] - coding_start_site
        )

        if mane_c_pos_change[0] > mane_c_pos_change[1]:
            mane_c_pos_change = mane_c_pos_change[1], mane_c_pos_change[0]

        return self._get_mane_c(mane_data, mane_c_pos_change,
                                cds_start_end)

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

    def _validate_references(self, ac, coding_start_site, start_pos, end_pos,
                             mane_transcript, expected_ref,
                             anno) -> StrictBool:
        """Return whether or not reference changes are the same.

        :param str ac: Query accession
        :param int coding_start_site: ac's coding start site
        :param int start_pos: Original start position change
        :param int end_pos: Origin end position change
        :param dict mane_transcript: Ensembl and RefSeq transcripts with
            corresponding position change
        :param str expected_ref: Reference at position given during input
        :param str anno: Annotation layer we are starting from.
            Must be either `p`, `c`, or `g`.
        :return: `True` if reference check passes. `False` otherwise.
        """
        if anno == 'c':
            start_pos += coding_start_site
            end_pos += coding_start_site

        if start_pos == end_pos:
            ref = self.seqrepo_access.sequence_at_position(ac, start_pos)
        else:
            ref = self.seqrepo_access.get_sequence(ac, start_pos, end_pos)
        if not ref:
            return False

        if mane_transcript:
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
                logger.info("Unable to validate reference for MANE Transcript")

            if expected_ref != mane_ref:
                logger.info(f"Expected ref, {expected_ref}, but got {mane_ref}"
                            f" on MANE accession, {mane_transcript['refseq']}")

        if expected_ref != ref:
            logger.warning(f"Expected ref, {expected_ref}, but got {ref} "
                           f"on accession, {ac}")
            return False

        return True

    def _validate_index(self, ac, pos, coding_start_site):
        """Validate that positions actually exist on accession

        :param str ac: Accession
        :param tuple[int, int] pos: Start position change, End position change
        :param int coding_start_site: coding start site for accession
        :return: `True` if positions exist on accession. `False` otherwise
        """
        len_of_seq = self.seqrepo_access.len_of_sequence(ac)
        return pos[0] + coding_start_site <= pos[1] + coding_start_site <= len_of_seq  # noqa: E501

    def get_longest_compatible_transcript(self, gene, start_pos, end_pos,
                                          start_annotation_layer,
                                          ref=None):
        """Get longest compatible transcript from a gene.
        Try GRCh38 first, then GRCh37.
        Transcript is compatible if it passes validation checks.

        :param str gene: Gene symbol
        :param int start_pos: Start position change
        :param int end_pos: End position change
        :param  str start_annotation_layer: Starting annotation layer.
            Must be either `p`, or `c`.
        :param str ref: Reference at position given during input
        :return:
        """
        anno = start_annotation_layer.lower()
        if anno not in ['p', 'c']:
            logger.warning(f"Annotation layer not supported: {anno}")
            return None

        if end_pos is None:
            end_pos = start_pos

        if anno == 'p':
            c_start_pos, c_end_pos = self._p_to_c_pos(start_pos, end_pos)
        else:
            c_start_pos, c_end_pos = start_pos, end_pos

        # Data Frame that contains transcripts associated to a gene
        df = self.uta.get_transcripts_from_gene(gene, c_start_pos, c_end_pos)
        nc_acs = list(df['alt_ac'].unique())
        nc_acs.sort(reverse=True)

        for nc_ac in nc_acs:
            # Most recent accession first
            tmp_df = df.loc[df['alt_ac'] == nc_ac]
            for index, row in tmp_df.iterrows():
                g = self._c_to_g(row['tx_ac'], (c_start_pos, c_end_pos))
                if not g:
                    return None

                # Validate references
                if ref:
                    if anno == 'p':
                        valid_references = self._validate_references(
                            row['pro_ac'], row['cds_start_i'], start_pos,
                            end_pos, {}, ref, 'p'
                        )
                    else:
                        valid_references = self._validate_references(
                            row['tx_ac'], row['cds_start_i'], c_start_pos,
                            c_end_pos, {}, ref, 'c'
                        )

                    if not valid_references:
                        continue

                if anno == 'p':
                    pos = start_pos, end_pos
                    ac = row['pro_ac']
                    coding_start_site = 0

                else:
                    pos = c_start_pos, c_end_pos
                    ac = row['tx_ac']
                    coding_start_site = row['cds_start_i']

                if not self._validate_index(
                    ac, pos, coding_start_site
                ):
                    logger.warning(f"{pos} are not valid positions on {ac}"
                                   f"with coding start site "
                                   f"{coding_start_site}")
                    continue

                return dict(
                    refseq=ac if ac.startswith('N') else None,
                    ensembl=ac if ac.startswith('E') else None,
                    pos=pos,
                    strand=g['strand'],
                    status='Longest Compatible Remaining'
                )
        return None

    def get_mane_transcript(self, ac, start_pos, end_pos,
                            start_annotation_layer, gene=None, ref=None,
                            normalize_endpoint=False) -> Optional[Dict]:
        """Return mane transcript.

        :param str ac: Accession
        :param int start_pos: Start position change
        :param int end_pos: End position change
        :param str start_annotation_layer: Starting annotation layer.
            Must be either `p`, `c`, or `g`.
        :param str gene: Gene symbol
        :param str ref: Reference at position given during input
        :param bool normalize_endpoint: `True` if normalize endpoint is being
            used. `False` otherwise. Only for normalize endpoint will we
            look up set of transcripts associated with a gene if MANE
            Transcript not found or validation checks not passed.
        :return: MANE transcript
        """
        anno = start_annotation_layer.lower()
        if end_pos is None:
            end_pos = start_pos

        if isinstance(start_pos, str) or isinstance(end_pos, str):
            try:
                start_pos = int(start_pos)
                end_pos = int(end_pos)
            except ValueError:
                logger.warning(f"{start_pos} and {end_pos} "
                               f"must be valid integers")
                return None

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
                if not mane:
                    continue

                valid_reading_frame = self._validate_reading_frames(
                    c_ac, c_pos[0], c_pos[1], mane
                )
                if not valid_reading_frame:
                    continue

                if anno == 'p':
                    mane = self._get_mane_p(current_mane_data, mane['pos'])

                if ref:
                    valid_references = self._validate_references(
                        ac, g['coding_start_site'], start_pos, end_pos,
                        mane, ref, anno
                    )
                    if not valid_references:
                        continue

                return mane
            if normalize_endpoint:
                if anno == 'p':
                    return self.get_longest_compatible_transcript(
                        g['gene'], start_pos, end_pos, 'p', ref
                    )
                else:
                    return self.get_longest_compatible_transcript(
                        g['gene'], c_pos[0], c_pos[1], 'c', ref
                    )
            else:
                return None
        elif anno == 'g':
            return self.g_to_mane_c(ac, start_pos, end_pos, gene=gene)
        else:
            logger.warning(f"Annotation layer not supported: {anno}")

    def g_to_grch38(self, ac, start_pos, end_pos) -> Optional[Dict]:
        """Return genomic coordinate on GRCh38 when not given gene context.

        :param str ac: Genomic accession
        :param int start_pos: Genomic start position change
        :param int end_pos: Genomic end position change
        :return: NC accession, start and end pos on GRCh38 assembly
        """
        if end_pos is None:
            end_pos = start_pos

        # Checking to see what chromosome and assembly we're on
        descr = self.uta.get_chr_assembly(ac)
        if not descr:
            # Already GRCh38 assembly
            if self._validate_index(ac, (start_pos, end_pos), 0):
                return dict(
                    ac=ac,
                    pos=(start_pos, end_pos)
                )
            else:
                return None
        chromosome, assembly = descr

        # Coordinate liftover
        lo = LiftOver(GRCH_TO_HG[assembly], 'hg38')
        liftover_start_i = self.uta.get_liftover(lo, chromosome, start_pos)
        if liftover_start_i is None:
            return None
        else:
            start_pos = liftover_start_i[1]

        if start_pos != end_pos:
            liftover_end_i = self.uta.get_liftover(lo, chromosome, end_pos)
            if liftover_end_i is None:
                return None
            else:
                end_pos = liftover_end_i[1]
        else:
            end_pos = start_pos

        newest_ac = self.uta.get_newest_assembly_ac(ac)
        if newest_ac:
            ac = newest_ac[0][0]
            if self._validate_index(ac, (start_pos, end_pos), 0):
                return dict(
                    ac=ac,
                    pos=(start_pos, end_pos)
                )

        return None

    def g_to_mane_c(self, ac, start_pos, end_pos, gene=None):
        """Return MANE Transcript on the c. coordinate.
        g->GRCh38->MANE c.
        If MANE c. cannot be found, we return the genomic coordinate on
            GRCh38

        :param str ac: Transcript accession on g. coordinate
        :param int start_pos: genomic change start position
        :param int end_pos: genomic change end position
        :param str gene: Gene symbol
        :return: MANE Transcripts with cDNA change on c. coordinate
        """
        if not self.uta.validate_genomic_ac(ac):
            logger.warning(f"Genomic accession does not exist: {ac}")
            return None

        if not gene:
            gene_symbol = self.uta.get_gene_from_ac(ac, start_pos, end_pos)
            if not gene_symbol:
                return None

            if len(gene_symbol) != 1:
                # Return GRCh38
                grch38 = self.g_to_grch38(ac, start_pos, end_pos)
                if not grch38:
                    return None

                return dict(
                    gene=None,
                    refseq=grch38['ac'],  # TODO: Is this always refseq
                    ensembl=None,
                    coding_start_site=None,
                    coding_end_site=None,
                    pos=grch38['pos'],
                    strand=None,
                    status='GRCh38'
                )

            # Exactly one gene found for an accession
            gene_symbol = gene_symbol[0][0]
        else:
            gene_symbol = gene

        mane_data =\
            self.mane_transcript_mappings.get_gene_mane_data(gene_symbol)
        if not mane_data:
            return None
        mane_data_len = len(mane_data)

        for i in range(mane_data_len):
            index = mane_data_len - i - 1
            current_mane_data = mane_data[index]

            mane_c_ac = current_mane_data['RefSeq_nuc']

            # Liftover to GRCh38
            grch38 = self.g_to_grch38(ac, start_pos, end_pos)
            mane_tx_genomic_data = None
            if grch38:
                # GRCh38 -> MANE C
                mane_tx_genomic_data = self.uta.get_mane_c_genomic_data(
                    mane_c_ac, None, grch38['pos'][0], grch38['pos'][1]
                )

            if not grch38 or not mane_tx_genomic_data:
                # GRCh38 did not work, so let's try original assembly
                mane_tx_genomic_data = self.uta.get_mane_c_genomic_data(
                    mane_c_ac, ac, start_pos, end_pos
                )
                if not mane_tx_genomic_data:
                    return None
                else:
                    logger.info("Not using most recent assembly")

            tx_pos_range = mane_tx_genomic_data['tx_pos_range']
            alt_pos_change = mane_tx_genomic_data['alt_pos_change']
            coding_start_site = mane_tx_genomic_data['coding_start_site']
            coding_end_site = mane_tx_genomic_data['coding_end_site']

            if mane_tx_genomic_data['strand'] == '-':
                alt_pos_change = (alt_pos_change[1] + 1, alt_pos_change[0] - 1)

            mane_c_pos_change = (
                tx_pos_range[0] + alt_pos_change[0] - coding_start_site,
                tx_pos_range[1] - alt_pos_change[1] - coding_start_site
            )

            if not self._validate_index(mane_c_ac, mane_c_pos_change,
                                        coding_start_site):
                logger.warning(f"{mane_c_pos_change} are not valid positions"
                               f" on {mane_c_ac}with coding start site "
                               f"{coding_start_site}")
                return None

            return self._get_mane_c(current_mane_data, mane_c_pos_change,
                                    (coding_start_site, coding_end_site))
