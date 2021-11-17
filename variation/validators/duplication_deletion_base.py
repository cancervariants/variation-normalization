"""The base class for Duplication and Deletion Validation."""
from typing import Optional, List, Dict
from variation.schemas.token_response_schema import Token
from variation.validators.validator import Validator
from variation.hgvs_dup_del_mode import HGVSDupDelMode
from variation.data_sources import SeqRepoAccess, TranscriptMappings, UTA
from variation.tokenizers import GeneSymbol
from variation.mane_transcript import MANETranscript
from ga4gh.vrs.dataproxy import SeqRepoDataProxy
from ga4gh.vrs.extras.translator import Translator
from gene.query import QueryHandler as GeneQueryHandler
import logging
from variation.schemas.normalize_response_schema\
    import HGVSDupDelMode as HGVSDupDelModeEnum

logger = logging.getLogger('variation')
logger.setLevel(logging.DEBUG)


class DuplicationDeletionBase(Validator):
    """The Deletion Validator Base class."""

    def __init__(self, seq_repo_access: SeqRepoAccess,
                 transcript_mappings: TranscriptMappings,
                 gene_symbol: GeneSymbol,
                 mane_transcript: MANETranscript,
                 uta: UTA, dp: SeqRepoDataProxy, tlr: Translator,
                 gene_normalizer: GeneQueryHandler):
        """Initialize the Deletion Base validator.

        :param SeqRepoAccess seq_repo_access: Access to SeqRepo data
        :param TranscriptMappings transcript_mappings: Access to transcript
            mappings
        :param GeneSymbol gene_symbol: Gene symbol tokenizer
        :param MANETranscript mane_transcript: Access MANE Transcript
            information
        :param UTA uta: Access to UTA queries
        :param GeneQueryHandler gene_normalizer: Access to gene-normalizer
        """
        super().__init__(
            seq_repo_access, transcript_mappings, gene_symbol, mane_transcript,
            uta, dp, tlr, gene_normalizer
        )
        self.hgvs_dup_del_mode = HGVSDupDelMode(seq_repo_access)

    def get_reference_sequence(self, ac, start, end, errors, cds_start=None)\
            -> Optional[str]:
        """Get deleted reference sequence.

        :param str ac: Accession
        :param int start: Start position
        :param int end: End position
        :param list errors: List of errors
        :param int cds_start: Coding start site
        :return: Reference sequence of nucleotides
        """
        if cds_start:
            start += cds_start
            if end is not None:
                end += cds_start

        if start and not end:
            ref_sequence = self.seqrepo_access.get_sequence(
                ac, start
            )
        elif start is not None and end is not None:
            ref_sequence = self.seqrepo_access.get_sequence(
                ac, start, end
            )
        else:
            ref_sequence = None

        if not ref_sequence:
            errors.append("Unable to get reference sequence.")
        return ref_sequence

    def check_reference_sequence(self, t, s, errors, cds_start=None) -> bool:
        """Check that reference sequence matches deleted sequence.

        :param str t: Accession
        :param Token s: Classification token
        :param list errors: List of errors
        :param int cds_start: Coding start site
        :return: `True` if ref_sequences matches deleted_sequence.
            `False` otherwise.
        """
        ref_sequence = self.get_reference_sequence(
            t, s.start_pos_del, s.end_pos_del, errors, cds_start=cds_start
        )

        if not errors and ref_sequence and s.deleted_sequence:
            if ref_sequence != s.deleted_sequence:
                errors.append(f"Expected deleted sequence {ref_sequence} "
                              f"but got {s.deleted_sequence}")

    def concise_description(self, transcript, token) -> str:
        """Return a HGVS description of the identified variation.

        :param str transcript: Transcript accession
        :param Token token: Classification token
        :return: HGVS expression
        """
        position = f"{token.start_pos_del}"
        if token.end_pos_del is not None:
            position += f"_{token.end_pos_del}"

        descr = f"{transcript}:{token.reference_sequence}.{position}del"
        if token.deleted_sequence:
            descr += f"{token.deleted_sequence}"
        return descr

    def add_normalized_genomic_dup_del(
            self, s: Token, t: str, start: int, end: int, gene: str,
            so_id: str, errors: List, hgvs_dup_del_mode: HGVSDupDelModeEnum,
            mane_data_found: Dict) -> None:
        """Add normalized genomic dup or del to mane data

        :param Token s: Classification token
        :param str t: Accession
        :param int start: Start position
        :param int end: ENd position
        :param str gene: Gene
        :param str so_id: Sequence ontology id
        :param List errors: List of errors
        :param HGVSDupDelModeEnum hgvs_dup_del_mode: Must be: `default`, `cnv`,
            `repeated_seq_expr`, `literal_seq_expr`.
            This parameter determines how to represent HGVS dup/del expressions
            as VRS objects.
        :param Dict mane_data_found: MANE Transcript information found
        """
        mane = self.mane_transcript.get_mane_transcript(
            t, start, end, s.reference_sequence, gene=gene,
            normalize_endpoint=True
        )

        if mane:
            s.reference_sequence = 'c'
            s.molecule_context = 'transcript'
            s.so_id = so_id

            allele = self.to_vrs_allele(
                mane['refseq'], mane['pos'][0], mane['pos'][1],
                s.reference_sequence, s.alt_type, errors,
                cds_start=mane['coding_start_site']
            )

            mane_variation = self.hgvs_dup_del_mode.interpret_variation(
                t, s.alt_type, allele, errors, hgvs_dup_del_mode)

            if mane_variation:
                self._add_dict_to_mane_data(
                    mane['refseq'], s, mane_variation,
                    mane_data_found, mane['status']
                )

    def validate_gene_or_accession_pos(self, t: str, pos_list: List,
                                       errors: List,
                                       gene: Optional[str] = None) -> None:
        """Validate positions on gene or accession.
        If not valid, add to list of errors

        :param str t: Accession
        :param List pos_list: List of positions to validate
        :param List errors: List of errors
        :param Optional[str] gene: Gene
        """
        if gene:
            len_pos_list = len(pos_list)
            pos1, pos2 = pos_list[0], pos_list[1]
            if len_pos_list == 2:
                pos3, pos4 = None, None
            elif len_pos_list == 4:
                pos3, pos4 = pos_list[2], pos_list[3]
            else:
                errors.append(f"Unexpected amount of positions:"
                              f" {len_pos_list}")
                return
            self._validate_gene_pos(
                gene, t, pos1, pos2, errors, pos3=pos3, pos4=pos4)
        else:
            for pos in pos_list:
                self._check_index(t, pos, errors)
