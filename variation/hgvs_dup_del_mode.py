"""Module for hgvs_dup_del_mode in normalize endpoint."""
import logging
from typing import Optional, Dict, Tuple, List
from variation.data_sources.seq_repo_access import SeqRepoAccess
from ga4gh.vrs import models
from ga4gh.core import ga4gh_identify
from variation.schemas.normalize_response_schema\
    import HGVSDupDelMode as HGVSDupDelModeEnum

logger = logging.getLogger('variation')
logger.setLevel(logging.DEBUG)


class HGVSDupDelMode:
    """Class for handling how to interpret HGVS duplications and deletions."""

    def __init__(self, seqrepo_access: SeqRepoAccess) -> None:
        """Initialize HGVS Dup Del Mode.

        :param SeqRepoAccess seqrepo_access: Access to seqrepo
        """
        self.seqrepo_access = seqrepo_access
        self.valid_modes = [mode.value for mode in
                            HGVSDupDelModeEnum.__members__.values()]

    def is_valid_mode(self, mode: str) -> bool:
        """Determine if mode is a valid input.

        :param str mode: Entered mode
        :return: `True` if valid mode. `False` otherwise.
        """
        hgvs_dup_del_mode = mode.strip().lower()
        return hgvs_dup_del_mode in self.valid_modes

    def default_mode(self, ac: str, alt_type: str, pos: Tuple[int, int],
                     del_or_dup: str, location: Dict,
                     chromosome: str = None,
                     allele: Dict = None) -> Optional[Dict]:
        """Use default characteristics to return a variation.
        If endpoints are ambiguous: cnv
            handling X chromosome, make cnv a definite range with base 1-2
            handling Y chromosome, base of 1
            handling anything else, base of 2
        elif len del or dup > 100bp:
            repeated_seq_expr with a derived_seq_expr subject
        else:
            literal_seq_expr (normalized LiteralSequenceExpression Allele)

        :param str ac: Accession
        :param str alt_type: Alteration type
        :param tuple pos: start_pos, end_pos
        :param str del_or_dup: Must be either `del` or `dup`
        :param dict location: Sequence Location object
        :param str chromosome: Chromosome
        :param dict allele: VRS Allele object represented as a dict
        :return: VRS Variation object represented as a dict
        """
        if 'uncertain' in alt_type or 'range' in alt_type:
            variation = self.cnv_mode(ac, del_or_dup,
                                      location, chromosome=chromosome)
        elif pos and (pos[1] - pos[0] > 100):
            variation = self.repeated_seq_expr_mode(alt_type, location)
        else:
            variation = self.literal_seq_expr_mode(allele, alt_type)
        return variation

    def cnv_mode(self, ac: str, del_or_dup: str, location: Dict,
                 chromosome: str = None) -> Optional[Dict]:
        """Return a VRS Copy Number Variation.

        :param str ac: Accession
        :param str del_or_dup: Must be either `del` or `dup`
        :param dict location: VRS SequenceLocation
        :param str chromosome: Chromosome
        :return: VRS Copy Number object represented as a dict
        """
        if chromosome is None:
            chromosome = self.seqrepo_access.ac_to_chromosome(ac)

        if chromosome is None:
            logger.warning(f"Unable to find chromosome on {ac}")
            return None

        if chromosome == 'X':
            copies = models.DefiniteRange(
                min=0 if del_or_dup == 'del' else 2,
                max=1 if del_or_dup == 'del' else 3
            )
        elif chromosome == 'Y':
            copies = models.Number(
                value=0 if del_or_dup == 'del' else 2
            )
        else:
            # Chr 1-22
            copies = models.Number(
                value=1 if del_or_dup == 'del' else 3
            )

        variation = models.CopyNumber(
            subject=models.DerivedSequenceExpression(
                location=location,
                reverse_complement=False
            ),
            copies=copies
        )
        return self._ga4gh_identify_variation(variation)

    def repeated_seq_expr_mode(self, alt_type: str,
                               location: Dict) -> Optional[Dict]:
        """Return a VRS Allele with a RepeatedSequenceExpression.
        The RepeatedSequenceExpression subject will be a
            DerivedSequenceExpression.

        :param str alt_type: Alteration type
        :param dict location: VRS SequenceLocation
        :return: VRS Allele object represented as a dict
        """
        if 'range' in alt_type:
            # Ranges should return an error
            return None

        if alt_type == 'duplication':
            count = models.Number(value=2)
        elif alt_type == 'deletion':
            count = models.Number(value=0)
        else:
            return None

        seq_expr = models.RepeatedSequenceExpression(
            seq_expr=models.DerivedSequenceExpression(
                location=location,
                reverse_complement=False
            ),
            count=count
        )

        variation = models.Allele(
            location=location,
            state=seq_expr
        )
        return self._ga4gh_identify_variation(variation)

    def literal_seq_expr_mode(self, allele: Dict,
                              alt_type: str) -> Optional[Dict]:
        """Return a VRS Allele with a normalized LiteralSequenceExpression.

        :param dict allele: normalized VRS Allele object represented as a dict
        :param str alt_type: Alteration type
        :return: VRS Allele object represented as a dict
        """
        if 'range' in alt_type or 'uncertain' in alt_type:
            return None

        variation = models.Allele(**allele) if allele else None
        return self._ga4gh_identify_variation(variation)

    def _ga4gh_identify_variation(self,
                                  variation: models.Variation) -> Optional[Dict]:  # noqa: E501
        """Return variation with GA4GH digest-based id.

        :param models.Variation variation: VRS variation object
        :return: VRS Variation with GA4GH digest-based id represented as a dict
        """
        if variation is None:
            return None
        else:
            variation._id = ga4gh_identify(variation)
            return variation.as_dict()

    def interpret_variation(self, ac: str, alt_type: str, allele: Dict,
                            errors: List,
                            hgvs_dup_del_mode: HGVSDupDelModeEnum,
                            pos: Optional[Tuple[int, int]] = None) -> Dict:
        """Interpret variation using HGVSDupDelMode

        :param str ac: Accession
        :param str alt_type: Alteration type
        :param dict allele: VRS Allele object
        :param List errors: List of errors
        :param HGVSDupDelModeEnum hgvs_dup_del_mode: Mode to use for
            interpreting HGVS duplications and deletions
        :param Optional[Tuple[int, int]] pos: Position changes
        :return: VRS Variation object
        """
        if 'deletion' in alt_type:
            del_or_dup = 'del'
        else:
            del_or_dup = 'dup'
        variation = None
        if allele is None:
            errors.append("Unable to get Allele")
        else:
            if hgvs_dup_del_mode == HGVSDupDelModeEnum.DEFAULT:
                variation = self.default_mode(
                    ac, alt_type, pos, del_or_dup,
                    allele['location'], allele=allele
                )
            elif hgvs_dup_del_mode == HGVSDupDelModeEnum.CNV:
                variation = self.cnv_mode(ac, del_or_dup, allele['location'])
            elif hgvs_dup_del_mode == HGVSDupDelModeEnum.REPEATED_SEQ_EXPR:
                variation = self.repeated_seq_expr_mode(
                    alt_type, allele['location']
                )
            elif hgvs_dup_del_mode == HGVSDupDelModeEnum.LITERAL_SEQ_EXPR:
                variation = self.literal_seq_expr_mode(allele, alt_type)
            if not variation:
                errors.append("Unable to get VRS Variation")
        return variation
