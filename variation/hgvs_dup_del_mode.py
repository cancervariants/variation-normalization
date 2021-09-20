"""Module for hgvs_dup_del_mode in normalize endpoint."""
import logging
from typing import Optional, Dict
from variation.data_sources.seq_repo_access import SeqRepoAccess
from ga4gh.vrs import models
from ga4gh.core import ga4gh_identify

logger = logging.getLogger('variation')
logger.setLevel(logging.DEBUG)


class HGVSDupDelMode:
    """Class for handling how to interpret HGVS duplications and deletions."""

    def __init__(self, seqrepo_access: SeqRepoAccess):
        """Initialize HGVS Dup Del Mode.

        :param SeqRepoAccess seqrepo_access: Access to seqrepo
        """
        self.seqrepo_access = seqrepo_access

    def _get_chr(self, ac) -> Optional[str]:
        """Get chromosome for accession.

        :param str ac: Accession
        :return: Chromosome
        """
        aliases = self.seqrepo_access.aliases(ac)
        return ([a.split(':')[-1] for a in aliases
                 if a.startswith('GRCh') and '.' not in a and 'chr' not in a] or [None])[0]  # noqa: E501

    def default_mode(self, ac, alt_type, pos, del_or_dup, location,
                     chromosome=None, allele=None) -> Optional[Dict]:
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
            # TODO: Check if ok to keep second condition
            variation = self.repeated_seq_expr_mode(alt_type, location)
        else:
            variation = self.literal_seq_expr_mode(allele)
        return variation

    def cnv_mode(self, ac, del_or_dup, location, chromosome=None)\
            -> Optional[Dict]:
        """Return a VRS Copy Number Variation.

        :param str ac: Accession
        :param str del_or_dup: Must be either `del` or `dup`
        :param dict location: VRS SequenceLocation
        :param str chromosome: Chromosome
        :return: VRS Copy Number object represented as a dict
        """
        if chromosome is None:
            chromosome = self._get_chr(ac)

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

    def repeated_seq_expr_mode(self, alt_type, location) -> Optional[Dict]:
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

    def literal_seq_expr_mode(self, allele) -> Optional[Dict]:
        """Return a VRS Allele with a normalized LiteralSequenceExpression.

        :param dict allele: normalized VRS Allele object represented as a dict
        :return: VRS Allele object represented as a dict
        """
        if allele:
            # TODO: I don't believe we'll have to do this once
            #  We are able to normalize LSE in vrs-python
            allele['state']['type'] = 'LiteralSequenceExpression'
            variation = models.Allele(**allele)
        else:
            variation = None
        return self._ga4gh_identify_variation(variation)

    def _ga4gh_identify_variation(self, variation) -> Optional[Dict]:
        """Return variation with GA4GH digest-based id.

        :param ga4gh.models.Variation variation: VRS variation object
        :return: VRS Variation with GA4GH digest-based id represented as a dict
        """
        if variation is None:
            return None
        else:
            variation._id = ga4gh_identify(variation)
            return variation.as_dict()
