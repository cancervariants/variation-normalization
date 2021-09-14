"""Module for hgvs_dup_del_mode in normalize endpoint."""
import logging
from typing import Optional, Dict
from variation.data_sources.seq_repo_access import SeqRepoAccess
from ga4gh.vrs.dataproxy import SeqRepoDataProxy
from ga4gh.vrs import models
from ga4gh.core import ga4gh_identify

logger = logging.getLogger('variation')
logger.setLevel(logging.DEBUG)


class HGVSDupDelMode:
    """HGVS Dup Del Mode class."""

    def __init__(self, seqrepo_access: SeqRepoAccess,
                 dp: SeqRepoDataProxy):
        """Initialize HGVS Dup Del Mode.

        :param SeqRepoAccess seqrepo_access: Access to seqrepo
        :param SeqRepoDataProxy dp: ga4gh dataproxy for seqrepo
        """
        self.seqrepo_access = seqrepo_access
        self.dp = dp

    def _get_chr(self, ac) -> Optional[str]:
        """Get chromosome for accession.

        :param str ac: Accession
        :return: Chromosome
        """
        aliases = self.seqrepo_access.aliases(ac)
        return ([a.split(':')[-1] for a in aliases
                 if a.startswith('GRCh') and '.' not in a and 'chr' not in a] or [None])[0]  # noqa: E501

    def default_mode(self, ac, alt_type, pos, del_or_dup, location,
                     chromosome=None, allele=None):
        """Use default characteristics to return a variation

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
        """
        variation = None
        if chromosome is None:
            chromosome = self._get_chr(ac)

        if chromosome is None:
            logger.warning(f"Unable to find chromosome on {ac}")
            return None

        if 'uncertain' in alt_type:
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
        elif pos and (pos[1] - pos[0] > 100):
            # TODO: RSE with a Derived Sequence Expr subject
            seq_expr = models.RepeatedSequenceExpression(
                seq_expr=models.DerivedSequenceExpression(
                    location=location,
                    reverse_complement=False
                ),
                count=models.Number(value=1)  # TODO: Change
            )

            # TODO: Is this a CNV or Allele?
            variation = models.Allele(
                location=self.seqrepo_access.seq_repo_client.translate_identifier(ac, 'ga4gh')[0],  # noqa: E501
                state=seq_expr
            )

        else:
            if allele:
                allele['state']['type'] = 'LiteralSequenceExpression'
                variation = models.Allele(**allele)
        if variation:
            variation = self._ga4gh_identify_variation(variation)
        return variation

    def cnv_mode(self, subject, copies):
        """CNV mode"""
        variation = models.CopyNumber(
            subject=subject,
            copies=copies
        )
        return self._ga4gh_identify_variation(variation)

    def repeated_seq_expr_mode(self):
        """RSE mode"""
        pass

    def literal_seq_expr_mode(self, location, state):
        """Literal seq expression mode"""
        variation = models.Allele(
            location=location,
            state=state
        )
        return self._ga4gh_identify_variation(variation)

    def _ga4gh_identify_variation(self, variation) -> Dict:
        """Return variation with GA4GH digest-based id.

        :param ga4gh.models.Variation variation: VRS variation object
        :return: VRS Variation with GA4GH digest-based id represented as a dict
        """
        variation._id = ga4gh_identify(variation)
        return variation.as_dict()
