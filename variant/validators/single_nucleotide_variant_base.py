"""The module for Single Nucleotide Variant Validation."""
from abc import abstractmethod
from .validator import Validator
from ga4gh.vrs import models
from ga4gh.core import ga4gh_identify
import logging

logger = logging.getLogger('variant')
logger.setLevel(logging.DEBUG)


class SingleNucleotideVariantBase(Validator):
    """The Single Nucleotide Variant Validator Base class."""

    def get_vrs_allele(self, sequence_id, s) -> dict:
        """Return VRS Allele object.

        :param str sequence_id: Sequence containing the sequence to be located
        :param Token s: A Classification token
        :return: A VRS Allele object as a dictionary
        """
        seq_location = models.SequenceLocation(
            sequence_id=sequence_id,
            interval=models.SimpleInterval(
                start=s.position - 1,
                end=s.position
            )
        )

        state = models.SequenceState(sequence=s.new_nucleotide)
        allele = models.Allele(location=seq_location, state=state)
        allele['_id'] = ga4gh_identify(allele)
        return allele.as_dict()

    @abstractmethod
    def get_hgvs_expr(self, classification, t, s, is_hgvs):
        """Return a HGVS expression.

        :param Classification classification: A classification for a list of
            tokens
        :param str t: Transcript retrieved from transcript mapping
        """
        raise NotImplementedError

    def get_allele_from_transcript(self, classification, t, s, errors):
        """Return allele from a given transcript.
        :param Classification s: Classification token
        :param str t: Transcript
        :param list errors: List of errors
        :return: Allele as a dictionary
        """
        hgvs_expr, _ = self.get_hgvs_expr(classification, t, s, False)
        return self.get_allele_from_hgvs(hgvs_expr, errors)

    def check_ref_nucleotide(self, ref_nuc, s, t, errors):
        """Assert that ref_nuc matches s.ref_nucleotide."""
        if ref_nuc != s.ref_nucleotide:
            errors.append(f'Needed to find {s.ref_nucleotide} at'
                          f' position {s.position} on {t}'
                          f' but found {ref_nuc}')

    def concise_description(self, transcript, token) -> str:
        """Return a description of the identified variant."""
        return f'{transcript} {token.ref_nucleotide}' \
               f'{token.position}{token.new_nucleotide}'
