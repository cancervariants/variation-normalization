"""The module for Single Nucleotide Variation Validation."""
from abc import abstractmethod
from .validator import Validator
from ga4gh.vrs import models
from ga4gh.core import ga4gh_identify
import logging

logger = logging.getLogger('variation')
logger.setLevel(logging.DEBUG)


class SingleNucleotideVariationBase(Validator):
    """The Single Nucleotide Variation Validator Base class."""

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
        allele = allele.as_dict()
        allele_seq_id = allele['location']['sequence_id']
        if allele_seq_id.startswith('ga4gh:GS.'):
            allele['location']['sequence_id'] = \
                allele_seq_id.replace('ga4gh:GS.', 'ga4gh:SQ.')
        return allele

    def silent_mutation_valid_invalid_results(self, classification_tokens,
                                              transcripts, classification,
                                              results, gene_tokens,
                                              normalize_endpoint) -> None:
        """Add validation result objects to a list of results for
        Silent Mutations.

        :param list classification_tokens: A list of Tokens
        :param list transcripts: A list of transcript strings
        :param Classification classification: A classification for a list of
            tokens
        :param list results: A list to store validation result objects
        :param list gene_tokens: List of GeneMatchTokens
        :param bool normalize_endpoint: `True` if normalize endpoint is being
            used. `False` otherwise.
        """
        valid_alleles = list()
        if 'HGVS' in classification.matching_tokens:
            is_hgvs = True
        else:
            is_hgvs = False

        mane_data = {
            'mane_select': dict(),
            'mane_plus_clinical': dict(),
            'longest_compatible_remaining': dict()
        }

        for s in classification_tokens:
            for t in transcripts:
                errors = list()

                if 'HGVS' in classification.matching_tokens:
                    # TODO: How to convert ENST_ to NM_ versioned
                    hgvs_expr = self.get_hgvs_expr(classification, t, s, True)
                else:
                    hgvs_expr = self.get_hgvs_expr(classification, t, s, False)
                t = hgvs_expr.split(':')[0]
                allele = None

                mane = self.mane_transcript.get_mane_transcript(
                    t, s.position, s.position, s.reference_sequence,
                    ref=s.ref_nucleotide, normalize_endpoint=normalize_endpoint
                )
                if mane:
                    if not gene_tokens:
                        gene_tokens.append(
                            self._gene_matcher.match(mane['gene'])
                        )

                    s.molecule_context = 'transcript'

                    mane_hgvs_expr =\
                        f"{mane['refseq']}:c.{mane['pos'][0]}="
                    self.add_mane_data(mane_hgvs_expr, mane, mane_data, s)

                try:
                    sequence_id = \
                        self.dp.translate_sequence_identifier(t, 'ga4gh')[0]
                except KeyError:
                    errors.append(f"GA4GH Data Proxy unable to translate "
                                  f"sequence identifier {t}")
                else:
                    s.new_nucleotide = \
                        self.seqrepo_access.sequence_at_position(t, s.position)
                    if s.new_nucleotide:
                        allele = self.get_vrs_allele(sequence_id, s)

                if allele:
                    len_of_seq = self.seqrepo_access.len_of_sequence(t)
                    if len_of_seq < s.position - 1:
                        errors.append('Sequence index error')

                self.add_validation_result(
                    allele, valid_alleles, results,
                    classification, s, t, gene_tokens, errors
                )

                if is_hgvs:
                    break
        self.add_mane_to_validation_results(
            mane_data, valid_alleles, results, classification, gene_tokens
        )

    @abstractmethod
    def get_hgvs_expr(self, classification, t, s, is_hgvs):
        """Return a HGVS expression.

        :param Classification classification: A classification for a list of
            tokens
        :param str t: Transcript retrieved from transcript mapping
        :param Token s: The classification token
        :param bool is_hgvs: Whether or not classification is HGVS token
        :return: hgvs expression
        """
        raise NotImplementedError

    def check_ref_nucleotide(self, ref_nuc, s, t, errors):
        """Assert that ref_nuc matches s.ref_nucleotide."""
        if ref_nuc != s.ref_nucleotide:
            errors.append(f'Needed to find {s.ref_nucleotide} at'
                          f' position {s.position} on {t}'
                          f' but found {ref_nuc}')

    def concise_description(self, transcript, token) -> str:
        """Return a HGVS description of the identified variation.

        :param str transcript: Transcript accession
        :param Token token: Classification token
        :return: HGVS expression
        """
        prefix = f'{transcript}:{token.reference_sequence}.{token.position}'
        if token.new_nucleotide == '=':
            change = "="
        else:
            change = f"{token.ref_nucleotide}>{token.new_nucleotide}"
        return prefix + change
