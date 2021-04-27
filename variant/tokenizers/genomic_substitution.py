"""A module for Genomic Substitution Tokenization."""
from variant.schemas.token_response_schema import GenomicSubstitutionToken
from .single_nucleotide_variant_base import SingleNucleotideVariantBase


class GenomicSubstitution(SingleNucleotideVariantBase):
    """Class for tokenizing SNV Substitution at the linear genomic
    reference sequence.
    """

    def return_token(self, params):
        """Return Genomic Substitution token."""
        if self.sub['reference_sequence'] == 'g' and \
                self.sub['ref_nucleotide'] is not None and \
                self.sub['new_nucleotide'] != '=':
            return GenomicSubstitutionToken(**params)
