"""A module for Protein Substitution Tokenization."""
import re
from typing import Optional
from .caches import AminoAcidCache
from .tokenizer import Tokenizer
from ..models import ProteinSubstitutionToken


class ProteinSubstitution(Tokenizer):
    """Class for tokenizing Protein Substitution."""

    def __init__(self, amino_acid_cache: AminoAcidCache) -> None:
        """Initialize the Protein Substitution Class.

        :param AminoAcidCache amino_acid_cache: Valid amino acid codes.
        """
        self.amino_acid_cache = amino_acid_cache
        self.splitter = re.compile(r'(\d+)')
        self.psub = {
            'amino_acid': None,
            'position': None,
            'new_amino_acid': None
        }

    def match(self, input_string) -> Optional[ProteinSubstitutionToken]:
        """Return a ProteinSubstitutionToken match if one exists.

        :param str input_string: The input string to match
        :return: A ProteinSubstitutionToken if a match exists. Otherwise, None.
        """
        if input_string is None:
            return None

        if ':p.' in input_string:
            psub_parts = self.splitter.split(input_string.split(':')[-1])
        else:
            psub_parts = self.splitter.split(input_string)

        self._get_psub(psub_parts)

        if all(self.psub.values()):
            if not self._is_valid_amino_acid(self.psub['amino_acid'],
                                             self.psub['new_amino_acid']):
                return None
            return ProteinSubstitutionToken(input_string,
                                            self.psub['amino_acid'],
                                            self.psub['new_amino_acid'],
                                            self.psub['position'])
        return None

    def _get_psub(self, psub_parts):
        """Return Protein Substitution tokens.

        :param list psub_parts: The split input string
        """
        # missense
        if len(psub_parts) == 3:
            if self._is_valid_amino_acid(psub_parts[0], psub_parts[2]):
                self._set_psub(psub_parts[0], psub_parts[1], psub_parts[2])
            elif 'p.' in psub_parts[0]:
                # Predicted consequences given in parentheses
                if '(' in psub_parts[0] and ')' in psub_parts[-1]:
                    self._set_psub(
                        psub_parts[0].split('p.')[-1].split('(')[-1],
                        psub_parts[1],
                        psub_parts[-1].split(')')[0]
                    )
                else:
                    self._set_psub(psub_parts[0].split('p.')[-1],
                                   psub_parts[1],
                                   psub_parts[-1])

    def _set_psub(self, amino_acid, position, new_amino_acid):
        self.psub['amino_acid'] = amino_acid
        self.psub['position'] = position
        self.psub['new_amino_acid'] = new_amino_acid

    def _is_valid_amino_acid(self, amino_acid, new_amino_acid):
        """Return whether or not amino acids are valid.

        :param str amino_acid: Reference amino acid
        :param str new_amino_acid: New amino acid
        """
        return (amino_acid in self.amino_acid_cache and
                new_amino_acid in self.amino_acid_cache)
