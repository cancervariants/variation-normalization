"""Module for translation."""
from abc import ABC, abstractmethod
from variation.schemas.ga4gh_vrs import Allele
from variation.schemas.validation_response_schema import ValidationResult
from variation.schemas.classification_response_schema import ClassificationType


class Translator(ABC):
    """The translation class."""

    @abstractmethod
    def can_translate(self, type: ClassificationType) -> bool:
        """Determine if it's possible to translate a classification."""
        raise NotImplementedError

    @abstractmethod
    def is_token_instance(self, token):
        """Check tthat the token is the correct instance."""
        raise NotImplementedError

    def translate(self, res: ValidationResult) -> Allele:
        """Translate to VRS representation."""
        instance_tokens = [t for t in res.classification.all_tokens if
                           self.is_token_instance(t)]
        len_instance_tokens = len(instance_tokens)

        if len_instance_tokens > 1:
            if len_instance_tokens == 2:
                if {t.token_type for t in instance_tokens} \
                        == {'PolypeptideTruncation'}:
                    tokens = sorted([t.token.lower() for t in instance_tokens],
                                    key=len)
                    if len(tokens[1]) > len(tokens[0]):
                        t = f"{tokens[0]} ({tokens[0].replace('ter', '*')})"
                        if t.lower() == tokens[1].lower():
                            return Allele(**res.allele)

            raise Exception(f'Should not have more than one '
                            f'{self.__class__.__name__} '
                            f'token if the result is valid')

        if not res.allele['location']:
            raise Exception('Cannot translate a variation with no location')

        return Allele(**res.allele)
