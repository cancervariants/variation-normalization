"""Module for translation."""
from abc import ABC, abstractmethod
from variant.schemas.ga4gh_vrs import Allele
from variant.schemas.validation_response_schema import ValidationResult
from variant.schemas.classification_response_schema import ClassificationType


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

        if len(instance_tokens) > 1:
            raise Exception(f'Should not have more than one '
                            f'{self.__class__.__name__} '
                            f'token if the result is valid')

        if not res.allele['location']:
            raise Exception('Cannot translate a variant with no location')

        return Allele(**res.allele)
