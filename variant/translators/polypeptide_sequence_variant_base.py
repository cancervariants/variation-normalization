"""Module for Polypeptide Sequence Variant Translation Base Class.."""
from abc import abstractmethod
from .translator import Translator
from variant.schemas.validation_response_schema import ValidationResult
from variant.schemas.classification_response_schema import ClassificationType
from variant.schemas.ga4gh_vrs import Allele


class PolypeptideSequenceVariantBase(Translator):
    """Polypeptide Sequence Variant Translation Base class."""

    def translate(self, res: ValidationResult) -> Allele:
        """Translate to VRS representation."""
        psub_tokens = [t for t in res.classification.all_tokens if
                       self.is_token_instance(t)]
        if len(psub_tokens) > 1:
            raise Exception(f'Should not have more than one '
                            f'{self.__class__.__name__} '
                            f'token if the result is valid')

        if not res.allele['location']:
            raise Exception('Cannot translate a variant with no location')

        return Allele(**res.allele)

    @abstractmethod
    def can_translate(self, type: ClassificationType) -> bool:
        """Check that the classification type matches."""
        raise NotImplementedError

    @abstractmethod
    def is_token_instance(self, token):
        """Check tthat the token is the correct instance."""
        raise NotImplementedError
