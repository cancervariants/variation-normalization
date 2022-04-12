"""Module for translation."""
from abc import ABC, abstractmethod
from typing import Dict, Optional

from ga4gh.vrsatile.pydantic.vrs_models import Allele, AbsoluteCopyNumber, \
    RelativeCopyNumber
from pydantic.error_wrappers import ValidationError

from variation.schemas.validation_response_schema import ValidationResult
from variation.schemas.classification_response_schema import ClassificationType
from variation.schemas.token_response_schema import Token


class Translator(ABC):
    """The translation class."""

    @abstractmethod
    def can_translate(self, type: ClassificationType) -> bool:
        """Determine if it's possible to translate a classification."""
        raise NotImplementedError

    @abstractmethod
    def is_token_instance(self, token: Token) -> bool:
        """Check that the token is the correct instance for a translator."""
        raise NotImplementedError

    def translate(self, res: ValidationResult) -> Optional[Dict]:
        """Translate to VRS Variation representation."""
        instance_tokens = [t for t in res.classification.all_tokens if
                           self.is_token_instance(t)]
        len_instance_tokens = len(instance_tokens)

        variation_type = res.variation["type"]
        if len_instance_tokens > 1:
            if len_instance_tokens == 2:
                if {t.token_type for t in instance_tokens} \
                        == {"PolypeptideTruncation"}:
                    tokens = sorted([t.token.lower() for t in instance_tokens],
                                    key=len)
                    if len(tokens[1]) > len(tokens[0]):
                        t = f"{tokens[0]} ({tokens[0].replace('ter', '*')})"
                        if t.lower() == tokens[1].lower():
                            if variation_type == "Allele":
                                try:
                                    Allele(**res.variation)
                                except ValidationError:
                                    return None
                                else:
                                    return res.variation

            raise Exception(f"Should not have more than one "
                            f"{self.__class__.__name__} "
                            f"token if the result is valid")

        if variation_type == "Allele":
            try:
                Allele(**res.variation)
            except ValidationError:
                variation = None
            else:
                variation = res.variation
        elif variation_type == "AbsoluteCopyNumber":
            try:
                AbsoluteCopyNumber(**res.variation)
            except ValidationError:
                variation = None
            else:
                variation = res.variation
        elif variation_type == "RelativeCopyNumber":
            try:
                RelativeCopyNumber(**res.variation)
            except ValidationError:
                variation = None
            else:
                variation = res.variation
        else:
            raise Exception(f"{variation_type} not supported in "
                            f"Variation Normalization")
        return variation
