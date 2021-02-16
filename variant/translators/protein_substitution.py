"""Module for Protein Substitution Translation."""
from .translator import Translator
from variant.schemas.validation_response_schema import ValidationResult
from variant.schemas.classification_response_schema import ClassificationType
from variant.schemas.token_response_schema import ProteinSubstitutionToken
from variant.schemas.ga4gh_vrs import Allele


class ProteinSubstitution(Translator):
    """The Protein Substitution Translator class."""

    def translate(self, res: ValidationResult) -> Allele:
        """Translate Protein Substitution to VRS representation."""
        psub_tokens = [t for t in res.classification.all_tokens if
                       isinstance(t, ProteinSubstitutionToken)]
        if len(psub_tokens) > 1:
            raise Exception('Should not have more than one ProteinSubstitution'
                            ' token if the result is valid')

        if not res.allele['location']:
            raise Exception('Cannot translate a variant with no location')

        return Allele(**res.allele)

    def can_translate(self, type: ClassificationType) -> bool:
        """Return if classification type is Protein Substitution."""
        return type == ClassificationType.PROTEIN_SUBSTITUTION
