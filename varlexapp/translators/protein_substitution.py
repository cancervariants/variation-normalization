"""Module for Protein Substitution Translation."""
from .translator import Translator
from ..models import ValidationResult
from ..models import ClassificationType
from ..models import VariantRepresentation
from ..models import ProteinSubstitutionToken


class ProteinSubstitution(Translator):
    """The Protein Substitution Translator class."""

    def translate(self, res: ValidationResult) -> VariantRepresentation:
        """Translate Protein Substitution to VRS representation."""
        psub_tokens = [t for t in res.classification.all_tokens if
                       isinstance(t, ProteinSubstitutionToken)]
        if len(psub_tokens) > 1:
            raise Exception('Should not have more than one ProteinSubstitution'
                            ' token if the result is valid')

        if not res.location:
            raise Exception('Cannot translate a variant with no location')

        state = res.location['state']
        _id = res.location['_id']
        vrs_type = res.location['type']

        for field in ['_id', 'state', 'type']:
            if field in res.location:
                del res.location[field]

        if 'location' in res.location:
            res.location = res.location['location']

        return VariantRepresentation(_id, res.location, state, vrs_type)

    def can_translate(self, type: ClassificationType) -> bool:
        """Return if classification type is Protein Substitution."""
        return type == ClassificationType.PROTEIN_SUBSTITUTION
