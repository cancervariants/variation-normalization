from .translator import Translator
from ..models import ValidationResult
from ..models import ClassificationType
from ..models import VariantRepresentation
from ..models import SequenceState
from ..models import Location
from ..models import SimpleInterval
from ..models import ProteinSubstitutionToken

class ProteinSubstitution(Translator):
    def translate(self, res: ValidationResult) -> VariantRepresentation:
        psub_tokens = [t for t in res.classification.all_tokens if isinstance(t, ProteinSubstitutionToken)]
        if len(psub_tokens) > 1:
            raise Exception('Should not have more than one ProteinSubstitution token if the result is valid')

        if not res.location:
            raise Exception('Cannot translate a variant with no location')

        state = SequenceState(psub_tokens[0].alt_protein)

        return VariantRepresentation(res.location, state)


    def can_translate(self, type: ClassificationType) -> bool:
        return type == ClassificationType.PROTEIN_SUBSTITUTION
