from typing import List

from .validator import Validator
from ..models import ValidationResult, ClassificationType, Classification, ProteinSubstitutionToken, LookupType, Location, SimpleInterval
from ..data_sources import TranscriptMappings, SeqRepoAccess

class ProteinSubstitution(Validator):
    def __init__(self, seq_repo_client, transcript_mappings: TranscriptMappings) -> None:
        self.transcript_mappings = transcript_mappings
        self.seq_repo_client = seq_repo_client

    def validate(self, classification: Classification) -> List[ValidationResult]:
        results = list()
        errors = list()
        gene_tokens = [t for t in classification.all_tokens if t.token_type == 'GeneSymbol']
        psub_tokens = [t for t in classification.all_tokens if isinstance(t, ProteinSubstitutionToken)]

        if len(gene_tokens) > 1:
            errors.append('More than one gene symbol found for a single protein substitution')

        if len(errors) > 0:
            return [ValidationResult(classification, False, 0, None, '', '', errors)]

        transcripts = self.transcript_mappings.protein_transcripts(gene_tokens[0].token, LookupType.GENE_SYMBOL)

        if not transcripts:
            errors.append(f'No transcripts found for gene symbol {gene_tokens[0].token}')
            return [ValidationResult(classification, False, 0, None, '', '', errors)]

        for s in psub_tokens:
            for t in transcripts:
                valid = True
                errors = list()
                ref_protein = self.seq_repo_client.protein_at_position(t, s.pos)
                interval = SimpleInterval(s.pos, s.pos + 1)
                location = Location(f'ensembl:{t}', interval)
                if ref_protein != s.ref_protein:
                    errors.append(f'Needed to find {s.ref_protein} at position {s.pos} on {t} but found {ref_protein}')
                    valid = False
                if valid:
                    results.append(ValidationResult(
                        classification,
                        True,
                        1,
                        location,
                        self.concise_description(t, s),
                        self.human_description(t, s),
                        []
                    ))
                else:
                    results.append(ValidationResult(
                        classification,
                        False,
                        1,
                        location,
                        self.concise_description(t, s),
                        self.human_description(t, s),
                        errors
                    ))
                errors = list()

        return results


    def validates_classification_type(self, classification_type: ClassificationType) -> bool:
        return classification_type == ClassificationType.PROTEIN_SUBSTITUTION


    def concise_description(self, transcript, psub_token: ProteinSubstitutionToken) -> str:
        return f'{transcript} {psub_token.ref_protein}{psub_token.pos}{psub_token.alt_protein}'

    def human_description(self, transcript, psub_token: ProteinSubstitutionToken) -> str:
        return f'A protein substitution from {psub_token.ref_protein} to {psub_token.alt_protein} at position {psub_token.pos} on transcript {transcript}'

