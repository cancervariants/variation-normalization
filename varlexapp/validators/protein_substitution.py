from typing import List

from .validator import Validator
from ..models import ValidationResult, ClassificationType, Classification, ProteinSubstitutionToken
from .data_sources import SeqRepoAccess

class ProteinSubstitution(Validator):
    def __init__(self, seq_repo_client: SeqRepoAccess) -> None:
        self.seq_repo_client = seq_repo_client

    def validate(self, classification: Classification) -> List[ValidationResult]:
        results = list()
        errors = list()
        gene_tokens = [t for t in classification.all_tokens if t.token_type == 'GeneSymbol']
        psub_tokens = [t for t in classification.all_tokens if isinstance(t, ProteinSubstitutionToken)]

        if len(gene_tokens) > 1:
            errors.append('More than one gene symbol found for a single protein substitution')

        if len(errors) > 0:
            return [ValidationResult(classification, False, 0, '', '', errors)]

        transcripts = self.seq_repo_client.transcripts_for_gene_symbol(gene_tokens[0].token)

        for s in psub_tokens:
            for t in transcripts:
                valid = True
                errors = list()
                tx_protein = self.seq_repo_client.protein_at_position(t, s.pos)
                if tx_protein != s.ref_protein:
                    errors.append(f'Needed to find {s.ref_protein} at position {s.pos} on {t["tx_ac"]} but found {tx_protein}')
                    valid = False
                if valid:
                    results.append(ValidationResult(
                        classification,
                        True,
                        1,
                        self.concise_description(t, s),
                        self.human_description(t, s),
                        []
                    ))
                else:
                    results.append(ValidationResult(
                        classification,
                        False,
                        1,
                        self.concise_description(t, s),
                        self.human_description(t, s),
                        errors
                    ))
                errors = list()

        return results


    def validates_classification_type(self, classification_type: ClassificationType) -> bool:
        return classification_type == ClassificationType.PROTEIN_SUBSTITUTION


    def concise_description(self, transcsript, psub_token: ProteinSubstitutionToken) -> str:
        return f'{transcsript["tx_ac"]} {psub_token.ref_protein}{psub_token.pos}{psub_token.alt_protein}'

    def human_description(self, transcsript, psub_token: ProteinSubstitutionToken) -> str:
        return f'A protein substitution from {psub_token.ref_protein} to {psub_token.alt_protein} at position {psub_token.pos} on transcript {transcsript["tx_ac"]}'

