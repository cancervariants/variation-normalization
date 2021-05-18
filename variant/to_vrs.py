"""Module for to VRS translation."""
from typing import Tuple, Optional, List
from variant.schemas.ga4gh_vrs import Allele
from variant.schemas.validation_response_schema import ValidationSummary
from variant.classifiers import Classify
from variant.tokenizers import Tokenize
from variant.validators import Validate
from variant.translators import Translate
from variant.data_sources import SeqRepoAccess, TranscriptMappings
from variant.tokenizers import GeneSymbol
from variant.tokenizers.caches import GeneSymbolCache, AminoAcidCache
from urllib.parse import unquote


class ToVRS:
    """The class for translating variant strings to VRS representations."""

    def __init__(self):
        """Initialize the ToVRS class."""
        self.tokenizer = Tokenize()
        self.classifier = Classify()
        self.seq_repo_access = SeqRepoAccess()
        self.transcript_mappings = TranscriptMappings()
        self.gene_symbol = GeneSymbol(GeneSymbolCache())
        self.amino_acid_cache = AminoAcidCache()
        self.validator = Validate(self.seq_repo_access,
                                  self.transcript_mappings, self.gene_symbol,
                                  self.amino_acid_cache)
        self.translator = Translate()

    def get_validations(self, q)\
            -> Tuple[ValidationSummary, Optional[List[str]]]:
        """Return validation results for a given variant.

        :param str q: Variant to get validation results for
        :return: ValidationSummary for the variant and list of warnings
        """
        warnings = list()
        tokens = self.tokenizer.perform(unquote(q.strip()), warnings)
        classifications = self.classifier.perform(tokens)
        validations = self.validator.perform(classifications, warnings)
        if not warnings:
            warnings = validations.warnings
        return validations, warnings

    def get_translations(self, validations, warnings)\
            -> Tuple[Optional[List[Allele]], Optional[List[str]]]:
        """Return a list translations from a ValidationSummary.

        :param ValidationSummary validations: Valid and Invalid results
        :param list warnings: List of warnings
        :return: A list of unique translations from valid results
        """
        translations = []
        if not warnings:
            warnings.append("Unable to validate variant")
        for valid_variant in validations.valid_results:
            result = self.translator.perform(valid_variant)
            if result not in translations:
                translations.append(result)
        return translations, warnings
