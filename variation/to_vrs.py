"""Module for to VRS translation."""
from typing import Tuple, Optional, List
from variation.schemas.ga4gh_vrs import Allele
from variation.schemas.validation_response_schema import ValidationSummary
from variation.classifiers import Classify
from variation.tokenizers import Tokenize
from variation.validators import Validate
from variation.translators import Translate
from variation.data_sources import SeqRepoAccess, TranscriptMappings, \
    UTA, MANETranscriptMappings
from variation.mane_transcript import MANETranscript
from variation.tokenizers import GeneSymbol
from variation.tokenizers.caches import GeneSymbolCache, AminoAcidCache
from urllib.parse import unquote


class ToVRS:
    """The class for translating variation strings to VRS representations."""

    def __init__(self):
        """Initialize the ToVRS class."""
        self.tokenizer = Tokenize()
        self.classifier = Classify()
        self.seq_repo_access = SeqRepoAccess()
        self.transcript_mappings = TranscriptMappings()
        self.gene_symbol = GeneSymbol(GeneSymbolCache())
        self.amino_acid_cache = AminoAcidCache()
        self.uta = UTA()
        self.mane_transcript_mappings = MANETranscriptMappings()
        self.mane_transcript = MANETranscript(
            self.seq_repo_access, self.transcript_mappings,
            self.mane_transcript_mappings, self.uta
        )
        self.validator = Validate(self.seq_repo_access,
                                  self.transcript_mappings, self.gene_symbol,
                                  self.mane_transcript, self.uta,
                                  self.amino_acid_cache)
        self.translator = Translate()

    def get_validations(self, q, normalize_endpoint=False)\
            -> Tuple[ValidationSummary, Optional[List[str]]]:
        """Return validation results for a given variation.

        :param str q: variation to get validation results for
        :param bool normalize_endpoint: `True` if normalize endpoint is being
            used. `False` otherwise.
        :return: ValidationSummary for the variation and list of warnings
        """
        warnings = list()
        tokens = self.tokenizer.perform(unquote(q.strip()), warnings)
        classifications = self.classifier.perform(tokens)
        validations = self.validator.perform(
            classifications, normalize_endpoint, warnings
        )
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
        for valid_variation in validations.valid_results:
            result = self.translator.perform(valid_variation)
            if result not in translations:
                translations.append(result)
        if not translations and not warnings:
            warnings.append("Unable to validate variation")
        return translations, warnings
