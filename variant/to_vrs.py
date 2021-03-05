"""Module for to VRS translation."""
from variant.classifiers import Classify
from variant.tokenizers import Tokenize
from variant.validators import Validate
from variant.translators import Translate
from variant.data_sources import SeqRepoAccess, TranscriptMappings
from variant.tokenizers import GeneSymbol
from variant.tokenizers.caches import GeneSymbolCache, AminoAcidCache


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

    def get_validations(self, q):
        """Return validation results for a given variant.

        :param str q: Variant to get validation results for
        :return: ValidationSummary for the variant
        """
        tokens = self.tokenizer.perform(q.strip())
        classifications = self.classifier.perform(tokens)
        validations = self.validator.perform(classifications)
        return validations

    def get_translations(self, validations):
        """Return a list translations from a ValidationSummary.

        :param ValidationSummary validations: Valid and Invalid results
        :return: A list of unique translations from valid results
        """
        translations = []
        for valid_variant in validations.valid_results:
            result = self.translator.perform(valid_variant)
            if result not in translations:
                translations.append(result)
        return translations
