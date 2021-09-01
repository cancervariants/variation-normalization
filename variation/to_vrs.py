"""Module for to VRS translation."""
from typing import Tuple, Optional, List, Union
from ga4gh.vrsatile.pydantic.vrs_model import Allele, Haplotype, CopyNumber,\
    VariationSet, Text
from variation.schemas.validation_response_schema import ValidationSummary
from variation.classifiers import Classify
from variation.tokenizers import Tokenize
from variation.validators import Validate
from variation.translators import Translate
from variation.data_sources import SeqRepoAccess, TranscriptMappings, \
    UTA, MANETranscriptMappings
from variation.mane_transcript import MANETranscript
from variation.tokenizers import GeneSymbol
from variation.tokenizers.caches import AminoAcidCache
from urllib.parse import unquote


class ToVRS:
    """The class for translating variation strings to VRS representations."""

    def __init__(self, tokenizer: Tokenize, classifier: Classify,
                 seqrepo_access: SeqRepoAccess,
                 transcript_mappings: TranscriptMappings,
                 gene_symbol: GeneSymbol, amino_acid_cache: AminoAcidCache,
                 uta: UTA, mane_transcript_mappings: MANETranscriptMappings,
                 mane_transcript: MANETranscript, validator: Validate,
                 translator: Translate):
        """Initialize the ToVRS class."""
        self.tokenizer = tokenizer
        self.classifier = classifier
        self.seq_repo_access = seqrepo_access
        self.transcript_mappings = transcript_mappings
        self.gene_symbol = gene_symbol
        self.amino_acid_cache = amino_acid_cache
        self.uta = uta
        self.mane_transcript_mappings = mane_transcript_mappings
        self.mane_transcript = mane_transcript
        self.validator = validator
        self.translator = translator

    def get_validations(self, q, normalize_endpoint=False)\
            -> Tuple[Optional[ValidationSummary], Optional[List[str]]]:
        """Return validation results for a given variation.

        :param str q: variation to get validation results for
        :param bool normalize_endpoint: `True` if normalize endpoint is being
            used. `False` otherwise.
        :return: ValidationSummary for the variation and list of warnings
        """
        warnings = list()
        if q is None:
            return None, ["No variation entered"]
        tokens = self.tokenizer.perform(unquote(q.strip()), warnings)
        classifications = self.classifier.perform(tokens)
        validations = self.validator.perform(
            classifications, normalize_endpoint, warnings
        )
        if not warnings:
            warnings = validations.warnings
        return validations, warnings

    def get_translations(self, validations, warnings)\
            -> Tuple[Optional[Union[List[Allele], List[CopyNumber],
                                    List[Text], List[Haplotype],
                                    List[VariationSet]]],
                     Optional[List[str]]]:
        """Return a list translations from a ValidationSummary.

        :param ValidationSummary validations: Valid and Invalid results
        :param list warnings: List of warnings
        :return: A list of unique translations from valid results
        """
        translations = []
        if validations is not None:
            for valid_variation in validations.valid_results:
                result = self.translator.perform(valid_variation)
                if result not in translations:
                    translations.append(result)
            if not translations and not warnings:
                warnings.append("Unable to validate variation")
        return translations, warnings
