"""Module for to VRS translation."""
from typing import Tuple, Optional, List, Union
from ga4gh.vrsatile.pydantic.vrs_models import Allele, Haplotype, CopyNumber,\
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
from gene.query import QueryHandler as GeneQueryHandler
from variation.schemas.normalize_response_schema\
    import HGVSDupDelMode as HGVSDupDelModeEnum


class ToVRS:
    """The class for translating variation strings to VRS representations."""

    def __init__(self, tokenizer: Tokenize, classifier: Classify,
                 seqrepo_access: SeqRepoAccess,
                 transcript_mappings: TranscriptMappings,
                 gene_symbol: GeneSymbol, amino_acid_cache: AminoAcidCache,
                 uta: UTA, mane_transcript_mappings: MANETranscriptMappings,
                 mane_transcript: MANETranscript, validator: Validate,
                 translator: Translate,
                 gene_normalizer: GeneQueryHandler) -> None:
        """Initialize the ToVRS class.

        :param Tokenize tokenizer: Tokenizer class for tokenizing
        :param Classify classifier: Classifier class for classifying tokens
        :param SeqRepoAccess seqrepo_access: Access to SeqRepo
        :param TranscriptMappings transcript_mappings: Transcript mappings
            data class
        :param GeneSymbol gene_symbol: Class for identifying gene symbols
        :param AminoAcidCache amino_acid_cache: Amino Acid data class
        :param UTA uta: UTA DB and queries
        :param MANETranscriptMappings mane_transcript_mappings: Class for
            getting mane transcript data from gene
        :param MANETranscript mane_transcript: Mane transcript data class
        :param Validate validator: Validator class for validating valid inputs
        :param Translate translator: Translating valid inputs
        :param GeneQueryHandler gene_normalizer: Gene normalizer access
        """
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
        self.gene_normalizer = gene_normalizer

    def get_validations(
            self, q: str, normalize_endpoint: bool = False,
            hgvs_dup_del_mode: HGVSDupDelModeEnum = HGVSDupDelModeEnum.DEFAULT
    ) -> Tuple[Optional[ValidationSummary], Optional[List[str]]]:
        """Return validation results for a given variation.

        :param str q: variation to get validation results for
        :param bool normalize_endpoint: `True` if normalize endpoint is being
            used. `False` otherwise.
        :param HGVSDupDelModeEnum hgvs_dup_del_mode: Must be: `default`, `cnv`,
            `repeated_seq_expr`, `literal_seq_expr`.
            This parameter determines how to interpret HGVS dup/del expressions
            in VRS.
        :return: ValidationSummary for the variation and list of warnings
        """
        warnings = list()
        if q is None:
            return None, ["No variation entered"]
        tokens = self.tokenizer.perform(unquote(q.strip()), warnings)
        classifications = self.classifier.perform(tokens)
        validations = self.validator.perform(
            classifications, normalize_endpoint, warnings,
            hgvs_dup_del_mode=hgvs_dup_del_mode
        )
        if not warnings:
            warnings = validations.warnings
        return validations, warnings

    def get_translations(self, validations: ValidationSummary,
                         warnings: List)\
            -> Tuple[Optional[Union[List[Allele], List[CopyNumber],
                                    List[Text], List[Haplotype],
                                    List[VariationSet]]],
                     Optional[List[str]]]:
        """Return a list translations from a ValidationSummary.

        :param ValidationSummary validations: Valid and Invalid results
        :param List warnings: List of warnings
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
