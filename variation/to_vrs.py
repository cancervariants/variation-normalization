"""Module for toVRS endpoint."""
from typing import Tuple, Optional, List, Union
from ga4gh.vrsatile.pydantic.vrs_models import Allele, Haplotype, CopyNumber,\
    VariationSet, Text

from variation.hgvs_dup_del_mode import HGVSDupDelMode
from variation.schemas.app_schemas import Endpoint
from variation.schemas.hgvs_to_copy_number_schema import CopyNumberType
from variation.schemas.token_response_schema import Nomenclature
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
                 gene_normalizer: GeneQueryHandler,
                 hgvs_dup_del_mode: HGVSDupDelMode) -> None:
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
        :param HGVSDupDelMode hgvs_dup_del_mode: Class for handling
            HGVS dup/del expressions
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
        self.hgvs_dup_del_mode = hgvs_dup_del_mode

    def get_validations(
            self, q: str, endpoint_name: Optional[Endpoint] = None,
            hgvs_dup_del_mode: Optional[Union[HGVSDupDelModeEnum, CopyNumberType]] = None,  # noqa: E501
            baseline_copies: Optional[int] = None
    ) -> Tuple[Optional[ValidationSummary], Optional[List[str]]]:
        """Return validation results for a given variation.

        :param str q: variation to get validation results for
        :param Optional[Endpoint] endpoint_name: Then name of the endpoint being used
        :param Optional[Union[HGVSDupDelModeEnum, CopyNumberType]] hgvs_dup_del_mode:
            This parameter determines how to interpret HGVS dup/del expressions
            in VRS.
        :param Optional[int] baseline_copies: Baseline copies number
        :return: ValidationSummary for the variation and list of warnings
        """
        warnings = list()
        if q is None:
            return None, ["No variation entered"]
        tokens = self.tokenizer.perform(unquote(q.strip()), warnings)

        # gnomad vcf should always be a literal seq expression (allele)
        nomenclature = {t.nomenclature for t in tokens}
        if Nomenclature.GNOMAD_VCF in nomenclature:
            hgvs_dup_del_mode = HGVSDupDelModeEnum.LITERAL_SEQ_EXPR
        else:
            if endpoint_name == Endpoint.NORMALIZE:
                if hgvs_dup_del_mode:
                    hgvs_dup_del_mode = hgvs_dup_del_mode.strip().lower()
                    if not self.hgvs_dup_del_mode.is_valid_dup_del_mode(hgvs_dup_del_mode):  # noqa: E501
                        warnings.append(
                            f"hgvs_dup_del_mode must be one of: "
                            f"{self.hgvs_dup_del_mode.valid_dup_del_modes}")
                        return None, warnings
                else:
                    hgvs_dup_del_mode = HGVSDupDelModeEnum.DEFAULT
            elif endpoint_name in [Endpoint.HGVS_TO_ABSOLUTE_CN,
                                   Endpoint.HGVS_TO_RELATIVE_CN]:
                if not hgvs_dup_del_mode:
                    warnings.append(f"hgvs_dup_del_mode must be either "
                                    f"{CopyNumberType.ABSOLUTE.value} or "
                                    f"{CopyNumberType.RELATIVE.value}")
                    return None, warnings
                if not self.hgvs_dup_del_mode.is_valid_copy_number_mode(hgvs_dup_del_mode):  # noqa: E501
                    warnings.append(f"hgvs_dup_del_mode must be one of "
                                    f"{self.hgvs_dup_del_mode.valid_copy_number_modes}")
            else:
                hgvs_dup_del_mode = HGVSDupDelModeEnum.DEFAULT

        classifications = self.classifier.perform(tokens)
        validations = self.validator.perform(
            classifications, endpoint_name=endpoint_name, warnings=warnings,
            hgvs_dup_del_mode=hgvs_dup_del_mode, baseline_copies=baseline_copies
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
