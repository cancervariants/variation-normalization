"""Module for to_vrs endpoint."""
from typing import Tuple, Optional, List, Union
from urllib.parse import unquote

from gene.query import QueryHandler as GeneQueryHandler
from ga4gh.vrsatile.pydantic.vrs_models import Allele, Haplotype, AbsoluteCopyNumber,\
    VariationSet, Text
from uta_tools.data_sources import SeqRepoAccess, TranscriptMappings, UTADatabase, \
    MANETranscriptMappings, MANETranscript

from variation.hgvs_dup_del_mode import HGVSDupDelMode
from variation.schemas.app_schemas import Endpoint
from variation.schemas.hgvs_to_copy_number_schema import VALID_CLASSIFICATION_TYPES,\
    RelativeCopyClass
from variation.schemas.token_response_schema import Nomenclature
from variation.schemas.validation_response_schema import ValidationSummary
from variation.classifiers import Classify
from variation.tokenizers import Tokenize
from variation.validators import Validate
from variation.translators import Translate
from variation.tokenizers import GeneSymbol
from variation.tokenizers.caches import AminoAcidCache
from variation.schemas.normalize_response_schema\
    import HGVSDupDelMode as HGVSDupDelModeEnum


class ToVRS:
    """The class for translating variation strings to VRS representations."""

    def __init__(self, tokenizer: Tokenize, classifier: Classify,
                 seqrepo_access: SeqRepoAccess,
                 transcript_mappings: TranscriptMappings,
                 gene_symbol: GeneSymbol, amino_acid_cache: AminoAcidCache,
                 uta: UTADatabase, mane_transcript_mappings: MANETranscriptMappings,
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
        :param UTADatabase uta: UTA DB and queries
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

    async def get_validations(
            self, q: str, endpoint_name: Optional[Endpoint] = None,
            hgvs_dup_del_mode: Optional[HGVSDupDelModeEnum] = None,  # noqa: E501
            baseline_copies: Optional[int] = None,
            relative_copy_class: Optional[RelativeCopyClass] = None,
            do_liftover: bool = False
    ) -> Tuple[Optional[ValidationSummary], Optional[List[str]]]:
        """Return validation results for a given variation.

        :param str q: variation to get validation results for
        :param Optional[Endpoint] endpoint_name: Then name of the endpoint being used
        :param Optional HGVSDupDelModeEnum hgvs_dup_del_mode: This parameter determines
            how to interpret HGVS dup/del expressions in VRS.
        :param Optional[int] baseline_copies: Baseline copies number
        :param Optional[RelativeCopyClass] relative_copy_class: The relative copy class
        :param bool do_liftover: Whether or not to liftover to GRCh38 assembly
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
            if endpoint_name in [Endpoint.NORMALIZE, Endpoint.HGVS_TO_ABSOLUTE_CN,
                                 Endpoint.HGVS_TO_RELATIVE_CN]:
                if hgvs_dup_del_mode:
                    hgvs_dup_del_mode = hgvs_dup_del_mode.strip().lower()
                    if endpoint_name == Endpoint.NORMALIZE:
                        if not self.hgvs_dup_del_mode.is_valid_dup_del_mode(hgvs_dup_del_mode):  # noqa: E501
                            warnings.append(
                                f"hgvs_dup_del_mode must be one of: "
                                f"{self.hgvs_dup_del_mode.valid_dup_del_modes}")
                            return None, warnings
                    else:
                        if not self.hgvs_dup_del_mode.is_valid_copy_number_mode(hgvs_dup_del_mode):  # noqa: E501
                            warnings.append(f"hgvs_dup_del_mode must be one of "
                                            f"{self.hgvs_dup_del_mode.valid_copy_number_modes}")  # noqa: E501
                            return None, warnings
                    if hgvs_dup_del_mode == HGVSDupDelModeEnum.ABSOLUTE_CNV:
                        if not baseline_copies:
                            warnings.append(f"{hgvs_dup_del_mode} mode "
                                            f"requires `baseline_copies`")
                            return None, warnings
                elif not hgvs_dup_del_mode and endpoint_name == Endpoint.NORMALIZE:
                    hgvs_dup_del_mode = HGVSDupDelModeEnum.DEFAULT
                else:
                    warnings.append(f"hgvs_dup_del_mode must be either "
                                    f"{HGVSDupDelModeEnum.ABSOLUTE_CNV.value} or "
                                    f"{HGVSDupDelModeEnum.RELATIVE_CNV.value}")
                    return None, warnings
            else:
                hgvs_dup_del_mode = HGVSDupDelModeEnum.DEFAULT

        classifications = self.classifier.perform(tokens)

        if endpoint_name in [Endpoint.HGVS_TO_ABSOLUTE_CN,
                             Endpoint.HGVS_TO_RELATIVE_CN]:
            tmp_classifications = []
            for c in classifications:
                conditions = (
                    c.classification_type in VALID_CLASSIFICATION_TYPES,
                    ("HGVS" in c.matching_tokens or "ReferenceSequence" in c.matching_tokens)  # noqa: E501
                )
                if all(conditions):
                    tmp_classifications.append(c)
            classifications = tmp_classifications
            if not classifications:
                warnings = [f"{q} is not a supported HGVS genomic "
                            f"duplication or deletion"]
                return None, warnings

        validations = await self.validator.perform(
            classifications, endpoint_name=endpoint_name, warnings=warnings,
            hgvs_dup_del_mode=hgvs_dup_del_mode, baseline_copies=baseline_copies,
            relative_copy_class=relative_copy_class, do_liftover=do_liftover
        )
        if not warnings:
            warnings = validations.warnings
        return validations, warnings

    def get_translations(self, validations: ValidationSummary,
                         warnings: List)\
            -> Tuple[Optional[Union[List[Allele], List[AbsoluteCopyNumber],
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
