"""Module for to_vrs endpoint."""
from typing import Tuple, Optional, List, Union, Dict
from urllib.parse import unquote
from datetime import datetime

from ga4gh.vrsatile.pydantic.vrs_models import Allele, Haplotype, AbsoluteCopyNumber,\
    VariationSet, Text
from ga4gh.vrs import models
from ga4gh.core import ga4gh_identify
from uta_tools.schemas import ResidueMode
from uta_tools.data_sources import SeqRepoAccess
from ga4gh.vrs.dataproxy import SeqRepoDataProxy

from variation.schemas.normalize_response_schema\
    import HGVSDupDelMode as HGVSDupDelModeEnum, ServiceMeta
from variation.hgvs_dup_del_mode import HGVSDupDelMode
from variation.schemas.app_schemas import Endpoint
from variation.schemas.hgvs_to_copy_number_schema import VALID_CLASSIFICATION_TYPES,\
    RelativeCopyClass
from variation.schemas.to_vrs_response_schema import ToVRSService
from variation.schemas.token_response_schema import Nomenclature
from variation.schemas.validation_response_schema import ValidationSummary
from variation.classifiers import Classify
from variation.tokenizers import Tokenize
from variation.validators import Validate
from variation.translators import Translate
from variation.vrs_representation import VRSRepresentation
from variation.version import __version__


class ToVRS(VRSRepresentation):
    """The class for translating variation strings to VRS representations."""

    def __init__(self, seqrepo_access: SeqRepoAccess, dp: SeqRepoDataProxy,
                 tokenizer: Tokenize, classifier: Classify, validator: Validate,
                 translator: Translate, hgvs_dup_del_mode: HGVSDupDelMode) -> None:
        """Initialize the ToVrsAndVrsatile class.

        :param SeqRepoAccess seqrepo_access: Access to SeqRepo via UTA Tools
        :param SeqRepoDataProxy dp: Access to SeqRepo via VRS Python
        :param Tokenize tokenizer: Tokenizer class for tokenizing
        :param Classify classifier: Classifier class for classifying tokens
        :param Validate validator: Validator class for validating valid inputs
        :param Translate translator: Translating valid inputs
        :param HGVSDupDelMode hgvs_dup_del_mode: Class for handling
            HGVS dup/del expressions
        """
        super().__init__(dp, seqrepo_access)
        self.tokenizer = tokenizer
        self.classifier = classifier
        self.validator = validator
        self.translator = translator
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

    def get_ref_allele_seq(self, allele: Dict,
                           identifier: str) -> Optional[str]:
        """Return ref allele seq for transcript.

        :param Dict allele: VRS Allele object
        :param str identifier: Identifier for allele
        :return: Ref seq allele
        """
        start = None
        end = None
        interval = allele["location"]["interval"]
        ival_type = interval["type"]
        if ival_type == "SequenceInterval":
            if interval["start"]["type"] == "Number":
                start = interval["start"]["value"]
                end = interval["end"]["value"]

                if start == end:
                    return None

        if start is None and end is None:
            return None

        ref, _ = self.seqrepo_access.get_reference_sequence(
            identifier, start, end, residue_mode=ResidueMode.INTER_RESIDUE)

        return ref

    async def to_vrs(self, q: str,
                     untranslatable_returns_text: bool = False) -> ToVRSService:
        """Return a VRS-like representation of all validated variations for a query.  # noqa: E501

        :param str q: The variation to translate
        :param bool untranslatable_returns_text: `True` return VRS Text Object when
            unable to translate or normalize query. `False` returns empty list when
            unable to translate or normalize query.
        :return: ToVRSService containing VRS variations and warnings
        """
        validations, warnings = await self.get_validations(q)
        translations, warnings = self.get_translations(validations, warnings)

        if not translations:
            if untranslatable_returns_text and q and q.strip():
                text = models.Text(definition=q, type="Text")
                text._id = ga4gh_identify(text)
                translations = [Text(**text.as_dict())]
            else:
                translations = []

        return ToVRSService(
            search_term=q,
            variations=translations,
            service_meta_=ServiceMeta(
                version=__version__,
                response_datetime=datetime.now()
            ),
            warnings=warnings
        )
