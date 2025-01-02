"""Module for to_vrs endpoint."""

import datetime
from urllib.parse import unquote

from cool_seq_tool.handlers import SeqRepoAccess
from ga4gh.vrs import models

from variation import __version__
from variation.classify import Classify
from variation.schemas.app_schemas import Endpoint
from variation.schemas.normalize_response_schema import (
    HGVSDupDelModeOption,
    ServiceMeta,
)
from variation.schemas.to_vrs_response_schema import ToVRSService
from variation.schemas.translation_response_schema import TranslationResult
from variation.schemas.validation_response_schema import ValidationResult
from variation.tokenize import Tokenize
from variation.translate import Translate
from variation.utils import get_vrs_loc_seq
from variation.validate import Validate
from variation.vrs_representation import VRSRepresentation


class ToVRS(VRSRepresentation):
    """The class for translating variation strings to VRS representations."""

    def __init__(
        self,
        seqrepo_access: SeqRepoAccess,
        tokenizer: Tokenize,
        classifier: Classify,
        validator: Validate,
        translator: Translate,
    ) -> None:
        """Initialize the ToVRS class.

        :param SeqRepoAccess seqrepo_access: Access to SeqRepo
        :param Tokenize tokenizer: Tokenizer class for tokenizing
        :param Classify classifier: Classifier class for classifying tokens
        :param Validate validator: Validator class for validating valid inputs
        :param Translate translator: Translating valid inputs
        """
        super().__init__(seqrepo_access)
        self.tokenizer = tokenizer
        self.classifier = classifier
        self.validator = validator
        self.translator = translator

    async def get_translations(
        self,
        valid_results: list[ValidationResult],
        warnings: list,
        endpoint_name: Endpoint | None = None,
        hgvs_dup_del_mode: HGVSDupDelModeOption = HGVSDupDelModeOption.DEFAULT,
        baseline_copies: int | None = None,
        copy_change: models.CopyChange | None = None,
        do_liftover: bool = False,
    ) -> tuple[list[TranslationResult], list[str]]:
        """Get translation results

        :param valid_results: List of valid results for a given input
        :param warnings: List of warnings
        :param endpoint_name: Name of endpoint that is being used
        :param hgvs_dup_del_mode: Mode to use for interpreting HGVS duplications and
            deletions
        :param baseline_copies: The baseline copies for a copy number count variation
        :param copy_change: The copy change for a copy number change variation
        :param do_liftover: Whether or not to liftover to GRC3h8 assembly
        :return: Tuple containing list of translations and list of warnings
        """
        translations = []
        for valid_result in valid_results:
            tr = await self.translator.perform(
                valid_result,
                warnings,
                endpoint_name=endpoint_name,
                hgvs_dup_del_mode=hgvs_dup_del_mode,
                baseline_copies=baseline_copies,
                copy_change=copy_change,
                do_liftover=do_liftover,
            )
            if tr and tr not in translations:
                translations.append(tr)

        if not translations and not warnings:
            warnings.append("Unable to translate variation")

        return translations, warnings

    def _get_vrs_variations(self, translations: list[TranslationResult]) -> list[dict]:
        """Get translated VRS Variations.

        This method will also add ``sequence`` to the variation's location

        :param translations: List of translation results
        :return: List of unique VRS Variations
        """
        variations = []
        _added_variation_ids = set()

        # Ensure only unique VRS variations are in the list of variations returned
        for tr in translations:
            if tr.vrs_variation["id"] not in _added_variation_ids:
                vrs_variation = tr.vrs_variation
                vrs_variation["location"]["sequence"] = get_vrs_loc_seq(
                    self.seqrepo_access,
                    tr.vrs_seq_loc_ac,
                    vrs_variation["location"]["start"],
                    vrs_variation["location"]["end"],
                )
                variations.append(vrs_variation)
                _added_variation_ids.add(vrs_variation["id"])
        return variations

    async def to_vrs(self, q: str) -> ToVRSService:
        """Return a VRS-like representation of all validated variations for a query.

        :param str q: The variation to translate (HGVS, gnomAD VCF, or free text) on
            GRCh37 or GRCh38 assembly
        :return: ToVRSService containing VRS variations and warnings
        """
        warnings = []
        variations = []
        params = {
            "search_term": q,
            "variations": variations,
            "service_meta_": ServiceMeta(
                version=__version__,
                response_datetime=datetime.datetime.now(tz=datetime.UTC),
            ),
            "warnings": warnings,
        }

        # Get tokens for input query
        tokens = self.tokenizer.perform(unquote(q.strip()), warnings)
        if warnings:
            params["warnings"] = warnings
            return ToVRSService(**params)

        # Get classification for list of tokens
        classification = self.classifier.perform(tokens)
        if not classification:
            params["warnings"] = [f"Unable to find classification for: {q}"]
            return ToVRSService(**params)

        # Get validation summary for classification
        validation_summary = await self.validator.perform(classification)
        if validation_summary.valid_results:
            # Get translated VRS representation for valid results
            translations, warnings = await self.get_translations(
                validation_summary.valid_results,
                warnings,
                endpoint_name=Endpoint.TO_VRS,
                hgvs_dup_del_mode=HGVSDupDelModeOption.DEFAULT,
                do_liftover=False,
            )
        else:
            translations = []
            warnings = validation_summary.warnings

        params["warnings"] = warnings
        params["variations"] = self._get_vrs_variations(translations)
        return ToVRSService(**params)
