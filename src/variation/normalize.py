"""Module for Variation Normalization."""

import datetime
from urllib.parse import unquote

from cool_seq_tool.handlers import SeqRepoAccess
from cool_seq_tool.sources import UtaDatabase
from ga4gh.vrs import models

from variation import __version__
from variation.classify import Classify
from variation.schemas.app_schemas import Endpoint
from variation.schemas.classification_response_schema import ClassificationType
from variation.schemas.normalize_response_schema import (
    HGVSDupDelModeOption,
    NormalizeService,
    ServiceMeta,
)
from variation.schemas.token_response_schema import GnomadVcfToken, Token
from variation.schemas.translation_response_schema import (
    AC_PRIORITY_LABELS,
    TranslationResult,
    VrsSeqLocAcStatus,
)
from variation.schemas.validation_response_schema import ValidationSummary
from variation.to_vrs import ToVRS
from variation.tokenize import Tokenize
from variation.translate import Translate
from variation.utils import get_vrs_loc_seq, update_warnings_for_no_resp
from variation.validate import Validate


class Normalize(ToVRS):
    """The Normalize class used to normalize a given variation."""

    def __init__(
        self,
        seqrepo_access: SeqRepoAccess,
        tokenizer: Tokenize,
        classifier: Classify,
        validator: Validate,
        translator: Translate,
        uta: UtaDatabase,
    ) -> None:
        """Initialize Normalize class.

        :param seqrepo_access: Access to SeqRepo
        :param tokenizer: Tokenizer class for tokenizing
        :param classifier: Classifier class for classifying tokens
        :param validator: Validator class for validating valid inputs
        :param translator: Translating valid inputs
        :param UtaDatabase uta: Access to db containing alignment data
        """
        super().__init__(
            seqrepo_access,
            tokenizer,
            classifier,
            validator,
            translator,
        )
        self.uta = uta

    @staticmethod
    def _get_priority_translation_result(
        translations: list[TranslationResult], ac_status: VrsSeqLocAcStatus
    ) -> TranslationResult | None:
        """Get prioritized translation result. Tries to find translation results with
        the same `vrs_seq_loc_ac_status` as `ac_status`. If more than one translation
        result is found, will try to find translation result where `og_ac` (original
        accession used to get the translation) is the same as `vrs_seq_loc_ac`. If not
        found, will sort the translations and return the first translation result in
        the list

        :param translations: List of translation results
        :param ac_status: Accession status to filter by
        :return: Prioritized translation result with `ac_status` if found. Else, `None`
        """
        preferred_translations = [
            t for t in translations if t.vrs_seq_loc_ac_status == ac_status
        ]
        len_preferred_translations = len(preferred_translations)

        # Need to handle cases where there are multiple translations.
        # Different `og_ac`'s can lead to different translation results.
        # We must be consistent in what we return in /normalize
        if len_preferred_translations > 1:
            preferred_translations.sort(
                key=lambda t: (t.og_ac.split(".")[0], int(t.og_ac.split(".")[1])),
                reverse=True,
            )
            og_ac_preferred_match = (
                [t for t in preferred_translations if t.og_ac == t.vrs_seq_loc_ac]
                or [None]
            )[0]

            # We'll first see if `og_ac` (starting ac) matches the `ac_status`
            # accession. If that doesn't match, we'll just sort the original
            # acs and return the first element. Later on, we'll want to figure
            # out a better way to do this.
            if og_ac_preferred_match:
                translation_result = og_ac_preferred_match
            else:
                translation_result = preferred_translations[0]
        elif len_preferred_translations == 1:
            translation_result = preferred_translations[0]
        else:
            translation_result = None

        return translation_result

    @staticmethod
    def get_hgvs_dup_del_mode(
        tokens: list[Token],
        hgvs_dup_del_mode: HGVSDupDelModeOption | None = None,
        baseline_copies: int | None = None,
    ) -> tuple[HGVSDupDelModeOption | None, str | None]:
        """Get option to use for hgvs dup del mode

        :param tokens: List of tokens found in an input query
        :param hgvs_dup_del_mode: The hgvs dup del mode option provided in the input
            query. Mode to use for interpreting HGVS duplications and deletions.
            gnomad vcf token will always set to `HGVSDupDelModeOption.ALLELE`.
        :param baseline_copies: The baseline copies provided in the input query.
            Required when `hgvs_dup_del_mode == HGVSDupDelModeOption.COPY_NUMBER_COUNT`.
        :return: Tuple containing the hgvs dup del mode option and warnings
        """
        warning = None
        if len(tokens) == 1 and isinstance(tokens[0], GnomadVcfToken):
            hgvs_dup_del_mode = HGVSDupDelModeOption.ALLELE
        else:
            if not hgvs_dup_del_mode:
                hgvs_dup_del_mode = HGVSDupDelModeOption.DEFAULT

            if (
                hgvs_dup_del_mode == HGVSDupDelModeOption.COPY_NUMBER_COUNT
                and not baseline_copies
            ):
                warning = f"{hgvs_dup_del_mode.value} mode requires `baseline_copies`"
                return None, warning

        return hgvs_dup_del_mode, warning

    def _get_location_seq(
        self,
        validation_summary: ValidationSummary,
        variation: dict,
        priority_translation_result: TranslationResult,
    ) -> str | None:
        """Get reference sequence for a Sequence Location

        Does not support:
        - Ambiguous genomic deletions or duplications
        - Amplifications
        - Variations that are not Allele or Copy Number

        :param validation_summary: Validation summary for classification containing
            valid and invalid results
        :param variation: VRS Variation object
        :param priority_translation_result: Prioritized translation result
        :return: Reference sequence for a sequence location if found
        """
        valid_result = validation_summary.valid_results[0]
        classification_type = valid_result.classification.classification_type
        if classification_type not in {
            ClassificationType.GENOMIC_DELETION_AMBIGUOUS,
            ClassificationType.GENOMIC_DUPLICATION_AMBIGUOUS,
            ClassificationType.AMPLIFICATION,
        } and variation["type"] in {"Allele", "CopyNumberChange", "CopyNumberCount"}:
            return get_vrs_loc_seq(
                self.seqrepo_access,
                priority_translation_result.vrs_seq_loc_ac,
                variation["location"]["start"],
                variation["location"]["end"],
            )
        return None

    async def normalize(
        self,
        q: str,
        hgvs_dup_del_mode: HGVSDupDelModeOption | None = HGVSDupDelModeOption.DEFAULT,
        baseline_copies: int | None = None,
        copy_change: models.CopyChange | None = None,
    ) -> NormalizeService:
        """Normalize a given variation.

        :param q: HGVS, gnomAD VCF or Free Text description on GRCh37 or GRCh38 assembly
        :param hgvs_dup_del_mode: This parameter determines how to interpret HGVS
            dup/del expressions in VRS.
        :param baseline_copies: Baseline copies for HGVS duplications and deletions
        :param copy_change: The copy change for HGVS duplications and deletions
            represented as Copy Number Change Variation.
        :return: NormalizeService with variation and warnings
        """
        label = q.strip()
        variation = None
        warnings = []
        params = {
            "variation_query": q,
            "variation": variation,
            "warnings": warnings,
            "service_meta_": ServiceMeta(
                version=__version__,
                response_datetime=datetime.datetime.now(tz=datetime.UTC),
            ),
        }

        # Get tokens for input query
        tokens = self.tokenizer.perform(unquote(q.strip()), warnings)
        if warnings:
            update_warnings_for_no_resp(label, warnings)
            params["warnings"] = warnings
            return NormalizeService(**params)

        # Get HGVS dup del mode option to use
        hgvs_dup_del_mode, warning = self.get_hgvs_dup_del_mode(
            tokens, hgvs_dup_del_mode=hgvs_dup_del_mode, baseline_copies=baseline_copies
        )
        if warning:
            warnings.append(warning)
            update_warnings_for_no_resp(label, warnings)
            params["warnings"] = warnings
            return NormalizeService(**params)

        # Get classification for list of tokens
        classification = self.classifier.perform(tokens)
        if not classification:
            warnings.append(f"Unable to find classification for: {q}")
            params["warnings"] = warnings
            return NormalizeService(**params)

        # Get validation summary for classification
        validation_summary = await self.validator.perform(classification)
        if not validation_summary:
            update_warnings_for_no_resp(label, validation_summary.warnings)
            params["warnings"] = warnings
            return NormalizeService(**params)

        variation = None
        if validation_summary.valid_results:
            # Get translated VRS representations for valid results
            translations, warnings = await self.get_translations(
                validation_summary.valid_results,
                warnings,
                endpoint_name=Endpoint.NORMALIZE,
                hgvs_dup_del_mode=hgvs_dup_del_mode,
                baseline_copies=baseline_copies,
                copy_change=copy_change,
                do_liftover=True,
            )
            if translations:
                # Get prioritized translation result so that output is always the same
                for ac_status in AC_PRIORITY_LABELS:
                    translation_result = self._get_priority_translation_result(
                        translations, ac_status
                    )
                    if translation_result:
                        if (
                            translation_result.vrs_seq_loc_ac_status
                            == VrsSeqLocAcStatus.NA
                        ):
                            classification_type = translation_result.validation_result.classification.classification_type.value
                            if classification_type.startswith(("protein", "cdna")):
                                # Only supports protein/cDNA at the moment
                                warnings.append("Unable to find MANE representation")
                        break

                try:
                    variation = translation_result.vrs_variation
                except AttributeError as e:
                    warnings.append(str(e))
                else:
                    variation["location"]["sequence"] = self._get_location_seq(
                        validation_summary, variation, translation_result
                    )

                if not variation:
                    update_warnings_for_no_resp(label, warnings)
            else:
                update_warnings_for_no_resp(label, warnings)
        else:
            # No valid results were found for input query
            update_warnings_for_no_resp(label, warnings)

        params["variation"] = variation
        params["warnings"] = warnings
        return NormalizeService(**params)
