"""Module for toVRS endpoint."""
from typing import Tuple, Optional, List, Union
from ga4gh.vrsatile.pydantic.vrs_models import Allele, Haplotype, CopyNumber,\
    VariationSet, Text
from variation.base import Base
from variation.schemas.token_response_schema import Nomenclature
from variation.schemas.validation_response_schema import ValidationSummary
from urllib.parse import unquote
from variation.schemas.normalize_response_schema\
    import HGVSDupDelMode as HGVSDupDelModeEnum


class ToVRS(Base):
    """The class for translating variation strings to VRS representations."""

    def get_validations(
            self, q: str, normalize_endpoint: bool = False,
            hgvs_dup_del_mode: Optional[HGVSDupDelModeEnum] = None
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

        # gnomad vcf should always be a literal seq expression (allele)
        nomenclature = {t.nomenclature for t in tokens}
        if Nomenclature.GNOMAD_VCF in nomenclature:
            hgvs_dup_del_mode = HGVSDupDelModeEnum.LITERAL_SEQ_EXPR
        else:
            if normalize_endpoint:
                if hgvs_dup_del_mode:
                    hgvs_dup_del_mode = hgvs_dup_del_mode.strip().lower()
                    if not self.hgvs_dup_del_mode.is_valid_mode(hgvs_dup_del_mode):  # noqa: E501
                        warnings.append(
                            f"hgvs_dup_del_mode must be one of: "
                            f"{self.hgvs_dup_del_mode.valid_modes}")
                        return None, warnings
                else:
                    hgvs_dup_del_mode = HGVSDupDelModeEnum.DEFAULT
            else:
                hgvs_dup_del_mode = HGVSDupDelModeEnum.DEFAULT

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
