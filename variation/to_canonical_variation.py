"""Module containing To Canonical Variation work"""
from typing import Optional, List, Tuple, Dict
import copy
import json
from datetime import datetime

import python_jsonschema_objects
from ga4gh.vrsatile.pydantic.vrs_models import Text
from ga4gh.vrsatile.pydantic.vrsatile_models import CanonicalVariation
from ga4gh.vrs import models
from ga4gh.core import ga4gh_identify, sha512t24u
from ga4gh.vrs.dataproxy import SeqRepoDataProxy
from ga4gh.vrs.extras.translator import Translator
from uta_tools.schemas import Assembly
from uta_tools.data_sources import UTADatabase, SeqRepoAccess

from variation.classifiers.classify import Classify
from variation.hgvs_dup_del_mode import HGVSDupDelMode
from variation.schemas.app_schemas import Endpoint
from variation.schemas.normalize_response_schema import ServiceMeta, \
    ToCanonicalVariationFmt, ToCanonicalVariationService
from variation.schemas.normalize_response_schema\
    import HGVSDupDelMode as HGVSDupDelModeEnum
from variation.schemas.hgvs_to_copy_number_schema import VALID_RELATIVE_COPY_CLASS,\
    RelativeCopyClass
from variation.to_vrs import ToVRS
from variation.tokenizers.tokenize import Tokenize
from variation.translators.translate import Translate
from variation.utils import get_mane_valid_result, no_variation_entered
from variation.validators.validate import Validate
from variation.version import __version__


class ToCanonicalVariation(ToVRS):
    """Class for translating to canonical variation"""

    def __init__(self, seqrepo_access: SeqRepoAccess, dp: SeqRepoDataProxy,
                 tokenizer: Tokenize, classifier: Classify, validator: Validate,
                 translator: Translate, hgvs_dup_del_mode: HGVSDupDelMode,
                 tlr: Translator, uta: UTADatabase) -> None:
        """Initialize the to canonical variation class

        :param SeqRepoAccess seqrepo_access: Access to SeqRepo via UTA Tools
        :param SeqRepoDataProxy dp: Access to SeqRepo via VRS Python
        :param Tokenize tokenizer: Tokenizer class for tokenizing
        :param Classify classifier: Classifier class for classifying tokens
        :param Validate validator: Validator class for validating valid inputs
        :param Translate translator: Translating valid inputs
        :param HGVSDupDelMode hgvs_dup_del_mode: Class for handling
            HGVS dup/del expressions
        :param Translator tlr: Class for translating nomenclatures to and from VRS
        :param UTADatabase uta: Access to UTA queries
        """
        super().__init__(seqrepo_access, dp, tokenizer, classifier, validator,
                         translator, hgvs_dup_del_mode)
        self.tlr = tlr
        self.uta = uta

    async def to_canonical_variation(
        self, q: str,
        fmt: ToCanonicalVariationFmt,
        complement: bool = False,
        do_liftover: bool = False,
        hgvs_dup_del_mode: Optional[HGVSDupDelModeEnum] = None,
        relative_copy_class: Optional[RelativeCopyClass] = None,
        baseline_copies: Optional[int] = None,
        untranslatable_returns_text: bool = False
    ) -> ToCanonicalVariationService:
        """Given query as ToCanonicalVariationFmt, return canonical variation

        :param str q: Query to translate to canonical variation
        :param ToCanonicalVariationFmt fmt: The representation for `q`
        :param bool complement: This field indicates that a categorical variation is
            defined to include (false) or exclude (true) variation concepts matching the
            categorical variation. This is equivalent to a logical NOT operation on the
            categorical variation properties.
        :param bool do_liftover: Whether or not to liftover to GRCh38 assembly.
        :param Optional[HGVSDupDelModeEnum] hgvs_dup_del_mode: Determines how to
            interpret HGVS dup/del expressions in VRS. Must be one of: `default`,
            `absolute_cnv`, `relative_cnv`, `repeated_seq_expr`, `literal_seq_expr`
        :param RelativeCopyClass relative_copy_class: The relative copy class
        :param Optional[int] baseline_copies: Baseline copies number
        :param bool untranslatable_returns_text: `True` return VRS Text Object when
            unable to translate or normalize query. `False` return `None` when
            unable to translate or normalize query.
        :return: ToCanonicalVariationService containing Canonical Variation and
            list of warnings
        """
        q = q.strip()
        variation = None
        warnings = list()
        if q:
            if fmt == ToCanonicalVariationFmt.SPDI:
                variation, warnings = await self.spdi_to_canonical_variation(
                    q, warnings, do_liftover=do_liftover)
            elif fmt == ToCanonicalVariationFmt.HGVS:
                variation, warnings = await self.hgvs_to_canonical_variation(
                    q, warnings, do_liftover=do_liftover,
                    hgvs_dup_del_mode=hgvs_dup_del_mode,
                    relative_copy_class=relative_copy_class,
                    baseline_copies=baseline_copies)
            else:
                warnings = [f"fmt, {fmt}, is not supported. "
                            f"Must be either `spdi` or `hgvs`"]

            if variation and not warnings:
                variation_type = variation["type"]
                if variation_type == "Allele":
                    variation["location"]["_id"] = ga4gh_identify(
                        models.SequenceLocation(**variation["location"]))
                elif variation_type in ["RelativeCopyNumber", "AbsoluteCopyNumber"]:
                    variation["subject"]["_id"] = ga4gh_identify(
                        models.SequenceLocation(**variation["subject"]))
                else:
                    warnings = [f"Variation type, {variation_type}, not supported"]

            if not variation and untranslatable_returns_text:
                text = models.Text(definition=q, type="Text")
                text._id = ga4gh_identify(text)
                variation = Text(**text.as_dict()).dict(by_alias=True)

            if variation:
                canonical_variation = {
                    "type": "CanonicalVariation",
                    "complement": complement,
                    "variation": variation
                }

                cpy_canonical_variation = copy.deepcopy(canonical_variation)
                cpy_canonical_variation["variation"] = canonical_variation["variation"]["_id"].split(".")[-1]  # noqa: E501
                serialized = json.dumps(
                    cpy_canonical_variation, sort_keys=True, separators=(",", ":"),
                    indent=None
                ).encode("utf-8")
                digest = sha512t24u(serialized)
                # VCC = variation categorical canonical
                canonical_variation["_id"] = f"ga4gh:VCC.{digest}"
                canonical_variation = CanonicalVariation(**canonical_variation)
            else:
                canonical_variation = None
        else:
            canonical_variation, warnings = no_variation_entered()

        return ToCanonicalVariationService(
            query=q,
            canonical_variation=canonical_variation,
            warnings=warnings,
            service_meta_=ServiceMeta(
                version=__version__,
                response_datetime=datetime.now()
            )
        )

    async def spdi_to_canonical_variation(
        self, q: str, warnings: List, do_liftover: bool = False
    ) -> Tuple[Optional[Dict], List]:
        """Given SPDI representation, return canonical variation

        :param str q: SPDI query
        :param List warnings: List of warnings
        :param bool do_liftover: Whether or not to liftover to GRCh38 assembly
        :return: Canonical Variation and warnings
        """
        variation = None
        match = self.tlr.spdi_re.match(q)
        if not match:
            warnings.append(f"{q} is not a valid SPDI expression")
            return variation, warnings

        spdi_parts = q.split(":")
        ac = spdi_parts[0]
        start_pos = int(spdi_parts[1])  # inter-residue, 0 based
        deleted_seq = spdi_parts[2]
        inserted_seq = spdi_parts[3]
        end_pos = start_pos + len(deleted_seq)

        if do_liftover:
            newest_assembly_acs = await self.uta.get_newest_assembly_ac(ac)
            if not newest_assembly_acs:
                warnings.append(f"Unable to get newest assemblies for {ac}")
                return variation, None
            new_ac = newest_assembly_acs[0][0]
            if new_ac != ac:
                ac = new_ac
                chromosome, warning = self.seqrepo_access.ac_to_chromosome(ac)
                if not chromosome:
                    warnings.append(warning)
                    return variation, warnings
                chromosome = f"chr{chromosome}"

                pos = start_pos + len(deleted_seq)
                liftover_resp = self.uta.get_liftover(chromosome, pos, Assembly.GRCH38)
                if not liftover_resp:
                    warnings.append(f"Position {pos} does not exist on "
                                    f"chromosome {chromosome}")
                    return variation, warnings
                start_pos = liftover_resp[1] - len(deleted_seq)
                end_pos = start_pos + len(deleted_seq)
                q = f"{ac}:{start_pos}:{deleted_seq}:{inserted_seq}"

        try:
            variation = self.tlr.translate_from(
                q, fmt=ToCanonicalVariationFmt.SPDI.value)
        except (ValueError, python_jsonschema_objects.validators.ValidationError) as e:  # noqa: E501
            warnings.append(f"vrs-python translator raised error: {e}")
        except KeyError as e:
            warnings.append(f"vrs-python translator raised error: "
                            f"seqrepo could not translate identifier {e}")
        else:
            # Validate SPDI
            try:
                sequence = self.seqrepo_access.seqrepo_client.fetch(
                    ac, start_pos, end=end_pos)
            except ValueError as e:
                warnings.append(str(e))
            else:
                if not sequence:
                    warnings.append(f"Position, {start_pos}, does not exist on {ac}")
                else:
                    if deleted_seq != sequence:
                        warnings.append(f"Expected to find reference sequence"
                                        f" {deleted_seq} but found {sequence} on {ac}")

            if warnings:
                variation = None
            else:
                variation = variation.as_dict()
        return variation, warnings

    async def hgvs_to_canonical_variation(
        self, q: str, warnings: List, do_liftover: bool = False,
        hgvs_dup_del_mode: Optional[HGVSDupDelModeEnum] = None,
        relative_copy_class: Optional[RelativeCopyClass] = None,
        baseline_copies: Optional[int] = None
    ) -> Tuple[Optional[Dict], List]:
        """Given HGVS representation, return canonical variation

        :param str q: HGVS query
        :param List warnings: List of warnings
        :param bool do_liftover: Whether or not to liftover to GRCh38 assembly
        :param Optional[HGVSDupDelModeEnum] hgvs_dup_del_mode: Determines how to
            interpret HGVS dup/del expressions in VRS. Must be one of: `default`,
            `absolute_cnv`, `relative_cnv`, `repeated_seq_expr`, `literal_seq_expr`
        :param RelativeCopyClass relative_copy_class: The relative copy class
        :param Optional[int] baseline_copies: Baseline copies number
        :return: Canonical Variation and warnings
        """
        variation = None
        match = self.tlr.hgvs_re.match(q)
        if not match:
            warnings.append(f"{q} is not a valid HGVS expression")
            return variation, warnings

        tokens = self.tokenizer.perform(q, warnings)
        classifications = self.classifier.perform(tokens)
        hgvs_classifications = list()
        for c in classifications:
            if "HGVS" in c.matching_tokens or "ReferenceSequence" in c.matching_tokens:
                hgvs_classifications.append(c)
        if not hgvs_classifications:
            warnings = [f"{q} is not a supported HGVS expression"]
            return variation, warnings

        if hgvs_dup_del_mode == HGVSDupDelModeEnum.RELATIVE_CNV:
            if relative_copy_class:
                if relative_copy_class.lower() not in VALID_RELATIVE_COPY_CLASS:
                    return None, [f"{relative_copy_class} is not a valid relative "
                                  f"copy class: {VALID_RELATIVE_COPY_CLASS}"]
        elif hgvs_dup_del_mode == HGVSDupDelModeEnum.ABSOLUTE_CNV:
            if not baseline_copies:
                return None, [f"{hgvs_dup_del_mode} requires `baseline_copies`"]
        elif not hgvs_dup_del_mode:
            hgvs_dup_del_mode = HGVSDupDelModeEnum.DEFAULT

        validations = await self.validator.perform(
            classifications, endpoint_name=Endpoint.TO_CANONICAL, warnings=warnings,
            hgvs_dup_del_mode=hgvs_dup_del_mode, do_liftover=do_liftover,
            relative_copy_class=relative_copy_class, baseline_copies=baseline_copies)
        if not warnings:
            warnings = validations.warnings

        if do_liftover:
            if len(validations.valid_results) > 0:
                valid_result = get_mane_valid_result(q, validations, warnings)
                if valid_result:
                    variation = valid_result.variation
        else:
            translations, warnings = self.get_translations(
                validations, warnings)
            if translations:
                variation = translations[0]
        return variation, warnings
