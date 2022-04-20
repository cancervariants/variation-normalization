"""Module for hgvs_dup_del_mode in normalize endpoint."""
import logging
from typing import Optional, Dict, Tuple, List
import copy
import json

from ga4gh.vrs import models
from ga4gh.core import ga4gh_identify, sha512t24u
from uta_tools.data_sources import SeqRepoAccess

from variation.schemas.hgvs_to_copy_number_schema import RelativeCopyClass
from variation.schemas.normalize_response_schema\
    import HGVSDupDelMode as HGVSDupDelModeEnum

logger = logging.getLogger("variation")
logger.setLevel(logging.DEBUG)


class HGVSDupDelMode:
    """Class for handling how to interpret HGVS duplications and deletions."""

    def __init__(self, seqrepo_access: SeqRepoAccess) -> None:
        """Initialize HGVS Dup Del Mode.

        :param SeqRepoAccess seqrepo_access: Access to seqrepo
        """
        self.seqrepo_access = seqrepo_access
        self.valid_dup_del_modes = {mode.value for mode in
                                    HGVSDupDelModeEnum.__members__.values()}
        self.valid_copy_number_modes = {mode.value for mode in
                                        HGVSDupDelModeEnum.__members__.values() if
                                        mode.value.endswith("_cnv")}

    def is_valid_dup_del_mode(self, mode: str) -> bool:
        """Determine if mode is a valid input for HGVS Dup Del Mode.

        :param str mode: Entered mode
        :return: `True` if valid mode. `False` otherwise.
        """
        hgvs_dup_del_mode = mode.strip().lower()
        return hgvs_dup_del_mode in self.valid_dup_del_modes

    def is_valid_copy_number_mode(self, mode: str) -> bool:
        """Determine if mode is a valid input for copy number mode

        :param str mode: Entered mode
        :return: `True` if valid mode. `False` otherwise.
        """
        copy_number_type_mode = mode.strip().lower()
        return copy_number_type_mode in self.valid_copy_number_modes

    def default_mode(
        self, alt_type: str, pos: Tuple[int, int], del_or_dup: str,
        location: Dict, allele: Dict = None, baseline_copies: Optional[int] = None,
        relative_copy_class: Optional[RelativeCopyClass] = None
    ) -> Optional[Dict]:
        """Use default characteristics to return a variation.
        If baseline_copies not provided and endpoints are ambiguous: relative_cnv
            if relative_copy_class not provided:
                relative_copy_class = `partial loss` if del, `low-level gain` if dup
        elif baseline_copies provided: absolute_cnv
            copies are baseline + 1 for dup, baseline - 1 for del
        elif len del or dup > 100bp (use outermost coordinates):
            repeated_seq_expr with a derived_seq_expr subject (Allele)
        else:
            literal_seq_expr (normalized LiteralSequenceExpression Allele)

        :param str alt_type: Alteration type
        :param tuple pos: start_pos, end_pos
        :param str del_or_dup: Must be either `del` or `dup`
        :param Dict location: Sequence Location object
        :param Dict allele: VRS Allele object represented as a dict
        :param Optional[int] baseline_copies: Baseline copies for Absolute Relative
            Copy Number Variation
        :param Optional[RelativeCopyClass] relative_copy_class: Relative copy class
            for Relative Copy Number Variation
        :return: VRS Variation object represented as a dict
        """
        variation = None
        if not baseline_copies and ("uncertain" in alt_type or "range" in alt_type):
            variation = self.relative_copy_number_mode(del_or_dup, location,
                                                       relative_copy_class)
        elif baseline_copies:
            variation = self.absolute_copy_number_mode(del_or_dup, location,
                                                       baseline_copies)
        elif pos and (pos[1] - pos[0] > 100):
            variation = self.repeated_seq_expr_mode(alt_type, location)
        else:
            variation = self.literal_seq_expr_mode(allele, alt_type)
        return variation

    def _ga4gh_identify_cnv(self, variation: Dict, is_abs: bool = True) -> Dict:
        """Add ga4gh digest to variation

        :param Dict variation: VRS Copy Number Variation
        :param bool is_abs: `True` if Absolute Copy Number.
            `False` if Relative Copy Number.
        :return: Variation with ga4gh digest identifiers
        """
        copy_variation = copy.deepcopy(variation)
        location_id = variation["subject"]["_id"].split(".")[-1]
        copy_variation["subject"] = location_id
        serialized = json.dumps(
            copy_variation, sort_keys=True, separators=(",", ":"), indent=None
        ).encode("utf-8")
        digest = sha512t24u(serialized)
        if is_abs:
            variation["_id"] = f"ga4gh:VAC.{digest}"
        else:
            variation["_id"] = f"ga4gh:VRC.{digest}"
        return variation

    def absolute_copy_number_mode(self, del_or_dup: str, location: Dict,
                                  baseline_copies: int) -> Optional[Dict]:
        """Return a VRS Copy Number Variation.

        :param str del_or_dup: Must be either `del` or `dup`
        :param Dict location: VRS SequenceLocation
        :param int baseline_copies: Baseline copies number
        :return: VRS Copy Number object represented as a dict
        """
        copies = models.Number(
            value=baseline_copies - 1 if del_or_dup == "del" else baseline_copies + 1,
            type="Number"
        )
        variation = {
            "type": "AbsoluteCopyNumber",
            "subject": location,
            "copies": copies.as_dict()
        }
        return self._ga4gh_identify_cnv(variation, is_abs=True)

    def relative_copy_number_mode(
        self, del_or_dup: str, location: Dict,
        relative_copy_class: Optional[RelativeCopyClass] = None
    ) -> Optional[Dict]:
        """Return relative copy number variation

        :param str del_or_dup: Must be either `del` or `dup`
        :param Dict location: VRS SequenceLocation
        :param Optional[RelativeCopyClass] relative_copy_class: The relative copy class
        :return: Relative copy number variation as a dict
        """
        if not relative_copy_class:
            if del_or_dup == "del":
                relative_copy_class = RelativeCopyClass.PARTIAL_LOSS.value
            else:
                relative_copy_class = RelativeCopyClass.LOW_LEVEL_GAIN.value
        variation = {
            "type": "RelativeCopyNumber",
            "subject": location,
            "relative_copy_class": relative_copy_class
        }
        return self._ga4gh_identify_cnv(variation, is_abs=False)

    def repeated_seq_expr_mode(self, alt_type: str,
                               location: Dict) -> Optional[Dict]:
        """Return a VRS Allele with a RepeatedSequenceExpression.
        The RepeatedSequenceExpression subject will be a
            DerivedSequenceExpression.

        :param str alt_type: Alteration type
        :param Dict location: VRS SequenceLocation
        :return: VRS Allele object represented as a dict
        """
        if "range" in alt_type:
            # Ranges should return an error
            return None

        if alt_type == "duplication":
            count = models.Number(value=2, type="Number")
        elif alt_type == "deletion":
            count = models.Number(value=0, type="Number")
        else:
            return None

        seq_expr = models.RepeatedSequenceExpression(
            seq_expr=models.DerivedSequenceExpression(
                location=location,
                reverse_complement=False,
                type="DerivedSequenceExpression"
            ),
            count=count,
            type="RepeatedSequenceExpression"
        )

        variation = models.Allele(
            location=location,
            state=seq_expr,
            type="Allele"
        )
        return self._ga4gh_identify_variation(variation)

    def literal_seq_expr_mode(self, allele: Dict,
                              alt_type: str) -> Optional[Dict]:
        """Return a VRS Allele with a normalized LiteralSequenceExpression.

        :param Dict allele: normalized VRS Allele object represented as a dict
        :param str alt_type: Alteration type
        :return: VRS Allele object represented as a dict
        """
        if "range" in alt_type or "uncertain" in alt_type:
            return None

        variation = models.Allele(**allele) if allele else None
        return self._ga4gh_identify_variation(variation)

    @staticmethod
    def _ga4gh_identify_variation(variation: models.Variation) -> Optional[Dict]:
        """Return variation with GA4GH digest-based id.

        :param models.Variation variation: VRS variation object
        :return: VRS Variation with GA4GH digest-based id represented as a dict
        """
        if variation is None:
            return None
        else:
            variation._id = ga4gh_identify(variation)
            return variation.as_dict()

    def interpret_variation(
        self, alt_type: str, allele: Dict, errors: List,
        hgvs_dup_del_mode: HGVSDupDelModeEnum, pos: Optional[Tuple[int, int]] = None,
        baseline_copies: Optional[int] = None,
        relative_copy_class: Optional[RelativeCopyClass] = None
    ) -> Dict:
        """Interpret variation using HGVSDupDelMode

        :param str alt_type: Alteration type
        :param Dict allele: VRS Allele object
        :param List errors: List of errors
        :param HGVSDupDelModeEnum hgvs_dup_del_mode: Mode to use for
            interpreting HGVS duplications and deletions
        :param Optional[Tuple[int, int]] pos: Position changes
        :param Optional[int] baseline_copies: Baseline copies number
        :param Optional[RelativeCopyClass] relative_copy_class: The relative copy class
        :return: VRS Variation object
        """
        if "deletion" in alt_type:
            del_or_dup = "del"
        else:
            del_or_dup = "dup"
        variation = None
        if allele is None:
            errors.append("Unable to get Allele")
        else:
            if hgvs_dup_del_mode == HGVSDupDelModeEnum.DEFAULT:
                variation = self.default_mode(
                    alt_type, pos, del_or_dup, allele["location"],
                    allele=allele, baseline_copies=baseline_copies,
                    relative_copy_class=relative_copy_class)
            elif hgvs_dup_del_mode == HGVSDupDelModeEnum.REPEATED_SEQ_EXPR:
                variation = self.repeated_seq_expr_mode(
                    alt_type, allele["location"]
                )
            elif hgvs_dup_del_mode == HGVSDupDelModeEnum.LITERAL_SEQ_EXPR:
                variation = self.literal_seq_expr_mode(allele, alt_type)
            elif hgvs_dup_del_mode == HGVSDupDelModeEnum.ABSOLUTE_CNV:
                if baseline_copies:
                    variation = self.absolute_copy_number_mode(
                        del_or_dup, allele["location"], baseline_copies)
                else:
                    errors.append("`baseline_copies` must be provided for Absolute"
                                  " Copy Number Variation")
            elif hgvs_dup_del_mode == HGVSDupDelModeEnum.RELATIVE_CNV:
                variation = self.relative_copy_number_mode(
                    del_or_dup, allele["location"],
                    relative_copy_class=relative_copy_class)
            if not variation:
                errors.append("Unable to get VRS Variation")
        return variation
