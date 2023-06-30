"""Module for hgvs_dup_del_mode in normalize endpoint."""
from typing import Optional, Dict, Tuple, List

from ga4gh.vrs import models, normalize
from ga4gh.core import ga4gh_identify
from cool_seq_tool.data_sources import SeqRepoAccess

from variation.schemas.hgvs_to_copy_number_schema import CopyChange
from variation.schemas.normalize_response_schema import (
    HGVSDupDelMode as HGVSDupDelModeEnum
)
from variation.schemas.token_response_schema import AltType, AMBIGUOUS_REGIONS


class HGVSDupDelMode:
    """Class for handling how to interpret HGVS duplications and deletions."""

    def __init__(self, seqrepo_access: SeqRepoAccess) -> None:
        """Initialize HGVS Dup Del Mode.

        :param SeqRepoAccess seqrepo_access: Access to seqrepo
        """
        self.seqrepo_access = seqrepo_access

    def default_mode(
        self, alt_type: AltType, pos: Tuple[int, int], del_or_dup: str,
        location: Dict, vrs_seq_loc_ac: str, baseline_copies: Optional[int] = None,
        copy_change: Optional[CopyChange] = None
    ) -> Optional[Dict]:
        """Use default characteristics to return a variation.
        If baseline_copies not provided and endpoints are ambiguous - copy_number_change
            if copy_change not provided:
                copy_change - `efo:0030067` (loss) if del, `efo:0030070` (gain) if dup
        elif baseline_copies provided: copy_number_count
            copies are baseline + 1 for dup, baseline - 1 for del
        elif len del or dup > 100bp (use outermost coordinates)
            repeated_seq_expr with a derived_seq_expr subject (Allele)
        else:
            literal_seq_expr (normalized LiteralSequenceExpression Allele)

        :param alt_type: Alteration type
        :param pos: start_pos, end_pos
        :param del_or_dup: Must be either `del` or `dup`
        :param location: Sequence Location object
        :param vrs_seq_loc_ac: Accession used in VRS Sequence Location
        :param baseline_copies: Baseline copies for Copy Number Count variation
        :param Optional[CopyChange] copy_change: copy change
            for Copy Number Change Variation
        :return: VRS Variation object represented as a dict
        """
        variation = None
        if not baseline_copies and alt_type in AMBIGUOUS_REGIONS:
            variation = self.copy_number_change_mode(del_or_dup, location, copy_change)
        elif baseline_copies:
            variation = self.copy_number_count_mode(del_or_dup, location,
                                                    baseline_copies)
        elif pos and (pos[1] - pos[0] > 100):
            variation = self.repeated_seq_expr_mode(alt_type, location)
        else:
            variation = self.literal_seq_expr_mode(location, alt_type, vrs_seq_loc_ac)
        return variation

    def copy_number_count_mode(self, del_or_dup: str, location: Dict,
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
            "type": "CopyNumberCount",
            "subject": location,
            "copies": copies.as_dict()
        }
        variation["subject"]["_id"] = ga4gh_identify(
            models.SequenceLocation(**location)
        )
        variation["_id"] = ga4gh_identify(models.CopyNumberCount(**variation))
        return variation

    def copy_number_change_mode(
        self, del_or_dup: str, location: Dict,
        copy_change: Optional[CopyChange] = None
    ) -> Optional[Dict]:
        """Return copy number change variation

        :param str del_or_dup: Must be either `del` or `dup`
        :param Dict location: VRS SequenceLocation
        :param Optional[CopyChange] copy_change: The copy change
        :return: Copy Number Change variation as a dict
        """
        if not copy_change:
            if del_or_dup == "del":
                copy_change = CopyChange.LOSS.value
            else:
                copy_change = CopyChange.GAIN.value
        variation = {
            "type": "CopyNumberChange",
            "subject": location,
            "copy_change": copy_change
        }
        variation["subject"]["_id"] = ga4gh_identify(
            models.SequenceLocation(**location)
        )
        variation["_id"] = ga4gh_identify(models.CopyNumberChange(**variation))
        return variation

    def repeated_seq_expr_mode(self, alt_type: AltType,
                               location: Dict) -> Optional[Dict]:
        """Return a VRS Allele with a RepeatedSequenceExpression.
        The RepeatedSequenceExpression subject will be a
            DerivedSequenceExpression.

        :param AltType alt_type: Alteration type
        :param Dict location: VRS SequenceLocation
        :return: VRS Allele object represented as a dict
        """
        if "range" in alt_type.value:
            # Ranges should return an error
            return None

        if alt_type == AltType.DUPLICATION:
            count = models.Number(value=2, type="Number")
        elif alt_type == AltType.DELETION:
            count = models.Number(value=0, type="Number")
        else:
            return None

        location["_id"] = ga4gh_identify(models.SequenceLocation(**location))

        seq_expr = models.RepeatedSequenceExpression(
            seq_expr=models.DerivedSequenceExpression(
                location=location,
                reverse_complement=False,
                type="DerivedSequenceExpression"
            ),
            count=count,
            type="RepeatedSequenceExpression"
        )

        allele = models.Allele(
            location=location,
            state=seq_expr,
            type="Allele"
        )

        try:
            allele = normalize(allele, self.seqrepo_access)
        except (KeyError, AttributeError):
            return None
        else:
            allele.state.seq_expr.location = allele.location
            allele.location._id = ga4gh_identify(allele.location)
            allele._id = ga4gh_identify(allele)
            return allele.as_dict()

    def literal_seq_expr_mode(
        self, location: Dict, alt_type: AltType, vrs_seq_loc_ac: str
    ) -> Optional[Dict]:
        """Return a VRS Allele with a normalized LiteralSequenceExpression.

        :param Dict location: VRS Location
        :param AltType alt_type: Alteration type
        :param vrs_seq_loc_ac: Accession used in VRS Sequence Location
        :return: VRS Allele object represented as a dict
        """
        if alt_type in AMBIGUOUS_REGIONS:
            return None

        # This will only be deletion or dup
        if alt_type == AltType.DELETION:
            state = ""
        else:
            # start is start - 1, end is end
            ival = location["interval"]
            ref, _ = self.seqrepo_access.get_reference_sequence(
                vrs_seq_loc_ac, ival["start"]["value"] + 1, ival["end"]["value"] + 1
            )

            if ref:
                state = ref + ref
            else:
                return None

        allele = models.Allele(**{
            "type": "Allele",
            "location": location,
            "state": models.LiteralSequenceExpression(
                sequence=state, type="LiteralSequenceExpression"
            ).as_dict()
        })
        try:
            allele = normalize(allele, self.seqrepo_access)
        except (KeyError, AttributeError):
            return None
        else:
            allele.location._id = ga4gh_identify(allele.location)
            allele._id = ga4gh_identify(allele)
            return allele.as_dict()

    def interpret_variation(
        self, alt_type: AltType, location: Dict, errors: List,
        hgvs_dup_del_mode: HGVSDupDelModeEnum, vrs_seq_loc_ac: str,
        pos: Optional[Tuple[int, int]] = None,
        baseline_copies: Optional[int] = None,
        copy_change: Optional[CopyChange] = None
    ) -> Dict:
        """Interpret variation using HGVSDupDelMode

        :param alt_type: Alteration type
        :param location: VRS Location object
        :param errors: List of errors
        :param hgvs_dup_del_mode: Mode to use for interpreting HGVS duplications and
            deletions
        :param vrs_seq_loc_ac: Accession used in VRS Sequence Location
        :param pos: Position changes
        :param baseline_copies: Baseline copies number
        :param copy_change: The copy change
        :return: VRS Variation object
        """
        if "deletion" in alt_type.value:
            del_or_dup = "del"
        else:
            del_or_dup = "dup"

        variation = None
        if hgvs_dup_del_mode == HGVSDupDelModeEnum.DEFAULT:
            variation = self.default_mode(
                alt_type, pos, del_or_dup, location, vrs_seq_loc_ac,
                baseline_copies=baseline_copies, copy_change=copy_change
            )
        elif hgvs_dup_del_mode == HGVSDupDelModeEnum.REPEATED_SEQ_EXPR:
            variation = self.repeated_seq_expr_mode(alt_type, location)
        elif hgvs_dup_del_mode == HGVSDupDelModeEnum.LITERAL_SEQ_EXPR:
            variation = self.literal_seq_expr_mode(location, alt_type, vrs_seq_loc_ac)
        elif hgvs_dup_del_mode == HGVSDupDelModeEnum.COPY_NUMBER_COUNT:
            if baseline_copies:
                variation = self.copy_number_count_mode(
                    del_or_dup, location, baseline_copies
                )
            else:
                errors.append(
                    "`baseline_copies` must be provided for Copy Number Count Variation"
                )
        elif hgvs_dup_del_mode == HGVSDupDelModeEnum.COPY_NUMBER_CHANGE:
            variation = self.copy_number_change_mode(
                del_or_dup, location, copy_change=copy_change
            )

        if not variation:
            errors.append("Unable to get VRS Variation")

        return variation
