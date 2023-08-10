"""Module for hgvs_dup_del_mode in normalize endpoint."""
from typing import Dict, List, Optional, Tuple, Union

from cool_seq_tool.data_sources import SeqRepoAccess
from ga4gh.core import ga4gh_identify
from ga4gh.vrs import models, normalize
from ga4gh.vrsatile.pydantic.vrs_models import CopyChange

from variation.schemas.normalize_response_schema import HGVSDupDelModeOption
from variation.schemas.token_response_schema import AMBIGUOUS_REGIONS, AltType

# Define deletion alt types
DELS = {AltType.DELETION_AMBIGUOUS, AltType.DELETION}


class HGVSDupDelMode:
    """Class for handling how to interpret HGVS duplications and deletions."""

    def __init__(self, seqrepo_access: SeqRepoAccess) -> None:
        """Initialize HGVS Dup Del Mode.

        :param SeqRepoAccess seqrepo_access: Access to seqrepo
        """
        self.seqrepo_access = seqrepo_access

    def default_mode(
        self,
        alt_type: Union[
            AltType.DELETION,
            AltType.DELETION_AMBIGUOUS,
            AltType.DUPLICATION,
            AltType.DUPLICATION_AMBIGUOUS,
        ],
        pos: Tuple[int, int],
        location: Dict,
        vrs_seq_loc_ac: str,
        baseline_copies: Optional[int] = None,
        copy_change: Optional[CopyChange] = None,
    ) -> Optional[Dict]:
        """Use default characteristics to return a variation.
        If baseline_copies not provided and endpoints are ambiguous - copy_number_change
            if copy_change not provided:
                copy_change - `efo:0030067` (loss) if del, `efo:0030070` (gain) if dup
        elif baseline_copies provided: copy_number_count
            copies are baseline + 1 for dup, baseline - 1 for del
        elif len dup > 100bp (use outermost coordinates)
            repeated_seq_expr with a derived_seq_expr subject (Allele)
        else:
            literal_seq_expr (normalized LiteralSequenceExpression Allele)

        :param alt_type: The type of alteration
        :param pos: start_pos, end_pos
        :param location: Sequence Location object
        :param vrs_seq_loc_ac: Accession used in VRS Sequence Location
        :param baseline_copies: Baseline copies for Copy Number Count variation
        :param copy_change: copy change for Copy Number Change Variation
        :return: VRS Variation object represented as a dict
        """
        variation = None
        if not baseline_copies and alt_type in AMBIGUOUS_REGIONS:
            variation = self.copy_number_change_mode(alt_type, location, copy_change)
        elif baseline_copies:
            variation = self.copy_number_count_mode(alt_type, location, baseline_copies)
        elif (alt_type not in DELS) and pos and (pos[1] - pos[0] > 100):
            variation = self.repeated_seq_expr_mode(alt_type, location)
        else:
            variation = self.literal_seq_expr_mode(location, alt_type, vrs_seq_loc_ac)
        return variation

    def copy_number_count_mode(
        self,
        alt_type: Union[
            AltType.DELETION,
            AltType.DELETION_AMBIGUOUS,
            AltType.DUPLICATION,
            AltType.DUPLICATION_AMBIGUOUS,
        ],
        location: Dict,
        baseline_copies: int,
    ) -> Optional[Dict]:
        """Return a VRS Copy Number Variation.

        :param alt_type: The type of alteration
        :param location: VRS SequenceLocation
        :param baseline_copies: Baseline copies number
        :return: VRS Copy Number object represented as a dict
        """
        copies = models.Number(
            value=baseline_copies - 1 if alt_type in DELS else baseline_copies + 1,
            type="Number",
        )
        variation = {
            "type": "CopyNumberCount",
            "subject": location,
            "copies": copies.as_dict(),
        }
        variation["subject"]["_id"] = ga4gh_identify(
            models.SequenceLocation(**location)
        )
        variation["_id"] = ga4gh_identify(models.CopyNumberCount(**variation))
        return variation

    def copy_number_change_mode(
        self,
        alt_type: Union[
            AltType.DELETION,
            AltType.DELETION_AMBIGUOUS,
            AltType.DUPLICATION,
            AltType.DUPLICATION_AMBIGUOUS,
        ],
        location: Dict,
        copy_change: Optional[CopyChange] = None,
    ) -> Optional[Dict]:
        """Return copy number change variation

        :param alt_type: The type of alteration
        :param Dict location: VRS SequenceLocation
        :param Optional[CopyChange] copy_change: The copy change
        :return: Copy Number Change variation as a dict
        """
        if not copy_change:
            if alt_type in DELS:
                copy_change = CopyChange.LOSS.value
            else:
                copy_change = CopyChange.GAIN.value
        variation = {
            "type": "CopyNumberChange",
            "subject": location,
            "copy_change": copy_change,
        }
        variation["subject"]["_id"] = ga4gh_identify(
            models.SequenceLocation(**location)
        )
        variation["_id"] = ga4gh_identify(models.CopyNumberChange(**variation))
        return variation

    def repeated_seq_expr_mode(
        self, alt_type: AltType, location: Dict
    ) -> Optional[Dict]:
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
                type="DerivedSequenceExpression",
            ),
            count=count,
            type="RepeatedSequenceExpression",
        )

        allele = models.Allele(location=location, state=seq_expr, type="Allele")

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
        self,
        location: Dict,
        alt_type: AltType,
        vrs_seq_loc_ac: str,
    ) -> Optional[Dict]:
        """Return a VRS Allele with a normalized LiteralSequenceExpression.

        :param Dict location: VRS Location
        :param AltType alt_type: Alteration type
        :param vrs_seq_loc_ac: Accession used in VRS Sequence Location
        :return: VRS Allele object represented as a dict
        """
        if alt_type in AMBIGUOUS_REGIONS:
            return None

        if alt_type == AltType.DUPLICATION:
            # start is start - 1, end is end
            ival = location["interval"]
            ref, _ = self.seqrepo_access.get_reference_sequence(
                vrs_seq_loc_ac, ival["start"]["value"] + 1, ival["end"]["value"] + 1
            )

            if ref:
                state = ref + ref
            else:
                return None
        else:
            state = ""

        allele = models.Allele(
            **{
                "type": "Allele",
                "location": location,
                "state": models.LiteralSequenceExpression(
                    sequence=state, type="LiteralSequenceExpression"
                ).as_dict(),
            }
        )

        try:
            allele = normalize(allele, self.seqrepo_access)
        except (KeyError, AttributeError):
            return None
        else:
            allele.location._id = ga4gh_identify(allele.location)
            allele._id = ga4gh_identify(allele)
            return allele.as_dict()

    def interpret_variation(
        self,
        alt_type: AltType,
        location: Dict,
        errors: List,
        hgvs_dup_del_mode: HGVSDupDelModeOption,
        vrs_seq_loc_ac: str,
        pos: Optional[Tuple[int, int]] = None,
        baseline_copies: Optional[int] = None,
        copy_change: Optional[CopyChange] = None,
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
        variation = None
        if hgvs_dup_del_mode == HGVSDupDelModeOption.DEFAULT:
            variation = self.default_mode(
                alt_type,
                pos,
                location,
                vrs_seq_loc_ac,
                baseline_copies=baseline_copies,
                copy_change=copy_change,
            )
        elif hgvs_dup_del_mode == HGVSDupDelModeOption.REPEATED_SEQ_EXPR:
            variation = self.repeated_seq_expr_mode(alt_type, location)
        elif hgvs_dup_del_mode == HGVSDupDelModeOption.LITERAL_SEQ_EXPR:
            variation = self.literal_seq_expr_mode(location, alt_type, vrs_seq_loc_ac)
        elif hgvs_dup_del_mode == HGVSDupDelModeOption.COPY_NUMBER_COUNT:
            if baseline_copies:
                variation = self.copy_number_count_mode(
                    alt_type, location, baseline_copies
                )
            else:
                errors.append(
                    "`baseline_copies` must be provided for Copy Number Count Variation"
                )
        elif hgvs_dup_del_mode == HGVSDupDelModeOption.COPY_NUMBER_CHANGE:
            variation = self.copy_number_change_mode(
                alt_type, location, copy_change=copy_change
            )

        if not variation:
            errors.append("Unable to get VRS Variation")

        return variation
