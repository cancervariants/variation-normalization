"""Module for hgvs_dup_del_mode in normalize endpoint."""
from typing import Dict, List, Optional, Union

from cool_seq_tool.handlers import SeqRepoAccess
from ga4gh.core import ga4gh_identify
from ga4gh.vrs import models, normalize

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
        location: Dict,
        vrs_seq_loc_ac: str,
        baseline_copies: Optional[int] = None,
        copy_change: Optional[models.CopyChange] = None,
        alt: Optional[str] = None,
    ) -> Optional[Dict]:
        """Use default characteristics to return a variation.
        If baseline_copies not provided and endpoints are ambiguous - copy_number_change
            if copy_change not provided:
                copy_change - `efo:0030067` (loss) if del, `efo:0030070` (gain) if dup
        elif baseline_copies provided: copy_number_count
            copies are baseline + 1 for dup, baseline - 1 for del
        else
            allele

        :param alt_type: The type of alteration
        :param location: Sequence Location object
        :param vrs_seq_loc_ac: Accession used in VRS Sequence Location
        :param baseline_copies: Baseline copies for Copy Number Count variation
        :param copy_change: copy change for Copy Number Change Variation
        :param alt: Alteration
        :return: VRS Variation object represented as a dict
        """
        variation = None
        if not baseline_copies and alt_type in AMBIGUOUS_REGIONS:
            variation = self.copy_number_change_mode(alt_type, location, copy_change)
        elif baseline_copies:
            variation = self.copy_number_count_mode(alt_type, location, baseline_copies)
        else:
            variation = self.allele_mode(location, alt_type, vrs_seq_loc_ac, alt)
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
    ) -> Dict:
        """Return a VRS Copy Number Variation.

        :param alt_type: The type of alteration
        :param location: VRS SequenceLocation
        :param baseline_copies: Baseline copies number
        :return: VRS Copy Number object represented as a dict
        """
        copies = baseline_copies - 1 if alt_type in DELS else baseline_copies + 1
        seq_loc = models.SequenceLocation(**location)
        seq_loc.id = ga4gh_identify(seq_loc)
        cn = models.CopyNumberCount(copies=copies, location=seq_loc)
        cn.id = ga4gh_identify(cn)
        return cn.model_dump(exclude_none=True)

    def copy_number_change_mode(
        self,
        alt_type: Union[
            AltType.DELETION,
            AltType.DELETION_AMBIGUOUS,
            AltType.DUPLICATION,
            AltType.DUPLICATION_AMBIGUOUS,
        ],
        location: Dict,
        copy_change: Optional[models.CopyChange] = None,
    ) -> Dict:
        """Return copy number change variation

        :param alt_type: The type of alteration
        :param location: VRS SequenceLocation
        :param copy_change: The copy change
        :return: Copy Number Change variation as a dict
        """
        if not copy_change:
            copy_change = (
                models.CopyChange.EFO_0030067
                if alt_type in DELS
                else models.CopyChange.EFO_0030070
            )

        seq_loc = models.SequenceLocation(**location)
        seq_loc.id = ga4gh_identify(seq_loc)
        cx = models.CopyNumberChange(location=seq_loc, copyChange=copy_change)
        cx.id = ga4gh_identify(cx)
        return cx.model_dump(exclude_none=True)

    def allele_mode(
        self,
        location: Dict,
        alt_type: AltType,
        vrs_seq_loc_ac: str,
        alt: str,
    ) -> Optional[Dict]:
        """Return a VRS Allele with a normalized LiteralSequenceExpression or
        ReferenceLengthExpression.

        :param location: VRS Location
        :param alt_type: Alteration type
        :param vrs_seq_loc_ac: Accession used in VRS Sequence Location
        :param alt: Alteration
        :return: VRS Allele object represented as a dict
        """
        if alt_type in AMBIGUOUS_REGIONS:
            return None

        if alt_type == AltType.DUPLICATION:
            # start is start - 1, end is end
            ref, _ = self.seqrepo_access.get_reference_sequence(
                vrs_seq_loc_ac,
                location["start"] + 1,
                location["end"] + 1,
            )

            if ref:
                state = ref + ref
            else:
                return None
        else:
            state = alt or ""

        allele = models.Allele(
            location=models.SequenceLocation(**location),
            state=models.LiteralSequenceExpression(sequence=state),
        )

        try:
            allele = normalize(allele, self.seqrepo_access)
        except (KeyError, AttributeError):
            return None
        else:
            allele.location.id = ga4gh_identify(allele.location)
            allele.id = ga4gh_identify(allele)
            return allele.model_dump(exclude_none=True)

    def interpret_variation(
        self,
        alt_type: AltType,
        location: Dict,
        errors: List,
        hgvs_dup_del_mode: HGVSDupDelModeOption,
        vrs_seq_loc_ac: str,
        baseline_copies: Optional[int] = None,
        copy_change: Optional[models.CopyChange] = None,
        alt: Optional[str] = None,
    ) -> Dict:
        """Interpret variation using HGVSDupDelMode

        :param alt_type: Alteration type
        :param location: VRS Location object
        :param errors: List of errors
        :param hgvs_dup_del_mode: Mode to use for interpreting HGVS duplications and
            deletions
        :param vrs_seq_loc_ac: Accession used in VRS Sequence Location
        :param baseline_copies: Baseline copies number
        :param copy_change: The copy change
        :param alt: The alteration
        :return: VRS Variation object
        """
        variation = None
        if hgvs_dup_del_mode == HGVSDupDelModeOption.DEFAULT:
            variation = self.default_mode(
                alt_type,
                location,
                vrs_seq_loc_ac,
                baseline_copies=baseline_copies,
                copy_change=copy_change,
                alt=alt,
            )
        elif hgvs_dup_del_mode == HGVSDupDelModeOption.ALLELE:
            variation = self.allele_mode(location, alt_type, vrs_seq_loc_ac, alt)
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
