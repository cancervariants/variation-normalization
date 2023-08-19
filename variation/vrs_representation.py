"""Module for generating VRS objects"""
from typing import Dict, List, Optional, Tuple, Union

from bioutils.accessions import coerce_namespace
from cool_seq_tool.data_sources import SeqRepoAccess
from cool_seq_tool.schemas import AnnotationLayer
from ga4gh.core import ga4gh_identify
from ga4gh.vrs import models, normalize
from ga4gh.vrsatile.pydantic.vrs_models import Allele
from pydantic import ValidationError

from variation.schemas.token_response_schema import (
    AMBIGUOUS_REGIONS,
    AltType,
)


class VRSRepresentation:
    """Class for representing VRS objects"""

    def __init__(self, seqrepo_access: SeqRepoAccess) -> None:
        """Initialize the VRSRepresentation class

        :param SeqRepoAccess seqrepo_access: Access to SeqRepo
        """
        self.seqrepo_access = seqrepo_access

    @staticmethod
    def get_start_end(
        coordinate: str, start: int, end: int, cds_start: int, errors: List
    ) -> Optional[Tuple[int, int]]:
        """Get start and end coordinates.

        :param str coordinate: Coordinate used. Must be either `p`, `c`, or `g`
        :param int start: Start position change
        :param int end: End position change
        :param int cds_start: Coding start site
        :param List errors: List of errors
        :return: Tuple[start, end]
        """
        try:
            start = int(start)
            if end is None:
                end = start
            end = int(end)
        except (ValueError, TypeError):
            errors.append("Start/End must be valid ints")
            return None

        if coordinate == "c":
            if cds_start:
                start += cds_start
                end += cds_start
        return start, end

    @staticmethod
    def get_start_indef_range(start: int) -> models.IndefiniteRange:
        """Return indefinite range given start coordinate

        :param int start: Start position (assumes 1-based)
        :return: Indefinite range model
        """
        return models.IndefiniteRange(
            value=start - 1, comparator="<=", type="IndefiniteRange"
        )

    @staticmethod
    def get_end_indef_range(end: int) -> models.IndefiniteRange:
        """Return indefinite range given end coordinate

        :param int end: End position (assumes 1-based)
        :return: Indefinite range model
        """
        return models.IndefiniteRange(
            value=end, comparator=">=", type="IndefiniteRange"
        )

    @staticmethod
    def get_start_end_definite_range(
        start0: int, start1: int, end0: int, end1: int
    ) -> Tuple[models.DefiniteRange, models.DefiniteRange]:
        """Return sequence start and end as definite ranges

        :param start0: Start left pos (assumes 1-based)
        :paramnt start1: Start right pos (assumes 1-based)
        :param end0: End left pos (assumes 1-based)
        :param end1: End right pos (assumes 1-based)
        :return: Start and end as definite ranges
        """
        start = models.DefiniteRange(
            min=start0 - 1, max=start1 - 1, type="DefiniteRange"
        )
        end = models.DefiniteRange(min=end0 + 1, max=end1 + 1, type="DefiniteRange")
        return start, end

    @staticmethod
    def get_sequence_loc(
        ac: str,
        start: Union[models.Number, models.DefiniteRange, models.IndefiniteRange],
        end: Union[models.Number, models.DefiniteRange, models.IndefiniteRange],
    ) -> models.Location:
        """Return VRS location

        :param ac: Accession
        :param start: start pos
        :param end: end pos
        :return: VRS Location model
        """
        return models.SequenceLocation(
            sequence_id=coerce_namespace(ac),
            start=start,
            end=end,
            type="SequenceLocation",
        )

    def vrs_allele(
        self,
        ac: str,
        start: Union[models.Number, models.DefiniteRange, models.IndefiniteRange],
        end: Union[models.Number, models.DefiniteRange, models.IndefiniteRange],
        sstate: Union[
            models.LiteralSequenceExpression,
            models.DerivedSequenceExpression,
            models.RepeatedSequenceExpression,
        ],
        alt_type: AltType,
        errors: List[str],
    ) -> Optional[Dict]:
        """Create a VRS Allele object.

        :param ac: Accession
        :param start: start pos
        :param end: end pos
        :param sstate: State
        :param alt_type: Type of alteration
        :param errors: List of errors
        :return: VRS Allele object represented as a Dict
        """
        try:
            location = self.get_sequence_loc(ac, start, end)
        except ValueError as e:
            errors.append(f"Unable to get sequence location: {e}")
            return None
        allele = models.Allele(location=location, state=sstate, type="Allele")
        # Ambiguous regions do not get normalized
        if alt_type not in AMBIGUOUS_REGIONS:
            try:
                allele = normalize(allele, self.seqrepo_access)
            except (KeyError, AttributeError) as e:
                errors.append(f"vrs-python unable to normalize allele: {e}")
                return None

        if not allele:
            errors.append("Unable to get allele")
            return None

        seq_id, w = self.seqrepo_access.translate_identifier(
            allele.location.sequence_id._value, "ga4gh"
        )
        if seq_id:
            seq_id = seq_id[0]
            allele.location.sequence_id = seq_id
            allele.location.id = ga4gh_identify(allele.location)
            allele.id = ga4gh_identify(allele)
            allele_dict = allele.as_dict()
            try:
                Allele(**allele_dict)
            except ValidationError as e:
                errors.append(str(e))
                return None
            else:
                return allele_dict
        else:
            errors.append(w)
            return None

    def to_vrs_allele(
        self,
        ac: str,
        start: int,
        end: int,
        coordinate: AnnotationLayer,
        alt_type: AltType,
        errors: List[str],
        cds_start: Optional[int] = None,
        alt: Optional[str] = None,
    ) -> Optional[Dict]:
        """Translate accession and position to VRS Allele Object.

        :param ac: Accession
        :param start: Start position change
        :param end: End position change
        :param coordinate: Coordinate used
        :param alt_type: Type of alteration
        :param errors: List of errors
        :param cds_start: Coding start site
        :param alt: Alteration
        :return: VRS Allele Object
        """
        coords = self.get_start_end(coordinate, start, end, cds_start, errors)
        if not coords:
            return None
        if coords[0] > coords[1]:
            new_end, new_start = coords
        else:
            new_start, new_end = coords

        # Right now, this follows HGVS conventions
        # This will change once we support other representations
        if alt_type == AltType.INSERTION:
            state = alt
            new_end = new_start
        elif alt_type in {
            AltType.SUBSTITUTION,
            AltType.STOP_GAIN,
            AltType.DELETION,
            AltType.DELINS,
            AltType.REFERENCE_AGREE,
            AltType.NONSENSE,
        }:
            if alt_type == AltType.REFERENCE_AGREE:
                state, _ = self.seqrepo_access.get_reference_sequence(ac, new_start)
                if state is None:
                    errors.append(
                        f"Unable to get sequence on {ac} from " f"{new_start}"
                    )
                    return None
            else:
                state = alt or ""

            if alt_type == AltType.SUBSTITUTION:
                # This accounts for MNVs
                new_end += len(state) - 1

            new_start -= 1
        elif alt_type == AltType.DUPLICATION:
            ref, _ = self.seqrepo_access.get_reference_sequence(
                ac, new_start, new_end + 1
            )
            if ref is not None:
                state = ref + ref
            else:
                errors.append(
                    f"Unable to get sequence on {ac} from {new_start} to {new_end + 1}"
                )
                return None
            new_start -= 1
        else:
            errors.append(f"alt_type not supported: {alt_type}")
            return None

        start_vo = models.Number(value=new_start, type="Number")
        end_vo = models.Number(value=new_end, type="Number")
        sstate = models.LiteralSequenceExpression(
            sequence=state, type="LiteralSequenceExpression"
        )
        return self.vrs_allele(ac, start_vo, end_vo, sstate, alt_type, errors)
