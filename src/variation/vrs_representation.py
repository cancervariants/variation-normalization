"""Module for generating VRS objects"""

from cool_seq_tool.handlers import SeqRepoAccess
from cool_seq_tool.schemas import AnnotationLayer, ResidueMode
from ga4gh.core import ga4gh_identify
from ga4gh.core.models import Extension
from ga4gh.vrs import models, normalize
from pydantic import ValidationError

from variation.schemas.token_response_schema import (
    AMBIGUOUS_REGIONS,
    AltType,
)
from variation.utils import get_refget_accession


class VRSRepresentation:
    """Class for representing VRS objects"""

    def __init__(self, seqrepo_access: SeqRepoAccess) -> None:
        """Initialize the VRSRepresentation class

        :param SeqRepoAccess seqrepo_access: Access to SeqRepo
        """
        self.seqrepo_access = seqrepo_access

    @staticmethod
    def get_start_end(
        coordinate: str, start: int, end: int, cds_start: int, errors: list
    ) -> tuple[int, int] | None:
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

        if coordinate == "c" and cds_start:
            start += cds_start
            end += cds_start
        return start, end

    @staticmethod
    def get_start_indef_range(start: int) -> models.Range:
        """Return indefinite range given start coordinate

        :param int start: Start position (assumes 1-based)
        :return: Range
        """
        return models.Range([None, start - 1])

    @staticmethod
    def get_end_indef_range(end: int) -> models.Range:
        """Return indefinite range given end coordinate

        :param int end: End position (assumes 1-based)
        :return: Range model
        """
        return models.Range([end, None])

    @staticmethod
    def get_sequence_loc(
        refget_accession: str,
        start: int | models.Range,
        end: int | models.Range,
    ) -> models.Location:
        """Return VRS location

        :param refget_accession: Refget accession (SQ.)
        :param start: start pos
        :param end: end pos
        :return: VRS Location model
        """
        return models.SequenceLocation(
            sequenceReference=models.SequenceReference(
                refgetAccession=refget_accession
            ),
            start=start,
            end=end,
        )

    def vrs_allele(
        self,
        ac: str,
        start: int | models.Range,
        end: int | models.Range,
        sstate: models.LiteralSequenceExpression | models.ReferenceLengthExpression,
        alt_type: AltType,
        errors: list[str],
        extensions: list[Extension] | None = None,
    ) -> dict | None:
        """Create a VRS Allele object.

        :param ac: Accession
        :param start: start pos
        :param end: end pos
        :param sstate: State
        :param alt_type: Type of alteration
        :param errors: List of errors
        :param extensions: List of extensions for variation
        :return: VRS Allele object represented as a Dict
        """
        refget_accession = get_refget_accession(self.seqrepo_access, ac, errors)
        if not refget_accession:
            return None

        try:
            location = self.get_sequence_loc(refget_accession, start, end)
        except ValueError as e:
            errors.append(f"Unable to get sequence location: {e}")
            return None
        allele = models.Allele(location=location, state=sstate, extensions=extensions)
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

        allele.location.id = ga4gh_identify(allele.location)
        allele.id = ga4gh_identify(allele)
        allele_dict = allele.model_dump(exclude_none=True)
        try:
            models.Allele(**allele_dict)
        except ValidationError as e:
            errors.append(str(e))
            return None
        else:
            return allele_dict

    def to_vrs_allele(
        self,
        ac: str,
        start: int,
        end: int,
        coordinate: AnnotationLayer,
        alt_type: AltType,
        errors: list[str],
        cds_start: int | None = None,
        alt: str | None = None,
        residue_mode: ResidueMode = ResidueMode.RESIDUE,
        extensions: list[Extension] | None = None,
    ) -> dict | None:
        """Translate accession and position to VRS Allele Object.

        :param ac: Accession
        :param start: Start position change
        :param end: End position change
        :param coordinate: Coordinate used
        :param alt_type: Type of alteration
        :param errors: List of errors
        :param cds_start: Coding start site
        :param alt: Alteration
        :param residue_mode: Residue mode for ``start`` and ``end`` positions
        :param extensions: List of extensions for variation
        :return: VRS Allele Object
        """
        coords = self.get_start_end(coordinate, start, end, cds_start, errors)
        if not coords:
            return None
        if coords[0] > coords[1]:
            new_end, new_start = coords
        else:
            new_start, new_end = coords

        if residue_mode == ResidueMode.RESIDUE:
            new_start -= 1
            residue_mode = ResidueMode.INTER_RESIDUE

        # Right now, this follows HGVS conventions
        # This will change once we support other representations
        if alt_type == AltType.INSERTION:
            state = alt
            new_start += 1
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
                state, _ = self.seqrepo_access.get_reference_sequence(
                    ac, start=new_start, end=new_end, residue_mode=residue_mode
                )
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

        elif alt_type == AltType.DUPLICATION:
            ref, _ = self.seqrepo_access.get_reference_sequence(
                ac, start=new_start, end=new_end, residue_mode=residue_mode
            )
            if ref is not None:
                state = ref + ref
            else:
                errors.append(
                    f"Unable to get sequence on {ac} from {new_start} to {new_end + 1}"
                )
                return None
        else:
            errors.append(f"alt_type not supported: {alt_type}")
            return None

        sstate = models.LiteralSequenceExpression(sequence=state)
        return self.vrs_allele(
            ac, new_start, new_end, sstate, alt_type, errors, extensions=extensions
        )
