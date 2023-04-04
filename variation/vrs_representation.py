"""Module for generating VRS objects"""
from typing import List, Optional, Tuple, Union, Dict

from ga4gh.vrsatile.pydantic.vrs_models import CURIE, RelativeCopyClass
from ga4gh.vrs import models, normalize
from ga4gh.core import ga4gh_identify
from cool_seq_tool.data_sources import SeqRepoAccess
from bioutils.accessions import coerce_namespace


class VRSRepresentation:
    """Class for representing VRS objects"""

    def __init__(self, seqrepo_access: SeqRepoAccess) -> None:
        """Initialize the VRSRepresentation class

        :param SeqRepoAccess seqrepo_access: Access to SeqRepo
        """
        self.seqrepo_access = seqrepo_access

    @staticmethod
    def get_start_end(
            coordinate: str, start: int, end: int, cds_start: int,
            errors: List) -> Optional[Tuple[int, int]]:
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
        return models.IndefiniteRange(value=start - 1, comparator="<=",
                                      type="IndefiniteRange")

    @staticmethod
    def get_end_indef_range(end: int) -> models.IndefiniteRange:
        """Return indefinite range given end coordinate

        :param int end: End position (assumes 1-based)
        :return: Indefinite range model
        """
        return models.IndefiniteRange(value=end, comparator=">=",
                                      type="IndefiniteRange")

    @staticmethod
    def get_start_end_definite_range(
        start1: int, start2: int, end1: int, end2: int
    ) -> Tuple[models.DefiniteRange, models.DefiniteRange]:
        """Return sequence start and end as definite ranges

        :param int start1: Start left pos (assumes 1-based)
        :param int start2: Start right pos (assumes 1-based)
        :param int end1: End left pos (assumes 1-based)
        :param int end2: End right pos (assumes 1-based)
        :return: Start and end as definite ranges
        """
        start = models.DefiniteRange(min=start1 - 1, max=start2 - 1,
                                     type="DefiniteRange")
        end = models.DefiniteRange(min=end1 + 1, max=end2 + 1, type="DefiniteRange")
        return start, end

    @staticmethod
    def get_sequence_loc(
        ac: str,
        start: Union[models.Number, models.DefiniteRange, models.IndefiniteRange],
        end: Union[models.Number, models.DefiniteRange, models.IndefiniteRange]
    ) -> models.Location:
        """Return VRS location

        :param str ac: Accession
        :param start: start pos
        :type start:
            models.Number or models.DefiniteRange or models.IndefiniteRange
        :param end: end pos
        :type end:
            models.Number or models.DefiniteRange or models.IndefiniteRange
        :return: VRS Location model
        """
        return models.SequenceLocation(
            sequence_id=coerce_namespace(ac),
            start=start,
            end=end,
            type="SequenceLocation")

    def vrs_allele(
        self,
        ac: str,
        start: Union[models.Number, models.DefiniteRange, models.IndefiniteRange],
        end: Union[models.Number, models.DefiniteRange, models.IndefiniteRange],
        sstate: Union[models.LiteralSequenceExpression,
                      models.DerivedSequenceExpression,
                      models.RepeatedSequenceExpression],
        alt_type: str,
        errors: List
    ) -> Optional[Dict]:
        """Create a VRS Allele object.

        :param str ac: Accession
        :param start: start pos
        :type start:
            models.Number or models.DefiniteRange or models.IndefiniteRange
        :param end: end pos
        :type end:
            models.Number or models.DefiniteRange or models.IndefiniteRange
        :param sstate: State
        :type sstate: models.LiteralSequenceExpression or
            models.DerivedSequenceExpression or
            models.RepeatedSequenceExpression
        :param str alt_type: Type of alteration
        :param List errors: List of errors
        :return: VRS Allele object represented as a Dict
        """
        try:
            location = self.get_sequence_loc(ac, start, end)
        except ValueError as e:
            errors.append(f"Unable to get sequence location: {e}")
            return None
        allele = models.Allele(location=location, state=sstate, type="Allele")
        # Ambiguous regions do not get normalized
        if alt_type not in ["uncertain_deletion", "uncertain_duplication",
                            "duplication_range", "deletion_range"]:
            try:
                allele = normalize(allele, self.seqrepo_access)
            except (KeyError, AttributeError) as e:
                errors.append(f"vrs-python unable to normalize allele: {e}")
                return None

        if not allele:
            errors.append("Unable to get allele")
            return None

        seq_id, w = self.seqrepo_access.translate_identifier(
            allele.location.sequence_id._value, "ga4gh")
        if seq_id:
            seq_id = seq_id[0]
            allele.location.sequence_id = seq_id
            allele.location.id = ga4gh_identify(allele.location)
            allele.id = ga4gh_identify(allele)
            return allele.as_dict()
        else:
            errors.append(w)
            return None

    def to_vrs_allele(
            self, ac: str, start: int, end: int, coordinate: str,
            alt_type: str, errors: List, cds_start: int = None,
            alt: str = None) -> Optional[Dict]:
        """Translate accession and position to VRS Allele Object.

        :param str ac: Accession
        :param int start: Start position change
        :param int end: End position change
        :param str coordinate: Coordinate used. Must be either `p`, `c`, or `g`
        :param str alt_type: Type of alteration
        :param List errors: List of errors
        :param int cds_start: Coding start site
        :param str alt: Alteration
        :return: VRS Allele Object
        """
        coords = self.get_start_end(coordinate, start, end, cds_start, errors)
        if not coords:
            return None
        if coords[0] > coords[1]:
            end, start = coords
        else:
            start, end = coords

        # Right now, this follows HGVS conventions
        # This will change once we support other representations
        if alt_type == "insertion":
            state = alt
            end = start
        elif alt_type in ["substitution", "deletion", "delins",
                          "silent_mutation", "nonsense"]:
            if alt_type == "silent_mutation":
                state, _ = self.seqrepo_access.get_reference_sequence(ac, start)
                if state is None:
                    errors.append(f"Unable to get sequence on {ac} from {start}")
                    return None
            else:
                state = alt or ""
            start -= 1
        elif alt_type == "duplication":
            ref, _ = self.seqrepo_access.get_reference_sequence(
                ac, start, end + 1)
            if ref is not None:
                state = ref + ref
            else:
                errors.append(f"Unable to get sequence on {ac} from "
                              f"{start} to {end + 1}")
                return None
            start -= 1
        else:
            errors.append(f"alt_type not supported: {alt_type}")
            return None

        start = models.Number(value=start, type="Number")
        end = models.Number(value=end, type="Number")
        sstate = models.LiteralSequenceExpression(sequence=state,
                                                  type="LiteralSequenceExpression")
        return self.vrs_allele(ac, start, end, sstate, alt_type, errors)

    def to_vrs_allele_ranges(
        self, ac: str, coordinate: str, alt_type: str, errors: List,
        start: Union[models.Number, models.DefiniteRange, models.IndefiniteRange],
        end: Union[models.Number, models.DefiniteRange, models.IndefiniteRange]
    ) -> Optional[Dict]:
        """Translate variation ranges to VRS Allele Object.

        :param str ac: Accession
        :param str coordinate: Coordinate used. Must be either `p`, `c`, or `g`
        :param str alt_type: Type of alteration
        :param List errors: List of errors
        :param start: start pos
        :type start:
            models.Number or models.DefiniteRange or models.IndefiniteRange
        :param end: end pos
        :type end:
            models.Number or models.DefiniteRange or models.IndefiniteRange
        :return: VRS Allele object
        """
        if coordinate == "c":
            # TODO: Once we add support for ranges on c. coord
            return None
        if alt_type in ["uncertain_deletion", "uncertain_duplication",
                        "duplication_range", "deletion_range"]:
            sstate = models.LiteralSequenceExpression(
                sequence="", type="LiteralSequenceExpression"
            )
        else:
            errors.append("No state")
            return None

        return self.vrs_allele(ac, start, end, sstate, alt_type, errors)

    def to_rel_cnv(
        self,
        location: Union[models.SequenceLocation, models.ChromosomeLocation, CURIE],
        relative_copy_class: RelativeCopyClass
    ) -> Dict:
        """Return VRS Relative Copy Number Variation

        :param location: Location for relative copy number
        :type location: models.SequenceLocation or models.ChromosomeLocation or CURIE
        :param RelativeCopyClass relative_copy_class: Relative copy class value
        :return: VRS Relative Copy Number Variation represented as a dictionary
        """
        location.id = ga4gh_identify(location)
        rcn = models.RelativeCopyNumber(
            location=location,
            relative_copy_class=relative_copy_class.value
        )
        rcn.id = ga4gh_identify(rcn)
        return rcn.as_dict()
