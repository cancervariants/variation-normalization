"""Module for generating VRS objects"""
from typing import List, Optional, Tuple, Union, Dict

from ga4gh.vrs.dataproxy import SeqRepoDataProxy
from ga4gh.vrs import models, normalize
from ga4gh.core import ga4gh_identify
from uta_tools.data_sources import SeqRepoAccess
from bioutils.accessions import coerce_namespace


class VRSRepresentation:
    """Class for representing VRS objects"""

    def __init__(self, dp: SeqRepoDataProxy, seqrepo_access: SeqRepoAccess) -> None:
        """Initialize the VRSRepresentation class

        :param SeqRepoAccess seqrepo_access: Access to SeqRepo via UTA Tools
        :param SeqRepoDataProxy dp: Access to SeqRepo via VRS Python
        """
        self.seqrepo_access = seqrepo_access
        self.dp = dp

    @staticmethod
    def get_ival_start_end(
            coordinate: str, start: int, end: int, cds_start: int,
            errors: List) -> Optional[Tuple[int, int]]:
        """Get ival_start and ival_end coordinates.

        :param str coordinate: Coordinate used. Must be either `p`, `c`, or `g`
        :param int start: Start position change
        :param int end: End position change
        :param int cds_start: Coding start site
        :param List errors: List of errors
        :return: Tuple[ival_start, ival_end]
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
    def get_ival_certain_range(start1: int, start2: int, end1: int,
                               end2: int) -> models.SequenceInterval:
        """Return sequence interval

        :param int start1: Start left pos (assumes 1-based)
        :param int start2: Start right pos (assumes 1-based)
        :param int end1: End left pos (assumes 1-based)
        :param int end2: End right pos (assumes 1-based)
        :return: Sequence Interval model
        """
        return models.SequenceInterval(
            start=models.DefiniteRange(min=start1 - 1, max=start2 - 1,
                                       type="DefiniteRange"),
            end=models.DefiniteRange(min=end1 + 1, max=end2 + 1,
                                     type="DefiniteRange"),
            type="SequenceInterval"
        )

    @staticmethod
    def get_sequence_loc(
            ac: str, interval: models.SequenceInterval) -> models.Location:
        """Return VRS location

        :param str ac: Accession
        :param models.SequenceInterval interval: VRS sequence interval
        :return: VRS Location model
        """
        return models.SequenceLocation(
            sequence_id=coerce_namespace(ac),
            interval=interval, type="SequenceLocation")

    def vrs_allele(self, ac: str, interval: models.SequenceInterval,
                   sstate: Union[models.LiteralSequenceExpression,
                                 models.DerivedSequenceExpression,
                                 models.RepeatedSequenceExpression],
                   alt_type: str, errors: List) -> Optional[Dict]:
        """Create a VRS Allele object.

        :param str ac: Accession
        :param SequenceInterval interval: Sequence Interval
        :param sstate: State
        :type sstate: models.LiteralSequenceExpression or
            models.DerivedSequenceExpression or
            models.RepeatedSequenceExpression
        :param str alt_type: Type of alteration
        :param List errors: List of errors
        :return: VRS Allele object represented as a Dict
        """
        try:
            location = self.get_sequence_loc(ac, interval)
        except ValueError as e:
            errors.append(f"Unable to get sequence location: {e}")
            return None
        allele = models.Allele(location=location, state=sstate, type="Allele")
        # Ambiguous regions do not get normalized
        if alt_type not in ["uncertain_deletion", "uncertain_duplication",
                            "duplication_range", "deletion_range"]:
            try:
                allele = normalize(allele, self.dp)
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
            allele.location._id = ga4gh_identify(allele.location)
            allele._id = ga4gh_identify(allele)
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
        ival_coords = self.get_ival_start_end(coordinate, start, end,
                                              cds_start, errors)
        if not ival_coords:
            return None
        if ival_coords[0] > ival_coords[1]:
            ival_end, ival_start = ival_coords
        else:
            ival_start, ival_end = ival_coords

        # Right now, this follows HGVS conventions
        # This will change once we support other representations
        if alt_type == "insertion":
            state = alt
            ival_end = ival_start
        elif alt_type in ["substitution", "deletion", "delins",
                          "silent_mutation", "nonsense"]:
            if alt_type == "silent_mutation":
                state, _ = self.seqrepo_access.get_reference_sequence(ac, ival_start)
                if state is None:
                    errors.append(f"Unable to get sequence on {ac} from "
                                  f"{ival_start}")
                    return None
            else:
                state = alt or ""
            ival_start -= 1
        elif alt_type == "duplication":
            ref, _ = self.seqrepo_access.get_reference_sequence(
                ac, ival_start, ival_end + 1)
            if ref is not None:
                state = ref + ref
            else:
                errors.append(f"Unable to get sequence on {ac} from "
                              f"{ival_start} to {ival_end + 1}")
                return None
            ival_start -= 1
        else:
            errors.append(f"alt_type not supported: {alt_type}")
            return None

        interval = models.SequenceInterval(
            start=models.Number(value=ival_start, type="Number"),
            end=models.Number(value=ival_end, type="Number"),
            type="SequenceInterval")
        sstate = models.LiteralSequenceExpression(sequence=state,
                                                  type="LiteralSequenceExpression")
        return self.vrs_allele(ac, interval, sstate, alt_type, errors)

    def to_vrs_allele_ranges(
            self, ac: str, coordinate: str, alt_type: str, errors: List,
            ival: models.SequenceInterval) -> Optional[Dict]:
        """Translate variation ranges to VRS Allele Object.

        :param str ac: Accession
        :param str coordinate: Coordinate used. Must be either `p`, `c`, or `g`
        :param str alt_type: Type of alteration
        :param List errors: List of errors
        :param models.SequenceInterval ival: Sequence Interval
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

        return self.vrs_allele(ac, ival, sstate, alt_type, errors)
