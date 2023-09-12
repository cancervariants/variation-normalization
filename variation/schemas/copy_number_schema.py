"""Module containing schemas for services"""
import re
from enum import Enum
from typing import Dict, Optional

from ga4gh.vrs import models
from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    StrictBool,
    StrictInt,
    StrictStr,
    model_validator,
)

from variation.schemas.normalize_response_schema import ServiceResponse
from variation.version import __version__


class ParsedPosType(str, Enum):
    """Define position type for parsed to cnv endpoints"""

    NUMBER = "number"
    DEFINITE_RANGE = "definite_range"
    INDEFINITE_RANGE = "indefinite_range"


class Comparator(str, Enum):
    """A range comparator."""

    LT_OR_EQUAL = "<="
    GT_OR_EQUAL = ">="


class ClinVarAssembly(str, Enum):
    """Define assemblies in ClinVar"""

    GRCH38 = "GRCh38"
    GRCH37 = "GRCh37"
    NCBI36 = "NCBI36"
    HG38 = "hg38"
    HG19 = "hg19"
    HG18 = "hg18"


def validate_parsed_fields(cls, v: Dict) -> Dict:
    """Validate base copy number query fields
    - `accession` or both `assembly` and `chromosome` must be provided
    - `start1` is required when `start_pos_type` is a definite
    range.
    - `end1` is required when `end_pos_type` is a Definite Range.
    - `start_pos_comparator` is required when `start_pos_type` is an Indefinite
        Range
    - `end_pos_comparator` is required when `end_pos_type` is an Indefinite Range
    - End positions must be greater than start positions
    """
    ac_assembly_chr_msg = "Must provide either `accession` or both `assembly` and `chromosome`"  # noqa: E501
    assembly = v.assembly
    chromosome = v.chromosome
    assembly_chr_set = assembly and chromosome
    assert v.accession or assembly_chr_set, ac_assembly_chr_msg  # noqa: E501

    if assembly_chr_set:
        pattern = r"^chr(X|Y|([1-9]|1[0-9]|2[0-2]))$"
        assert re.match(
            pattern, chromosome
        ), f"`chromosome`, {chromosome}, does not match r'{pattern}'"  # noqa: E501

    start0 = v.start0
    start1 = v.start1
    if v.start_pos_type == ParsedPosType.DEFINITE_RANGE:
        assert start1 is not None, "`start1` is required for definite ranges"
        assert start1 > start0, "`start0` must be less than `start1`"
    elif v.start_pos_type == ParsedPosType.INDEFINITE_RANGE:
        assert (
            v.start_pos_comparator
        ), "`start_pos_comparator` is required for indefinite ranges"  # noqa: E501

    end0 = v.end0
    end1 = v.end1
    if v.end_pos_type == ParsedPosType.DEFINITE_RANGE:
        assert end1 is not None, "`end1` is required for definite ranges"
        assert end1 > end0, "`end0` must be less than `end1`"
    elif v.end_pos_type == ParsedPosType.INDEFINITE_RANGE:
        assert (
            v.end_pos_comparator
        ), "`end_pos_comparator` is required for indefinite ranges"  # noqa: E501

    err_msg = "end positions must be greater than start"
    if start1 is None:
        assert end0 > start0, err_msg
    else:
        assert end0 > start1, err_msg


class ParsedToCopyNumberQuery(BaseModel):
    """Define base model for parsed to copy number queries"""

    assembly: Optional[ClinVarAssembly] = Field(
        default=None,
        description=(
            "Assembly. Ignored, along with `chromosome`, if `accession` is " "provided."
        ),
    )
    chromosome: Optional[StrictStr] = Field(
        default=None,
        description=(
            "Chromosome. Must contain `chr` prefix, i.e. 'chr7'. Must provide "
            "when `assembly` is provided."
        ),
    )
    accession: Optional[StrictStr] = Field(
        default=None,
        description=(
            "Genomic RefSeq accession. If `accession` is provided, will "
            "ignore `assembly` and `chromosome`. If `accession` is not "
            "provided, must provide both `assembly` and `chromosome`."
        ),
    )
    start0: StrictInt = Field(
        description=(
            "Start position (residue coords). If `start_pos_type` is a "
            "Definite Range, this will be the min start position."
        ),
    )
    end0: StrictInt = Field(
        description=(
            "End position (residue coords). If `end_pos_type` is a definite "
            "range, this will be the min end position."
        ),
    )
    start_pos_comparator: Optional[Comparator] = Field(
        default=None,
        description=(
            "Must provide when `start_pos_type` is an Indefinite Range. "
            "Indicates which direction the range is indefinite. To represent "
            "(#_?), set to '<='. To represent (?_#), set to '>='."
        ),
    )
    end_pos_comparator: Optional[Comparator] = Field(
        default=None,
        description=(
            "Must provide when `end_pos_type` is an Indefinite Range. "
            "Indicates which direction the range is indefinite. To represent "
            "(#_?), set to '<='. To represent (?_#), set to '>='."
        ),
    )
    start_pos_type: ParsedPosType = Field(
        default=ParsedPosType.NUMBER,
        description="The type of the start value in the VRS SequenceLocation",
    )
    end_pos_type: ParsedPosType = Field(
        default=ParsedPosType.NUMBER,
        description="Type of the end value in the VRS SequenceLocation",
    )
    start1: Optional[StrictInt] = Field(
        default=None,
        description=(
            "Only provided when `start_pos_type` is a Definite Range, this "
            "will be the max start position."
        ),
    )
    end1: Optional[StrictInt] = Field(
        default=None,
        description=(
            "Only provided when `end_pos_type` is a Definite Range, this "
            "will be the max end position."
        ),
    )
    do_liftover: StrictBool = Field(
        default=False, description="Whether or not to liftover to GRCh38 assembly"
    )


class ParsedToCnVarQuery(ParsedToCopyNumberQuery):
    """Define query for parsed to copy number count variation endpoint"""

    copies0: StrictInt = Field(
        description=(
            "Number of copies. When `copies_type` is a Number or Indefinite "
            "Range, this will be the `value` for copies. When `copies_type` "
            "is an Definite Range, this will be the `min` copies."
        ),
    )
    copies1: Optional[StrictInt] = Field(
        default=None,
        description=(
            "Must provide when `copies_type` is a Definite Range. This will "
            "be the `max` copies."
        ),
    )
    copies_type: ParsedPosType = Field(
        default=ParsedPosType.NUMBER,
        description="Type for the `copies` in the `location`",
    )
    copies_comparator: Optional[Comparator] = Field(
        default=None,
        description=(
            "Must provide when `copies_type` is an Indefinite Range. "
            "Indicates which direction the range is indefinite."
        ),
    )

    @model_validator(mode="after")
    def validate_fields(cls, v: Dict) -> Dict:
        """Validate fields.

        - `copies1` should exist when `copies_type == ParsedPosType.DEFINITE_RANGE`
        - `copies_comparator` should exist when
            `copies_type == ParsedPosType.INDEFINITE_RANGE`
        """
        validate_parsed_fields(cls, v)
        copies1 = v.copies1
        copies_type = v.copies_type
        copies_comparator = v.copies_comparator

        if copies_type == ParsedPosType.DEFINITE_RANGE:
            assert (
                copies1
            ), "`copies1` must be provided for `copies_type == ParsedPosType.DEFINITE_RANGE`"  # noqa: E501
        elif copies_type == ParsedPosType.INDEFINITE_RANGE:
            assert (
                copies_comparator
            ), "`copies_comparator` must be provided for `copies_type == ParsedPosType.INDEFINITE_RANGE`"  # noqa: E501

        return v

    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "assembly": "GRCh37",
                "chromosome": "chr1",
                "accession": None,
                "start0": 143134063,
                "end0": 143284670,
                "copies0": 3,
                "copies1": None,
                "copies_comparator": None,
                "copies_type": "number",
                "start_pos_comparator": "<=",
                "end_pos_comparator": ">=",
                "start_pos_type": "indefinite_range",
                "end_pos_type": "indefinite_range",
                "start1": None,
                "end1": None,
                "do_liftover": False,
            }
        }
    )


class ParsedToCnVarService(ServiceResponse):
    """A response for translating parsed components to Copy Number Count"""

    copy_number_count: Optional[models.CopyNumberCount] = None

    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "copy_number_count": {
                    "id": "ga4gh:CN.Qrs0TaGCcJiibMvhcML6BTSCVtX95FBl",
                    "type": "CopyNumberCount",
                    "location": {
                        "id": "ga4gh:SL.g6xj5oKF99OysSxcfHyGYbh8NFNn2r61",
                        "type": "SequenceLocation",
                        "sequenceReference": {
                            "type": "SequenceReference",
                            "refgetAccession": "SQ.S_KjnFVz-FE7M0W6yoaUDgYxLPc1jyWU",
                        },
                        "start": [None, 143134062],
                        "end": [143284670, None],
                    },
                    "copies": 3,
                },
                "service_meta_": {
                    "name": "variation-normalizer",
                    "version": __version__,
                    "response_datetime": "2022-01-26T22:23:41.821673",
                    "url": "https://github.com/cancervariants/variation-normalization",
                },
            }
        }
    )


class ParsedToCxVarQuery(ParsedToCopyNumberQuery):
    """Define query for parsed to copy number change variation endpoint"""

    copy_change: models.CopyChange

    @model_validator(mode="after")
    def validate_fields(cls, v: Dict) -> Dict:
        """Validate fields"""
        validate_parsed_fields(cls, v)
        return v

    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "assembly": "GRCh38",
                "chromosome": "chrY",
                "accession": None,
                "start0": 10001,
                "end0": 1223133,
                "copy_change": "efo:0030069",
                "start_pos_type": "number",
                "end_pos_type": "number",
                "start1": None,
                "end1": None,
                "do_liftover": False,
            }
        }
    )


class ParsedToCxVarService(ServiceResponse):
    """A response for translating parsed components to Copy Number Change"""

    copy_number_change: Optional[models.CopyNumberChange] = None

    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "copy_number_change": {
                    "type": "CopyNumberChange",
                    "id": "ga4gh:CX.BTNwndSs3RylLhtL9Y45GePsVX35eeTT",
                    "location": {
                        "type": "SequenceLocation",
                        "id": "ga4gh:SL.Pu3oAKHColJSZ3zY_Xu5MeezINaTFlNq",
                        "sequenceReference": {
                            "type": "SequenceReference",
                            "refgetAccession": "SQ.8_liLu1aycC0tPQPFmUaGXJLDs5SbPZ5",
                        },
                        "start": 10000,
                        "end": 1223133,
                    },
                    "copyChange": "efo:0030069",
                },
                "service_meta_": {
                    "name": "variation-normalizer",
                    "version": __version__,
                    "response_datetime": "2022-01-26T22:23:41.821673",
                    "url": "https://github.com/cancervariants/variation-normalization",
                },
            }
        }
    )


class AmplificationToCxVarQuery(BaseModel):
    """Define query for amplification to copy number change variation endpoint"""

    gene: str
    sequence_id: Optional[str] = None
    start: Optional[int] = None
    end: Optional[int] = None
    sequence_location: Optional[models.SequenceLocation] = None


class AmplificationToCxVarService(ServiceResponse):
    """A response for translating Amplification queries to Copy Number Change"""

    query: Optional[AmplificationToCxVarQuery] = None
    amplification_label: Optional[str]
    copy_number_change: Optional[models.CopyNumberChange]

    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "query": {
                    "gene": "braf",
                    "sequence_id": None,
                    "start": None,
                    "end": None,
                    "sequence_location": None,
                },
                "amplification_label": "BRAF Amplification",
                "copy_number_change": {
                    "id": "ga4gh:CX.89PECTeQjhhXnNW9yg24DheWOQMgmKk2",
                    "type": "CopyNumberChange",
                    "location": {
                        "id": "ga4gh:SL.uNBZoxhjhohl24VlIut-JxPJAGfJ7EQE",
                        "type": "SequenceLocation",
                        "sequenceReference": {
                            "type": "SequenceReference",
                            "refgetAccession": "SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
                        },
                        "start": 140713327,
                        "end": 140924929,
                    },
                    "copyChange": "efo:0030072",
                },
                "service_meta_": {
                    "version": __version__,
                    "response_datetime": "2022-09-29T15:08:18.696882",
                    "name": "variation-normalizer",
                    "url": "https://github.com/cancervariants/variation-normalization",
                },
            }
        }
    )
