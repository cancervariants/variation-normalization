"""Module containing schemas for services"""

import re
from enum import Enum

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

from variation import __version__
from variation.schemas.normalize_response_schema import ServiceResponse


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


def validate_parsed_fields(cls, v: dict) -> dict:  # noqa: ARG001
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
    ac_assembly_chr_msg = (
        "Must provide either `accession` or both `assembly` and `chromosome`"
    )
    assembly = v.assembly
    chromosome = v.chromosome
    assembly_chr_set = assembly and chromosome
    assert v.accession or assembly_chr_set, ac_assembly_chr_msg

    if assembly_chr_set:
        pattern = r"^chr(X|Y|([1-9]|1[0-9]|2[0-2]))$"
        assert re.match(
            pattern, chromosome
        ), f"`chromosome`, {chromosome}, does not match r'{pattern}'"

    start0 = v.start0
    start1 = v.start1
    if v.start_pos_type == ParsedPosType.DEFINITE_RANGE:
        assert start1 is not None, "`start1` is required for definite ranges"
        assert start1 > start0, "`start0` must be less than `start1`"
    elif v.start_pos_type == ParsedPosType.INDEFINITE_RANGE:
        assert (
            v.start_pos_comparator
        ), "`start_pos_comparator` is required for indefinite ranges"

    end0 = v.end0
    end1 = v.end1
    if v.end_pos_type == ParsedPosType.DEFINITE_RANGE:
        assert end1 is not None, "`end1` is required for definite ranges"
        assert end1 > end0, "`end0` must be less than `end1`"
    elif v.end_pos_type == ParsedPosType.INDEFINITE_RANGE:
        assert (
            v.end_pos_comparator
        ), "`end_pos_comparator` is required for indefinite ranges"

    err_msg = "end positions must be greater than start"
    if start1 is None:
        assert end0 > start0, err_msg
    else:
        assert end0 > start1, err_msg


class ParsedToCopyNumberQuery(BaseModel):
    """Define base model for parsed to copy number queries"""

    assembly: ClinVarAssembly | None = Field(
        default=None,
        description=(
            "Assembly. Ignored, along with `chromosome`, if `accession` is " "provided."
        ),
    )
    chromosome: StrictStr | None = Field(
        default=None,
        description=(
            "Chromosome. Must contain `chr` prefix, i.e. 'chr7'. Must provide "
            "when `assembly` is provided."
        ),
    )
    accession: StrictStr | None = Field(
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
    start_pos_comparator: Comparator | None = Field(
        default=None,
        description=(
            "Must provide when `start_pos_type` is an Indefinite Range. "
            "Indicates which direction the range is indefinite. To represent "
            "(#_?), set to '<='. To represent (?_#), set to '>='."
        ),
    )
    end_pos_comparator: Comparator | None = Field(
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
    start1: StrictInt | None = Field(
        default=None,
        description=(
            "Only provided when `start_pos_type` is a Definite Range, this "
            "will be the max start position."
        ),
    )
    end1: StrictInt | None = Field(
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
    copies1: StrictInt | None = Field(
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
    copies_comparator: Comparator | None = Field(
        default=None,
        description=(
            "Must provide when `copies_type` is an Indefinite Range. "
            "Indicates which direction the range is indefinite."
        ),
    )

    @model_validator(mode="after")
    def validate_fields(cls, v: dict) -> dict:
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
            assert copies1, "`copies1` must be provided for `copies_type == ParsedPosType.DEFINITE_RANGE`"
        elif copies_type == ParsedPosType.INDEFINITE_RANGE:
            assert copies_comparator, "`copies_comparator` must be provided for `copies_type == ParsedPosType.INDEFINITE_RANGE`"

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

    copy_number_count: models.CopyNumberCount | None = None

    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "copy_number_count": {
                    "id": "ga4gh:CN.pbVk38-x5YGW7yhEtaBnWYjrzcb25L16",
                    "digest": "pbVk38-x5YGW7yhEtaBnWYjrzcb25L16",
                    "type": "CopyNumberCount",
                    "location": {
                        "id": "ga4gh:SL.6jZXELPqf5JDeN4CpOGde8foTUkHi1jy",
                        "digest": "6jZXELPqf5JDeN4CpOGde8foTUkHi1jy",
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
    def validate_fields(cls, v: dict) -> dict:
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

    copy_number_change: models.CopyNumberChange | None = None

    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "copy_number_change": {
                    "type": "CopyNumberChange",
                    "id": "ga4gh:CX.5kaJC-7Jj851bfJ6EipsHV413feg1T4T",
                    "digest": "5kaJC-7Jj851bfJ6EipsHV413feg1T4T",
                    "location": {
                        "type": "SequenceLocation",
                        "id": "ga4gh:SL.Iz_azSFTEulx7tCluLgGhE1n0hTLUocb",
                        "digest": "Iz_azSFTEulx7tCluLgGhE1n0hTLUocb",
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
    sequence_id: str | None = None
    start: int | None = None
    end: int | None = None
    sequence_location: models.SequenceLocation | None = None


class AmplificationToCxVarService(ServiceResponse):
    """A response for translating Amplification queries to Copy Number Change"""

    query: AmplificationToCxVarQuery | None = None
    amplification_label: str | None
    copy_number_change: models.CopyNumberChange | None

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
                    "id": "ga4gh:CX._UsXDMCLtPwsVKiNByhbwfS569K1wLWW",
                    "digest": "_UsXDMCLtPwsVKiNByhbwfS569K1wLWW",
                    "type": "CopyNumberChange",
                    "location": {
                        "id": "ga4gh:SL.0nPwKHYNnTmJ06G-gSmz8BEhB_NTp-0B",
                        "digest": "0nPwKHYNnTmJ06G-gSmz8BEhB_NTp-0B",
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
