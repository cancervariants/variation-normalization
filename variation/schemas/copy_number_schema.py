"""Module containing schemas for services"""
import re
from enum import Enum
from typing import Any, Dict, Literal, Optional, Type, Union

from ga4gh.vrsatile.pydantic.vrs_models import (
    Comparator,
    CopyChange,
    CopyNumberChange,
    CopyNumberCount,
    SequenceLocation,
    Text,
    VRSTypes,
)
from pydantic import BaseModel, Field, StrictBool, StrictInt, StrictStr, root_validator
from pydantic.main import ModelMetaclass

from variation.schemas.normalize_response_schema import ServiceResponse
from variation.version import __version__


class ClinVarAssembly(str, Enum):
    """Define assemblies in ClinVar"""

    GRCH38 = "GRCh38"
    GRCH37 = "GRCh37"
    NCBI36 = "NCBI36"
    HG38 = "hg38"
    HG19 = "hg19"
    HG18 = "hg18"


def validate_parsed_fields(cls: ModelMetaclass, v: Dict) -> Dict:
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
    assembly = v.get("assembly")
    chromosome = v.get("chromosome")
    assembly_chr_set = assembly and chromosome
    assert v.get("accession") or assembly_chr_set, ac_assembly_chr_msg

    if assembly_chr_set:
        pattern = r"^chr(X|Y|([1-9]|1[0-9]|2[0-2]))$"
        assert re.match(
            pattern, chromosome
        ), f"`chromosome`, {chromosome}, does not match r'{pattern}'"

    start0 = v["start0"]
    start1 = v.get("start1")
    if v["start_pos_type"] == VRSTypes.DEFINITE_RANGE:
        assert start1 is not None, "`start1` is required for definite ranges"
        assert start1 > start0, "`start0` must be less than `start1`"
    elif v["start_pos_type"] == VRSTypes.INDEFINITE_RANGE:
        assert v.get(
            "start_pos_comparator"
        ), "`start_pos_comparator` is required for indefinite ranges"

    end0 = v["end0"]
    end1 = v.get("end1")
    if v["end_pos_type"] == VRSTypes.DEFINITE_RANGE:
        assert end1 is not None, "`end1` is required for definite ranges"
        assert end1 > end0, "`end0` must be less than `end1`"
    elif v["end_pos_type"] == VRSTypes.INDEFINITE_RANGE:
        assert v.get(
            "end_pos_comparator"
        ), "`end_pos_comparator` is required for indefinite ranges"

    err_msg = "end positions must be greater than start"
    if start1 is None:
        assert end0 > start0, err_msg
    else:
        assert end0 > start1, err_msg


class ParsedToCopyNumberQuery(BaseModel):
    """Define base model for parsed to copy number queries"""

    assembly: Optional[ClinVarAssembly] = Field(
        None,
        description=(
            "Assembly. Ignored, along with `chromosome`, if `accession` is " "provided."
        ),
    )
    chromosome: Optional[StrictStr] = Field(
        None,
        description=(
            "Chromosome. Must contain `chr` prefix, i.e. 'chr7'. Must provide "
            "when `assembly` is provided."
        ),
    )
    accession: Optional[StrictStr] = Field(
        None,
        description=(
            "Genomic RefSeq accession. If `accession` is provided, will "
            "ignore `assembly` and `chromosome`. If `accession` is not "
            "provided, must provide both `assembly` and `chromosome`."
        ),
    )
    start0: StrictInt = Field(
        ...,
        description=(
            "Start position (residue coords). If `start_pos_type` is a "
            "Definite Range, this will be the min start position."
        ),
    )
    end0: StrictInt = Field(
        ...,
        description=(
            "End position (residue coords). If `end_pos_type` is a definite "
            "range, this will be the min end position."
        ),
    )
    start_pos_comparator: Optional[Comparator] = Field(
        None,
        description=(
            "Must provide when `start_pos_type` is an Indefinite Range. "
            "Indicates which direction the range is indefinite. To represent "
            "(#_?), set to '<='. To represent (?_#), set to '>='."
        ),
    )
    end_pos_comparator: Optional[Comparator] = Field(
        None,
        description=(
            "Must provide when `end_pos_type` is an Indefinite Range. "
            "Indicates which direction the range is indefinite. To represent "
            "(#_?), set to '<='. To represent (?_#), set to '>='."
        ),
    )
    start_pos_type: Literal[
        VRSTypes.NUMBER, VRSTypes.DEFINITE_RANGE, VRSTypes.INDEFINITE_RANGE
    ] = Field(
        VRSTypes.NUMBER,
        description="The type of the start value in the VRS SequenceLocation",
    )
    end_pos_type: Literal[
        VRSTypes.NUMBER, VRSTypes.DEFINITE_RANGE, VRSTypes.INDEFINITE_RANGE
    ] = Field(
        VRSTypes.NUMBER, description="Type of the end value in the VRS SequenceLocation"
    )
    start1: Optional[StrictInt] = Field(
        None,
        description=(
            "Only provided when `start_pos_type` is a Definite Range, this "
            "will be the max start position."
        ),
    )
    end1: Optional[StrictInt] = Field(
        None,
        description=(
            "Only provided when `end_pos_type` is a Definite Range, this "
            "will be the max end position."
        ),
    )
    do_liftover: StrictBool = Field(
        False, description="Whether or not to liftover to GRCh38 assembly"
    )
    untranslatable_returns_text: StrictBool = Field(
        False,
        description=(
            "When set to `True`, return VRS Text Object when unable to "
            "translate or normalize query. When set to `False`, return `None` "
            "when unable to translate or normalize query."
        ),
    )


class ParsedToCnVarQuery(ParsedToCopyNumberQuery):
    """Define query for parsed to copy number count variation endpoint"""

    copies0: StrictInt = Field(
        ...,
        description=(
            "Number of copies. When `copies_type` is a Number or Indefinite "
            "Range, this will be the `value` for copies. When `copies_type` "
            "is an Definite Range, this will be the `min` copies."
        ),
    )
    copies1: Optional[StrictInt] = Field(
        None,
        description=(
            "Must provide when `copies_type` is a Definite Range. This will "
            "be the `max` copies."
        ),
    )
    copies_type: Literal[
        VRSTypes.NUMBER, VRSTypes.DEFINITE_RANGE, VRSTypes.INDEFINITE_RANGE
    ] = Field(VRSTypes.NUMBER, description="Type for the `copies` in the `subject`")
    copies_comparator: Optional[Comparator] = Field(
        None,
        description=(
            "Must provide when `copies_type` is an Indefinite Range. "
            "Indicates which direction the range is indefinite."
        ),
    )

    @root_validator(pre=False, skip_on_failure=True)
    def validate_fields(cls: ModelMetaclass, v: Dict) -> Dict:
        """Validate fields.

        - `copies1` should exist when `copies_type == VRSTypes.DEFINITE_RANGE`
        - `copies_comparator` should exist when
            `copies_type == VRSTypes.INDEFINITE_RANGE`
        """
        validate_parsed_fields(cls, v)
        copies1 = v.get("copies1")
        copies_type = v.get("copies_type")
        copies_comparator = v.get("copies_comparator")

        if copies_type == VRSTypes.DEFINITE_RANGE:
            assert (
                copies1
            ), "`copies1` must be provided for `copies_type == DefiniteRange`"
        elif copies_type == VRSTypes.INDEFINITE_RANGE:
            assert (
                copies_comparator
            ), "`copies_comparator` must be provided for `copies_type == IndefiniteRange`"  # noqa: E501

        return v

    class Config:
        """Configure model."""

        @staticmethod
        def schema_extra(
            schema: Dict[str, Any], model: Type["ParsedToCnVarQuery"]
        ) -> None:
            """Configure OpenAPI schema."""
            if "title" in schema.keys():
                schema.pop("title", None)
            for prop in schema.get("properties", {}).values():
                prop.pop("title", None)
            schema["example"] = {
                "assembly": "GRCh37",
                "chromosome": "chr1",
                "accession": None,
                "start0": 143134063,
                "end0": 143284670,
                "copies0": 3,
                "copies1": None,
                "copies_comparator": None,
                "copies_type": "Number",
                "start_pos_comparator": "<=",
                "end_pos_comparator": ">=",
                "start_pos_type": "IndefiniteRange",
                "end_pos_type": "IndefiniteRange",
                "start1": None,
                "end1": None,
                "do_liftover": False,
                "untranslatable_returns_text": False,
            }


class ParsedToCnVarService(ServiceResponse):
    """A response for translating parsed components to Copy Number Count"""

    copy_number_count: Optional[Union[Text, CopyNumberCount]]

    class Config:
        """Configure model."""

        @staticmethod
        def schema_extra(
            schema: Dict[str, Any], model: Type["ParsedToCnVarService"]
        ) -> None:
            """Configure OpenAPI schema."""
            if "title" in schema.keys():
                schema.pop("title", None)
            for prop in schema.get("properties", {}).values():
                prop.pop("title", None)
            schema["example"] = {
                "copy_number_count": {
                    "_id": "ga4gh:CN.N6C9rWBjrNuiIhJkPxdPlRKvSGKoFynr",
                    "type": "CopyNumberCount",
                    "subject": {
                        "_id": "ga4gh:VSL.JTsxd9PiPZaIPL9Tl3ss78GYYnDeogvf",
                        "type": "SequenceLocation",
                        "sequence_id": "ga4gh:SQ.S_KjnFVz-FE7M0W6yoaUDgYxLPc1jyWU",
                        "interval": {
                            "start": {
                                "type": "IndefiniteRange",
                                "value": 143134062,
                                "comparator": "<=",
                            },
                            "end": {
                                "type": "IndefiniteRange",
                                "value": 143284670,
                                "comparator": ">=",
                            },
                        },
                    },
                    "copies": {"type": "Number", "value": 3},
                },
                "service_meta_": {
                    "name": "variation-normalizer",
                    "version": "0.2.17",
                    "response_datetime": "2022-01-26T22:23:41.821673",
                    "url": "https://github.com/cancervariants/variation-normalization",
                },
            }


class ParsedToCxVarQuery(ParsedToCopyNumberQuery):
    """Define query for parsed to copy number change variation endpoint"""

    copy_change: CopyChange

    @root_validator(pre=False, skip_on_failure=True)
    def validate_fields(cls: ModelMetaclass, v: Dict) -> Dict:
        """Validate fields"""
        validate_parsed_fields(cls, v)
        return v

    class Config:
        """Configure model."""

        @staticmethod
        def schema_extra(
            schema: Dict[str, Any], model: Type["ParsedToCxVarQuery"]
        ) -> None:
            """Configure OpenAPI schema."""
            if "title" in schema.keys():
                schema.pop("title", None)
            for prop in schema.get("properties", {}).values():
                prop.pop("title", None)
            schema["example"] = {
                "assembly": "GRCh38",
                "chromosome": "chrY",
                "accession": None,
                "start0": 10001,
                "end0": 1223133,
                "copy_change": "efo:0030069",
                "start_pos_type": "Number",
                "end_pos_type": "Number",
                "start1": None,
                "end1": None,
                "do_liftover": False,
                "untranslatable_returns_text": False,
            }


class ParsedToCxVarService(ServiceResponse):
    """A response for translating parsed components to Copy Number Change"""

    copy_number_change: Optional[Union[Text, CopyNumberChange]]

    class Config:
        """Configure model."""

        @staticmethod
        def schema_extra(
            schema: Dict[str, Any], model: Type["ParsedToCxVarService"]
        ) -> None:
            """Configure OpenAPI schema."""
            if "title" in schema.keys():
                schema.pop("title", None)
            for prop in schema.get("properties", {}).values():
                prop.pop("title", None)
            schema["example"] = {
                "copy_number_change": {
                    "type": "CopyNumberChange",
                    "_id": "ga4gh:CX.KYAQwf8-DQu23LsDbFHP0BiRzrCmu46x",
                    "subject": {
                        "type": "SequenceLocation",
                        "_id": "ga4gh:VSL.djv1Oq_qNjDialZakqQGoUAJvohREPBL",
                        "sequence_id": "ga4gh:SQ.8_liLu1aycC0tPQPFmUaGXJLDs5SbPZ5",
                        "interval": {
                            "start": {"type": "Number", "value": 10000},
                            "end": {"type": "Number", "value": 1223133},
                        },
                    },
                    "copy_change": "efo:0030069",
                },
                "service_meta_": {
                    "name": "variation-normalizer",
                    "version": __version__,
                    "response_datetime": "2022-01-26T22:23:41.821673",
                    "url": "https://github.com/cancervariants/variation-normalization",
                },
            }


class AmplificationToCxVarQuery(BaseModel):
    """Define query for amplification to copy number change variation endpoint"""

    gene: str
    sequence_id: Optional[str]
    start: Optional[int]
    end: Optional[int]
    sequence_location: Optional[SequenceLocation]


class AmplificationToCxVarService(ServiceResponse):
    """A response for translating Amplification queries to Copy Number Change"""

    query: Optional[AmplificationToCxVarQuery] = None
    amplification_label: Optional[str]
    copy_number_change: Optional[Union[Text, CopyNumberChange]]

    class Config:
        """Configure model."""

        @staticmethod
        def schema_extra(
            schema: Dict[str, Any], model: Type["AmplificationToCxVarService"]
        ) -> None:
            """Configure OpenAPI schema."""
            if "title" in schema.keys():
                schema.pop("title", None)
            for prop in schema.get("properties", {}).values():
                prop.pop("title", None)
            schema["example"] = {
                "query": {
                    "gene": "braf",
                    "sequence_id": None,
                    "start": None,
                    "end": None,
                    "sequence_location": None,
                },
                "amplification_label": "BRAF Amplification",
                "copy_number_change": {
                    "_id": "ga4gh:CX.TZBOQe5xFojvFJ1XjQQD0633rStHtGUs",
                    "type": "CopyNumberChange",
                    "subject": {
                        "_id": "ga4gh:VSL.xZU3kL8F6t2ca6WH_26CWKfNW9-owhR4",
                        "type": "SequenceLocation",
                        "sequence_id": "ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
                        "interval": {
                            "start": {"type": "Number", "value": 140713327},
                            "end": {"type": "Number", "value": 140924929},
                        },
                    },
                    "copy_change": "efo:0030072",
                },
                "service_meta_": {
                    "version": "0.7.dev0",
                    "response_datetime": "2022-09-29T15:08:18.696882",
                    "name": "variation-normalizer",
                    "url": "https://github.com/cancervariants/variation-normalization",
                },
            }
