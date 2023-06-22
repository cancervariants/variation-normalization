"""Module containing schemas for services"""
from enum import Enum
from typing import Optional, Union, Dict, Any, Type, Literal

from pydantic import BaseModel, StrictStr, root_validator, StrictInt
from pydantic.main import ModelMetaclass
from ga4gh.vrsatile.pydantic.vrs_models import CopyNumberCount, Text, \
    SequenceLocation, CopyNumberChange, CopyChange, VRSTypes

from variation.version import __version__
from variation.schemas.normalize_response_schema import ServiceResponse


class ClinVarAssembly(str, Enum):
    """Define assemblies in ClinVar"""

    GRCH38 = "GRCh38"
    GRCH37 = "GRCh37"
    NCBI36 = "NCBI36"
    HG38 = "hg38"
    HG19 = "hg19"
    HG18 = "hg18"


class ParsedToCopyNumberQuery(BaseModel):
    """Define base model for parsed to copy number queries"""

    assembly: Optional[ClinVarAssembly] = None
    chr: Optional[StrictStr] = None
    accession: Optional[StrictStr] = None
    start0: StrictInt
    end0: StrictInt
    start_pos_type: Literal[
        VRSTypes.NUMBER, VRSTypes.DEFINITE_RANGE, VRSTypes.INDEFINITE_RANGE
    ]
    end_pos_type: Literal[
        VRSTypes.NUMBER, VRSTypes.DEFINITE_RANGE, VRSTypes.INDEFINITE_RANGE
    ]
    start1: Optional[StrictInt]
    end1: Optional[StrictInt]

    @root_validator(pre=False, skip_on_failure=True)
    def validate_fields(cls: ModelMetaclass, v: Dict) -> Dict:
        """Validate fields.
        - `accession` or both `assembly` and `chr` must be provided
        - `start1` is required when `start_pos_type` is a definite
        range.
        - `end1` is required when `end_pos_type` is a definite range.
        - End positions must be greater than start positions
        """
        ac_assembly_chr_msg = "Must provide either `accession` or both `assembly` and `chr`"  # noqa: E501
        assert v.get("accession") or (v.get("assembly") and v.get("chr")), ac_assembly_chr_msg  # noqa: E501

        start0 = v["start0"]
        start1 = v.get("start1")
        if v["start_pos_type"] == VRSTypes.DEFINITE_RANGE:
            assert start1 is not None, "`start1` is required for definite ranges"
            assert start1 > start0, "`start0` must be less than `start1`"

        end0 = v["end0"]
        end1 = v.get("end1")
        if v["end_pos_type"] == VRSTypes.DEFINITE_RANGE:
            assert end1 is not None, "`end1` is required for definite ranges"
            assert end1 > end0, "`end0` must be less than `end1`"

        err_msg = "end positions must be greater than start"
        if start1 is None:
            assert end0 > start0, err_msg
        else:
            assert end0 > start1, err_msg

        return v


class ParsedToCnVarQuery(ParsedToCopyNumberQuery):
    """Define query for parsed to copy number count variation endpoint"""

    total_copies: StrictInt


class ParsedToCnVarService(ServiceResponse):
    """A response for translating parsed components to Copy Number Count"""

    query: Optional[ParsedToCnVarQuery] = None
    copy_number_count: Optional[Union[Text, CopyNumberCount]]

    class Config:
        """Configure model."""

        @staticmethod
        def schema_extra(schema: Dict[str, Any],
                         model: Type["ParsedToCnVarService"]) -> None:
            """Configure OpenAPI schema."""
            if "title" in schema.keys():
                schema.pop("title", None)
            for prop in schema.get("properties", {}).values():
                prop.pop("title", None)
            schema["example"] = {
                "query": {
                    "assembly": "GRCh37",
                    "chr": "1",
                    "accession": None,
                    "start0": 143134063,
                    "end0": 143284670,
                    "total_copies": 3,
                    "start_pos_type": "IndefiniteRange",
                    "end_pos_type": "IndefiniteRange",
                    "start1": None,
                    "end1": None
                },
                "copy_number_count": {
                    "id": "ga4gh:CN._IYaKE4CoDa01tkcgOuqPhnYbZ5RuPcj",
                    "type": "CopyNumberCount",
                    "subject": {
                        "id": "ga4gh:SL.RIgksXkT_kWCJv3poK4WQ9PK5_YSRBuh",
                        "type": "SequenceLocation",
                        "sequence_id": "ga4gh:SQ.S_KjnFVz-FE7M0W6yoaUDgYxLPc1jyWU",
                        "start": {
                            "type": "IndefiniteRange",
                            "value": 143134062,
                            "comparator": "<="
                        },
                        "end": {
                            "type": "IndefiniteRange",
                            "value": 143284670,
                            "comparator": ">="
                        }
                    },
                    "copies": {"type": "Number", "value": 3}
                },
                "service_meta_": {
                    "name": "variation-normalizer",
                    "version": "0.2.17",
                    "response_datetime": "2022-01-26T22:23:41.821673",
                    "url": "https://github.com/cancervariants/variation-normalization"
                }
            }


class ParsedToCxVarQuery(ParsedToCopyNumberQuery):
    """Define query for parsed to copy number change variation endpoint"""

    copy_change: CopyChange


class ParsedToCxVarService(ServiceResponse):
    """A response for translating parsed components to Copy Number Change"""

    query: Optional[ParsedToCxVarQuery] = None
    copy_number_change: Optional[Union[Text, CopyNumberChange]]

    class Config:
        """Configure model."""

        @staticmethod
        def schema_extra(schema: Dict[str, Any],
                         model: Type["ParsedToCxVarService"]) -> None:
            """Configure OpenAPI schema."""
            if "title" in schema.keys():
                schema.pop("title", None)
            for prop in schema.get("properties", {}).values():
                prop.pop("title", None)
            schema["example"] = {
                "query": {
                    "assembly": "GRCh38",
                    "chr": "chrY",
                    "accession": None,
                    "start0": 10001,
                    "end0": 1223133,
                    "copy_change": "efo:0030069",
                    "start_pos_type": "Number",
                    "end_pos_type": "Number",
                    "start1": None,
                    "end1": None
                },
                "copy_number_change": {
                    "type": "CopyNumberChange",
                    "id": "ga4gh:CX.UirzxujWnAIklYHh4VxSnFglfDROHYv6",
                    "subject": {
                        "type": "SequenceLocation",
                        "id": "ga4gh:SL.x075Sp6tCfGZcpHHmJ1e5oUdAW0CvN0X",
                        "sequence_id": "ga4gh:SQ.8_liLu1aycC0tPQPFmUaGXJLDs5SbPZ5",
                        "start": {
                            "type": "Number",
                            "value": 10000
                        },
                        "end": {
                            "type": "Number",
                            "value": 1223133
                        }
                    },
                    "copy_change": "efo:0030069"
                },
                "service_meta_": {
                    "name": "variation-normalizer",
                    "version": __version__,
                    "response_datetime": "2022-01-26T22:23:41.821673",
                    "url": "https://github.com/cancervariants/variation-normalization"
                }
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
        def schema_extra(schema: Dict[str, Any],
                         model: Type["AmplificationToCxVarService"]) -> None:
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
                    "sequence_location": None
                },
                "amplification_label": "BRAF Amplification",
                "copy_number_change": {
                    "id": "ga4gh:CX.1RJp1zW60x2t4Exc4965_a3CvYFtsL4q",
                    "type": "CopyNumberChange",
                    "subject": {
                        "id": "ga4gh:SL.po-AExwyqkstDx3JWYn6plIlxn5eojv4",
                        "type": "SequenceLocation",
                        "sequence_id": "ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
                        "start": {"type": "Number", "value": 140713327},
                        "end": {"type": "Number", "value": 140924929}
                    },
                    "copy_change": "efo:0030072"
                },
                "service_meta_": {
                    "version": "0.7.dev0",
                    "response_datetime": "2022-09-29T15:08:18.696882",
                    "name": "variation-normalizer",
                    "url": "https://github.com/cancervariants/variation-normalization"
                }
            }
