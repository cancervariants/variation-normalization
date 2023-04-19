"""Module containing schemas for services"""
from enum import Enum
from typing import Optional, Union, Dict, Any, Type

from pydantic import BaseModel, StrictStr
from ga4gh.vrsatile.pydantic.vrs_models import CopyNumberCount, Text, \
    SequenceLocation, CopyNumberChange

from variation.schemas.normalize_response_schema import ServiceResponse


class ClinVarAssembly(str, Enum):
    """Define assemblies in ClinVar"""

    GRCH38 = "GRCh38"
    GRCH37 = "GRCh37"
    NCBI36 = "NCBI36"
    HG38 = "hg38"
    HG19 = "hg19"
    HG18 = "hg18"


class ParsedToCnVarQuery(BaseModel):
    """Define query for parsed to copy number count variation endpoint"""

    assembly: Optional[ClinVarAssembly] = None
    chr: Optional[StrictStr] = None
    accession: Optional[StrictStr] = None
    start: int
    end: int
    total_copies: int


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
                    "start": 143134063,
                    "end": 143284670,
                    "total_copies": 3
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
