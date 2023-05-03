"""Module containing schemas used in HGVS To Copy Number endpoints"""
from typing import Optional, Type, Any, Dict, Union

from ga4gh.vrsatile.pydantic.vrs_models import CopyChange, CopyNumberCount, \
    CopyNumberChange, Text
from pydantic import StrictStr

from variation.schemas.classification_response_schema import ClassificationType
from variation.schemas.normalize_response_schema import ServiceResponse


VALID_CLASSIFICATION_TYPES = [
    ClassificationType.GENOMIC_DUPLICATION,
    ClassificationType.GENOMIC_DELETION,
    ClassificationType.GENOMIC_DELETION_RANGE,
    ClassificationType.GENOMIC_UNCERTAIN_DELETION
]

VALID_COPY_CHANGE = [rcc.value for rcc in CopyChange.__members__.values()]


class HgvsToCopyNumberCountService(ServiceResponse):
    """A response for translating HGVS to copy number count."""

    hgvs_expr: StrictStr
    copy_number_count: Optional[Union[CopyNumberCount, Text]]

    class Config:
        """Configure model."""

        @staticmethod
        def schema_extra(schema: Dict[str, Any],
                         model: Type["HgvsToCopyNumberCountService"]) -> None:
            """Configure OpenAPI schema."""
            if "title" in schema.keys():
                schema.pop("title", None)
            for prop in schema.get("properties", {}).values():
                prop.pop("title", None)
            schema["example"] = {
                "hgvs_expr": "NC_000003.12:g.49531262dup",
                "copy_number_count": {
                    "_id": "ga4gh:CN.wIUwSQ9MQdv-2dsoDo-RjI97iK3Mn5m6",
                    "type": "CopyNumberCount",
                    "subject": {
                        "_id": "ga4gh:VSL.G_J9WrfooiONRgjbmGPuCBYbBYFQnYOg",
                        "type": "SequenceLocation",
                        "sequence_id": "ga4gh:SQ.Zu7h9AggXxhTaGVsy7h_EZSChSZGcmgX",
                        "interval": {
                            "type": "SequenceInterval",
                            "start": {"type": "Number", "value": 49531260},
                            "end": {"type": "Number", "value": 49531262}
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


class HgvsToCopyNumberChangeService(ServiceResponse):
    """A response for translating HGVS to copy number change."""

    hgvs_expr: StrictStr
    copy_number_change: Optional[Union[CopyNumberChange, Text]]

    class Config:
        """Configure model."""

        @staticmethod
        def schema_extra(schema: Dict[str, Any],
                         model: Type["HgvsToCopyNumberChangeService"]) -> None:
            """Configure OpenAPI schema."""
            if "title" in schema.keys():
                schema.pop("title", None)
            for prop in schema.get("properties", {}).values():
                prop.pop("title", None)
            schema["example"] = {
                "hgvs_expr": "NC_000003.12:g.49531262dup",
                "copy_number_change": {
                    "_id": "ga4gh:CX.hGuvyiJmDtx4-MRjsLja0fb_DqOE2chN",
                    "type": "CopyNumberChange",
                    "subject": {
                        "_id": "ga4gh:VSL.G_J9WrfooiONRgjbmGPuCBYbBYFQnYOg",
                        "type": "SequenceLocation",
                        "sequence_id": "ga4gh:SQ.Zu7h9AggXxhTaGVsy7h_EZSChSZGcmgX",
                        "interval": {
                            "type": "SequenceInterval",
                            "start": {"type": "Number", "value": 49531260},
                            "end": {"type": "Number", "value": 49531262}
                        }
                    },
                    "copy_change": "efo:0030069"
                },
                "service_meta_": {
                    "name": "variation-normalizer",
                    "version": "0.2.17",
                    "response_datetime": "2022-01-26T22:23:41.821673",
                    "url": "https://github.com/cancervariants/variation-normalization"
                }
            }
