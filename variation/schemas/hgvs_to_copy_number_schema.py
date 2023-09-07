"""Module containing schemas used in HGVS To Copy Number endpoints"""
from typing import Any, Dict, Optional, Type, Union

from ga4gh.vrsatile.pydantic.vrs_models import CopyNumberChange, CopyNumberCount, Text
from pydantic import StrictStr

from variation.schemas.normalize_response_schema import ServiceResponse


class HgvsToCopyNumberCountService(ServiceResponse):
    """A response for translating HGVS to copy number count."""

    hgvs_expr: StrictStr
    copy_number_count: Optional[Union[CopyNumberCount, Text]]

    class Config:
        """Configure model."""

        @staticmethod
        def schema_extra(
            schema: Dict[str, Any], model: Type["HgvsToCopyNumberCountService"]
        ) -> None:
            """Configure OpenAPI schema."""
            if "title" in schema.keys():
                schema.pop("title", None)
            for prop in schema.get("properties", {}).values():
                prop.pop("title", None)
            schema["example"] = {
                "hgvs_expr": "NC_000003.12:g.49531262dup",
                "copy_number_count": {
                    "_id": "ga4gh:CN.lxbM1jOtrVgrwy_SHSSd3o2QkCRRswyf",
                    "type": "CopyNumberCount",
                    "subject": {
                        "_id": "ga4gh:VSL.0dgeuVKngTm5HWjNjcZ9PO-fnbNmKmBv",
                        "type": "SequenceLocation",
                        "sequence_id": "ga4gh:SQ.Zu7h9AggXxhTaGVsy7h_EZSChSZGcmgX",
                        "interval": {
                            "type": "SequenceInterval",
                            "start": {"type": "Number", "value": 49531261},
                            "end": {"type": "Number", "value": 49531262},
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


class HgvsToCopyNumberChangeService(ServiceResponse):
    """A response for translating HGVS to copy number change."""

    hgvs_expr: StrictStr
    copy_number_change: Optional[Union[CopyNumberChange, Text]]

    class Config:
        """Configure model."""

        @staticmethod
        def schema_extra(
            schema: Dict[str, Any], model: Type["HgvsToCopyNumberChangeService"]
        ) -> None:
            """Configure OpenAPI schema."""
            if "title" in schema.keys():
                schema.pop("title", None)
            for prop in schema.get("properties", {}).values():
                prop.pop("title", None)
            schema["example"] = {
                "hgvs_expr": "NC_000003.12:g.49531262dup",
                "copy_number_change": {
                    "_id": "ga4gh:CX.49PTi3fDMxTdYRLp-svfrHrHc_pIAWT3",
                    "type": "CopyNumberChange",
                    "subject": {
                        "_id": "ga4gh:VSL.0dgeuVKngTm5HWjNjcZ9PO-fnbNmKmBv",
                        "type": "SequenceLocation",
                        "sequence_id": "ga4gh:SQ.Zu7h9AggXxhTaGVsy7h_EZSChSZGcmgX",
                        "interval": {
                            "type": "SequenceInterval",
                            "start": {"type": "Number", "value": 49531261},
                            "end": {"type": "Number", "value": 49531262},
                        },
                    },
                    "copy_change": "efo:0030069",
                },
                "service_meta_": {
                    "name": "variation-normalizer",
                    "version": "0.2.17",
                    "response_datetime": "2022-01-26T22:23:41.821673",
                    "url": "https://github.com/cancervariants/variation-normalization",
                },
            }
