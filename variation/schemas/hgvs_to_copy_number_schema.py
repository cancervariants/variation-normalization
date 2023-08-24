"""Module containing schemas used in HGVS To Copy Number endpoints"""
from typing import Optional

from ga4gh.vrs import models
from pydantic import ConfigDict, StrictStr

from variation.schemas.normalize_response_schema import ServiceResponse
from variation.version import __version__


class HgvsToCopyNumberCountService(ServiceResponse):
    """A response for translating HGVS to copy number count."""

    hgvs_expr: StrictStr
    copy_number_count: Optional[models.CopyNumberCount] = None

    # FIXME: This might change due to vrs-python normalize
    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "hgvs_expr": "NC_000003.12:g.49531262dup",
                "copy_number_count": {
                    "id": "ga4gh:CN.5Ft6N62SSDslmYFdyaurhWYfSlsfidUl",
                    "type": "CopyNumberCount",
                    "subject": {
                        "id": "ga4gh:SL.LC-4oEU03VAX1bJD4pQii8tEYl9UbiBw",
                        "type": "SequenceLocation",
                        "sequence": "ga4gh:SQ.Zu7h9AggXxhTaGVsy7h_EZSChSZGcmgX",
                        "start": 49531260,
                        "end": 49531262,
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


class HgvsToCopyNumberChangeService(ServiceResponse):
    """A response for translating HGVS to copy number change."""

    hgvs_expr: StrictStr
    copy_number_change: Optional[models.CopyNumberChange] = None

    # FIXME: This might change due to vrs-python normalize
    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "hgvs_expr": "NC_000003.12:g.49531262dup",
                "copy_number_change": {
                    "id": "ga4gh:CX.rLoNInslL4_U4MEze31LHZ6hqMBeRo52",
                    "type": "CopyNumberChange",
                    "subject": {
                        "id": "ga4gh:SL.LC-4oEU03VAX1bJD4pQii8tEYl9UbiBw",
                        "type": "SequenceLocation",
                        "sequence": "ga4gh:SQ.Zu7h9AggXxhTaGVsy7h_EZSChSZGcmgX",
                        "start": 49531260,
                        "end": 49531262,
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
