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

    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "hgvs_expr": "NC_000003.12:g.49531262dup",
                "copy_number_count": {
                    "id": "ga4gh:CN.07iM14yvZ80N_AiaM7G_V4f1pCkmFYz4",
                    "type": "CopyNumberCount",
                    "location": {
                        "id": "ga4gh:SL.y4-cVA2VxMCDxb9gV2oFrzC386yrEVqh",
                        "type": "SequenceLocation",
                        "sequenceReference": {
                            "type": "SequenceReference",
                            "refgetAccession": "SQ.Zu7h9AggXxhTaGVsy7h_EZSChSZGcmgX",
                        },
                        "start": 49531261,
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

    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "hgvs_expr": "NC_000003.12:g.49531262dup",
                "copy_number_change": {
                    "id": "ga4gh:CX.d8BWSLNKN0K4n8ySG0jWPCr4cJIqEf5g",
                    "type": "CopyNumberChange",
                    "location": {
                        "id": "ga4gh:SL.y4-cVA2VxMCDxb9gV2oFrzC386yrEVqh",
                        "type": "SequenceLocation",
                        "sequenceReference": {
                            "type": "SequenceReference",
                            "refgetAccession": "SQ.Zu7h9AggXxhTaGVsy7h_EZSChSZGcmgX",
                        },
                        "start": 49531261,
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
