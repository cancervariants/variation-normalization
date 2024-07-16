"""Module containing schemas used in HGVS To Copy Number endpoints"""

from ga4gh.vrs import models
from pydantic import ConfigDict, StrictStr

from variation import __version__
from variation.schemas.normalize_response_schema import ServiceResponse


class HgvsToCopyNumberCountService(ServiceResponse):
    """A response for translating HGVS to copy number count."""

    hgvs_expr: StrictStr
    copy_number_count: models.CopyNumberCount | None = None

    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "hgvs_expr": "NC_000003.12:g.49531262dup",
                "copy_number_count": {
                    "id": "ga4gh:CN.gF1l6Zh6aY3vy_TR7rrat6FTmwiwIukY",
                    "digest": "gF1l6Zh6aY3vy_TR7rrat6FTmwiwIukY",
                    "type": "CopyNumberCount",
                    "location": {
                        "id": "ga4gh:SL.2vbgFGHGB0QGODwgZNi05fWbROkkjf04",
                        "digest": "2vbgFGHGB0QGODwgZNi05fWbROkkjf04",
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
    copy_number_change: models.CopyNumberChange | None = None

    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "hgvs_expr": "NC_000003.12:g.49531262dup",
                "copy_number_change": {
                    "id": "ga4gh:CX.Zzws_y4cnoooQ7WXjg2B3nKIyFWXzOg3",
                    "digest": "Zzws_y4cnoooQ7WXjg2B3nKIyFWXzOg3",
                    "type": "CopyNumberChange",
                    "location": {
                        "id": "ga4gh:SL.2vbgFGHGB0QGODwgZNi05fWbROkkjf04",
                        "digest": "2vbgFGHGB0QGODwgZNi05fWbROkkjf04",
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
