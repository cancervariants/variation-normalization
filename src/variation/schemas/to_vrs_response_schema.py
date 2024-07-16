"""Module for to_vrs endpoint response schema."""

from ga4gh.vrs import models
from pydantic import BaseModel, ConfigDict, StrictStr

from variation import __version__
from variation.schemas.normalize_response_schema import ServiceMeta


class ToVRSService(BaseModel):
    """Define model for translation response."""

    search_term: StrictStr
    warnings: list[StrictStr] = []
    variations: (
        list[models.Allele]
        | list[models.CopyNumberCount]
        | list[models.CopyNumberChange]
    ) = []
    service_meta_: ServiceMeta

    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "search_term": "BRAF V600E",
                "warnings": [],
                "variations": [
                    {
                        "id": "ga4gh:VA.GGJOybg6mckctDlXxLY1kHZQ6dbB0U75",
                        "digest": "GGJOybg6mckctDlXxLY1kHZQ6dbB0U75",
                        "location": {
                            "id": "ga4gh:SL.0Y2ZW1zB9mf0qGesJx1kDYtwfB-67Gnf",
                            "digest": "0Y2ZW1zB9mf0qGesJx1kDYtwfB-67Gnf",
                            "end": 600,
                            "start": 599,
                            "sequenceReference": {
                                "type": "SequenceReference",
                                "refgetAccession": "SQ.ZJwurRo2HLY018wghYjDKSfIlEH0Y8At",
                            },
                            "type": "SequenceLocation",
                        },
                        "state": {"sequence": "E", "type": "LiteralSequenceExpression"},
                        "type": "Allele",
                    },
                    {
                        "id": "ga4gh:VA.j4XnsLZcdzDIYa5pvvXM7t1wn9OITr0L",
                        "digest": "j4XnsLZcdzDIYa5pvvXM7t1wn9OITr0L",
                        "location": {
                            "id": "ga4gh:SL.t-3DrWALhgLdXHsupI-e-M00aL3HgK3y",
                            "digest": "t-3DrWALhgLdXHsupI-e-M00aL3HgK3y",
                            "end": 600,
                            "start": 599,
                            "sequenceReference": {
                                "type": "SequenceReference",
                                "refgetAccession": "SQ.cQvw4UsHHRRlogxbWCB8W-mKD4AraM9y",
                            },
                            "type": "SequenceLocation",
                        },
                        "state": {"sequence": "E", "type": "LiteralSequenceExpression"},
                        "type": "Allele",
                    },
                    {
                        "id": "ga4gh:VA.jSy0uhhLefGH3396djPYcJeSVyDvRYGc",
                        "digest": "jSy0uhhLefGH3396djPYcJeSVyDvRYGc",
                        "location": {
                            "id": "ga4gh:SL.CowRxWqyJfqVwlfs5YdswnrnZDQL5QCi",
                            "digest": "CowRxWqyJfqVwlfs5YdswnrnZDQL5QCi",
                            "start": 599,
                            "end": 600,
                            "sequenceReference": {
                                "type": "SequenceReference",
                                "refgetAccession": "SQ.lKdPZpuT-VNvRuKDjsUItNgutfWYgWQd",
                            },
                            "type": "SequenceLocation",
                        },
                        "state": {"sequence": "E", "type": "LiteralSequenceExpression"},
                        "type": "Allele",
                    },
                    {
                        "id": "ga4gh:VA.T1wKHKsXNF6gm7BJJ9kkyCgYrdzeG5Eh",
                        "digest": "T1wKHKsXNF6gm7BJJ9kkyCgYrdzeG5Eh",
                        "location": {
                            "id": "ga4gh:SL.5_mc1jrXsfhj-Cp8iFEEtvgHWXRtPerb",
                            "digest": "5_mc1jrXsfhj-Cp8iFEEtvgHWXRtPerb",
                            "start": 599,
                            "end": 600,
                            "sequenceReference": {
                                "type": "SequenceReference",
                                "refgetAccession": "SQ.0Q-SgJX1V3seUUIu3qVUtEa55CQsGmEU",
                            },
                            "type": "SequenceLocation",
                        },
                        "state": {"sequence": "E", "type": "LiteralSequenceExpression"},
                        "type": "Allele",
                    },
                ],
                "service_meta_": {
                    "name": "variation-normalizer",
                    "version": __version__,
                    "response_datetime": "2023-08-24T09:05:03.622667",
                    "url": "https://github.com/cancervariants/variation-normalization",
                },
            }
        }
    )
