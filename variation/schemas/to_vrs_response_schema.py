"""Module for to_vrs endpoint response schema."""
from typing import List, Union

from ga4gh.vrs import models
from pydantic import BaseModel, ConfigDict, StrictStr

from variation.schemas.normalize_response_schema import ServiceMeta
from variation.version import __version__


class ToVRSService(BaseModel):
    """Define model for translation response."""

    search_term: StrictStr
    warnings: List[StrictStr] = []
    variations: Union[
        List[models.Allele],
        List[models.CopyNumberCount],
        List[models.CopyNumberChange],
    ] = []
    service_meta_: ServiceMeta

    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "search_term": "BRAF V600E",
                "warnings": [],
                "variations": [
                    {
                        "id": "ga4gh:VA.kY8XpXZ4ueOQoPbYxx5nvj0zGXvsmZg7",
                        "type": "Allele",
                        "location": {
                            "id": "ga4gh:SL.lekHz-yBPRgkcIWIYoYum5xNkdQ-H9D-",
                            "type": "SequenceLocation",
                            "sequenceReference": {
                                "type": "SequenceReference",
                                "refgetAccession": "SQ.0Q-SgJX1V3seUUIu3qVUtEa55CQsGmEU",
                            },
                            "start": 599,
                            "end": 600,
                        },
                        "state": {"type": "LiteralSequenceExpression", "sequence": "E"},
                    },
                    {
                        "id": "ga4gh:VA.s79REj1tC9sM5mhsc1m_yuXcJF3SH3su",
                        "type": "Allele",
                        "location": {
                            "id": "ga4gh:SL.pFtHSvBcQqv2jZYMg5BkXwGNBRpL_kvQ",
                            "type": "SequenceLocation",
                            "sequenceReference": {
                                "type": "SequenceReference",
                                "refgetAccession": "SQ.ZJwurRo2HLY018wghYjDKSfIlEH0Y8At",
                            },
                            "start": 599,
                            "end": 600,
                        },
                        "state": {"type": "LiteralSequenceExpression", "sequence": "E"},
                    },
                    {
                        "id": "ga4gh:VA.tJ7kARp1YlxgcUesm7DXSj_SDXjXg3-u",
                        "type": "Allele",
                        "location": {
                            "id": "ga4gh:SL.pqBCJoUGtDNLKib8F5ke3cyqqrkZzfud",
                            "type": "SequenceLocation",
                            "sequenceReference": {
                                "type": "SequenceReference",
                                "refgetAccession": "SQ.cQvw4UsHHRRlogxbWCB8W-mKD4AraM9y",
                            },
                            "start": 599,
                            "end": 600,
                        },
                        "state": {"type": "LiteralSequenceExpression", "sequence": "E"},
                    },
                    {
                        "id": "ga4gh:VA.KpFlKus-2m141vUNdFlOC1ejisAo4CTv",
                        "type": "Allele",
                        "location": {
                            "id": "ga4gh:SL.mtXl9xjWA8YmN-3LyCeJtvZtNznAIhWo",
                            "type": "SequenceLocation",
                            "sequenceReference": {
                                "type": "SequenceReference",
                                "refgetAccession": "SQ.lKdPZpuT-VNvRuKDjsUItNgutfWYgWQd",
                            },
                            "start": 599,
                            "end": 600,
                        },
                        "state": {"type": "LiteralSequenceExpression", "sequence": "E"},
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
