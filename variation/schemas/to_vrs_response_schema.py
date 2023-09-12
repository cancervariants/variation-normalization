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
                        "id": "ga4gh:VA.PJu8CCaVzEyqXMAEcMNegyDWyvT_jzNn",
                        "location": {
                            "id": "ga4gh:SL.EpHaD2ygDuPMvyURI9L4yetEwF3W0G7G",
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
                        "id": "ga4gh:VA.4XBXAxSAk-WyAu5H0S1-plrk_SCTW1PO",
                        "location": {
                            "id": "ga4gh:SL.ZA1XNKhCT_7m2UtmnYb8ZYOVS4eplMEK",
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
                        "id": "ga4gh:VA.c-oRhbu7nDrBrSW2fPbFlDM15V6jiaho",
                        "location": {
                            "id": "ga4gh:SL.gkevJbLNOScKXhxhzOZXiG3hW8zeyo-q",
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
                        "id": "ga4gh:VA.3ex0cvKXjHbq8NLuitOAfVwSPzqZUFrR",
                        "location": {
                            "id": "ga4gh:SL.Q4MXez2kHFPQqGJKLP8quVHAskuCrOAA",
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
