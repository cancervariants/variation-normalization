"""Module for to_vrs endpoint response schema."""
from typing import List, Optional, Union

from ga4gh.vrsatile.pydantic.vrs_models import (
    Allele,
    CopyNumberChange,
    CopyNumberCount,
    Haplotype,
    Text,
    VariationSet,
)
from pydantic import BaseModel, ConfigDict, StrictStr

from variation.schemas.normalize_response_schema import ServiceMeta


class ToVRSService(BaseModel):
    """Define model for translation response."""

    search_term: StrictStr
    warnings: Optional[List[StrictStr]] = None
    variations: Optional[
        Union[
            List[Allele],
            List[Text],
            List[Haplotype],
            List[CopyNumberCount],
            List[CopyNumberChange],
            List[VariationSet],
        ]
    ] = None
    service_meta_: ServiceMeta

    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "search_term": "BRAF V600E",
                "variations": [
                    {
                        "_id": "ga4gh:VA.7ys8TiDzrk04O3Upd63__rOBCEhv3P5d",
                        "type": "Allele",
                        "location": {
                            "_id": "ga4gh:VSL.Vxqx2bv42rWeu08Eg7JpkdQkMCNLskoz",
                            "type": "SequenceLocation",
                            "sequence_id": "ga4gh:SQ.ZJwurRo2HLY018wghYjDKSfIlEH0Y8At",
                            "interval": {
                                "type": "SequenceInterval",
                                "start": {"type": "Number", "value": 599},
                                "end": {"type": "Number", "value": 600},
                            },
                        },
                        "state": {"type": "LiteralSequenceExpression", "sequence": "E"},
                    },
                    {
                        "_id": "ga4gh:VA.vimwyw0pFTwatfFhi3rhhb153ARWsPrW",
                        "type": "Allele",
                        "location": {
                            "_id": "ga4gh:VSL.FVmsWpfSOA3B2ryq0k995oHMuSGiFvMa",
                            "type": "SequenceLocation",
                            "sequence_id": "ga4gh:SQ.lKdPZpuT-VNvRuKDjsUItNgutfWYgWQd",
                            "interval": {
                                "type": "SequenceInterval",
                                "start": {"type": "Number", "value": 599},
                                "end": {"type": "Number", "value": 600},
                            },
                        },
                        "state": {"type": "LiteralSequenceExpression", "sequence": "E"},
                    },
                    {
                        "_id": "ga4gh:VA.FzlrH5feNcQ3S9GayMU9EF008j-8Pbz5",
                        "type": "Allele",
                        "location": {
                            "_id": "ga4gh:VSL.QDLST2nKpPWwIArdO57L2VIWPNZ0DiN3",
                            "type": "SequenceLocation",
                            "sequence_id": "ga4gh:SQ.0Q-SgJX1V3seUUIu3qVUtEa55CQsGmEU",
                            "interval": {
                                "type": "SequenceInterval",
                                "start": {"type": "Number", "value": 599},
                                "end": {"type": "Number", "value": 600},
                            },
                        },
                        "state": {"type": "LiteralSequenceExpression", "sequence": "E"},
                    },
                    {
                        "_id": "ga4gh:VA.ZDdoQdURgO2Daj2NxLj4pcDnjiiAsfbO",
                        "type": "Allele",
                        "location": {
                            "_id": "ga4gh:VSL.2cHIgn7iLKk4x9z3zLkSTTFMV0e48DR4",
                            "type": "SequenceLocation",
                            "sequence_id": "ga4gh:SQ.cQvw4UsHHRRlogxbWCB8W-mKD4AraM9y",
                            "interval": {
                                "type": "SequenceInterval",
                                "start": {"type": "Number", "value": 599},
                                "end": {"type": "Number", "value": 600},
                            },
                        },
                        "state": {"type": "LiteralSequenceExpression", "sequence": "E"},
                    },
                ],
                "service_meta_": {
                    "name": "variation-normalizer",
                    "version": "0.2.17",
                    "response_datetime": "2022-01-26T22:23:41.821673",
                    "url": "https://github.com/cancervariants/variation-normalization",
                },
            }
        }
    )
