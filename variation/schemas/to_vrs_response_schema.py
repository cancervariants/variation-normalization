"""Module for to_vrs endpoint response schema."""
from typing import List, Dict, Type, Any, Optional, Union

from pydantic import BaseModel
from pydantic.types import StrictStr
from ga4gh.vrsatile.pydantic.vrs_models import Allele, Text, Haplotype, \
    CopyNumberCount, VariationSet, CopyNumberChange

from variation.schemas.normalize_response_schema import ServiceMeta


class ToVRSService(BaseModel):
    """Define model for translation response."""

    search_term: StrictStr
    warnings: Optional[List[StrictStr]]
    variations: Optional[Union[List[Allele], List[Text], List[Haplotype],
                               List[CopyNumberCount], List[CopyNumberChange],
                               List[VariationSet]]]
    service_meta_: ServiceMeta

    class Config:
        """Configure model."""

        @staticmethod
        def schema_extra(schema: Dict[str, Any],
                         model: Type["ToVRSService"]) -> None:
            """Configure OpenAPI schema."""
            if "title" in schema.keys():
                schema.pop("title", None)
            for prop in schema.get("properties", {}).values():
                prop.pop("title", None)
            schema["example"] = {
                "search_term": "BRAF V600E",
                "variations": [
                    {
                        "id": "ga4gh:VA.Lado-KeGuXLnOX2-j920r67p9Z-3TCzb",
                        "type": "Allele",
                        "location": {
                            "id": "ga4gh:SL.x5j9BVC93gA1W7YrmMzFIP1Iyk4ipaal",
                            "type": "SequenceLocation",
                            "sequence_id": "ga4gh:SQ.ZJwurRo2HLY018wghYjDKSfIlEH0Y8At",
                            "start": {"type": "Number", "value": 599},
                            "end": {"type": "Number", "value": 600}
                        },
                        "state": {
                            "type": "LiteralSequenceExpression",
                            "sequence": "E"
                        }
                    },
                    {
                        "id": "ga4gh:VA.YZ6lbo8JEu_JdGCQYZ2_emQUBZjCxA4z",
                        "type": "Allele",
                        "location": {
                            "id": "ga4gh:SL.cxWCch3FsZlYpW3xpPknj4WltsZedn4B",
                            "type": "SequenceLocation",
                            "sequence_id": "ga4gh:SQ.lKdPZpuT-VNvRuKDjsUItNgutfWYgWQd",
                            "start": {"type": "Number", "value": 599},
                            "end": {"type": "Number", "value": 600}
                        },
                        "state": {
                            "type": "LiteralSequenceExpression",
                            "sequence": "E"
                        }
                    },
                    {
                        "id": "ga4gh:VA.Ej4bIR7tJHGvVHrWkYTC0Jv7AlGfSnyc",
                        "type": "Allele",
                        "location": {
                            "id": "ga4gh:SL.3ZYFXkQ-98YRuh0PQ1E42W2Ex5J6DMHS",
                            "type": "SequenceLocation",
                            "sequence_id": "ga4gh:SQ.0Q-SgJX1V3seUUIu3qVUtEa55CQsGmEU",
                            "start": {"type": "Number", "value": 599},
                            "end": {"type": "Number", "value": 600}
                        },
                        "state": {
                            "type": "LiteralSequenceExpression",
                            "sequence": "E"
                        }
                    },
                    {
                        "id": "ga4gh:VA.h313H4CQh6pogbbSJ3H5pI1cPoh9YMm_",
                        "type": "Allele",
                        "location": {
                            "id": "ga4gh:SL.xfBTztcmMstx8jrrdgPiE_BUoLHLFMMS",
                            "type": "SequenceLocation",
                            "sequence_id": "ga4gh:SQ.cQvw4UsHHRRlogxbWCB8W-mKD4AraM9y",
                            "start": {
                                "type": "Number",
                                "value": 599
                            },
                            "end": {
                                "type": "Number",
                                "value": 600
                            }
                        },
                        "state": {
                            "type": "LiteralSequenceExpression",
                            "sequence": "E"
                        }
                    }
                ],
                "service_meta_": {
                    "name": "variation-normalizer",
                    "version": "0.2.17",
                    "response_datetime": "2022-01-26T22:23:41.821673",
                    "url": "https://github.com/cancervariants/variation-normalization"  # noqa: E501
                }

            }
