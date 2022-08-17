"""Module for to_vrs endpoint response schema."""
from typing import List, Dict, Type, Any, Optional, Union

from pydantic import BaseModel
from pydantic.types import StrictStr
from ga4gh.vrsatile.pydantic.vrs_models import Allele, Text, Haplotype, \
    AbsoluteCopyNumber, VariationSet

from variation.schemas.normalize_response_schema import ServiceMeta


class ToVRSService(BaseModel):
    """Define model for translation response."""

    search_term: StrictStr
    warnings: Optional[List[StrictStr]]
    variations: Optional[Union[List[Allele], List[Text], List[Haplotype],
                               List[AbsoluteCopyNumber], List[VariationSet]]] = []
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
                        "_id": "ga4gh:VA.7ys8TiDzrk04O3Upd63__rOBCEhv3P5d",
                        "type": "Allele",
                        "location": {
                            "_id": "ga4gh:VSL.Vxqx2bv42rWeu08Eg7JpkdQkMCNLskoz",  # noqa: E501
                            "type": "SequenceLocation",
                            "sequence_id": "ga4gh:SQ.ZJwurRo2HLY018wghYjDKSfIlEH0Y8At",  # noqa: E501
                            "interval": {
                                "type": "SequenceInterval",
                                "start": {
                                    "type": "Number",
                                    "value": 599
                                },
                                "end": {
                                    "type": "Number",
                                    "value": 600
                                }
                            }
                        },
                        "state": {
                            "type": "LiteralSequenceExpression",
                            "sequence": "E"
                        }
                    },
                    {
                        "_id": "ga4gh:VA.vimwyw0pFTwatfFhi3rhhb153ARWsPrW",
                        "type": "Allele",
                        "location": {
                            "_id": "ga4gh:VSL.FVmsWpfSOA3B2ryq0k995oHMuSGiFvMa",  # noqa: E501
                            "type": "SequenceLocation",
                            "sequence_id": "ga4gh:SQ.lKdPZpuT-VNvRuKDjsUItNgutfWYgWQd",  # noqa: E501
                            "interval": {
                                "type": "SequenceInterval",
                                "start": {
                                    "type": "Number",
                                    "value": 599
                                },
                                "end": {
                                    "type": "Number",
                                    "value": 600
                                }
                            }
                        },
                        "state": {
                            "type": "LiteralSequenceExpression",
                            "sequence": "E"
                        }
                    },
                    {
                        "_id": "ga4gh:VA.FzlrH5feNcQ3S9GayMU9EF008j-8Pbz5",
                        "type": "Allele",
                        "location": {
                            "_id": "ga4gh:VSL.QDLST2nKpPWwIArdO57L2VIWPNZ0DiN3",  # noqa: E501
                            "type": "SequenceLocation",
                            "sequence_id": "ga4gh:SQ.0Q-SgJX1V3seUUIu3qVUtEa55CQsGmEU",  # noqa: E501
                            "interval": {
                                "type": "SequenceInterval",
                                "start": {
                                    "type": "Number",
                                    "value": 599
                                },
                                "end": {
                                    "type": "Number",
                                    "value": 600
                                }
                            }
                        },
                        "state": {
                            "type": "LiteralSequenceExpression",
                            "sequence": "E"
                        }
                    },
                    {
                        "_id": "ga4gh:VA.ZDdoQdURgO2Daj2NxLj4pcDnjiiAsfbO",
                        "type": "Allele",
                        "location": {
                            "_id": "ga4gh:VSL.2cHIgn7iLKk4x9z3zLkSTTFMV0e48DR4",  # noqa: E501
                            "type": "SequenceLocation",
                            "sequence_id": "ga4gh:SQ.cQvw4UsHHRRlogxbWCB8W-mKD4AraM9y",  # noqa: E501
                            "interval": {
                                "type": "SequenceInterval",
                                "start": {
                                    "type": "Number",
                                    "value": 599
                                },
                                "end": {
                                    "type": "Number",
                                    "value": 600
                                }
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
