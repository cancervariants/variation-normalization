"""Module for Translation Response Schema."""
from pydantic import BaseModel
from typing import List, Dict, Type, Any, Optional
from variant.schemas.ga4gh_vrs import Allele
from variant.schemas.normalize_response_schema import ServiceMeta


class ToVRSService(BaseModel):
    """Define model for translation response."""

    search_term: str
    warnings: Optional[List[str]]
    variants: List[Allele]
    service_meta_: ServiceMeta

    class Config:
        """Configure model."""

        @staticmethod
        def schema_extra(schema: Dict[str, Any],
                         model: Type['ToVRSService']) -> None:
            """Configure OpenAPI schema."""
            if 'title' in schema.keys():
                schema.pop('title', None)
            for prop in schema.get('properties', {}).values():
                prop.pop('title', None)
            schema['example'] = {
                "search_term": "BRAF V600E",
                "variants": [
                    {
                        "_id": "ga4gh:VA.u6sKlz0mMQvARmrlnt0Aksz6EbSkmL8z",
                        "location": {
                            "interval": {
                                "start": 599,
                                "end": 600,
                                "type": "SimpleInterval"
                            },
                            "sequence_id": "ga4gh:SQ.ZJwurRo2HLY018wghYjDKSfIlEH0Y8At",  # noqa: E501
                            "type": "SequenceLocation"
                        },
                        "state": {
                            "sequence": "E",
                            "type": "SequenceState"
                        },
                        "type": "Allele"
                    }
                ],
                "service_meta_": {
                    'name': 'variant-normalizer',
                    'version': '0.1.0',
                    'response_datetime': '2021-04-05T16:44:15.367831',
                    'url': 'https://github.com/cancervariants/variant-normalization'  # noqa: E501
                }

            }
