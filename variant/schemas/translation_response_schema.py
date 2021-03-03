"""Module for Translation Response Schema."""
from pydantic import BaseModel
from typing import List, Dict, Type, Any
from variant.schemas.ga4gh_vrs import Allele


class TranslationResponseSchema(BaseModel):
    """Define model for translation response."""

    search_term: str
    variants: List[Allele]

    class Config:
        """Configure model."""

        @staticmethod
        def schema_extra(schema: Dict[str, Any],
                         model: Type['TranslationResponseSchema']) -> None:
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
                ]
            }
