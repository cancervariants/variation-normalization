"""Module for toVRS endpoint response schema."""
from pydantic import BaseModel
from typing import List, Dict, Type, Any, Optional, Union
from pydantic.types import StrictStr
from ga4gh.vrsatile.pydantic.vrs_model import Allele, Text, Haplotype, \
    CopyNumber, VariationSet
from variation.schemas.normalize_response_schema import ServiceMeta


class ToVRSService(BaseModel):
    """Define model for translation response."""

    search_term: StrictStr
    warnings: Optional[List[StrictStr]]
    variations: Optional[Union[List[Allele], List[Text], List[Haplotype],
                               List[CopyNumber], List[VariationSet]]] = None
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
                "variations": [
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
                    'name': 'variation-normalizer',
                    'version': '0.1.0',
                    'response_datetime': '2021-04-05T16:44:15.367831',
                    'url': 'https://github.com/cancervariants/variation-normalization'  # noqa: E501
                }

            }
