"""Module for normalize endpoint response schema."""
from pydantic import BaseModel
from variant.schemas.ga4gh_vod import VariationDescriptor
from typing import List, Optional, Dict, Any, Type
from datetime import datetime


class ServiceMeta(BaseModel):
    """Metadata regarding the variant-normalization service."""

    name = 'variant-normalizer'
    version: str
    response_datetime: datetime
    url = 'https://github.com/cancervariants/variant-normalization'

    class Config:
        """Configure schema example."""

        @staticmethod
        def schema_extra(schema: Dict[str, Any],
                         model: Type['ServiceMeta']) -> None:
            """Configure OpenAPI schema"""
            if 'title' in schema.keys():
                schema.pop('title', None)
            for prop in schema.get('properties', {}).values():
                prop.pop('title', None)
            schema['example'] = {
                'name': 'variant-normalizer',
                'version': '0.1.0',
                'response_datetime': '2021-04-05T16:44:15.367831',
                'url': 'https://github.com/cancervariants/variant-normalization'  # noqa: E501
            }


class NormalizeService(BaseModel):
    """A response to normalizing a variant to a single GA4GH Value Object Descriptor."""  # noqa: E501

    variant_query: str
    warnings: Optional[List[str]]
    variation_descriptor: Optional[VariationDescriptor]
    service_meta_: ServiceMeta

    class Config:
        """Configure model."""

        @staticmethod
        def schema_extra(schema: Dict[str, Any],
                         model: Type['NormalizeService']) -> None:
            """Configure OpenAPI schema."""
            if 'title' in schema.keys():
                schema.pop('title', None)
            for prop in schema.get('properties', {}).values():
                prop.pop('title', None)
            schema['example'] = {
                "variant_query": "BRAF V600E",
                "variation_descriptor": {
                    "id": "normalize.variant:BRAF%20V600E",
                    "type": "VariationDescriptor",
                    "label": "NP_001361187.1:p.Val640Glu",
                    "value_id": "ga4gh:VA.9dA0egRAIfVFDL1sdU1VP7HsBcG0-DtE",
                    "value": {
                        "location": {
                            "interval": {
                                "end": 640,
                                "start": 639,
                                "type": "SimpleInterval"
                            },
                            "sequence_id": "ga4gh:SQ.WaAJ_cXXn9YpMNfhcq9lnzIvaB9ALawo",  # noqa: E501
                            "type": "SequenceLocation"
                        },
                        "state": {
                            "sequence": "E",
                            "type": "SequenceState"
                        },
                        "type": "Allele"
                    },
                    "molecule_context": "protein",
                    "structural_type": "SO:0001606",
                    "ref_allele_seq": "V",
                    "gene_context": {
                        "id": "normalize.gene:BRAF",
                        "type": "GeneDescriptor",
                        "label": "BRAF",
                        "value": {
                            "id": "hgnc:1097",
                            "type": "Gene"
                        },
                        "xrefs": [
                            "ncbigene:673",
                            "ensembl:ENSG00000157764"
                        ],
                        "alternate_labels": [
                            "B-Raf proto-oncogene, serine/threonine kinase",
                            "BRAF1"
                        ],
                        "extensions": [
                            {
                                "name": "symbol_status",
                                "value": "approved",
                                "type": "Extension"
                            },
                            {
                                "name": "associated_with",
                                "value": [
                                    "vega:OTTHUMG00000157457",
                                    "ucsc:uc003vwc.5",
                                    "ccds:CCDS5863",
                                    "ccds:CCDS87555",
                                    "uniprot:P15056",
                                    "pubmed:2284096",
                                    "pubmed:1565476",
                                    "cosmic:BRAF",
                                    "omim:164757",
                                    "orphanet:119066",
                                    "iuphar:1943",
                                    "ena.embl:M95712",
                                    "refseq:NM_004333"
                                ],
                                "type": "Extension"
                            },
                            {
                                "name": "chromosome_location",
                                "value": {
                                    "_id": "ga4gh:VCL.O6yCQ1cnThOrTfK9YUgMlTfM6HTqbrKw",  # noqa: E501
                                    "type": "ChromosomeLocation",
                                    "species_id": "taxonomy:9606",
                                    "chr": "7",
                                    "interval": {
                                        "end": "q34",
                                        "start": "q34",
                                        "type": "CytobandInterval"
                                    }
                                },
                                "type": "Extension"
                            }
                        ]
                    }
                },
                "service_meta_": {
                    'name': 'variant-normalizer',
                    'version': '0.1.0',
                    'response_datetime': '2021-04-05T16:44:15.367831',
                    'url': 'https://github.com/cancervariants/variant-normalization'  # noqa: E501
                }
            }
