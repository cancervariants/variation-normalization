"""Module for GA4GH Value Object Descriptor."""
from pydantic import BaseModel
from typing import List, Optional, Dict, Any, Type, Union
from enum import Enum


class Extension(BaseModel):
    """Extend descriptions with other attributes unique to a content provider. -GA4GH"""  # noqa: E501

    type = 'Extension'
    name: str
    value: Union[str, dict, List[str]]

    @staticmethod
    def schema_extra(schema: Dict[str, Any],
                     model: Type['Extension']) -> None:
        """Configure OpenAPI schema."""
        if 'title' in schema.keys():
            schema.pop('title', None)
        for prop in schema.get('properties', {}).values():
            prop.pop('title', None)
        schema['example'] = {
            'type': 'Extension',
            'name': 'symbol_status',
            'value': 'approved'
        }


class ValueObjectDescriptor(BaseModel):
    """GA4GH Value Object Descriptor."""

    id: str
    type: str
    label: Optional[str]
    description: Optional[str]
    value_id: Optional[str]
    value: Optional[dict]
    xrefs: Optional[List[str]]
    alternate_labels: Optional[List[str]]
    extensions: Optional[List[Extension]]


class MoleculeContext(str, Enum):
    """Define constraints for types of molecule context."""

    GENOMIC = 'genomic'
    TRANSCRIPT = 'transcript'
    PROTEIN = 'protein'


class Expression(BaseModel):
    """Enable descriptions based on a specified nomenclature or syntax for representing an object. - GA4GH"""  # noqa: E501

    type = 'Expression'
    syntax: str
    value: str
    version: Optional[str]

    class Config:
        """Configure model."""

        @staticmethod
        def schema_extra(schema: Dict[str, Any],
                         model: Type['Expression']) -> None:
            """Configure OpenAPI schema."""
            if 'title' in schema.keys():
                schema.pop('title', None)
            for prop in schema.get('properties', {}).values():
                prop.pop('title', None)
            schema['example'] = {
                'type': 'Expression',
                'syntax': 'hgvs:protein',
                'value': 'NP_005219.2:p.Leu858Arg'
            }


class GeneDescriptor(ValueObjectDescriptor):
    """Reference GA4GH Gene Value Objects."""

    type = 'GeneDescriptor'

    class Config:
        """Configure model."""

        @staticmethod
        def schema_extra(schema: Dict[str, Any],
                         model: Type['GeneDescriptor']) -> None:
            """Configure OpenAPI schema."""
            if 'title' in schema.keys():
                schema.pop('title', None)
            for prop in schema.get('properties', {}).values():
                prop.pop('title', None)
            schema['example'] = {
                "id": "normalize.gene:BRAF",
                "type": "GeneDescriptor",
                "label": "BRAF",
                "value": {
                    "id": "hgnc:1097",
                    "type": "Gene"
                },
                "xrefs": ["ncbigene:673", "ensembl:ENSG00000157764"],
                "alternate_labels": [
                    "B-Raf proto-oncogene, serine/threonine kinase",
                    "BRAF1"
                ],
                "extensions": [
                    {
                        "type": "Extension",
                        "name": "symbol_status",
                        "value": "approved"
                    },
                    {
                        "type": "Extension",
                        "name": "associated_with",
                        "value": [
                            "vega:OTTHUMG00000157457", "ucsc:uc003vwc.5",
                            "ccds:CCDS5863", "ccds:CCDS87555",
                            "uniprot:P15056",
                            "pubmed:2284096", "pubmed:1565476", "cosmic:BRAF",
                            "omim:164757", "orphanet:119066", "iuphar:1943",
                            "ena.embl:M95712", "refseq:NM_004333"
                        ]
                    },
                    {
                        "type": "Extension",
                        "name": "chromosome_location",
                        "value": {
                            "_id":
                                "ga4gh:VCL.O6yCQ1cnThOrTfK9YUgMlTfM6HTqbrKw",
                            "type": "ChromosomeLocation",
                            "species_id": "taxonomy:9606",
                            "chr": "7",
                            "interval": {
                                "end": "q34",
                                "start": "q34",
                                "type": "CytobandInterval"
                            }
                        }
                    }
                ]
            }


class Gene(BaseModel):
    """GA4GH Gene Value Object."""

    id: str
    type = "Gene"


class VariationDescriptor(ValueObjectDescriptor):
    """Reference GA4GH Variation Value Objects."""

    type = 'VariationDescriptor'
    molecule_context: Optional[MoleculeContext]
    structural_type: Optional[str]
    expressions: Optional[Expression]
    ref_allele_seq: Optional[str]
    gene_context: Optional[GeneDescriptor]

    class Config:
        """Configure model."""

        @staticmethod
        def schema_extra(schema: Dict[str, Any],
                         model: Type['VariationDescriptor']) -> None:
            """Configure OpenAPI schema."""
            if 'title' in schema.keys():
                schema.pop('title', None)
            for prop in schema.get('properties', {}).values():
                prop.pop('title', None)
            schema['example'] = {
                "id": "normalize.variation:BRAF%20v600e",
                "type": "VariationDescriptor",
                "value_id": "ga4gh:VA.u6sKlz0mMQvARmrlnt0Aksz6EbSkmL8z",
                "value": {
                    "location": {
                        "interval": {
                            "end": 600,
                            "start": 599,
                            "type": "SimpleInterval"
                        },
                        "sequence_id":
                            "ga4gh:SQ.ZJwurRo2HLY018wghYjDKSfIlEH0Y8At",
                        "type": "SequenceLocation"
                    },
                    "state": {
                        "sequence": "E",
                        "type": "SequenceState"
                    },
                    "type": "Allele"
                },
                "label": "BRAF V600E",
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
                            "type": "Extension",
                            "name": "symbol_status",
                            "value": "approved"
                        },
                        {
                            "type": "Extension",
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
                            ]
                        },
                        {
                            "type": "Extension",
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
                            }
                        }
                    ]
                }
            }
