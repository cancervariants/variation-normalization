"""Module for modeling VRSATILE objects."""
from pydantic import BaseModel, validator, root_validator
from pydantic.types import StrictStr
from typing import List, Optional, Dict, Any, Type, Union
from enum import Enum
from gene.schemas import GeneDescriptor, Extension, check_curie
from variation.schemas.ga4gh_vrs import Allele, Text
import re


class MoleculeContext(str, Enum):
    """Define constraints for types of molecule context."""

    GENOMIC = 'genomic'
    TRANSCRIPT = 'transcript'
    PROTEIN = 'protein'


class Expression(BaseModel):
    """Enable descriptions based on a specified nomenclature or syntax for representing an object. - GA4GH"""  # noqa: E501

    type = 'Expression'
    syntax: StrictStr
    value: StrictStr
    version: Optional[StrictStr]

    _validate_syntax = validator('syntax', allow_reuse=True)(check_curie)

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


class VariationDescriptor(BaseModel):
    """GA4GH Value Object Descriptor."""

    id: StrictStr
    type = 'VariationDescriptor'
    value_id: Optional[StrictStr]
    value: Optional[Union[Allele, Text]]
    label: Optional[StrictStr]
    extensions: Optional[List[Extension]]
    molecule_context: Optional[MoleculeContext]
    structural_type: Optional[StrictStr]
    expressions: Optional[Expression]
    gene_context: Optional[GeneDescriptor]
    vrs_ref_allele_seq: Optional[StrictStr]

    _validate_id = validator('id', allow_reuse=True)(check_curie)
    _validate_value_id = validator('value_id', allow_reuse=True)(check_curie)
    _validate_structural_type = \
        validator('structural_type', allow_reuse=True)(check_curie)

    @root_validator(pre=True)
    def check_value_or_value_id_present(cls, values):
        """Check that at least one of {`value`, `value_id`} is provided."""
        msg = 'Must give values for either `value`, `value_id`, or both'
        value, value_id = values.get('value'), values.get('value_id')
        assert value or value_id, msg
        return values

    @validator('vrs_ref_allele_seq')
    def check_vrs_ref_allele_seq(cls, v):
        """Validate vrs_ref_allele_seq"""
        if v is not None:
            assert bool(re.match(r"^[A-Z*\-]*$", v))
        return v

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
                "vrs_ref_allele_seq": "V",
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
