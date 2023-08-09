"""Module for normalize endpoint response schema."""
from datetime import datetime
from enum import Enum
from typing import Any, Dict, List, Optional, Type

from ga4gh.vrsatile.pydantic.vrsatile_models import VariationDescriptor
from pydantic import BaseModel, root_validator
from pydantic.types import StrictStr


class HGVSDupDelModeOption(str, Enum):
    """Define options for HGVSDupDelMode.
    This mode determines how to interpret HGVS dup/del.
    """

    DEFAULT = "default"
    COPY_NUMBER_COUNT = "copy_number_count"
    COPY_NUMBER_CHANGE = "copy_number_change"
    REPEATED_SEQ_EXPR = "repeated_seq_expr"  # VRS Allele
    LITERAL_SEQ_EXPR = "literal_seq_expr"  # VRS Allele


class ServiceMeta(BaseModel):
    """Metadata regarding the variation-normalization service."""

    name = "variation-normalizer"
    version: StrictStr
    response_datetime: datetime
    url = "https://github.com/cancervariants/variation-normalization"

    class Config:
        """Configure schema example."""

        @staticmethod
        def schema_extra(schema: Dict[str, Any], model: Type["ServiceMeta"]) -> None:
            """Configure OpenAPI schema"""
            if "title" in schema.keys():
                schema.pop("title", None)
            for prop in schema.get("properties", {}).values():
                prop.pop("title", None)
            schema["example"] = {
                "name": "variation-normalizer",
                "version": "0.1.0",
                "response_datetime": "2021-04-05T16:44:15.367831",
                "url": "https://github.com/cancervariants/variation-normalization",
            }


class ServiceResponse(BaseModel):
    """Base response model for services"""

    warnings: Optional[List[StrictStr]] = []
    service_meta_: ServiceMeta

    @root_validator(pre=False)
    def unique_warnings(cls, values):
        """Ensure unique warnings"""
        values["warnings"] = list(set(values["warnings"]))
        return values


class NormalizeService(ServiceResponse):
    """A response to normalizing a variation to a single GA4GH Value Object
    Descriptor.
    """

    variation_query: StrictStr
    variation_descriptor: Optional[VariationDescriptor]

    class Config:
        """Configure model."""

        @staticmethod
        def schema_extra(
            schema: Dict[str, Any], model: Type["NormalizeService"]
        ) -> None:
            """Configure OpenAPI schema."""
            if "title" in schema.keys():
                schema.pop("title", None)
            for prop in schema.get("properties", {}).values():
                prop.pop("title", None)
            schema["example"] = {
                "variation_query": "BRAF V600E",
                "variation_descriptor": {
                    "id": "normalize.variation:BRAF%20V600E",
                    "label": "BRAF V600E",
                    "type": "VariationDescriptor",
                    "variation_id": "ga4gh:VA.8JkgnqIgYqufNl-OV_hpRG_aWF9UFQCE",
                    "variation": {
                        "_id": "ga4gh:VA.8JkgnqIgYqufNl-OV_hpRG_aWF9UFQCE",
                        "location": {
                            "_id": "ga4gh:VSL.AqrQ-EkAvTrXOFn70_8i3dXF5shBBZ5i",
                            "interval": {
                                "end": {"value": 640, "type": "Number"},
                                "start": {"value": 639, "type": "Number"},
                                "type": "SequenceInterval",
                            },
                            "sequence_id": "ga4gh:SQ.WaAJ_cXXn9YpMNfhcq9lnzIvaB9ALawo",
                            "type": "SequenceLocation",
                        },
                        "state": {"sequence": "E", "type": "LiteralSequenceExpression"},
                        "type": "Allele",
                    },
                    "molecule_context": "protein",
                    "vrs_ref_allele_seq": "V",
                    "gene_context": {
                        "id": "normalize.gene:BRAF",
                        "type": "GeneDescriptor",
                        "label": "BRAF",
                        "gene_id": "hgnc:1097",
                        "xrefs": ["ncbigene:673", "ensembl:ENSG00000157764"],
                        "alternate_labels": [
                            "BRAF1",
                            "RAFB1",
                            "B-raf",
                            "B-RAF1",
                            "NS7",
                            "BRAF-1",
                        ],
                        "extensions": [
                            {
                                "type": "Extension",
                                "name": "symbol_status",
                                "value": "approved",
                            },
                            {
                                "type": "Extension",
                                "name": "approved_name",
                                "value": "B-Raf proto-oncogene, serine/threonine kinase",  # noqa: E501
                            },
                            {
                                "type": "Extension",
                                "name": "associated_with",
                                "value": [
                                    "ucsc:uc003vwc.5",
                                    "pubmed:1565476",
                                    "omim:164757",
                                    "vega:OTTHUMG00000157457",
                                    "ccds:CCDS5863",
                                    "iuphar:1943",
                                    "ccds:CCDS87555",
                                    "orphanet:119066",
                                    "refseq:NM_004333",
                                    "ena.embl:M95712",
                                    "pubmed:2284096",
                                    "uniprot:P15056",
                                    "cosmic:BRAF",
                                ],
                            },
                            {
                                "type": "Extension",
                                "name": "hgnc_locations",
                                "value": [
                                    {
                                        "_id": "ga4gh:VCL.O6yCQ1cnThOrTfK9YUgMlTfM6HTqbrKw",  # noqa: E501
                                        "type": "ChromosomeLocation",
                                        "species_id": "taxonomy:9606",
                                        "chr": "7",
                                        "interval": {
                                            "end": "q34",
                                            "start": "q34",
                                            "type": "CytobandInterval",
                                        },
                                    }
                                ],
                            },
                            {
                                "type": "Extension",
                                "name": "ensembl_locations",
                                "value": [
                                    {
                                        "_id": "ga4gh:VSL.amNWL6i7F2nbSZAf2QLTRTujxuDrd0pR",  # noqa: E501
                                        "type": "SequenceLocation",
                                        "sequence_id": "ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",  # noqa: E501
                                        "interval": {
                                            "start": {
                                                "type": "Number",
                                                "value": 140719326,
                                            },
                                            "end": {
                                                "type": "Number",
                                                "value": 140924929,
                                            },
                                            "type": "SequenceInterval",
                                        },
                                    }
                                ],
                            },
                            {
                                "type": "Extension",
                                "name": "ncbi_locations",
                                "value": [
                                    {
                                        "_id": "ga4gh:VCL.O6yCQ1cnThOrTfK9YUgMlTfM6HTqbrKw",  # noqa: E501
                                        "type": "ChromosomeLocation",
                                        "species_id": "taxonomy:9606",
                                        "chr": "7",
                                        "interval": {
                                            "end": "q34",
                                            "start": "q34",
                                            "type": "CytobandInterval",
                                        },
                                    },
                                    {
                                        "_id": "ga4gh:VSL.xZU3kL8F6t2ca6WH_26CWKfNW9-owhR4",  # noqa: E501
                                        "type": "SequenceLocation",
                                        "sequence_id": "ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",  # noqa: E501
                                        "interval": {
                                            "start": {
                                                "type": "Number",
                                                "value": 140713327,
                                            },
                                            "end": {
                                                "type": "Number",
                                                "value": 140924929,
                                            },
                                            "type": "SequenceInterval",
                                        },
                                    },
                                ],
                            },
                            {
                                "type": "Extension",
                                "name": "hgnc_locus_type",
                                "value": "gene with protein product",
                            },
                            {
                                "type": "Extension",
                                "name": "ncbi_gene_type",
                                "value": "protein-coding",
                            },
                            {
                                "type": "Extension",
                                "name": "ensembl_biotype",
                                "value": "protein_coding",
                            },
                        ],
                    },
                },
                "service_meta_": {
                    "name": "variation-normalizer",
                    "version": "0.2.17",
                    "response_datetime": "2022-01-26T22:23:41.821673",
                    "url": "https://github.com/cancervariants/variation-normalization",
                },
            }


class TranslateIdentifierService(ServiceResponse):
    """A response to translating identifiers."""

    identifier_query: StrictStr
    aliases: List[StrictStr]

    class Config:
        """Configure model."""

        @staticmethod
        def schema_extra(
            schema: Dict[str, Any], model: Type["TranslateIdentifierService"]
        ) -> None:
            """Configure OpenAPI schema."""
            if "title" in schema.keys():
                schema.pop("title", None)
            for prop in schema.get("properties", {}).values():
                prop.pop("title", None)
            schema["example"] = {
                "identifier_query": "NP_004324.2",
                "warnings": None,
                "aliases": [
                    "Ensembl:ENSP00000288602.6",
                    "ensembl:ENSP00000288602.6",
                    "Ensembl:ENSP00000493543.1",
                    "ensembl:ENSP00000493543.1",
                    "MD5:74c9b69323bd112084c1b5b385e7e6c5",
                    "NCBI:NP_004324.2",
                    "refseq:NP_004324.2",
                    "SEGUID:sfzILpNpX8UFB/vgH9LOKLpl/+g",
                    "SHA1:b1fcc82e93695fc50507fbe01fd2ce28ba65ffe8",
                    "VMC:GS_cQvw4UsHHRRlogxbWCB8W-mKD4AraM9y",
                    "sha512t24u:cQvw4UsHHRRlogxbWCB8W-mKD4AraM9y",
                    "ga4gh:SQ.cQvw4UsHHRRlogxbWCB8W-mKD4AraM9y",
                ],
                "service_meta_": {
                    "name": "variation-normalizer",
                    "version": "0.2.14",
                    "response_datetime": "2021-11-18T14:10:53.909158",
                    "url": "https://github.com/cancervariants/variation-normalization",
                },
            }
