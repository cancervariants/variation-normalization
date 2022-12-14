"""Module for normalize endpoint response schema."""
from enum import Enum
from typing import List, Optional, Dict, Any, Type
from datetime import datetime

from pydantic import BaseModel
from pydantic.types import StrictStr
from ga4gh.vrsatile.pydantic.vrsatile_models import VariationDescriptor, \
    CanonicalVariation


class HGVSDupDelMode(str, Enum):
    """Define options for HGVSDupDelMode.
    This mode determines how to interpret HGVS dup/del.
    """

    DEFAULT = "default"
    ABSOLUTE_CNV = "absolute_cnv"
    RELATIVE_CNV = "relative_cnv"
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
        def schema_extra(schema: Dict[str, Any],
                         model: Type["ServiceMeta"]) -> None:
            """Configure OpenAPI schema"""
            if "title" in schema.keys():
                schema.pop("title", None)
            for prop in schema.get("properties", {}).values():
                prop.pop("title", None)
            schema["example"] = {
                "name": "variation-normalizer",
                "version": "0.1.0",
                "response_datetime": "2021-04-05T16:44:15.367831",
                "url": "https://github.com/cancervariants/variation-normalization"
            }


class ServiceResponse(BaseModel):
    """Base response model for services"""

    warnings: Optional[List[StrictStr]]
    service_meta_: ServiceMeta


class NormalizeService(ServiceResponse):
    """A response to normalizing a variation to a single GA4GH Value Object Descriptor."""  # noqa: E501

    variation_query: StrictStr
    variation_descriptor: Optional[VariationDescriptor]

    class Config:
        """Configure model."""

        @staticmethod
        def schema_extra(schema: Dict[str, Any],
                         model: Type["NormalizeService"]) -> None:
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
                    "variation": {
                        "id": "ga4gh:VA.h313H4CQh6pogbbSJ3H5pI1cPoh9YMm_",
                        "location": {
                            "id": "ga4gh:SL.xfBTztcmMstx8jrrdgPiE_BUoLHLFMMS",
                            "end": {"value": 600, "type": "Number"},
                            "start": {"value": 599, "type": "Number"},
                            "sequence_id": "ga4gh:SQ.cQvw4UsHHRRlogxbWCB8W-mKD4AraM9y",
                            "type": "SequenceLocation"
                        },
                        "state": {
                            "sequence": "E",
                            "type": "LiteralSequenceExpression"
                        },
                        "type": "Allele"
                    },
                    "molecule_context": "protein",
                    "structural_type": "SO:0001606",
                    "vrs_ref_allele_seq": "V",
                    "gene_context": {
                        "id": "normalize.gene:BRAF",
                        "type": "GeneDescriptor",
                        "label": "BRAF",
                        "gene": "hgnc:1097",
                        "xrefs": [
                            "ncbigene:673",
                            "ensembl:ENSG00000157764"
                        ],
                        "alternate_labels": [
                            "BRAF1",
                            "RAFB1",
                            "B-raf",
                            "B-RAF1",
                            "NS7",
                            "BRAF-1"
                        ],
                        "extensions": [
                            {
                                "name": "symbol_status",
                                "value": "approved",
                                "type": "Extension"
                            },
                            {
                                "type": "Extension",
                                "name": "approved_name",
                                "value": "B-Raf proto-oncogene, serine/threonine kinase"
                            },
                            {
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
                                    "cosmic:BRAF"
                                ],
                                "type": "Extension"
                            },
                            {
                                "name": "chromosome_location",
                                "value": {
                                    "id": "ga4gh:CL.ZZZYpOwuW1BLLJXc_Dm4eVZ5E0smVYCc",
                                    "type": "ChromosomeLocation",
                                    "species_id": "taxonomy:9606",
                                    "chr": "7",
                                    "end": "q34",
                                    "start": "q34",
                                },
                                "type": "Extension"
                            },
                            {
                                "type": "Extension",
                                "name": "hgnc_locus_type",
                                "value": "gene with protein product"
                            },
                            {
                                "type": "Extension",
                                "name": "ncbi_gene_type",
                                "value": "protein-coding"
                            },
                            {
                                "type": "Extension",
                                "name": "ensembl_biotype",
                                "value": "protein_coding"
                            }
                        ]
                    }
                },
                "service_meta_": {
                    "name": "variation-normalizer",
                    "version": "0.2.17",
                    "response_datetime": "2022-01-26T22:23:41.821673",
                    "url": "https://github.com/cancervariants/variation-normalization"
                }
            }


class TranslateIdentifierService(ServiceResponse):
    """A response to translating identifiers."""

    identifier_query: StrictStr
    aliases: List[StrictStr]

    class Config:
        """Configure model."""

        @staticmethod
        def schema_extra(schema: Dict[str, Any],
                         model: Type["TranslateIdentifierService"]) -> None:
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
                    "ga4gh:SQ.cQvw4UsHHRRlogxbWCB8W-mKD4AraM9y"
                ],
                "service_meta_": {
                    "name": "variation-normalizer",
                    "version": "0.2.14",
                    "response_datetime": "2021-11-18T14:10:53.909158",
                    "url": "https://github.com/cancervariants/variation-normalization"
                }
            }


class ToCanonicalVariationFmt(str, Enum):
    """Define formats for to_canonical endpoint"""

    HGVS = "hgvs"
    SPDI = "spdi"


class ToCanonicalVariationService(ServiceResponse):
    """A response model for the to canonical variation service"""

    query: str
    canonical_variation: Optional[CanonicalVariation]

    class Config:
        """Configure model."""

        @staticmethod
        def schema_extra(
                schema: Dict[str, Any],
                model: Type["ToCanonicalVariationService"]) -> None:
            """Configure OpenAPI schema."""
            if "title" in schema.keys():
                schema.pop("title", None)
            for prop in schema.get("properties", {}).values():
                prop.pop("title", None)
            schema["example"] = {
                "query": "NC_000007.14:140753335:A:T",
                "warnings": [],
                "canonical_variation": {
                    "id": "ga4gh:CAN.dP6z4p7SoGJFmlFQcjOQo2d1mXuo1QiY",
                    "type": "CanonicalVariation",
                    "canonical_context": {
                        "id": "ga4gh:VA.3xFHF399HGbG1JUf5uwcj3oWVKZJ70oX",
                        "type": "Allele",
                        "location": {
                            "id": "ga4gh:SL.WBLxdkoypnRME6b8tJtlOWqZKU1ruqY1",
                            "type": "SequenceLocation",
                            "sequence_id": "ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
                            "start": {"type": "Number", "value": 140753335},
                            "end": {"type": "Number", "value": 140753336}
                        },
                        "state": {
                            "type": "LiteralSequenceExpression",
                            "sequence": "T"
                        }
                    }
                },
                "service_meta_": {
                    "version": "0.2.20",
                    "response_datetime": "2022-02-20T17:16:19.415675",
                    "name": "variation-normalizer",
                    "url": "https://github.com/cancervariants/variation-normalization"
                }
            }
