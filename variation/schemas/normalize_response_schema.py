"""Module for normalize endpoint response schema."""
from datetime import datetime
from enum import Enum
from typing import Any, Dict, List, Optional, Type, Union

from ga4gh.vrs import models
from pydantic import BaseModel, model_validator
from pydantic.types import StrictStr


class HGVSDupDelModeOption(str, Enum):
    """Define options for HGVSDupDelMode.
    This mode determines how to interpret HGVS dup/del.
    """

    DEFAULT = "default"
    COPY_NUMBER_COUNT = "copy_number_count"
    COPY_NUMBER_CHANGE = "copy_number_change"
    REFERENCE_LEN_EXPR = "ref_len_expr"  # VRS Allele
    LITERAL_SEQ_EXPR = "literal_seq_expr"  # VRS Allele


class ServiceMeta(BaseModel):
    """Metadata regarding the variation-normalization service."""

    name: StrictStr = "variation-normalizer"
    version: StrictStr
    response_datetime: datetime
    url: StrictStr = "https://github.com/cancervariants/variation-normalization"

    class Config:
        """Configure schema example."""

        @staticmethod
        def json_schema_extra(
            schema: Dict[str, Any], model: Type["ServiceMeta"]
        ) -> None:
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

    @model_validator(mode="before")
    def unique_warnings(cls, values):
        """Ensure unique warnings"""
        values["warnings"] = list(set(values["warnings"]))
        return values


class NormalizeService(ServiceResponse):
    """A response to normalizing a variation to a single GA4GH Value Object
    Descriptor.
    """

    variation_query: StrictStr
    variation: Optional[
        Union[models.Allele, models.CopyNumberCount, models.CopyNumberChange]
    ] = None

    class Config:
        """Configure model."""

        @staticmethod
        def json_schema_extra(
            schema: Dict[str, Any], model: Type["NormalizeService"]
        ) -> None:
            """Configure OpenAPI schema."""
            if "title" in schema.keys():
                schema.pop("title", None)
            for prop in schema.get("properties", {}).values():
                prop.pop("title", None)
            schema["example"] = {
                "variation_query": "BRAF V600E",
                "variation": {
                    "id": "ga4gh:VA.h313H4CQh6pogbbSJ3H5pI1cPoh9YMm_",
                    "location": {
                        "id": "ga4gh:SL.xfBTztcmMstx8jrrdgPiE_BUoLHLFMMS",
                        "end": {"value": 600, "type": "Number"},
                        "start": {"value": 599, "type": "Number"},
                        "sequence_id": "ga4gh:SQ.cQvw4UsHHRRlogxbWCB8W-mKD4AraM9y",
                        "type": "SequenceLocation",
                    },
                    "state": {"sequence": "E", "type": "LiteralSequenceExpression"},
                    "type": "Allele",
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
        def json_schema_extra(
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
