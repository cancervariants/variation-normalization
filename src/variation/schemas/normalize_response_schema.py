"""Module for normalize endpoint response schema."""
from datetime import datetime
from enum import Enum
from typing import List, Literal, Optional, Union

from ga4gh.vrs import models
from pydantic import BaseModel, ConfigDict, StrictStr, model_validator

from variation.version import __version__


class HGVSDupDelModeOption(str, Enum):
    """Define options for HGVSDupDelMode.
    This mode determines how to interpret HGVS dup/del.
    """

    DEFAULT = "default"
    COPY_NUMBER_COUNT = "copy_number_count"
    COPY_NUMBER_CHANGE = "copy_number_change"
    ALLELE = "allele"


class ServiceMeta(BaseModel):
    """Metadata regarding the variation-normalization service."""

    name: Literal["variation-normalizer"] = "variation-normalizer"
    version: StrictStr
    response_datetime: datetime
    url: Literal[
        "https://github.com/cancervariants/variation-normalization"
    ] = "https://github.com/cancervariants/variation-normalization"

    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "name": "variation-normalizer",
                "version": __version__,
                "response_datetime": "2021-04-05T16:44:15.367831",
                "url": "https://github.com/cancervariants/variation-normalization",
            }
        }
    )


class ServiceResponse(BaseModel):
    """Base response model for services"""

    warnings: List[StrictStr] = []
    service_meta_: ServiceMeta

    @model_validator(mode="after")
    def unique_warnings(cls, v):
        """Ensure unique warnings"""
        v.warnings = list(set(v.warnings))
        return v

    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "warnings": [],
                "service_meta_": {
                    "name": "variation-normalizer",
                    "version": __version__,
                    "response_datetime": "2021-04-05T16:44:15.367831",
                    "url": "https://github.com/cancervariants/variation-normalization",
                },
            }
        }
    )


class NormalizeService(ServiceResponse):
    """A response to normalizing a variation to a single GA4GH VRS Variation"""

    variation_query: StrictStr
    variation: Optional[
        Union[models.Allele, models.CopyNumberCount, models.CopyNumberChange]
    ] = None

    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "variation_query": "BRAF V600E",
                "variation": {
                    "id": "ga4gh:VA.j4XnsLZcdzDIYa5pvvXM7t1wn9OITr0L",
                    "digest": "j4XnsLZcdzDIYa5pvvXM7t1wn9OITr0L",
                    "location": {
                        "id": "ga4gh:SL.t-3DrWALhgLdXHsupI-e-M00aL3HgK3y",
                        "digest": "t-3DrWALhgLdXHsupI-e-M00aL3HgK3y",
                        "end": 600,
                        "start": 599,
                        "sequenceReference": {
                            "type": "SequenceReference",
                            "refgetAccession": "SQ.cQvw4UsHHRRlogxbWCB8W-mKD4AraM9y",
                        },
                        "type": "SequenceLocation",
                    },
                    "state": {"sequence": "E", "type": "LiteralSequenceExpression"},
                    "type": "Allele",
                },
                "service_meta_": {
                    "name": "variation-normalizer",
                    "version": __version__,
                    "response_datetime": "2022-01-26T22:23:41.821673",
                    "url": "https://github.com/cancervariants/variation-normalization",
                },
            }
        }
    )


class TranslateIdentifierService(ServiceResponse):
    """A response to translating identifiers."""

    identifier_query: StrictStr
    aliases: List[StrictStr] = []

    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "identifier_query": "NP_004324.2",
                "warnings": [],
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
                    "version": __version__,
                    "response_datetime": "2021-11-18T14:10:53.909158",
                    "url": "https://github.com/cancervariants/variation-normalization",
                },
            }
        }
    )
