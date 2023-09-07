"""Module containing schemas for services"""
from enum import Enum
from typing import Any, Dict, Type

from cool_seq_tool.schemas import ToCdnaService as ToCdna
from cool_seq_tool.schemas import ToGenomicService as ToGenomic

from variation.schemas.normalize_response_schema import ServiceMeta


class ClinVarAssembly(str, Enum):
    """Define assemblies in ClinVar"""

    GRCH38 = "GRCh38"
    GRCH37 = "GRCh37"
    NCBI36 = "NCBI36"
    HG38 = "hg38"
    HG19 = "hg19"
    HG18 = "hg18"


class ToCdnaService(ToCdna):
    """Service model response for protein -> cDNA"""

    service_meta: ServiceMeta

    class Config:
        """Configure model."""

        @staticmethod
        def schema_extra(schema: Dict[str, Any], model: Type["ToCdnaService"]) -> None:
            """Configure OpenAPI schema."""
            if "title" in schema.keys():
                schema.pop("title", None)
            for prop in schema.get("properties", {}).values():
                prop.pop("title", None)
            schema["example"] = {
                "c_data": {
                    "c_ac": "NM_004333.6",
                    "c_start_pos": 1797,
                    "c_end_pos": 1800,
                    "cds_start": 226,
                    "residue_mode": "inter-residue",
                },
                "warnings": list(),
                "service_meta": {
                    "version": "0.5.4",
                    "response_datetime": "2022-09-29T15:08:18.696882",
                    "name": "variation-normalizer",
                    "url": "https://github.com/cancervariants/variation-normalization",
                },
            }


class ToGenomicService(ToGenomic):
    """Model response for genomic representation"""

    service_meta: ServiceMeta

    class Config:
        """Configure model."""

        @staticmethod
        def schema_extra(schema: Dict[str, Any], model: Type["ToCdnaService"]) -> None:
            """Configure OpenAPI schema."""
            if "title" in schema.keys():
                schema.pop("title", None)
            for prop in schema.get("properties", {}).values():
                prop.pop("title", None)
            schema["example"] = {
                "g_data": {
                    "g_ac": "NC_000007.13",
                    "g_start_pos": 140453134,
                    "g_end_pos": 140453137,
                    "residue_mode": "inter-residue",
                },
                "warnings": list(),
                "service_meta": {
                    "version": "0.5.4",
                    "response_datetime": "2022-09-29T15:08:18.696882",
                    "name": "variation-normalizer",
                    "url": "https://github.com/cancervariants/variation-normalization",
                },
            }
