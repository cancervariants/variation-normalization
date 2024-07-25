"""Module containing schemas for services"""

from enum import Enum
from typing import Literal

from cool_seq_tool.schemas import ResidueMode
from pydantic import BaseModel, ConfigDict, StrictInt, StrictStr

from variation import __version__
from variation.schemas.normalize_response_schema import ServiceMeta


class ClinVarAssembly(str, Enum):
    """Define assemblies in ClinVar"""

    GRCH38 = "GRCh38"
    GRCH37 = "GRCh37"
    NCBI36 = "NCBI36"
    HG38 = "hg38"
    HG19 = "hg19"
    HG18 = "hg18"


class CdnaRepresentation(BaseModel, extra="forbid"):
    """Model response for cDNA representation"""

    c_ac: StrictStr
    c_start_pos: StrictInt
    c_end_pos: StrictInt
    cds_start: StrictInt
    residue_mode: Literal["inter-residue"] = ResidueMode.INTER_RESIDUE.value

    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "c_ac": "NM_004333.6",
                "c_start_pos": 1797,
                "c_end_pos": 1800,
                "cds_start": 226,
                "residue_mode": ResidueMode.INTER_RESIDUE.value,
            }
        }
    )


class ToCdnaService(BaseModel, extra="forbid"):
    """Service model response for protein -> cDNA"""

    c_data: CdnaRepresentation | None = None
    warnings: list[StrictStr] = []
    service_meta: ServiceMeta

    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "c_data": {
                    "c_ac": "NM_004333.6",
                    "c_start_pos": 1797,
                    "c_end_pos": 1800,
                    "cds_start": 226,
                    "residue_mode": "inter-residue",
                },
                "warnings": [],
                "service_meta": {
                    "version": __version__,
                    "response_datetime": "2022-09-29T15:08:18.696882",
                    "name": "variation-normalizer",
                    "url": "https://github.com/cancervariants/variation-normalization",
                },
            }
        }
    )


class GenomicRepresentation(BaseModel, extra="forbid"):
    """Model response for genomic representation"""

    g_ac: StrictStr
    g_start_pos: StrictInt
    g_end_pos: StrictInt
    residue_mode: Literal["inter-residue"] = ResidueMode.INTER_RESIDUE.value

    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "g_ac": "NC_000007.13",
                "g_start_pos": 140453134,
                "g_end_pos": 140453137,
                "residue_mode": ResidueMode.INTER_RESIDUE.value,
            }
        }
    )


class ToGenomicService(BaseModel, extra="forbid"):
    """Service model response for cDNA -> genomic"""

    g_data: GenomicRepresentation | None = None
    warnings: list[StrictStr] = []
    service_meta: ServiceMeta

    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "g_data": {
                    "g_ac": "NC_000007.13",
                    "g_start_pos": 140453134,
                    "g_end_pos": 140453137,
                    "residue_mode": "inter-residue",
                },
                "warnings": [],
                "service_meta": {
                    "version": __version__,
                    "response_datetime": "2022-09-29T15:08:18.696882",
                    "name": "variation-normalizer",
                    "url": "https://github.com/cancervariants/variation-normalization",
                },
            }
        }
    )
