"""Module for vrs-python translator endpoint response schema"""
from enum import Enum
from typing import Optional, List, Union

from ga4gh.vrsatile.pydantic.vrs_models import Allele
from pydantic import BaseModel
from pydantic.types import StrictStr

from variation.schemas.normalize_response_schema import ServiceMeta


class VrsPythonMeta(BaseModel):
    """Metadata regarding vrs-python dependency"""

    name = "vrs-python"
    version: str
    url = "https://github.com/ga4gh/vrs-python"


class TranslateFromFormat(str, Enum):
    """Enums for formats that vrs-python can translate from"""

    HGVS = "hgvs"
    BEACON = "beacon"
    GNOMAD = "gnomad"
    SPDI = "spdi"


class TranslateToFormat(str, Enum):
    """Enums for formats that vrs-python can translate to"""

    HGVS = "hgvs"
    SPDI = "spdi"


class TranslateToQuery(BaseModel):
    """Query fields for Translate To Service"""

    variation: Allele
    fmt: TranslateToFormat


class TranslateFromQuery(BaseModel):
    """Query fields for Translate From Service"""

    variation: StrictStr
    fmt: Optional[TranslateFromFormat] = None


class TranslateService(BaseModel):
    """Response schema for vrs-python translator endpoints"""

    query: Union[TranslateFromQuery, TranslateToQuery]
    warnings: Optional[List[StrictStr]]
    service_meta_: ServiceMeta
    vrs_python_meta_: VrsPythonMeta


class TranslateFromService(TranslateService):
    """Response schema for vrs-python translate from endpoint"""

    variation: Optional[Allele] = None


class TranslateToService(TranslateService):
    """Response schema for vrs-python translate to endpoint"""

    variations: List[StrictStr]


class ErrorResponse(BaseModel):
    """Response for Translate services when request body is invalid"""

    errors: List[StrictStr]
