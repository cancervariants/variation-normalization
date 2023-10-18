"""Module for vrs-python translator endpoint response schema"""
from enum import Enum
from typing import List, Literal, Optional, Union

from ga4gh.vrsatile.pydantic.vrs_models import Allele
from pydantic import BaseModel, ConfigDict, StrictStr

from variation.schemas.normalize_response_schema import ServiceMeta


class VrsPythonMeta(BaseModel):
    """Metadata regarding vrs-python dependency"""

    name: Literal["vrs-python"] = "vrs-python"
    version: StrictStr
    url: Literal[
        "https://github.com/ga4gh/vrs-python"
    ] = "https://github.com/ga4gh/vrs-python"


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

    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "variation": {
                    "_id": "ga4gh:VA.jCQx4yBcU6u6u1RcT9Tp0PjhaQ6ynicY",
                    "type": "Allele",
                    "location": {
                        "_id": "ga4gh:VSL.8cQ9y-2J75mg5ioJvbqtgwiskdUV4zuO",
                        "type": "SequenceLocation",
                        "sequence_id": "ga4gh:SQ.IW78mgV5Cqf6M24hy52hPjyyo5tCCd86",
                        "interval": {
                            "type": "SequenceInterval",
                            "start": {"type": "Number", "value": 140453135},
                            "end": {"type": "Number", "value": 140453136},
                        },
                    },
                    "state": {"type": "LiteralSequenceExpression", "sequence": "T"},
                },
                "fmt": "hgvs",
            }
        }
    )


class TranslateToHGVSQuery(BaseModel):
    """Query fields for Translate To HGVS Service"""

    variation: Allele
    namespace: Optional[str] = None

    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "variation": {
                    "_id": "ga4gh:VA.jCQx4yBcU6u6u1RcT9Tp0PjhaQ6ynicY",
                    "type": "Allele",
                    "location": {
                        "_id": "ga4gh:VSL.8cQ9y-2J75mg5ioJvbqtgwiskdUV4zuO",
                        "type": "SequenceLocation",
                        "sequence_id": "ga4gh:SQ.IW78mgV5Cqf6M24hy52hPjyyo5tCCd86",
                        "interval": {
                            "type": "SequenceInterval",
                            "start": {"type": "Number", "value": 140453135},
                            "end": {"type": "Number", "value": 140453136},
                        },
                    },
                    "state": {"type": "LiteralSequenceExpression", "sequence": "T"},
                },
                "namespace": "refseq",
            }
        }
    )


class TranslateFromQuery(BaseModel):
    """Query fields for Translate From Service"""

    variation: StrictStr
    fmt: Optional[TranslateFromFormat] = None


class TranslateService(BaseModel):
    """Response schema for vrs-python translator endpoints"""

    query: Union[TranslateFromQuery, TranslateToQuery, TranslateToHGVSQuery]
    warnings: Optional[List[StrictStr]] = None
    service_meta_: ServiceMeta
    vrs_python_meta_: VrsPythonMeta


class TranslateFromService(TranslateService):
    """Response schema for vrs-python translate from endpoint"""

    variation: Optional[Allele] = None


class TranslateToService(TranslateService):
    """Response schema for vrs-python translate to endpoint"""

    variations: List[StrictStr]
