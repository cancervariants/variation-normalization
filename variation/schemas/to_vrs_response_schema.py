"""Module for to_vrs endpoint response schema."""
from typing import List, Union

from ga4gh.vrs import models
from pydantic import BaseModel, ConfigDict, StrictStr

from variation.schemas.normalize_response_schema import ServiceMeta
from variation.version import __version__


class ToVRSService(BaseModel):
    """Define model for translation response."""

    search_term: StrictStr
    warnings: List[StrictStr] = []
    variations: Union[
        List[models.Allele],
        List[models.CopyNumberCount],
        List[models.CopyNumberChange],
    ] = []
    service_meta_: ServiceMeta

    # FIXME:
    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "search_term": "BRAF V600E",
                "warnings": [],
                "variations": [
                    {
                        "id": "ga4gh:VA.JKi4CDlXNAIfxn-VxPTFSRvamGqBBXbi",
                        "type": "Allele",
                        "location": {
                            "id": "ga4gh:SL.pnvvgXpxPek5MU0dChggjPZC18wxcvTY",
                            "type": "SequenceLocation",
                            "sequence": "ga4gh:SQ.ZJwurRo2HLY018wghYjDKSfIlEH0Y8At",
                            "start": 599,
                            "end": 600,
                        },
                        "state": {
                            "type": "ReferenceLengthExpression",
                            "length": 1,
                            "repeatSubunitLength": 0,
                        },
                    },
                    {
                        "id": "ga4gh:VA.dpSeCSWRpr5-nR_IWbtMzFigbCb1_pn1",
                        "type": "Allele",
                        "location": {
                            "id": "ga4gh:SL.ko4RJfU-2fvZrbCDpo6i-Ljcfi59TcQI",
                            "type": "SequenceLocation",
                            "sequence": "ga4gh:SQ.cQvw4UsHHRRlogxbWCB8W-mKD4AraM9y",
                            "start": 599,
                            "end": 600,
                        },
                        "state": {
                            "type": "ReferenceLengthExpression",
                            "length": 1,
                            "repeatSubunitLength": 0,
                        },
                    },
                    {
                        "id": "ga4gh:VA.27lvBV4CP0pNPw_nM-fpJOprookDojhs",
                        "type": "Allele",
                        "location": {
                            "id": "ga4gh:SL.VtMGH5nBT3MtxEM46Q3Dflw0UhkdXbeX",
                            "type": "SequenceLocation",
                            "sequence": "ga4gh:SQ.0Q-SgJX1V3seUUIu3qVUtEa55CQsGmEU",
                            "start": 599,
                            "end": 600,
                        },
                        "state": {
                            "type": "ReferenceLengthExpression",
                            "length": 1,
                            "repeatSubunitLength": 0,
                        },
                    },
                    {
                        "id": "ga4gh:VA.WB9SKyyLApRCl0yteR_j1daCfh4QcNJN",
                        "type": "Allele",
                        "location": {
                            "id": "ga4gh:SL.2fANtMDFqO2luhgUkVEP9-XAxPCLWRHv",
                            "type": "SequenceLocation",
                            "sequence": "ga4gh:SQ.lKdPZpuT-VNvRuKDjsUItNgutfWYgWQd",
                            "start": 599,
                            "end": 600,
                        },
                        "state": {
                            "type": "ReferenceLengthExpression",
                            "length": 1,
                            "repeatSubunitLength": 0,
                        },
                    },
                ],
                "service_meta_": {
                    "name": "variation-normalizer",
                    "version": __version__,
                    "response_datetime": "2023-08-24T09:05:03.622667",
                    "url": "https://github.com/cancervariants/variation-normalization",
                },
            }
        }
    )
