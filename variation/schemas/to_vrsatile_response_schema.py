"""Module for toVRSATILE endpoint response schema."""
from pydantic import BaseModel
from typing import List, Dict, Type, Any, Optional, Union
from pydantic.types import StrictStr
from ga4gh.vrsatile.pydantic.vrsatile_models import VariationDescriptor
from variation.schemas.normalize_response_schema import ServiceMeta


class ToVRSATILEService(BaseModel):
    """Define model for translation response."""

    search_term: StrictStr
    warnings: Optional[List[StrictStr]]
    VariationDescriptors: Optional[List[VariationDescriptor]]
    service_meta_: ServiceMeta

    class Config:
        """Configure model."""

        @staticmethod
        def schema_extra(schema: Dict[str, Any],
                         model: Type["ToVRSATILEService"]) -> None:
            """Configure OpenAPI schema."""
            if "title" in schema.keys():
                schema.pop("title", None)
            for prop in schema.get("properties", {}).values():
                prop.pop("title", None)
            schema["example"] = {
                "search_term": "BRAF V600E",
                "variationdescriptors": [
                    {   
                        "id": "normalize.variation:BRAF%20V600E",
                        "type": "VariationDescriptor",
                        "extensions": [
                            {
                                "type": "Extension",
                                "name": "possible accessions",
                                "value": ["NM_001354609.2"]
                            },
                            {
                                "type": "Extension",
                                "name": "human reference genome assembly",
                                "value": [
                                    "GRCh37"
                                ]
                            },
                            {
                                "type": "Extension",
                                "name": "mane status",
                                "value": [
                                    "N/A"
                                ]
                            }
                        ],
                        "variation_id": "ga4gh:VA.7ys8TiDzrk04O3Upd63__rOBCEhv3P5d",
                        "variation": {
                            "id": "ga4gh:VA.7ys8TiDzrk04O3Upd63__rOBCEhv3P5d",
                            "type": "Allele",
                            "location": {
                                "id": "ga4gh:VSL.Vxqx2bv42rWeu08Eg7JpkdQkMCNLskoz",
                                "type": "SequenceLocation",
                                "sequence_id": "ga4gh:SQ.ZJwurRo2HLY018wghYjDKSfIlEH0Y8At",
                                "interval": {
                                    "type": "SequenceInterval",
                                    "start": {
                                        "type": "Number",
                                        "value": 599
                                    },
                                    "end": {
                                        "type": "Number",
                                        "value": 600
                                    }
                                }
                            },
                            "state": {
                                "type": "LiteralSequenceExpression",
                                "sequence": "E"
                            }
                        },
                        "molecule_context": "protein",
                        "structural_type": "SO:0001606",
                        "gene_context": {
                            "id": "normalize.gene:BRAF",
                            "type": "GeneDescriptor",
                            "label": "BRAF",
                            "xrefs": [
                                "ensembl:ENSG00000157764",
                                "ncbigene:673"
                            ],
                            "alternate_labels": [
                                "RAFB1",
                                "B-RAF1",
                                "BRAF1",
                                "B-raf",
                                "NS7"
                            ], 
                            "extensions": [
                                {
                                    "type": "Extension",
                                    "name": "symbol_status",
                                    "value": "approved"
                                },
                                {
                                    "type": "Extension",
                                    "name": "approved_name",
                                    "value": "B-Raf proto-oncogene, serine/threonine kinase"
                                },
                                {
                                    "type": "Extension",
                                    "name": "chromosome_location",
                                    "value": {
                                        "species_id": "taxonomy:9606",
                                        "interval": {
                                            "type": "CytobandInterval",
                                            "start": "q34",
                                            "end": "q34"
                                        },
                                        "_id": "ga4gh:VCL.O6yCQ1cnThOrTfK9YUgMlTfM6HTqbrKw",
                                        "type": "ChromosomeLocation",
                                        "chr": "7"
                                    }
                                }, 
                                {
                                    "type": "Extension",
                                    "name": "associated_with",
                                    "value": [
                                        "pubmed:1565476",
                                        "pubmed:2284096",
                                        "orphanet:119066",
                                        "ucsc:uc003vwc.5",
                                        "refseq:NM_004333",
                                        "vega:OTTHUMG00000157457",
                                        "iuphar:1943",
                                        "ccds:CCDS5863",
                                        "cosmic:BRAF",
                                        "ena.embl:M95712",
                                        "uniprot:P15056",
                                        "ccds:CCDS87555",
                                        "omim:164757"
                                    ]
                                }
                            ],
                            "gene_id": "hgnc:1097"
                        }, 
                        "vrs_ref_allele_seq": "V"
                    },
                    {
                        "id": "normalize.variation:BRAF%20V600E",
                        "type": "VariationDescriptor",
                        "extensions": [
                            {
                                "type": "Extension",
                                "name": "possible accessions",
                                "value": [
                                    "NM_001378468.1"
                                ]
                            },
                            {
                                "type": "Extension",
                                "name": "human reference genome assembly",
                                "value": [
                                    "GRCh37"
                                ]
                            },
                            {
                                "type": "Extension",
                                "name": "mane status",
                                "value": [
                                    "N/A"
                                ]
                            }
                        ],
                        "variation_id": "ga4gh:VA.FzlrH5feNcQ3S9GayMU9EF008j-8Pbz5",
                        "variation": {
                            "id": "ga4gh:VA.FzlrH5feNcQ3S9GayMU9EF008j-8Pbz5",
                            "type": "Allele",
                            "location": {
                                "id": "ga4gh:VSL.QDLST2nKpPWwIArdO57L2VIWPNZ0DiN3",
                                "type": "SequenceLocation",
                                "sequence_id": "ga4gh:SQ.0Q-SgJX1V3seUUIu3qVUtEa55CQsGmEU",
                                "interval": {
                                    "type": "SequenceInterval",
                                    "start": {
                                        "type": "Number",
                                        "value": 599
                                    },
                                    "end": {
                                        "type": "Number",
                                        "value": 600
                                    }
                                }
                            },
                            "state": {
                                "type": "LiteralSequenceExpression",
                                "sequence": "E"
                            }
                        },
                        "molecule_context": "protein",
                        "structural_type": "SO:0001606",
                        "gene_context": {
                            "id": "normalize.gene:BRAF",
                            "type": "GeneDescriptor",
                            "label": "BRAF",
                            "xrefs": [
                                "ensembl:ENSG00000157764",
                                "ncbigene:673"
                            ],
                            "alternate_labels": [
                                "RAFB1",
                                "B-RAF1",
                                "BRAF1",
                                "B-raf",
                                "NS7"
                            ],
                            "extensions": [
                                {
                                    "type": "Extension",
                                    "name": "symbol_status",
                                    "value": "approved"
                                },
                                {
                                    "type": "Extension",
                                    "name": "approved_name",
                                    "value": "B-Raf proto-oncogene, serine/threonine kinase"
                                },
                                {
                                    "type": "Extension",
                                    "name": "chromosome_location",
                                    "value": {
                                        "species_id": "taxonomy:9606",
                                        "interval": {
                                            "type": "CytobandInterval",
                                            "start": "q34",
                                            "end": "q34"
                                        },
                                        "_id": "ga4gh:VCL.O6yCQ1cnThOrTfK9YUgMlTfM6HTqbrKw",
                                        "type": "ChromosomeLocation",
                                        "chr": "7"
                                    }
                                },
                                {
                                    "type": "Extension",
                                    "name": "associated_with",
                                    "value": [
                                        "pubmed:1565476",
                                        "pubmed:2284096",
                                        "orphanet:119066",
                                        "ucsc:uc003vwc.5",
                                        "refseq:NM_004333",
                                        "vega:OTTHUMG00000157457",
                                        "iuphar:1943",
                                        "ccds:CCDS5863",
                                        "cosmic:BRAF",
                                        "ena.embl:M95712",
                                        "uniprot:P15056",
                                        "ccds:CCDS87555",
                                        "omim:164757"
                                    ]
                                }
                            ],
                            "gene_id": "hgnc:1097"
                        },
                        "vrs_ref_allele_seq": "V"
                    },
                    {
                        "id": "normalize.variation:BRAF%20V600E",
                        "type": "VariationDescriptor",
                        "extensions": [
                            {
                                "type": "Extension",
                                "name": "possible accessions",
                                "value": [
                                    "NM_004333.6"
                                ]
                            },
                            {
                                "type": "Extension",
                                "name": "human reference genome assembly",
                                "value": [
                                    "GRCh37"
                                ]
                            },
                            {
                                "type": "Extension",
                                "name": "mane status",
                                "value": [
                                    "N/A"
                                ]
                            }
                        ],
                        "variation_id": "ga4gh:VA.ZDdoQdURgO2Daj2NxLj4pcDnjiiAsfbO",
                        "variation": {
                            "id": "ga4gh:VA.ZDdoQdURgO2Daj2NxLj4pcDnjiiAsfbO",
                            "type": "Allele",
                            "location": {
                                "id": "ga4gh:VSL.2cHIgn7iLKk4x9z3zLkSTTFMV0e48DR4",
                                "type": "SequenceLocation",
                                "sequence_id": "ga4gh:SQ.cQvw4UsHHRRlogxbWCB8W-mKD4AraM9y",
                                "interval": {
                                    "type": "SequenceInterval",
                                    "start": {
                                        "type": "Number",
                                        "value": 599
                                    },
                                    "end": {
                                        "type": "Number",
                                        "value": 600
                                    }
                                }
                            },
                            "state": {
                                "type": "LiteralSequenceExpression",
                                "sequence": "E"
                            }
                        },
                        "molecule_context": "protein",
                        "structural_type": "SO:0001606",
                        "gene_context": {
                            "id": "normalize.gene:BRAF",
                            "type": "GeneDescriptor",
                            "label": "BRAF",
                            "xrefs": [
                                "ensembl:ENSG00000157764",
                                "ncbigene:673"
                            ],
                            "alternate_labels": [
                                "RAFB1",
                                "B-RAF1",
                                "BRAF1",
                                "B-raf",
                                "NS7"
                            ],
                            "extensions": [
                                {
                                    "type": "Extension",
                                    "name": "symbol_status",
                                    "value": "approved"
                                },
                                {
                                    "type": "Extension",
                                    "name": "approved_name",
                                    "value": "B-Raf proto-oncogene, serine/threonine kinase"
                                },
                                {
                                    "type": "Extension",
                                    "name": "chromosome_location",
                                    "value": {
                                        "species_id": "taxonomy:9606",
                                        "interval": {
                                            "type": "CytobandInterval",
                                            "start": "q34",
                                            "end": "q34"
                                        },
                                        "_id": "ga4gh:VCL.O6yCQ1cnThOrTfK9YUgMlTfM6HTqbrKw",
                                        "type": "ChromosomeLocation",
                                        "chr": "7"
                                    }
                                },
                                {
                                    "type": "Extension",
                                    "name": "associated_with",
                                    "value": [
                                        "pubmed:1565476",
                                        "pubmed:2284096",
                                        "orphanet:119066",
                                        "ucsc:uc003vwc.5",
                                        "refseq:NM_004333",
                                        "vega:OTTHUMG00000157457",
                                        "iuphar:1943",
                                        "ccds:CCDS5863",
                                        "cosmic:BRAF",
                                        "ena.embl:M95712",
                                        "uniprot:P15056",
                                        "ccds:CCDS87555",
                                        "omim:164757"
                                    ]
                                }
                            ],
                            "gene_id": "hgnc:1097"
                        },
                        "vrs_ref_allele_seq": "V"
                    },
                    {
                        "id": "normalize.variation:BRAF%20V600E",
                        "type": "VariationDescriptor",
                        "extensions": [
                            {
                                "type": "Extension",
                                "name": "possible accessions",
                                "value": [
                                    "NM_001378474.1"
                                ]
                            },
                            {
                                "type": "Extension",
                                "name": "human reference genome assembly",
                                "value": [
                                    "GRCh37"
                                ]
                            },
                            {
                                "type": "Extension",
                                "name": "mane status",
                                "value": [
                                    "N/A"
                                ]
                            }
                        ],
                        "variation_id": "ga4gh:VA.vimwyw0pFTwatfFhi3rhhb153ARWsPrW",
                        "variation": {
                            "id": "ga4gh:VA.vimwyw0pFTwatfFhi3rhhb153ARWsPrW",
                            "type": "Allele",
                            "location": {
                                "id": "ga4gh:VSL.FVmsWpfSOA3B2ryq0k995oHMuSGiFvMa",
                                "type": "SequenceLocation",
                                "sequence_id": "ga4gh:SQ.lKdPZpuT-VNvRuKDjsUItNgutfWYgWQd",
                                "interval": {
                                    "type": "SequenceInterval",
                                    "start": {
                                        "type": "Number",
                                        "value": 599
                                    },
                                    "end": {
                                        "type": "Number",
                                        "value": 600
                                    }
                                }
                            },
                            "state": {
                                "type": "LiteralSequenceExpression",
                                "sequence": "E"
                            }
                        },
                        "molecule_context": "protein",
                        "structural_type": "SO:0001606",
                        "gene_context": {
                            "id": "normalize.gene:BRAF",
                            "type": "GeneDescriptor",
                            "label": "BRAF",
                            "xrefs": [
                                "ensembl:ENSG00000157764",
                                "ncbigene:673"
                            ],
                            "alternate_labels": [
                                "RAFB1",
                                "B-RAF1",
                                "BRAF1",
                                "B-raf",
                                "NS7"
                            ],
                            "extensions": [
                                {
                                    "type": "Extension",
                                    "name": "symbol_status",
                                    "value": "approved"
                                },
                                {
                                    "type": "Extension",
                                    "name": "approved_name",
                                    "value": "B-Raf proto-oncogene, serine/threonine kinase"
                                },
                                {
                                    "type": "Extension",
                                    "name": "chromosome_location",
                                    "value": {
                                        "species_id": "taxonomy:9606",
                                        "interval": {
                                            "type": "CytobandInterval",
                                            "start": "q34",
                                            "end": "q34"
                                        },
                                        "_id": "ga4gh:VCL.O6yCQ1cnThOrTfK9YUgMlTfM6HTqbrKw",
                                        "type": "ChromosomeLocation",
                                        "chr": "7"
                                    }
                                },
                                {
                                    "type": "Extension",
                                    "name": "associated_with",
                                    "value": [
                                        "pubmed:1565476",
                                        "pubmed:2284096",
                                        "orphanet:119066",
                                        "ucsc:uc003vwc.5",
                                        "refseq:NM_004333",
                                        "vega:OTTHUMG00000157457",
                                        "iuphar:1943",
                                        "ccds:CCDS5863",
                                        "cosmic:BRAF",
                                        "ena.embl:M95712",
                                        "uniprot:P15056",
                                        "ccds:CCDS87555",
                                        "omim:164757"
                                    ]
                                }
                            ],
                            "gene_id": "hgnc:1097"
                        },
                        "vrs_ref_allele_seq": "V"
                    },
                    {
                        "id": "normalize.variation:BRAF%20V600E",
                        "type": "VariationDescriptor",
                        "extensions": [
                            {
                                "type": "Extension",
                                "name": "possible accessions",
                                "value": [
                                    "NM_001354609.2",
                                    "NM_001378468.1",
                                    "NM_004333.6",
                                    "NM_001378474.1",
                                    "ENST00000288602.6"
                                ]
                            },
                            {
                                "type": "Extension",
                                "name": "human reference genome assembly",
                                "value": [
                                    "GRCh38"
                                ]
                            },
                            {
                                "type": "Extension",
                                "name": "mane status",
                                "value": [
                                    "mane transcript"
                                ]
                            }
                        ],
                        "variation_id": "ga4gh:VA.8JkgnqIgYqufNl-OV_hpRG_aWF9UFQCE",
                        "variation": {
                            "id": "ga4gh:VA.8JkgnqIgYqufNl-OV_hpRG_aWF9UFQCE",
                            "type": "Allele",
                            "location": {
                                "id": "ga4gh:VSL.AqrQ-EkAvTrXOFn70_8i3dXF5shBBZ5i",
                                "type": "SequenceLocation",
                                "sequence_id": "ga4gh:SQ.WaAJ_cXXn9YpMNfhcq9lnzIvaB9ALawo",
                                "interval": {
                                    "type": "SequenceInterval",
                                    "start": {
                                        "type": "Number",
                                        "value": 639
                                    },
                                    "end": {
                                        "type": "Number",
                                        "value": 640
                                    }
                                }
                            },
                            "state": {
                                "type": "LiteralSequenceExpression",
                                "sequence": "E"
                            }
                        },
                        "molecule_context": "protein",
                        "structural_type": "SO:0001606",
                        "gene_context": {
                            "id": "normalize.gene:BRAF",
                            "type": "GeneDescriptor",
                            "label": "BRAF",
                            "xrefs": [
                                "ensembl:ENSG00000157764",
                                "ncbigene:673"
                            ],
                            "alternate_labels": [
                                "RAFB1",
                                "B-RAF1",
                                "BRAF1",
                                "B-raf",
                                "NS7"
                            ],
                            "extensions": [
                                {
                                    "type": "Extension",
                                    "name": "symbol_status",
                                    "value": "approved"
                                },
                                {
                                    "type": "Extension",
                                    "name": "approved_name",
                                    "value": "B-Raf proto-oncogene, serine/threonine kinase"
                                },
                                {
                                    "type": "Extension",
                                    "name": "chromosome_location",
                                    "value": {
                                        "species_id": "taxonomy:9606",
                                        "interval": {
                                            "type": "CytobandInterval",
                                            "start": "q34",
                                            "end": "q34"
                                        },
                                        "_id": "ga4gh:VCL.O6yCQ1cnThOrTfK9YUgMlTfM6HTqbrKw",
                                        "type": "ChromosomeLocation",
                                        "chr": "7"
                                    }
                                },
                                {
                                    "type": "Extension",
                                    "name": "associated_with",
                                    "value": [
                                        "pubmed:1565476",
                                        "pubmed:2284096",
                                        "orphanet:119066",
                                        "ucsc:uc003vwc.5",
                                        "refseq:NM_004333",
                                        "vega:OTTHUMG00000157457",
                                        "iuphar:1943",
                                        "ccds:CCDS5863",
                                        "cosmic:BRAF",
                                        "ena.embl:M95712",
                                        "uniprot:P15056",
                                        "ccds:CCDS87555",
                                        "omim:164757"
                                    ]
                                }
                            ],
                            "gene_id": "hgnc:1097"
                        },
                        "vrs_ref_allele_seq": "V"
                    }
                ],
                "service_meta_": {
                    "name": "variation-normalizer",
                    "version": "0.2.13",
                    "response_datetime": "2021-11-18T14:10:53.909158",
                    "url": "https://github.com/cancervariants/variation-normalization"  # noqa: E501
                }

            }
