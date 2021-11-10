"""Module for testing HGVS Dup Del mode."""
import pytest
from variation.query import QueryHandler
from ga4gh.vrsatile.pydantic.vrsatile_model import VariationDescriptor
from tests.conftest import assertion_checks


@pytest.fixture(scope="module")
def test_normalize():
    """Build normalize test fixture."""
    class TestNormalize:

        def __init__(self):
            self.query_handler = QueryHandler()
            self.warnings = []

        def normalize(self, q, hgvs_dup_del_mode="default"):
            resp = self.query_handler.normalize(
                q, hgvs_dup_del_mode=hgvs_dup_del_mode)
            self.warnings = \
                self.query_handler.normalize_handler.warnings
            return resp
    return TestNormalize()


@pytest.fixture(scope='module')
def genomic_dup1():
    """Create test fixture containing params for genomic dup VD."""
    params = {
        "id": "normalize.variation:NC_000003.12%3Ag.49531262dup",
        "type": "VariationDescriptor",
        "variation_id": "",
        "variation": dict(),
        "molecule_context": "genomic",
        "structural_type": "SO:1000035",
        "vrs_ref_allele_seq": "GG"
    }
    return params


@pytest.fixture(scope='module')
def genomic_dup1_default(genomic_dup1):
    """Create a test fixture for genomic dup default and LSE."""
    genomic_dup1["variation_id"] = "ga4gh:VA.GEhb3s5YP5xr5wqrVLT-1h9xryBdj9Ia"
    genomic_dup1["variation"] = {
        "type": "Allele",
        "location": {
            "sequence_id": "ga4gh:SQ.Zu7h9AggXxhTaGVsy7h_EZSChSZGcmgX",
            "interval": {
                "type": "SimpleInterval",
                "start": 49531260,
                "end": 49531262,
            },
            "type": "SequenceLocation",
        },
        "state": {
            "type": "LiteralSequenceExpression",
            "sequence": "GGG"
        }
    }
    return VariationDescriptor(**genomic_dup1)


@pytest.fixture(scope='module')
def genomic_dup1_cnv(genomic_dup1):
    """Create a test fixture for genomic dup CNV."""
    genomic_dup1["variation"] = {
        "type": "CopyNumber",
        "subject": {
            "location": {
                "sequence_id": "ga4gh:SQ.Zu7h9AggXxhTaGVsy7h_EZSChSZGcmgX",
                "interval": {
                    "type": "SimpleInterval",
                    "start": 49531260,
                    "end": 49531262,
                },
                "type": "SequenceLocation",
            },
            "reverse_complement": False,
            "type": "DerivedSequenceExpression"
        },
        "copies": {
            "type": "Number",
            "value": 3
        }
    }
    genomic_dup1["variation_id"] = "ga4gh:VCN.iYQ_ZTSxXApMw-nMXBgzqToGbdSWx929"
    return VariationDescriptor(**genomic_dup1)


@pytest.fixture(scope='module')
def genomic_dup1_rse(genomic_dup1):
    """Create a test fixture for genomic dup RSE."""
    seq_loc = {
        "sequence_id": "ga4gh:SQ.Zu7h9AggXxhTaGVsy7h_EZSChSZGcmgX",
        "interval": {
            "type": "SimpleInterval",
            "start": 49531260,
            "end": 49531262,
        },
        "type": "SequenceLocation",
    }
    genomic_dup1["variation"] = {
        "type": "Allele",
        "location": seq_loc,
        "state": {
            "type": "RepeatedSequenceExpression",
            "seq_expr": {
                "type": "DerivedSequenceExpression",
                "location": seq_loc,
                "reverse_complement": False
            },
            "count": {
                "type": "Number",
                "value": 2
            }
        }
    }
    genomic_dup1["variation_id"] = "ga4gh:VA.FqgQtoHT1WVj5dY52YtcPyC_Iu86A4Al"
    return VariationDescriptor(**genomic_dup1)


@pytest.fixture(scope='module')
def genomic_dup2():
    """Create test fixture containing params for genomic dup VD."""
    params = {
        "id": "normalize.variation:NC_000016.10%3Ag.2087938_2087948dup",
        "type": "VariationDescriptor",
        "variation_id": "",
        "variation": dict(),
        "molecule_context": "genomic",
        "structural_type": "SO:1000035",
        "vrs_ref_allele_seq": "AAAGGTAGGGC"
    }
    return params


@pytest.fixture(scope='module')
def genomic_dup2_default(genomic_dup2):
    """Create a test fixture for genomic dup default and LSE."""
    genomic_dup2["variation_id"] = "ga4gh:VA.AFhoQKzjhj87w-Yx0Mfz1yRU31Wul9oi"
    genomic_dup2["variation"] = {
        "type": "Allele",
        "location": {
            "sequence_id": "ga4gh:SQ.yC_0RBj3fgBlvgyAuycbzdubtLxq-rE0",
            "interval": {
                "type": "SimpleInterval",
                "start": 2087937,
                "end": 2087948,
            },
            "type": "SequenceLocation",
        },
        "state": {
            "type": "LiteralSequenceExpression",
            "sequence": "AAAGGTAGGGCAAAGGTAGGGC"
        }
    }
    return VariationDescriptor(**genomic_dup2)


@pytest.fixture(scope='module')
def genomic_dup2_cnv(genomic_dup2):
    """Create a test fixture for genomic dup CNV."""
    genomic_dup2["variation"] = {
        "type": "CopyNumber",
        "subject": {
            "location": {
                "sequence_id": "ga4gh:SQ.yC_0RBj3fgBlvgyAuycbzdubtLxq-rE0",
                "interval": {
                    "type": "SimpleInterval",
                    "start": 2087937,
                    "end": 2087948,
                },
                "type": "SequenceLocation",
            },
            "reverse_complement": False,
            "type": "DerivedSequenceExpression"
        },
        "copies": {
            "type": "Number",
            "value": 3
        }
    }
    genomic_dup2["variation_id"] = "ga4gh:VCN.A7ctORrxXrtyqPbqEm5_UgFr_sv8IaYm"
    return VariationDescriptor(**genomic_dup2)


@pytest.fixture(scope='module')
def genomic_dup2_rse(genomic_dup2):
    """Create a test fixture for genomic dup RSE."""
    seq_loc = {
        "sequence_id": "ga4gh:SQ.yC_0RBj3fgBlvgyAuycbzdubtLxq-rE0",
        "interval": {
            "type": "SimpleInterval",
            "start": 2087937,
            "end": 2087948,
        },
        "type": "SequenceLocation",
    }
    genomic_dup2["variation"] = {
        "type": "Allele",
        "location": seq_loc,
        "state": {
            "type": "RepeatedSequenceExpression",
            "seq_expr": {
                "type": "DerivedSequenceExpression",
                "location": seq_loc,
                "reverse_complement": False
            },
            "count": {
                "type": "Number",
                "value": 2
            }
        }
    }
    genomic_dup2["variation_id"] = "ga4gh:VA.0wcKT02yJBPhYf5Uepdq_uyI2s4PGQYF"
    return VariationDescriptor(**genomic_dup2)


@pytest.fixture(scope='module')
def genomic_dup3():
    """Create test fixture containing params for genomic dup VD."""
    params = {
        "id": "normalize.variation:NC_000023.11%3Ag.%2831060227_31100351%29_%2833274278_33417151%29dup",  # noqa: E501
        "type": "VariationDescriptor",
        "variation_id": "",
        "variation": dict(),
        "molecule_context": "genomic",
        "structural_type": "SO:0001742",
        "vrs_ref_allele_seq": None
    }
    return params


@pytest.fixture(scope='module')
def genomic_dup3_default(genomic_dup3):
    """Create a test fixture for genomic dup default and cnv."""
    genomic_dup3["variation"] = {
        "type": "CopyNumber",
        "subject": {
            "location": {
                "sequence_id": "ga4gh:SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
                "interval": {
                    "type": "SequenceInterval",
                    "start": {
                        "min": 31060226,
                        "max": 31100350,
                        "type": "DefiniteRange"
                    },
                    "end": {
                        "min": 33274279,
                        "max": 33417152,
                        "type": "DefiniteRange"
                    }
                },
                "type": "SequenceLocation",
            },
            "reverse_complement": False,
            "type": "DerivedSequenceExpression"
        },
        "copies": {
            "type": "DefiniteRange",
            "min": 2,
            "max": 3
        }
    }
    genomic_dup3["variation_id"] = "ga4gh:VCN.IgQATuKrM_J5MDHm2VemKThFOkzz-7AZ"
    return VariationDescriptor(**genomic_dup3)


@pytest.fixture(scope='module')
def genomic_dup3_rse_lse(genomic_dup3):
    """Create test fixture for genomic dup rse and lse."""
    params = {
        "id": genomic_dup3["id"],
        "type": genomic_dup3["type"],
        "variation": {
            "type": "Text",
            "definition": "NC_000023.11:g.(31060227_31100351)_(33274278_33417151)dup"  # noqa: E501
        }
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope='module')
def genomic_dup4():
    """Create test fixture containing params for genomic dup VD."""
    params = {
        "id": "normalize.variation:NC_000020.11%3Ag.%28%3F_30417576%29_%2831394018_%3F%29dup",  # noqa: E501
        "type": "VariationDescriptor",
        "variation_id": "",
        "variation": dict(),
        "molecule_context": "genomic",
        "structural_type": "SO:0001742",
        "vrs_ref_allele_seq": None
    }
    return params


@pytest.fixture(scope='module')
def genomic_dup4_default(genomic_dup4):
    """Create a test fixture for genomic dup default and cnv."""
    genomic_dup4["variation"] = {
        "type": "CopyNumber",
        "subject": {
            "location": {
                "sequence_id": "ga4gh:SQ.-A1QmD_MatoqxvgVxBLZTONHz9-c7nQo",
                "interval": {
                    "type": "SequenceInterval",
                    "start": {
                        "value": 30417575,
                        "comparator": "<=",
                        "type": "IndefiniteRange"
                    },
                    "end": {
                        "value": 31394018,
                        "comparator": ">=",
                        "type": "IndefiniteRange"
                    }
                },
                "type": "SequenceLocation",
            },
            "reverse_complement": False,
            "type": "DerivedSequenceExpression"
        },
        "copies": {
            "type": "Number",
            "value": 3
        }
    }
    genomic_dup4["variation_id"] = "ga4gh:VCN.3rvfUmiIb4hSxVQhXKOonuOY6Q3xTkKx"
    return VariationDescriptor(**genomic_dup4)


@pytest.fixture(scope='module')
def genomic_dup4_rse_lse(genomic_dup4):
    """Create test fixture for genomic dup rse and lse."""
    params = {
        "id": genomic_dup4["id"],
        "type": genomic_dup4["type"],
        "variation": {
            "type": "Text",
            "definition": "NC_000020.11:g.(?_30417576)_(31394018_?)dup"
        }
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope='module')
def genomic_dup5():
    """Create test fixture containing params for genomic dup VD."""
    params = {
        "id": "normalize.variation:NC_000023.11%3Ag.%28%3F_154021812%29_154092209dup",  # noqa: E501
        "type": "VariationDescriptor",
        "variation_id": "",
        "variation": dict(),
        "molecule_context": "genomic",
        "structural_type": "SO:0001742",
        "vrs_ref_allele_seq": None
    }
    return params


@pytest.fixture(scope='module')
def genomic_dup5_default(genomic_dup5):
    """Create a test fixture for genomic dup default and cnv."""
    genomic_dup5["variation"] = {
        "type": "CopyNumber",
        "subject": {
            "location": {
                "sequence_id": "ga4gh:SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
                "interval": {
                    "type": "SequenceInterval",
                    "start": {
                        "value": 154021811,
                        "comparator": "<=",
                        "type": "IndefiniteRange"
                    },
                    "end": {
                        "value": 154092209,
                        "type": "Number"
                    }
                },
                "type": "SequenceLocation",
            },
            "reverse_complement": False,
            "type": "DerivedSequenceExpression"
        },
        "copies": {
            "type": "DefiniteRange",
            "min": 2,
            "max": 3
        }
    }
    genomic_dup5["variation_id"] = "ga4gh:VCN.eLAZZ-ht1h2dTtZqzhO9TVhBdFufv67-"
    return VariationDescriptor(**genomic_dup5)


@pytest.fixture(scope='module')
def genomic_dup5_rse_lse(genomic_dup5):
    """Create test fixture for genomic dup rse and lse."""
    params = {
        "id": genomic_dup5["id"],
        "type": genomic_dup5["type"],
        "variation": {
            "type": "Text",
            "definition": "NC_000023.11:g.(?_154021812)_154092209dup"
        }
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope='module')
def genomic_dup6():
    """Create test fixture containing params for genomic dup VD."""
    params = {
        "id": "normalize.variation:NC_000023.11%3Ag.154021812_%28154092209_%3F%29dup",  # noqa: E501
        "type": "VariationDescriptor",
        "variation_id": "",
        "variation": dict(),
        "molecule_context": "genomic",
        "structural_type": "SO:0001742",
        "vrs_ref_allele_seq": None
    }
    return params


@pytest.fixture(scope='module')
def genomic_dup6_default(genomic_dup6):
    """Create a test fixture for genomic dup default and cnv."""
    genomic_dup6["variation"] = {
        "type": "CopyNumber",
        "subject": {
            "location": {
                "sequence_id": "ga4gh:SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
                "interval": {
                    "type": "SequenceInterval",
                    "start": {
                        "value": 154021811,
                        "type": "Number"
                    },
                    "end": {
                        "value": 154092209,
                        "comparator": ">=",
                        "type": "IndefiniteRange"
                    }
                },
                "type": "SequenceLocation",
            },
            "reverse_complement": False,
            "type": "DerivedSequenceExpression"
        },
        "copies": {
            "type": "DefiniteRange",
            "min": 2,
            "max": 3
        }
    }
    genomic_dup6["variation_id"] = "ga4gh:VCN.Rekk_MmUQ777V76S51x7nZGjh4U3LkLy"
    return VariationDescriptor(**genomic_dup6)


@pytest.fixture(scope='module')
def genomic_dup6_rse_lse(genomic_dup6):
    """Create test fixture for genomic dup rse and lse."""
    params = {
        "id": genomic_dup6["id"],
        "type": genomic_dup6["type"],
        "variation": {
            "type": "Text",
            "definition": "NC_000023.11:g.154021812_(154092209_?)dup"
        }
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope='module')
def genomic_del1():
    """Create test fixture containing params for genomic del VD."""
    params = {
        "id": "normalize.variation:NC_000003.12%3Ag.10149811del",
        "type": "VariationDescriptor",
        "variation_id": "",
        "variation": dict(),
        "molecule_context": "genomic",
        "structural_type": "SO:0000159",
        "vrs_ref_allele_seq": "T"
    }
    return params


@pytest.fixture(scope='module')
def genomic_del1_default(genomic_del1):
    """Create a test fixture for genomic del default and LSE."""
    genomic_del1["variation_id"] = "ga4gh:VA.yahcBIB2mO-Vlp3UGssoASAzYy1WY-ad"
    genomic_del1["variation"] = {
        "type": "Allele",
        "location": {
            "sequence_id": "ga4gh:SQ.Zu7h9AggXxhTaGVsy7h_EZSChSZGcmgX",
            "interval": {
                "type": "SimpleInterval",
                "start": 10149810,
                "end": 10149811,
            },
            "type": "SequenceLocation",
        },
        "state": {
            "type": "LiteralSequenceExpression",
            "sequence": ""
        }
    }
    return VariationDescriptor(**genomic_del1)


@pytest.fixture(scope='module')
def genomic_del1_cnv(genomic_del1):
    """Create a test fixture for genomic del CNV."""
    genomic_del1["variation"] = {
        "type": "CopyNumber",
        "subject": {
            "location": {
                "sequence_id": "ga4gh:SQ.Zu7h9AggXxhTaGVsy7h_EZSChSZGcmgX",
                "interval": {
                    "type": "SimpleInterval",
                    "start": 10149810,
                    "end": 10149811,
                },
                "type": "SequenceLocation",
            },
            "reverse_complement": False,
            "type": "DerivedSequenceExpression"
        },
        "copies": {
            "type": "Number",
            "value": 1
        }
    }
    genomic_del1["variation_id"] = "ga4gh:VCN.RcEZ9NrDwmYgL45gZpUABrd4fR6O7zM0"
    return VariationDescriptor(**genomic_del1)


@pytest.fixture(scope='module')
def genomic_del1_rse(genomic_del1):
    """Create a test fixture for genomic del RSE."""
    seq_loc = {
        "sequence_id": "ga4gh:SQ.Zu7h9AggXxhTaGVsy7h_EZSChSZGcmgX",
        "interval": {
            "type": "SimpleInterval",
            "start": 10149810,
            "end": 10149811,
        },
        "type": "SequenceLocation",
    }
    genomic_del1["variation"] = {
        "type": "Allele",
        "location": seq_loc,
        "state": {
            "type": "RepeatedSequenceExpression",
            "seq_expr": {
                "type": "DerivedSequenceExpression",
                "location": seq_loc,
                "reverse_complement": False
            },
            "count": {
                "type": "Number",
                "value": 0
            }
        }
    }
    genomic_del1["variation_id"] = "ga4gh:VA.0zNnXFB8EaWWu7CtbzIbbiedNiovw8CC"
    return VariationDescriptor(**genomic_del1)


@pytest.fixture(scope='module')
def genomic_del2():
    """Create test fixture containing params for genomic del VD."""
    params = {
        "id": "normalize.variation:NC_000003.12%3Ag.10146595_10146613del",
        "type": "VariationDescriptor",
        "variation_id": "",
        "variation": dict(),
        "molecule_context": "genomic",
        "structural_type": "SO:0000159",
        "vrs_ref_allele_seq": "ATGTTGACGGACAGCCTAT"
    }
    return params


@pytest.fixture(scope='module')
def genomic_del2_default(genomic_del2):
    """Create a test fixture for genomic del default and LSE."""
    genomic_del2["variation_id"] = "ga4gh:VA.VGGqV6SjfLH3AfPMA8322Chq_Y4QNp5_"
    genomic_del2["variation"] = {
        "type": "Allele",
        "location": {
            "sequence_id": "ga4gh:SQ.Zu7h9AggXxhTaGVsy7h_EZSChSZGcmgX",
            "interval": {
                "type": "SimpleInterval",
                "start": 10146594,
                "end": 10146613,
            },
            "type": "SequenceLocation",
        },
        "state": {
            "type": "LiteralSequenceExpression",
            "sequence": ""
        }
    }
    return VariationDescriptor(**genomic_del2)


@pytest.fixture(scope='module')
def genomic_del2_cnv(genomic_del2):
    """Create a test fixture for genomic del CNV."""
    genomic_del2["variation"] = {
        "type": "CopyNumber",
        "subject": {
            "location": {
                "sequence_id": "ga4gh:SQ.Zu7h9AggXxhTaGVsy7h_EZSChSZGcmgX",
                "interval": {
                    "type": "SimpleInterval",
                    "start": 10146594,
                    "end": 10146613,
                },
                "type": "SequenceLocation",
            },
            "reverse_complement": False,
            "type": "DerivedSequenceExpression"
        },
        "copies": {
            "type": "Number",
            "value": 1
        }
    }
    genomic_del2["variation_id"] = "ga4gh:VCN.xWlKFaMA2EcgK0yrgib96JJRZUgJoolg"
    return VariationDescriptor(**genomic_del2)


@pytest.fixture(scope='module')
def genomic_del2_rse(genomic_del2):
    """Create a test fixture for genomic del RSE."""
    seq_loc = {
        "sequence_id": "ga4gh:SQ.Zu7h9AggXxhTaGVsy7h_EZSChSZGcmgX",
        "interval": {
            "type": "SimpleInterval",
            "start": 10146594,
            "end": 10146613,
        },
        "type": "SequenceLocation",
    }
    genomic_del2["variation"] = {
        "type": "Allele",
        "location": seq_loc,
        "state": {
            "type": "RepeatedSequenceExpression",
            "seq_expr": {
                "type": "DerivedSequenceExpression",
                "location": seq_loc,
                "reverse_complement": False
            },
            "count": {
                "type": "Number",
                "value": 0
            }
        }
    }
    genomic_del2["variation_id"] = "ga4gh:VA.gmPJJDG_Yhrji3Ie6-WtvI3LvnNWtAN1"
    return VariationDescriptor(**genomic_del2)


@pytest.fixture(scope='module')
def genomic_del3():
    """Create test fixture containing params for genomic del VD."""
    params = {
        "id": "normalize.variation:NC_000023.11%3Ag.%2831060227_31100351%29_%2833274278_33417151%29del",  # noqa: E501
        "type": "VariationDescriptor",
        "variation_id": "",
        "variation": dict(),
        "molecule_context": "genomic",
        "structural_type": "SO:0001743",
        "vrs_ref_allele_seq": None
    }
    return params


@pytest.fixture(scope='module')
def genomic_del3_default(genomic_del3):
    """Create a test fixture for genomic del default and cnv."""
    genomic_del3["variation"] = {
        "type": "CopyNumber",
        "subject": {
            "location": {
                "sequence_id": "ga4gh:SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
                "interval": {
                    "type": "SequenceInterval",
                    "start": {
                        "min": 31060226,
                        "max": 31100350,
                        "type": "DefiniteRange"
                    },
                    "end": {
                        "min": 33274279,
                        "max": 33417152,
                        "type": "DefiniteRange"
                    }
                },
                "type": "SequenceLocation",
            },
            "reverse_complement": False,
            "type": "DerivedSequenceExpression"
        },
        "copies": {
            "type": "DefiniteRange",
            "min": 0,
            "max": 1
        }
    }
    genomic_del3["variation_id"] = "ga4gh:VCN.9h2LkajTwHBdXYMRyrD9HkYwU9d7fIBr"
    return VariationDescriptor(**genomic_del3)


@pytest.fixture(scope='module')
def genomic_del3_rse_lse(genomic_del3):
    """Create test fixture for genomic del rse and lse."""
    params = {
        "id": genomic_del3["id"],
        "type": genomic_del3["type"],
        "variation": {
            "type": "Text",
            "definition": "NC_000023.11:g.(31060227_31100351)_(33274278_33417151)del"  # noqa: E501
        }
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope='module')
def genomic_del4():
    """Create test fixture containing params for genomic del VD."""
    params = {
        "id": "normalize.variation:NC_000023.11%3Ag.%28%3F_31120496%29_%2833339477_%3F%29del",  # noqa: E501
        "type": "VariationDescriptor",
        "variation_id": "",
        "variation": dict(),
        "molecule_context": "genomic",
        "structural_type": "SO:0001743",
        "vrs_ref_allele_seq": None
    }
    return params


@pytest.fixture(scope='module')
def genomic_del4_default(genomic_del4):
    """Create a test fixture for genomic del default and cnv."""
    genomic_del4["variation"] = {
        "type": "CopyNumber",
        "subject": {
            "location": {
                "sequence_id": "ga4gh:SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
                "interval": {
                    "type": "SequenceInterval",
                    "start": {
                        "value": 31120495,
                        "comparator": "<=",
                        "type": "IndefiniteRange"
                    },
                    "end": {
                        "value": 33339477,
                        "comparator": ">=",
                        "type": "IndefiniteRange"
                    }
                },
                "type": "SequenceLocation",
            },
            "reverse_complement": False,
            "type": "DerivedSequenceExpression"
        },
        "copies": {
            "type": "DefiniteRange",
            "min": 0,
            "max": 1
        }
    }
    genomic_del4["variation_id"] = "ga4gh:VCN.yQJnQz12MXlZGWx6BuzccVGrCCic_tMk"
    return VariationDescriptor(**genomic_del4)


@pytest.fixture(scope='module')
def genomic_del4_rse_lse(genomic_del4):
    """Create test fixture for genomic dup rse and lse."""
    params = {
        "id": genomic_del4["id"],
        "type": genomic_del4["type"],
        "variation": {
            "type": "Text",
            "definition": "NC_000023.11:g.(?_31120496)_(33339477_?)del"
        }
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope='module')
def genomic_del5():
    """Create test fixture containing params for genomic del VD."""
    params = {
        "id": "normalize.variation:NC_000023.11%3Ag.%28%3F_18575354%29_18653629del",  # noqa: E501
        "type": "VariationDescriptor",
        "variation_id": "",
        "variation": dict(),
        "molecule_context": "genomic",
        "structural_type": "SO:0001743",
        "vrs_ref_allele_seq": None
    }
    return params


@pytest.fixture(scope='module')
def genomic_del5_default(genomic_del5):
    """Create a test fixture for genomic del default and cnv."""
    genomic_del5["variation"] = {
        "type": "CopyNumber",
        "subject": {
            "location": {
                "sequence_id": "ga4gh:SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
                "interval": {
                    "type": "SequenceInterval",
                    "start": {
                        "value": 18575353,
                        "comparator": "<=",
                        "type": "IndefiniteRange"
                    },
                    "end": {
                        "value": 18653629,
                        "type": "Number"
                    }
                },
                "type": "SequenceLocation",
            },
            "reverse_complement": False,
            "type": "DerivedSequenceExpression"
        },
        "copies": {
            "type": "DefiniteRange",
            "min": 0,
            "max": 1
        }
    }
    genomic_del5["variation_id"] = "ga4gh:VCN._RIw5UC5bZeLeHnBLYAow7Ml-lv2nKJW"
    return VariationDescriptor(**genomic_del5)


@pytest.fixture(scope='module')
def genomic_del5_rse_lse(genomic_del5):
    """Create test fixture for genomic del rse and lse."""
    params = {
        "id": genomic_del5["id"],
        "type": genomic_del5["type"],
        "variation": {
            "type": "Text",
            "definition": "NC_000023.11:g.(?_18575354)_18653629del"
        }
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope='module')
def genomic_del6():
    """Create test fixture containing params for genomic del VD."""
    params = {
        "id": "normalize.variation:NC_000006.12%3Ag.133462764_%28133464858_%3F%29del",  # noqa: E501
        "type": "VariationDescriptor",
        "variation_id": "",
        "variation": dict(),
        "molecule_context": "genomic",
        "structural_type": "SO:0001743",
        "vrs_ref_allele_seq": None
    }
    return params


@pytest.fixture(scope='module')
def genomic_del6_default(genomic_del6):
    """Create a test fixture for genomic del default and cnv."""
    genomic_del6["variation"] = {
        "type": "CopyNumber",
        "subject": {
            "location": {
                "sequence_id": "ga4gh:SQ.0iKlIQk2oZLoeOG9P1riRU6hvL5Ux8TV",
                "interval": {
                    "type": "SequenceInterval",
                    "start": {
                        "value": 133462763,
                        "type": "Number"
                    },
                    "end": {
                        "value": 133464858,
                        "comparator": ">=",
                        "type": "IndefiniteRange"
                    }
                },
                "type": "SequenceLocation",
            },
            "reverse_complement": False,
            "type": "DerivedSequenceExpression"
        },
        "copies": {
            "type": "Number",
            "value": 1
        }
    }
    genomic_del6["variation_id"] = "ga4gh:VCN.F3U6Rmov1WO2mhmRHWumJb-YALOMkeeI"
    return VariationDescriptor(**genomic_del6)


@pytest.fixture(scope='module')
def genomic_del6_rse_lse(genomic_del6):
    """Create test fixture for genomic del rse and lse."""
    params = {
        "id": genomic_del6["id"],
        "type": genomic_del6["type"],
        "variation": {
            "type": "Text",
            "definition": "NC_000006.12:g.133462764_(133464858_?)del"
        }
    }
    return VariationDescriptor(**params)


def test_genomic_dup1(test_normalize, genomic_dup1_default,
                      genomic_dup1_cnv, genomic_dup1_rse):
    """Test that genomic duplication works correctly."""
    q = "NC_000003.12:g.49531262dup"  # 38
    resp = test_normalize.normalize(q, "default")
    assertion_checks(resp, genomic_dup1_default)

    resp = test_normalize.normalize(q, "default")
    assertion_checks(resp, genomic_dup1_default)

    resp = test_normalize.normalize(q, "cnv")
    assertion_checks(resp, genomic_dup1_cnv)

    resp = test_normalize.normalize(q, "repeated_seq_expr")
    assertion_checks(resp, genomic_dup1_rse)

    resp = test_normalize.normalize(q, "literal_seq_expr")
    assertion_checks(resp, genomic_dup1_default)

    q = "NC_000003.11:g.49568695dup"  # 37
    resp = test_normalize.normalize(q, "default")
    assertion_checks(resp, genomic_dup1_default, ignore_id=True)

    resp = test_normalize.normalize(q, "default")
    assertion_checks(resp, genomic_dup1_default, ignore_id=True)

    resp = test_normalize.normalize(q, "cnv")
    assertion_checks(resp, genomic_dup1_cnv, ignore_id=True)

    resp = test_normalize.normalize(q, "repeated_seq_expr")
    assertion_checks(resp, genomic_dup1_rse, ignore_id=True)

    resp = test_normalize.normalize(q, "literal_seq_expr")
    assertion_checks(resp, genomic_dup1_default, ignore_id=True)


def test_genomic_dup2(test_normalize, genomic_dup2_default, genomic_dup2_cnv,
                      genomic_dup2_rse):
    """Test that genomic duplication works correctly."""
    q = "NC_000016.10:g.2087938_2087948dup"  # 38
    resp = test_normalize.normalize(q, "default")
    assertion_checks(resp, genomic_dup2_default)

    resp = test_normalize.normalize(q, "cnv")
    assertion_checks(resp, genomic_dup2_cnv)

    resp = test_normalize.normalize(q, "repeated_seq_expr")
    assertion_checks(resp, genomic_dup2_rse)

    resp = test_normalize.normalize(q, "literal_seq_expr")
    assertion_checks(resp, genomic_dup2_default)

    q = "NC_000016.9:g.2137939_2137949dup"  # 37
    resp = test_normalize.normalize(q, "default")
    assertion_checks(resp, genomic_dup2_default, ignore_id=True)

    resp = test_normalize.normalize(q, "cnv")
    assertion_checks(resp, genomic_dup2_cnv, ignore_id=True)

    resp = test_normalize.normalize(q, "repeated_seq_expr")
    assertion_checks(resp, genomic_dup2_rse, ignore_id=True)

    resp = test_normalize.normalize(q, "literal_seq_expr")
    assertion_checks(resp, genomic_dup2_default, ignore_id=True)


def test_genomic_dup3(test_normalize, genomic_dup3_default,
                      genomic_dup3_rse_lse):
    """Test that genomic duplication works correctly."""
    q = "NC_000023.11:g.(31060227_31100351)_(33274278_33417151)dup"  # 38
    resp = test_normalize.normalize(q, "default")
    assertion_checks(resp, genomic_dup3_default)

    resp = test_normalize.normalize(q, "cnv")
    assertion_checks(resp, genomic_dup3_default)

    resp = test_normalize.normalize(q, "repeated_seq_expr")
    assertion_checks(resp, genomic_dup3_rse_lse)

    resp = test_normalize.normalize(q, "literal_seq_expr")
    assertion_checks(resp, genomic_dup3_rse_lse)

    q = "NC_000023.10:g.(31078344_31118468)_(33292395_33435268)dup"  # 37
    resp = test_normalize.normalize(q, "default")
    assertion_checks(resp, genomic_dup3_default, ignore_id=True)

    resp = test_normalize.normalize(q, "cnv")
    assertion_checks(resp, genomic_dup3_default, ignore_id=True)

    genomic_dup3_rse_lse.variation.definition = q
    resp = test_normalize.normalize(q, "repeated_seq_expr")
    assertion_checks(resp, genomic_dup3_rse_lse, ignore_id=True)

    resp = test_normalize.normalize(q, "literal_seq_expr")
    assertion_checks(resp, genomic_dup3_rse_lse, ignore_id=True)


def test_genomic_dup4(test_normalize, genomic_dup4_default,
                      genomic_dup4_rse_lse):
    """Test that genomic duplication works correctly."""
    q = "NC_000020.11:g.(?_30417576)_(31394018_?)dup"  # 38
    resp = test_normalize.normalize(q, "default")
    assertion_checks(resp, genomic_dup4_default)

    resp = test_normalize.normalize(q, "cnv")
    assertion_checks(resp, genomic_dup4_default)

    resp = test_normalize.normalize(q, "repeated_seq_expr")
    assertion_checks(resp, genomic_dup4_rse_lse)

    resp = test_normalize.normalize(q, "literal_seq_expr")
    assertion_checks(resp, genomic_dup4_rse_lse)

    q = "NC_000020.10:g.(?_29652252)_(29981821_?)dup"  # 37
    resp = test_normalize.normalize(q, "default")
    assertion_checks(resp, genomic_dup4_default, ignore_id=True)

    resp = test_normalize.normalize(q, "cnv")
    assertion_checks(resp, genomic_dup4_default, ignore_id=True)

    genomic_dup4_rse_lse.variation.definition = q
    resp = test_normalize.normalize(q, "repeated_seq_expr")
    assertion_checks(resp, genomic_dup4_rse_lse, ignore_id=True)

    resp = test_normalize.normalize(q, "literal_seq_expr")
    assertion_checks(resp, genomic_dup4_rse_lse, ignore_id=True)


def test_genomic_dup5(test_normalize, genomic_dup5_default,
                      genomic_dup5_rse_lse):
    """Test that genomic duplication works correctly."""
    q = "NC_000023.11:g.(?_154021812)_154092209dup"  # 38
    resp = test_normalize.normalize(q, "default")
    assertion_checks(resp, genomic_dup5_default)

    resp = test_normalize.normalize(q, "cnv")
    assertion_checks(resp, genomic_dup5_default)

    resp = test_normalize.normalize(q, "repeated_seq_expr")
    assertion_checks(resp, genomic_dup5_rse_lse)

    resp = test_normalize.normalize(q, "literal_seq_expr")
    assertion_checks(resp, genomic_dup5_rse_lse)

    q = "NC_000023.10:g.(?_153287263)_153357667dup"  # 38
    resp = test_normalize.normalize(q, "default")
    assertion_checks(resp, genomic_dup5_default, ignore_id=True)

    resp = test_normalize.normalize(q, "cnv")
    assertion_checks(resp, genomic_dup5_default, ignore_id=True)

    genomic_dup5_rse_lse.variation.definition = q
    resp = test_normalize.normalize(q, "repeated_seq_expr")
    assertion_checks(resp, genomic_dup5_rse_lse, ignore_id=True)

    resp = test_normalize.normalize(q, "literal_seq_expr")
    assertion_checks(resp, genomic_dup5_rse_lse, ignore_id=True)


def test_genomic_dup6(test_normalize, genomic_dup6_default,
                      genomic_dup6_rse_lse):
    """Test that genomic duplication works correctly."""
    q = "NC_000023.11:g.154021812_(154092209_?)dup"  # 38
    resp = test_normalize.normalize(q, "default")
    assertion_checks(resp, genomic_dup6_default)

    resp = test_normalize.normalize(q, "cnv")
    assertion_checks(resp, genomic_dup6_default)

    resp = test_normalize.normalize(q, "repeated_seq_expr")
    assertion_checks(resp, genomic_dup6_rse_lse)

    resp = test_normalize.normalize(q, "literal_seq_expr")
    assertion_checks(resp, genomic_dup6_rse_lse)

    q = "NC_000023.10:g.153287263_(153357667_?)dup"  # 37
    resp = test_normalize.normalize(q, "default")
    assertion_checks(resp, genomic_dup6_default, ignore_id=True)

    resp = test_normalize.normalize(q, "cnv")
    assertion_checks(resp, genomic_dup6_default, ignore_id=True)

    genomic_dup6_rse_lse.variation.definition = q
    resp = test_normalize.normalize(q, "repeated_seq_expr")
    assertion_checks(resp, genomic_dup6_rse_lse, ignore_id=True)

    resp = test_normalize.normalize(q, "literal_seq_expr")
    assertion_checks(resp, genomic_dup6_rse_lse, ignore_id=True)


def test_genomic_del1(test_normalize, genomic_del1_default, genomic_del1_cnv,
                      genomic_del1_rse):
    """Test that genomic deletion works correctly."""
    q = "NC_000003.12:g.10149811del"  # 38
    resp = test_normalize.normalize(q, "default")
    assertion_checks(resp, genomic_del1_default)

    resp = test_normalize.normalize(q, "cnv")
    assertion_checks(resp, genomic_del1_cnv)

    resp = test_normalize.normalize(q, "repeated_seq_expr")
    assertion_checks(resp, genomic_del1_rse)

    resp = test_normalize.normalize(q, "literal_seq_expr")
    assertion_checks(resp, genomic_del1_default)

    q = "NC_000003.11:g.10191495del"  # 37
    resp = test_normalize.normalize(q, "default")
    assertion_checks(resp, genomic_del1_default, ignore_id=True)

    resp = test_normalize.normalize(q, "cnv")
    assertion_checks(resp, genomic_del1_cnv, ignore_id=True)

    resp = test_normalize.normalize(q, "repeated_seq_expr")
    assertion_checks(resp, genomic_del1_rse, ignore_id=True)

    resp = test_normalize.normalize(q, "literal_seq_expr")
    assertion_checks(resp, genomic_del1_default, ignore_id=True)


def test_genomic_del2(test_normalize, genomic_del2_default, genomic_del2_cnv,
                      genomic_del2_rse):
    """Test that genomic deletion works correctly."""
    q = "NC_000003.12:g.10146595_10146613del"  # 38
    resp = test_normalize.normalize(q, "default")
    assertion_checks(resp, genomic_del2_default)

    resp = test_normalize.normalize(q, "cnv")
    assertion_checks(resp, genomic_del2_cnv)

    resp = test_normalize.normalize(q, "repeated_seq_expr")
    assertion_checks(resp, genomic_del2_rse)

    resp = test_normalize.normalize(q, "literal_seq_expr")
    assertion_checks(resp, genomic_del2_default)

    q = "NC_000003.11:g.10188279_10188297del"  # 37
    resp = test_normalize.normalize(q, "default")
    assertion_checks(resp, genomic_del2_default, ignore_id=True)

    resp = test_normalize.normalize(q, "cnv")
    assertion_checks(resp, genomic_del2_cnv, ignore_id=True)

    resp = test_normalize.normalize(q, "repeated_seq_expr")
    assertion_checks(resp, genomic_del2_rse, ignore_id=True)

    resp = test_normalize.normalize(q, "literal_seq_expr")
    assertion_checks(resp, genomic_del2_default, ignore_id=True)


def test_genomic_del3(test_normalize, genomic_del3_default,
                      genomic_del3_rse_lse):
    """Test that genomic deletion works correctly."""
    q = "NC_000023.11:g.(31060227_31100351)_(33274278_33417151)del"  # 38
    resp = test_normalize.normalize(q, "default")
    assertion_checks(resp, genomic_del3_default)

    resp = test_normalize.normalize(q, "cnv")
    assertion_checks(resp, genomic_del3_default)

    resp = test_normalize.normalize(q, "repeated_seq_expr")
    assertion_checks(resp, genomic_del3_rse_lse)

    resp = test_normalize.normalize(q, "literal_seq_expr")
    assertion_checks(resp, genomic_del3_rse_lse)

    q = "NC_000023.10:g.(31078344_31118468)_(33292395_33435268)del"  # 37
    resp = test_normalize.normalize(q, "default")
    assertion_checks(resp, genomic_del3_default, ignore_id=True)

    resp = test_normalize.normalize(q, "cnv")
    assertion_checks(resp, genomic_del3_default, ignore_id=True)

    genomic_del3_rse_lse.variation.definition = q
    resp = test_normalize.normalize(q, "repeated_seq_expr")
    assertion_checks(resp, genomic_del3_rse_lse, ignore_id=True)

    resp = test_normalize.normalize(q, "literal_seq_expr")
    assertion_checks(resp, genomic_del3_rse_lse, ignore_id=True)


def test_genomic_del4(test_normalize, genomic_del4_default,
                      genomic_del4_rse_lse):
    """Test that genomic deletion works correctly."""
    q = "NC_000023.11:g.(?_31120496)_(33339477_?)del"  # 38
    resp = test_normalize.normalize(q, "default")
    assertion_checks(resp, genomic_del4_default)

    resp = test_normalize.normalize(q, "cnv")
    assertion_checks(resp, genomic_del4_default)

    resp = test_normalize.normalize(q, "repeated_seq_expr")
    assertion_checks(resp, genomic_del4_rse_lse)

    resp = test_normalize.normalize(q, "literal_seq_expr")
    assertion_checks(resp, genomic_del4_rse_lse)

    q = "NC_000023.10:g.(?_31138613)_(33357594_?)del"  # 37
    resp = test_normalize.normalize(q, "default")
    assertion_checks(resp, genomic_del4_default, ignore_id=True)

    resp = test_normalize.normalize(q, "cnv")
    assertion_checks(resp, genomic_del4_default, ignore_id=True)

    genomic_del4_rse_lse.variation.definition = q
    resp = test_normalize.normalize(q, "repeated_seq_expr")
    assertion_checks(resp, genomic_del4_rse_lse, ignore_id=True)

    resp = test_normalize.normalize(q, "literal_seq_expr")
    assertion_checks(resp, genomic_del4_rse_lse, ignore_id=True)


def test_genomic_del5(test_normalize, genomic_del5_default,
                      genomic_del5_rse_lse):
    """Test that genomic deletion works correctly."""
    q = "NC_000023.11:g.(?_18575354)_18653629del"  # 38
    resp = test_normalize.normalize(q, "default")
    assertion_checks(resp, genomic_del5_default)

    resp = test_normalize.normalize(q, "cnv")
    assertion_checks(resp, genomic_del5_default)

    resp = test_normalize.normalize(q, "repeated_seq_expr")
    assertion_checks(resp, genomic_del5_rse_lse)

    resp = test_normalize.normalize(q, "literal_seq_expr")
    assertion_checks(resp, genomic_del5_rse_lse)

    q = "NC_000023.10:g.(?_18593474)_18671749del"  # 37
    resp = test_normalize.normalize(q, "default")
    assertion_checks(resp, genomic_del5_default, ignore_id=True)

    resp = test_normalize.normalize(q, "cnv")
    assertion_checks(resp, genomic_del5_default, ignore_id=True)

    genomic_del5_rse_lse.variation.definition = q
    resp = test_normalize.normalize(q, "repeated_seq_expr")
    assertion_checks(resp, genomic_del5_rse_lse, ignore_id=True)

    resp = test_normalize.normalize(q, "literal_seq_expr")
    assertion_checks(resp, genomic_del5_rse_lse, ignore_id=True)


def test_genomic_del6(test_normalize, genomic_del6_default,
                      genomic_del6_rse_lse):
    """Test that genomic deletion works correctly."""
    q = "NC_000006.12:g.133462764_(133464858_?)del"  # 38
    resp = test_normalize.normalize(q, "default")
    assertion_checks(resp, genomic_del6_default)

    resp = test_normalize.normalize(q, "cnv")
    assertion_checks(resp, genomic_del6_default)

    resp = test_normalize.normalize(q, "repeated_seq_expr")
    assertion_checks(resp, genomic_del6_rse_lse)

    resp = test_normalize.normalize(q, "literal_seq_expr")
    assertion_checks(resp, genomic_del6_rse_lse)

    q = "NC_000006.11:g.133783902_(133785996_?)del"  # 37
    resp = test_normalize.normalize(q, "default")
    assertion_checks(resp, genomic_del6_default, ignore_id=True)

    resp = test_normalize.normalize(q, "cnv")
    assertion_checks(resp, genomic_del6_default, ignore_id=True)

    genomic_del6_rse_lse.variation.definition = q
    resp = test_normalize.normalize(q, "repeated_seq_expr")
    assertion_checks(resp, genomic_del6_rse_lse, ignore_id=True)

    resp = test_normalize.normalize(q, "literal_seq_expr")
    assertion_checks(resp, genomic_del6_rse_lse, ignore_id=True)


def test_parameters(test_normalize):
    """Check that valid and invalid parameters work as intended."""
    q = "NC_000003.12:g.49531262dup"
    warnings = ["hgvs_dup_del_mode must be one of: ['default', 'cnv', "
                "'repeated_seq_expr', 'literal_seq_expr']"]
    resp = test_normalize.normalize(q, "copy_number")
    assert resp is None
    assert test_normalize.warnings == warnings

    resp = test_normalize.normalize(q, "repeated_seq_exprs")
    assert resp is None
    assert test_normalize.warnings == warnings

    warnings = ["hgvs_dup_del_mode cannot be None"]
    resp = test_normalize.normalize(q, '')
    assert resp is None
    assert test_normalize.warnings == warnings

    resp = test_normalize.normalize(q, None)
    assert resp is None
    assert test_normalize.warnings == warnings

    resp = test_normalize.normalize(q, " CnV ")
    assert resp is not None
    assert test_normalize.warnings == []
