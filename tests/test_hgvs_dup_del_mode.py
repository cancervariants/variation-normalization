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

        def normalize(self, q, hgvs_dup_del_mode="default"):
            return self.query_handler.normalize(
                q, hgvs_dup_del_mode=hgvs_dup_del_mode)
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


def test_genomic_dup1(test_normalize, genomic_dup1_default,
                      genomic_dup1_cnv, genomic_dup1_rse):
    """Test that genomic duplication works correctly."""
    q = "NC_000003.12:g.49531262dup"
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


def test_genomic_dup2(test_normalize, genomic_dup2_default, genomic_dup2_cnv,
                      genomic_dup2_rse):
    """Test that genomic duplication works correctly."""
    q = "NC_000016.10:g.2087938_2087948dup"
    resp = test_normalize.normalize(q, "default")
    assertion_checks(resp, genomic_dup2_default)

    resp = test_normalize.normalize(q, "cnv")
    assertion_checks(resp, genomic_dup2_cnv)

    resp = test_normalize.normalize(q, "repeated_seq_expr")
    assertion_checks(resp, genomic_dup2_rse)

    resp = test_normalize.normalize(q, "literal_seq_expr")
    assertion_checks(resp, genomic_dup2_default)


def test_genomic_dup3(test_normalize, genomic_dup3_default,
                      genomic_dup3_rse_lse):
    """Test that genomic duplication works correctly."""
    q = "NC_000023.11:g.(31060227_31100351)_(33274278_33417151)dup"
    resp = test_normalize.normalize(q, "default")
    assertion_checks(resp, genomic_dup3_default)

    resp = test_normalize.normalize(q, "cnv")
    assertion_checks(resp, genomic_dup3_default)

    resp = test_normalize.normalize(q, "repeated_seq_expr")
    assertion_checks(resp, genomic_dup3_rse_lse)

    resp = test_normalize.normalize(q, "literal_seq_expr")
    assertion_checks(resp, genomic_dup3_rse_lse)


def test_genomic_dup4(test_normalize, genomic_dup4_default,
                      genomic_dup4_rse_lse):
    """Test that genomic duplication works correctly."""
    q = "NC_000020.11:g.(?_30417576)_(31394018_?)dup"
    resp = test_normalize.normalize(q, "default")
    assertion_checks(resp, genomic_dup4_default)

    resp = test_normalize.normalize(q, "cnv")
    assertion_checks(resp, genomic_dup4_default)

    resp = test_normalize.normalize(q, "repeated_seq_expr")
    assertion_checks(resp, genomic_dup4_rse_lse)

    resp = test_normalize.normalize(q, "literal_seq_expr")
    assertion_checks(resp, genomic_dup4_rse_lse)


# def test_genomic_dup5(test_normalize):
#     """Test that genomic duplication works correctly."""
#     q = "NC_000023.11:g.(?_154021812)_154092209dup"
#
#
# def test_genomic_dup6(test_normalize):
#     """Test that genomic duplication works correctly."""
#     q = "NC_000023.11:g.154021812_(154092209_?)dup"
