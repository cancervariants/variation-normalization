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


def test_genomic_duplication(test_normalize, genomic_dup1_default,
                             genomic_dup1_cnv, genomic_dup1_rse):
    """Test that genomic duplication works correctly."""
    resp = test_normalize.normalize("NC_000003.12:g.49531262dup", "default")
    assertion_checks(resp, genomic_dup1_default)

    resp = test_normalize.normalize("NC_000003.12:g.49531262dup", "cnv")
    assertion_checks(resp, genomic_dup1_cnv)

    resp = test_normalize.normalize("NC_000003.12:g.49531262dup",
                                    "repeated_seq_expr")
    assertion_checks(resp, genomic_dup1_rse)

    resp = test_normalize.normalize("NC_000003.12:g.49531262dup",
                                    "literal_seq_expr")
    assertion_checks(resp, genomic_dup1_default)
