"""Module for testing the hgvs to absolute and relative copy number endpoints"""
import copy

import pytest


@pytest.fixture(scope="module")
def test_handler(test_query_handler):
    """Create test fixture for to copy number handler"""
    return test_query_handler.to_copy_number_handler


@pytest.fixture(scope="module")
def genomic_dup1_rel_38(genomic_dup1_seq_loc):
    """Create test fixture relative copy number variation"""
    return {
        "type": "RelativeCopyNumber",
        "_id": "ga4gh:VRC.4ClST-L7j5CLv0buuOS5p3Yvp8GCByjF",
        "subject": genomic_dup1_seq_loc,
        "relative_copy_class": "complete loss"
    }


@pytest.fixture(scope="module")
def genomic_dup1_37_loc():
    """Create test fixture GRCh37 duplication subject"""
    return {
        "_id": "ga4gh:VSL.CXcLL6RUPkro3dLXN0miGEzlzPYiqw2q",
        "sequence_id": "ga4gh:SQ.VNBualIltAyi2AI_uXcKU7M9XUOuA7MS",
        "interval": {
            "type": "SequenceInterval",
            "start": {"value": 49568693, "type": "Number"},
            "end": {"value": 49568695, "type": "Number"},
        },
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_dup1_abs_37(genomic_dup1_37_loc):
    """Create test fixture absolute copy number variation"""
    return {
        "type": "AbsoluteCopyNumber",
        "_id": "ga4gh:VAC.hGeC0s0RZSX1FEb7ggN6KhAjZUyVF5mU",
        "subject": genomic_dup1_37_loc,
        "copies": {"type": "Number", "value": 3}
    }


@pytest.fixture(scope="module")
def genomic_dup1_rel_37(genomic_dup1_37_loc):
    """Create test fixture relative copy number variation"""
    return {
        "type": "RelativeCopyNumber",
        "_id": "ga4gh:VRC.OqSxkdmM95Dq-V-ZqMe7wXFXO9T0-8gP",
        "subject": genomic_dup1_37_loc,
        "relative_copy_class": "complete loss"
    }


@pytest.fixture(scope="module")
def genomic_dup2_rel_38(genomic_dup2_seq_loc):
    """Create test fixture relative copy number variation"""
    return {
        "type": "RelativeCopyNumber",
        "_id": "ga4gh:VRC.WLuxOcKidzwkLQ4Z5AyH7bfZiSg-srFw",
        "subject": genomic_dup2_seq_loc,
        "relative_copy_class": "partial loss"
    }


@pytest.fixture(scope="module")
def genomic_dup2_37_loc():
    """Create test fixture GRCh37 duplication subject"""
    return {
        "_id": "ga4gh:VSL.VY7qo3C2Nl-bpRtjOQ6MODJ2FfGUJRV1",
        "sequence_id": "ga4gh:SQ.W6wLoIFOn4G7cjopxPxYNk2lcEqhLQFb",
        "interval": {
            "type": "SequenceInterval",
            "start": {"value": 2137938, "type": "Number"},
            "end": {"value": 2137949, "type": "Number"},
        },
        "type": "SequenceLocation"
    }


@pytest.fixture(scope="module")
def genomic_dup2_abs_37(genomic_dup2_37_loc):
    """Create test fixture absolute copy number variation"""
    return {
        "type": "AbsoluteCopyNumber",
        "_id": "ga4gh:VAC.L1_P5Tf41A-b29DN9Jg5T-c33iAhW1A8",
        "subject": genomic_dup2_37_loc,
        "copies": {"type": "Number", "value": 3}
    }


@pytest.fixture(scope="module")
def genomic_dup2_rel_37(genomic_dup2_37_loc):
    """Create test fixture relative copy number variation"""
    return {
        "type": "RelativeCopyNumber",
        "_id": "ga4gh:VRC.zkTN-KUrBo_Bc3l5Wfczpb0lOktb77qE",
        "subject": genomic_dup2_37_loc,
        "relative_copy_class": "partial loss"
    }


@pytest.fixture(scope="module")
def genomic_dup3_abs_38(genomic_del3_dup3_loc):
    """Create test fixture absolute copy number variation"""
    return {
        "type": "AbsoluteCopyNumber",
        "_id": "ga4gh:VAC.cQATJ6a1uGwXOHu-advv8lRsMgjNLKul",
        "subject": genomic_del3_dup3_loc,
        "copies": {"type": "Number", "value": 2}
    }


@pytest.fixture(scope="module")
def genomic_dup3_rel_38(genomic_del3_dup3_loc):
    """Create test fixture relative copy number variation"""
    return {
        "type": "RelativeCopyNumber",
        "_id": "ga4gh:VRC.dAPJs31x9cI8bL3yJP-8BMOKGMMcCCbH",
        "subject": genomic_del3_dup3_loc,
        "relative_copy_class": "high-level gain"
    }


@pytest.fixture(scope="module")
def genomic_dup3_37_loc():
    """Create test fixture GRCh37 duplication subject"""
    return {
        "_id": "ga4gh:VSL.hDvTkhN210ZKyZOhilS7FxaKUCJgj6JC",
        "sequence_id": "ga4gh:SQ.v7noePfnNpK8ghYXEqZ9NukMXW7YeNsm",
        "interval": {
            "type": "SequenceInterval",
            "start": {
                "min": 31078343,
                "max": 31118467,
                "type": "DefiniteRange"
            },
            "end": {
                "min": 33292396,
                "max": 33435269,
                "type": "DefiniteRange"
            }
        },
        "type": "SequenceLocation"
    }


@pytest.fixture(scope="module")
def genomic_dup3_abs_37(genomic_dup3_37_loc):
    """Create test fixture absolute copy number variation"""
    return {
        "type": "AbsoluteCopyNumber",
        "_id": "ga4gh:VAC.Pv9I4Dqk69w-tX0axaikVqid-pozxU74",
        "subject": genomic_dup3_37_loc,
        "copies": {"type": "Number", "value": 2}
    }


@pytest.fixture(scope="module")
def genomic_dup3_rel_37(genomic_dup3_37_loc):
    """Create test fixture relative copy number variation"""
    return {
        "type": "RelativeCopyNumber",
        "_id": "ga4gh:VRC.4nV0gzcAm_tPMbfkrqtezaA61uZUEWzt",
        "subject": genomic_dup3_37_loc,
        "relative_copy_class": "high-level gain"
    }


@pytest.fixture(scope="module")
def genomic_dup4_abs_38(genoimc_dup4_loc):
    """Create test fixture absolute copy number variation"""
    return {
        "type": "AbsoluteCopyNumber",
        "_id": "ga4gh:VAC.h594XLS8a4VA6j-ghLaghqXmof8hmF5z",
        "subject": genoimc_dup4_loc,
        "copies": {"type": "Number", "value": 3}
    }


@pytest.fixture(scope="module")
def genomic_dup4_rel_38(genoimc_dup4_loc):
    """Create test fixture relative copy number variation"""
    return {
        "type": "RelativeCopyNumber",
        "_id": "ga4gh:VRC.99KSUw5LDOqDK2ETM81En2fCYEFCfIiL",
        "subject": genoimc_dup4_loc,
        "relative_copy_class": "complete loss"
    }


@pytest.fixture(scope="module")
def genomic_dup4_37_loc():
    """Create test fixture GRCh37 duplication subject"""
    return {
        "_id": "ga4gh:VSL.ypfdWazLg2-FkZhRH0EQDGwH5IqlAyLh",
        "sequence_id": "ga4gh:SQ.iy_UbUrvECxFRX5LPTH_KPojdlT7BKsf",
        "interval": {
            "type": "SequenceInterval",
            "start": {
                "value": 29652251,
                "comparator": "<=",
                "type": "IndefiniteRange"
            },
            "end": {
                "value": 29981821,
                "comparator": ">=",
                "type": "IndefiniteRange"
            }
        },
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_dup4_abs_37(genomic_dup4_37_loc):
    """Create test fixture relative copy number variation"""
    return {
        "type": "AbsoluteCopyNumber",
        "_id": "ga4gh:VAC.CYSRgw_prwhAhTZGM9blaEvmDjj952Uf",
        "subject": genomic_dup4_37_loc,
        "copies": {"type": "Number", "value": 3}
    }


@pytest.fixture(scope="module")
def genomic_dup4_rel_37(genomic_dup4_37_loc):
    """Create test fixture relative copy number variation"""
    return {
        "type": "RelativeCopyNumber",
        "_id": "ga4gh:VRC.X4mJ6dPachx66ssyvWpLraCAYe04YPO4",
        "subject": genomic_dup4_37_loc,
        "relative_copy_class": "complete loss"
    }


@pytest.fixture(scope="module")
def genomic_dup5_abs_38(genomic_dup5_loc):
    """Create test fixture absolute copy number variation"""
    return {
        "type": "AbsoluteCopyNumber",
        "_id": "ga4gh:VAC.BUEI9XPTvjBvNUoREsXRsm8THNuR5Fe7",
        "subject": genomic_dup5_loc,
        "copies": {"type": "Number", "value": 4}
    }


@pytest.fixture(scope="module")
def genomic_dup5_rel_38(genomic_dup5_loc):
    """Create test fixture relative copy number variation"""
    return {
        "type": "RelativeCopyNumber",
        "_id": "ga4gh:VRC.vy8SSVFuaeZTkUCCv6izNCkF0zgbBG7G",
        "subject": genomic_dup5_loc,
        "relative_copy_class": "partial loss"
    }


@pytest.fixture(scope="module")
def genomic_dup5_37_loc():
    """Create test fixture GRCh37 duplication subject"""
    return {
        "_id": "ga4gh:VSL.5pVsde3FHXnOaO0jONcisMHAfcjqmXiD",
        "sequence_id": "ga4gh:SQ.v7noePfnNpK8ghYXEqZ9NukMXW7YeNsm",
        "interval": {
            "type": "SequenceInterval",
            "start": {
                "value": 153287262,
                "comparator": "<=",
                "type": "IndefiniteRange"
            },
            "end": {
                "value": 153357667,
                "type": "Number"
            }
        },
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_dup5_abs_37(genomic_dup5_37_loc):
    """Create test fixture absolute copy number variation"""
    return {
        "type": "AbsoluteCopyNumber",
        "_id": "ga4gh:VAC.1-pcrgINIDzVXrTgs7xshzBQVlhQ_dX8",
        "subject": genomic_dup5_37_loc,
        "copies": {"type": "Number", "value": 4}
    }


@pytest.fixture(scope="module")
def genomic_dup5_rel_37(genomic_dup5_37_loc):
    """Create test fixture relative copy number variation"""
    return {
        "type": "RelativeCopyNumber",
        "_id": "ga4gh:VRC.v0u8jwoIYvE03wW7jfzw9-O-8N-LQp_N",
        "subject": genomic_dup5_37_loc,
        "relative_copy_class": "partial loss"
    }


@pytest.fixture(scope="module")
def genomic_dup6_abs_38(genoimc_dup6_loc):
    """Create test fixture absolute copy number variation"""
    return {
        "type": "AbsoluteCopyNumber",
        "_id": "ga4gh:VAC.ZkgR6TD7VypzVrLAYFnSb-D7DXp62Yfn",
        "subject": genoimc_dup6_loc,
        "copies": {"type": "Number", "value": 2}
    }


@pytest.fixture(scope="module")
def genomic_dup6_rel_38(genoimc_dup6_loc):
    """Create test fixture relative copy number variation"""
    return {
        "type": "RelativeCopyNumber",
        "_id": "ga4gh:VRC.dUyDSJDl72tWZKG2AcQEe83mIW3RCcQE",
        "subject": genoimc_dup6_loc,
        "relative_copy_class": "copy neutral"
    }


@pytest.fixture(scope="module")
def genomic_dup6_37_loc():
    """Create test fixture GRCh37 duplication subject"""
    return {
        "_id": "ga4gh:VSL.4FthtlBwSSEa3O0CUPw-DnKYz4Jnfa3n",
        "sequence_id": "ga4gh:SQ.v7noePfnNpK8ghYXEqZ9NukMXW7YeNsm",
        "interval": {
            "type": "SequenceInterval",
            "start": {
                "value": 153287262,
                "type": "Number"
            },
            "end": {
                "value": 153357667,
                "comparator": ">=",
                "type": "IndefiniteRange"
            }
        },
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_dup6_abs_37(genomic_dup6_37_loc):
    """Create test fixture absolute copy number variation"""
    return {
        "type": "AbsoluteCopyNumber",
        "_id": "ga4gh:VAC.8u4HiQWnHhwGZtIBF_DEYQCNdaFKOHuN",
        "subject": genomic_dup6_37_loc,
        "copies": {"type": "Number", "value": 2}
    }


@pytest.fixture(scope="module")
def genomic_dup6_rel_37(genomic_dup6_37_loc):
    """Create test fixture relative copy number variation"""
    return {
        "type": "RelativeCopyNumber",
        "_id": "ga4gh:VRC._Dy2mGmY7n0gYqFcnZqmgTCdlcChFUnC",
        "subject": genomic_dup6_37_loc,
        "relative_copy_class": "copy neutral"
    }


@pytest.fixture(scope="module")
def genomic_del1_rel_38(genomic_del1_seq_loc):
    """Create test fixture relative copy number variation"""
    return {
        "type": "RelativeCopyNumber",
        "_id": "ga4gh:VRC.FHImCphKfSBeobf9HO6qpu_Bm5U9VfHz",
        "subject": genomic_del1_seq_loc,
        "relative_copy_class": "copy neutral"
    }


@pytest.fixture(scope="module")
def genomic_del1_37_loc():
    """Create test fixture GRCh37 deletion subject"""
    return {
        "_id": "ga4gh:VSL.IpNE_OIaSJSeiiJTEwrxyiRsgpUay8_n",
        "sequence_id": "ga4gh:SQ.VNBualIltAyi2AI_uXcKU7M9XUOuA7MS",
        "interval": {
            "type": "SequenceInterval",
            "start": {"value": 10191494, "type": "Number"},
            "end": {"value": 10191495, "type": "Number"},
        },
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_del1_abs_37(genomic_del1_37_loc):
    """Create test fixture absolute copy number variation"""
    return {
        "type": "AbsoluteCopyNumber",
        "_id": "ga4gh:VAC.u4Dcmgh8rV8WxRxdrd4W6qA0bTCkS8KP",
        "subject": genomic_del1_37_loc,
        "copies": {"type": "Number", "value": 1}
    }


@pytest.fixture(scope="module")
def genomic_del1_rel_37(genomic_del1_37_loc):
    """Create test fixture relative copy number variation"""
    return {
        "type": "RelativeCopyNumber",
        "_id": "ga4gh:VRC.2_BewkjhJtvj8J814SE0RVxp-OdYxJSo",
        "subject": genomic_del1_37_loc,
        "relative_copy_class": "copy neutral"
    }


@pytest.fixture(scope="module")
def genomic_del2_rel_38(genomic_del2_seq_loc):
    """Create test fixture relative copy number variation"""
    return {
        "type": "RelativeCopyNumber",
        "_id": "ga4gh:VRC.xyQ8SLzsBZmeo5M3Aoii-3iED0CeB4m5",
        "subject": genomic_del2_seq_loc,
        "relative_copy_class": "low-level gain"
    }


@pytest.fixture(scope="module")
def genomic_del2_37_loc():
    """Create test fixture GRCh37 deletion subject"""
    return {
        "_id": "ga4gh:VSL.B-feFlnwavUnci2jFAqG8WOQrqFAZJ6N",
        "sequence_id": "ga4gh:SQ.VNBualIltAyi2AI_uXcKU7M9XUOuA7MS",
        "interval": {
            "type": "SequenceInterval",
            "start": {"value": 10188278, "type": "Number"},
            "end": {"value": 10188297, "type": "Number"},
        },
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_del2_abs_37(genomic_del2_37_loc):
    """Create test fixture absolute copy number variation"""
    return {
        "type": "AbsoluteCopyNumber",
        "_id": "ga4gh:VAC.wWy6pg25aZzAuDwEAcO18yD0QuSArYlE",
        "subject": genomic_del2_37_loc,
        "copies": {"type": "Number", "value": 1}
    }


@pytest.fixture(scope="module")
def genomic_del2_rel_37(genomic_del2_37_loc):
    """Create test fixture relative copy number variation"""
    return {
        "type": "RelativeCopyNumber",
        "_id": "ga4gh:VRC._iD9fTzI_KBtCOlI_6izwiUzjyoEZOft",
        "subject": genomic_del2_37_loc,
        "relative_copy_class": "low-level gain"
    }


@pytest.fixture(scope="module")
def genomic_del3_abs_38(genomic_del3_dup3_loc):
    """Create test fixture absolute copy number variation"""
    return {
        "type": "AbsoluteCopyNumber",
        "_id": "ga4gh:VAC.cQATJ6a1uGwXOHu-advv8lRsMgjNLKul",
        "subject": genomic_del3_dup3_loc,
        "copies": {"type": "Number", "value": 2}
    }


@pytest.fixture(scope="module")
def genomic_del3_rel_38(genomic_del3_dup3_loc):
    """Create test fixture relative copy number variation"""
    return {
        "type": "RelativeCopyNumber",
        "_id": "ga4gh:VRC.14ULmT3TSo6-du1BQehVkX14OH_92j4F",
        "subject": genomic_del3_dup3_loc,
        "relative_copy_class": "complete loss"
    }


@pytest.fixture(scope="module")
def genomic_del3_37_loc():
    """Create test fixture GRCh37 deletion subject"""
    return {
        "_id": "ga4gh:VSL.hDvTkhN210ZKyZOhilS7FxaKUCJgj6JC",
        "sequence_id": "ga4gh:SQ.v7noePfnNpK8ghYXEqZ9NukMXW7YeNsm",
        "interval": {
            "type": "SequenceInterval",
            "start": {
                "min": 31078343,
                "max": 31118467,
                "type": "DefiniteRange"
            },
            "end": {
                "min": 33292396,
                "max": 33435269,
                "type": "DefiniteRange"
            }
        },
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_del3_abs_37(genomic_del3_37_loc):
    """Create test fixture absolute copy number variation"""
    return {
        "type": "AbsoluteCopyNumber",
        "_id": "ga4gh:VAC.Pv9I4Dqk69w-tX0axaikVqid-pozxU74",
        "subject": genomic_del3_37_loc,
        "copies": {"type": "Number", "value": 2}
    }


@pytest.fixture(scope="module")
def genomic_del3_rel_37(genomic_del3_37_loc):
    """Create test fixture relative copy number variation"""
    return {
        "type": "RelativeCopyNumber",
        "_id": "ga4gh:VRC.tA7GSAeHf_MMDlSbDNGwnmrIfPxDN3et",
        "subject": genomic_del3_37_loc,
        "relative_copy_class": "complete loss"
    }


@pytest.fixture(scope="module")
def genomic_del4_abs_38(genomic_del4_seq_loc):
    """Create test fixture absolute copy number variation"""
    return {
        "type": "AbsoluteCopyNumber",
        "_id": "ga4gh:VAC.XCyM9ayMTrarSAMc00sHCmsPcCV8ymIN",
        "subject": genomic_del4_seq_loc,
        "copies": {"type": "Number", "value": 4}
    }


@pytest.fixture(scope="module")
def genomic_del4_rel_38(genomic_del4_seq_loc):
    """Create test fixture relative copy number variation"""
    return {
        "type": "RelativeCopyNumber",
        "_id": "ga4gh:VRC.dtwRjvChZ6LuyDXyWTnVGQJidyfJsQfe",
        "subject": genomic_del4_seq_loc,
        "relative_copy_class": "partial loss"
    }


@pytest.fixture(scope="module")
def genomic_del4_37_loc():
    """Create test fixture GRCh37 deletion subject"""
    return {
        "_id": "ga4gh:VSL.7Chaj1X9NH2G9sSK1diUKgBEZO4pHqr8",
        "sequence_id": "ga4gh:SQ.v7noePfnNpK8ghYXEqZ9NukMXW7YeNsm",
        "interval": {
            "type": "SequenceInterval",
            "start": {
                "value": 31138612,
                "comparator": "<=",
                "type": "IndefiniteRange"
            },
            "end": {
                "value": 33357594,
                "comparator": ">=",
                "type": "IndefiniteRange"
            }
        },
        "type": "SequenceLocation"
    }


@pytest.fixture(scope="module")
def genomic_del4_abs_37(genomic_del4_37_loc):
    """Create test fixture absolute copy number variation"""
    return {
        "type": "AbsoluteCopyNumber",
        "_id": "ga4gh:VAC.gXM6rRlCid3C1DmUGT2XynmGXDvt80P6",
        "subject": genomic_del4_37_loc,
        "copies": {"type": "Number", "value": 4}
    }


@pytest.fixture(scope="module")
def genomic_del4_rel_37(genomic_del4_37_loc):
    """Create test fixture relative copy number variation"""
    return {
        "type": "RelativeCopyNumber",
        "_id": "ga4gh:VRC.-s62STiIcHBm1bo6iJw_Gjtnu-gea_x_",
        "subject": genomic_del4_37_loc,
        "relative_copy_class": "partial loss"
    }


@pytest.fixture(scope="module")
def genomic_del5_abs_38(genomic_del5_seq_loc):
    """Create test fixture absolute copy number variation"""
    return {
        "type": "AbsoluteCopyNumber",
        "_id": "ga4gh:VAC.9A8BspUwfxSsceIScGOBtivMMASDsaid",
        "subject": genomic_del5_seq_loc,
        "copies": {"type": "Number", "value": 2}
    }


@pytest.fixture(scope="module")
def genomic_del5_rel_38(genomic_del5_seq_loc):
    """Create test fixture relative copy number variation"""
    return {
        "type": "RelativeCopyNumber",
        "_id": "ga4gh:VRC.rs2ifklHL-urqO6p4wUUNIh60n_JbCsv",
        "subject": genomic_del5_seq_loc,
        "relative_copy_class": "copy neutral"
    }


@pytest.fixture(scope="module")
def genomic_del5_37_loc():
    """Create test fixture GRCh37 deletion subject"""
    return {
        "_id": "ga4gh:VSL.BXq1jg19bLLthS20Ux3MZJaXEeAAyuW9",
        "sequence_id": "ga4gh:SQ.v7noePfnNpK8ghYXEqZ9NukMXW7YeNsm",
        "interval": {
            "type": "SequenceInterval",
            "start": {
                "value": 18593473,
                "comparator": "<=",
                "type": "IndefiniteRange"
            },
            "end": {
                "value": 18671749,
                "type": "Number"
            }
        },
        "type": "SequenceLocation"
    }


@pytest.fixture(scope="module")
def genomic_del5_abs_37(genomic_del5_37_loc):
    """Create test fixture absolute copy number variation"""
    return {
        "type": "AbsoluteCopyNumber",
        "_id": "ga4gh:VAC.9zhlg8DRSsE87N5SngYrMDWXStzp_WOX",
        "subject": genomic_del5_37_loc,
        "copies": {"type": "Number", "value": 2}
    }


@pytest.fixture(scope="module")
def genomic_del5_rel_37(genomic_del5_37_loc):
    """Create test fixture relative copy number variation"""
    return {
        "type": "RelativeCopyNumber",
        "_id": "ga4gh:VRC.bRyIyKsXFopd8L5vhyZ55xucloQRdxQS",
        "subject": genomic_del5_37_loc,
        "relative_copy_class": "copy neutral"
    }


@pytest.fixture(scope="module")
def genomic_del6_abs_38(genomic_del6_seq_loc):
    """Create test fixture absolute copy number variation"""
    return {
        "type": "AbsoluteCopyNumber",
        "_id": "ga4gh:VAC.6RkHgDOiRMZKMKgI6rmG9C3T6WuMhcex",
        "subject": genomic_del6_seq_loc,
        "copies": {"type": "Number", "value": 1}
    }


@pytest.fixture(scope="module")
def genomic_del6_rel_38(genomic_del6_seq_loc):
    """Create test fixture relative copy number variation"""
    return {
        "type": "RelativeCopyNumber",
        "_id": "ga4gh:VRC.5sY9YBXF0L685VHdNB1ZVuD3PCI3AjDf",
        "subject": genomic_del6_seq_loc,
        "relative_copy_class": "low-level gain"
    }


@pytest.fixture(scope="module")
def genomic_del6_37_loc():
    """Create test fixture GRCh37 deletion subject"""
    return {
        "_id": "ga4gh:VSL.FtOUv03wexaJWg_U7DISUtrM5WPdX60d",
        "sequence_id": "ga4gh:SQ.KqaUhJMW3CDjhoVtBetdEKT1n6hM-7Ek",
        "interval": {
            "type": "SequenceInterval",
            "start": {
                "value": 133783901,
                "type": "Number"
            },
            "end": {
                "value": 133785996,
                "comparator": ">=",
                "type": "IndefiniteRange"
            }
        },
        "type": "SequenceLocation"
    }


@pytest.fixture(scope="module")
def genomic_del6_abs_37(genomic_del6_37_loc):
    """Create test fixture absolute copy number variation"""
    return {
        "type": "AbsoluteCopyNumber",
        "_id": "ga4gh:VAC.60XjT6dzYKX8rn6ocG4AVAxCoUFfdjI6",
        "subject": genomic_del6_37_loc,
        "copies": {"type": "Number", "value": 1}
    }


@pytest.fixture(scope="module")
def genomic_del6_rel_37(genomic_del6_37_loc):
    """Create test fixture relative copy number variation"""
    return {
        "type": "RelativeCopyNumber",
        "_id": "ga4gh:VRC.wg_fX9QY-izwXjdP3Xqsjj69df-hyFob",
        "subject": genomic_del6_37_loc,
        "relative_copy_class": "low-level gain"
    }


@pytest.mark.asyncio
async def test_genomic_dup1_absolute_cnv(test_handler, genomic_dup1_38_vac,
                                         genomic_dup1_abs_37):
    """Test that genomic duplication works correctly"""
    q = "NC_000003.12:g.49531262dup"  # 38
    resp = await test_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=2, do_liftover=False, )
    assert resp.absolute_copy_number.dict(by_alias=True) == genomic_dup1_38_vac

    q = "NC_000003.11:g.49568695dup"  # 37
    resp = await test_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=2, do_liftover=False)
    assert resp.absolute_copy_number.dict(by_alias=True) == genomic_dup1_abs_37

    resp = await test_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=2, do_liftover=True)
    assert resp.absolute_copy_number.dict(by_alias=True) == genomic_dup1_38_vac

    resp = await test_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=1, do_liftover=True)
    expected = copy.deepcopy(genomic_dup1_38_vac)
    expected["copies"]["value"] = 2
    expected["_id"] = "ga4gh:VAC.WssYmdwNjxfDZNsN0fIMt0dZvGuj5FiL"
    assert resp.absolute_copy_number.dict(by_alias=True) == expected


@pytest.mark.asyncio
async def test_genomic_dup1_relative_cnv(test_handler, genomic_dup1_rel_38,
                                         genomic_dup1_rel_37):
    """Test that genomic duplication works correctly"""
    q = "NC_000003.12:g.49531262dup"  # 38
    resp = await test_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="complete loss", do_liftover=False)
    assert resp.relative_copy_number.dict(by_alias=True) == genomic_dup1_rel_38

    q = "NC_000003.11:g.49568695dup"  # 37
    resp = await test_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="complete loss", do_liftover=False)
    assert resp.relative_copy_number.dict(by_alias=True) == genomic_dup1_rel_37

    resp = await test_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="complete loss", do_liftover=True)
    assert resp.relative_copy_number.dict(by_alias=True) == genomic_dup1_rel_38


@pytest.mark.asyncio
async def test_genomic_dup2_absolute_cnv(test_handler, genomic_dup2_38_vac,
                                         genomic_dup2_abs_37):
    """Test that genomic duplication works correctly"""
    q = "NC_000016.10:g.2087938_2087948dup"  # 38
    resp = await test_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=2, do_liftover=False)
    assert resp.absolute_copy_number.dict(by_alias=True) == genomic_dup2_38_vac

    q = "NC_000016.9:g.2137939_2137949dup"  # 37
    resp = await test_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=2, do_liftover=False)
    assert resp.absolute_copy_number.dict(by_alias=True) == genomic_dup2_abs_37

    resp = await test_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=2, do_liftover=True)
    assert resp.absolute_copy_number.dict(by_alias=True) == genomic_dup2_38_vac

    resp = await test_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=1, do_liftover=True)
    expected = copy.deepcopy(genomic_dup2_38_vac)
    expected["copies"]["value"] = 2
    expected["_id"] = "ga4gh:VAC.dJw2a6Ft5QLcMr55dW7AJhZAFIC4AOOZ"
    assert resp.absolute_copy_number.dict(by_alias=True) == expected


@pytest.mark.asyncio
async def test_genomic_dup2_relative_cnv(test_handler, genomic_dup2_rel_38,
                                         genomic_dup2_rel_37):
    """Test that genomic duplication works correctly"""
    q = "NC_000016.10:g.2087938_2087948dup"  # 38
    resp = await test_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="partial loss", do_liftover=False)
    assert resp.relative_copy_number.dict(by_alias=True) == genomic_dup2_rel_38

    q = "NC_000016.9:g.2137939_2137949dup"  # 37
    resp = await test_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="partial loss", do_liftover=False)
    assert resp.relative_copy_number.dict(by_alias=True) == genomic_dup2_rel_37

    resp = await test_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="partial loss", do_liftover=True)
    assert resp.relative_copy_number.dict(by_alias=True) == genomic_dup2_rel_38


@pytest.mark.asyncio
async def test_genomic_dup3_absolute_cnv(test_handler, genomic_dup3_abs_38,
                                         genomic_dup3_abs_37):
    """Test that genomic duplication works correctly"""
    q = "NC_000023.11:g.(31060227_31100351)_(33274278_33417151)dup"  # 38
    resp = await test_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=1, do_liftover=False)
    assert resp.absolute_copy_number.dict(by_alias=True) == genomic_dup3_abs_38

    q = "NC_000023.10:g.(31078344_31118468)_(33292395_33435268)dup"  # 37
    resp = await test_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=1, do_liftover=False)
    assert resp.absolute_copy_number.dict(by_alias=True) == genomic_dup3_abs_37

    resp = await test_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=1, do_liftover=True)
    assert resp.absolute_copy_number.dict(by_alias=True) == genomic_dup3_abs_38

    resp = await test_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=2, do_liftover=True)
    expected = copy.deepcopy(genomic_dup3_abs_38)
    expected["copies"] = {"value": 3, "type": "Number"}
    expected["_id"] = "ga4gh:VAC.F95CRBT6lq24zbMfzDGJRlIPr2YTzCP6"
    assert resp.absolute_copy_number.dict(by_alias=True) == expected


@pytest.mark.asyncio
async def test_genomic_dup3_relative_cnv(test_handler, genomic_dup3_rel_38,
                                         genomic_dup3_rel_37):
    """Test that genomic duplication works correctly"""
    q = "NC_000023.11:g.(31060227_31100351)_(33274278_33417151)dup"  # 38
    resp = await test_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="high-level gain", do_liftover=False)
    assert resp.relative_copy_number.dict(by_alias=True) == genomic_dup3_rel_38

    q = "NC_000023.10:g.(31078344_31118468)_(33292395_33435268)dup"  # 37
    resp = await test_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="high-level gain", do_liftover=False)
    assert resp.relative_copy_number.dict(by_alias=True) == genomic_dup3_rel_37

    resp = await test_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="high-level gain", do_liftover=True)
    assert resp.relative_copy_number.dict(by_alias=True) == genomic_dup3_rel_38


@pytest.mark.asyncio
async def test_genomic_dup4_absolute_cnv(test_handler, genomic_dup4_abs_38,
                                         genomic_dup4_abs_37):
    """Test that genomic duplication works correctly"""
    q = "NC_000020.11:g.(?_30417576)_(31394018_?)dup"  # 38
    resp = await test_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=2, do_liftover=False)
    assert resp.absolute_copy_number.dict(by_alias=True) == genomic_dup4_abs_38

    q = "NC_000020.10:g.(?_29652252)_(29981821_?)dup"  # 37
    resp = await test_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=2, do_liftover=False)
    assert resp.absolute_copy_number.dict(by_alias=True) == genomic_dup4_abs_37

    resp = await test_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=2, do_liftover=True)
    assert resp.absolute_copy_number.dict(by_alias=True) == genomic_dup4_abs_38

    resp = await test_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=1, do_liftover=True)
    expected = copy.deepcopy(genomic_dup4_abs_38)
    expected["copies"]["value"] = 2
    expected["_id"] = "ga4gh:VAC.4DWr_DnOSASkH3D_oQfwFSYpfbjic8h2"
    assert resp.absolute_copy_number.dict(by_alias=True) == expected


@pytest.mark.asyncio
async def test_genomic_dup4_relative_cnv(test_handler, genomic_dup4_rel_38,
                                         genomic_dup4_rel_37):
    """Test that genomic duplication works correctly"""
    q = "NC_000020.11:g.(?_30417576)_(31394018_?)dup"  # 38
    resp = await test_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="complete loss", do_liftover=False)
    assert resp.relative_copy_number.dict(by_alias=True) == genomic_dup4_rel_38

    q = "NC_000020.10:g.(?_29652252)_(29981821_?)dup"  # 37
    resp = await test_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="complete loss", do_liftover=False)
    assert resp.relative_copy_number.dict(by_alias=True) == genomic_dup4_rel_37

    resp = await test_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="complete loss", do_liftover=True)
    assert resp.relative_copy_number.dict(by_alias=True) == genomic_dup4_rel_38


@pytest.mark.asyncio
async def test_genomic_dup5_absolute_cnv(test_handler, genomic_dup5_abs_38,
                                         genomic_dup5_abs_37):
    """Test that genomic duplication works correctly"""
    q = "NC_000023.11:g.(?_154021812)_154092209dup"  # 38
    resp = await test_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=3, do_liftover=False)
    assert resp.absolute_copy_number.dict(by_alias=True) == genomic_dup5_abs_38

    q = "NC_000023.10:g.(?_153287263)_153357667dup"  # 37
    resp = await test_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=3, do_liftover=False)
    assert resp.absolute_copy_number.dict(by_alias=True) == genomic_dup5_abs_37

    resp = await test_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=3, do_liftover=True)
    assert resp.absolute_copy_number.dict(by_alias=True) == genomic_dup5_abs_38

    resp = await test_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=4, do_liftover=True)
    expected = copy.deepcopy(genomic_dup5_abs_38)
    expected["copies"] = {"value": 5, "type": "Number"}
    expected["_id"] = "ga4gh:VAC.8jpu_iCkX1viYFmmOhxjr_54-1Ra2htN"
    assert resp.absolute_copy_number.dict(by_alias=True) == expected


@pytest.mark.asyncio
async def test_genomic_dup5_relative_cnv(test_handler, genomic_dup5_rel_38,
                                         genomic_dup5_rel_37):
    """Test that genomic duplication works correctly"""
    q = "NC_000023.11:g.(?_154021812)_154092209dup"  # 38
    resp = await test_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="partial loss", do_liftover=False)
    assert resp.relative_copy_number.dict(by_alias=True) == genomic_dup5_rel_38

    q = "NC_000023.10:g.(?_153287263)_153357667dup"  # 37
    resp = await test_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="partial loss", do_liftover=False)
    assert resp.relative_copy_number.dict(by_alias=True) == genomic_dup5_rel_37

    resp = await test_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="partial loss", do_liftover=True)
    assert resp.relative_copy_number.dict(by_alias=True) == genomic_dup5_rel_38


@pytest.mark.asyncio
async def test_genomic_dup6_absolute_cnv(test_handler, genomic_dup6_abs_38,
                                         genomic_dup6_abs_37):
    """Test that genomic duplication works correctly"""
    q = "NC_000023.11:g.154021812_(154092209_?)dup"  # 38
    resp = await test_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=1, do_liftover=False)
    assert resp.absolute_copy_number.dict(by_alias=True) == genomic_dup6_abs_38

    q = "NC_000023.10:g.153287263_(153357667_?)dup"  # 37
    resp = await test_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=1, do_liftover=False)
    assert resp.absolute_copy_number.dict(by_alias=True) == genomic_dup6_abs_37

    resp = await test_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=1, do_liftover=True)
    assert resp.absolute_copy_number.dict(by_alias=True) == genomic_dup6_abs_38

    resp = await test_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=2, do_liftover=True)
    expected = copy.deepcopy(genomic_dup6_abs_38)
    expected["copies"] = {"value": 3, "type": "Number"}
    expected["_id"] = "ga4gh:VAC.duHtKr9ejX_YCRUqUfn-x_I0F9LmXYD8"
    assert resp.absolute_copy_number.dict(by_alias=True) == expected


@pytest.mark.asyncio
async def test_genomic_dup6_relative_cnv(test_handler, genomic_dup6_rel_38,
                                         genomic_dup6_rel_37):
    """Test that genomic duplication works correctly"""
    q = "NC_000023.11:g.154021812_(154092209_?)dup"  # 38
    resp = await test_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="copy neutral", do_liftover=False)
    assert resp.relative_copy_number.dict(by_alias=True) == genomic_dup6_rel_38

    q = "NC_000023.10:g.153287263_(153357667_?)dup"  # 37
    resp = await test_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="copy neutral", do_liftover=False)
    assert resp.relative_copy_number.dict(by_alias=True) == genomic_dup6_rel_37

    resp = await test_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="copy neutral", do_liftover=True)
    assert resp.relative_copy_number.dict(by_alias=True) == genomic_dup6_rel_38


@pytest.mark.asyncio
async def test_genomic_del1_absolute_cnv(test_handler, genomic_del1_38_vac,
                                         genomic_del1_abs_37):
    """Test that genomic deletion works correctly"""
    q = "NC_000003.12:g.10149811del"  # 38
    resp = await test_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=2, do_liftover=False)
    assert resp.absolute_copy_number.dict(by_alias=True) == genomic_del1_38_vac

    q = "NC_000003.11:g.10191495del"  # 37
    resp = await test_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=2, do_liftover=False)
    assert resp.absolute_copy_number.dict(by_alias=True) == genomic_del1_abs_37

    resp = await test_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=2, do_liftover=True)
    assert resp.absolute_copy_number.dict(by_alias=True) == genomic_del1_38_vac

    resp = await test_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=3, do_liftover=True)
    expected = copy.deepcopy(genomic_del1_38_vac)
    expected["copies"]["value"] = 2
    expected["_id"] = "ga4gh:VAC.GPMXa216PdJGGjbwm4fGMxz-Ixwx7k8a"
    assert resp.absolute_copy_number.dict(by_alias=True) == expected


@pytest.mark.asyncio
async def test_genomic_del1_relative_cnv(test_handler, genomic_del1_rel_38,
                                         genomic_del1_rel_37):
    """Test that genomic deletion works correctly"""
    q = "NC_000003.12:g.10149811del"  # 38
    resp = await test_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="copy neutral", do_liftover=False)
    assert resp.relative_copy_number.dict(by_alias=True) == genomic_del1_rel_38

    q = "NC_000003.11:g.10191495del"  # 37
    resp = await test_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="copy neutral", do_liftover=False)
    assert resp.relative_copy_number.dict(by_alias=True) == genomic_del1_rel_37

    resp = await test_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="copy neutral", do_liftover=True)
    assert resp.relative_copy_number.dict(by_alias=True) == genomic_del1_rel_38


@pytest.mark.asyncio
async def test_genomic_del2_absolute_cnv(test_handler, genomic_del2_38_vac,
                                         genomic_del2_abs_37):
    """Test that genomic deletion works correctly"""
    q = "NC_000003.12:g.10146595_10146613del"  # 38
    resp = await test_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=2, do_liftover=False)
    assert resp.absolute_copy_number.dict(by_alias=True) == genomic_del2_38_vac

    q = "NC_000003.11:g.10188279_10188297del"  # 37
    resp = await test_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=2, do_liftover=False)
    assert resp.absolute_copy_number.dict(by_alias=True) == genomic_del2_abs_37

    resp = await test_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=2, do_liftover=True)
    assert resp.absolute_copy_number.dict(by_alias=True) == genomic_del2_38_vac

    resp = await test_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=4, do_liftover=True)
    expected = copy.deepcopy(genomic_del2_38_vac)
    expected["copies"]["value"] = 3
    expected["_id"] = "ga4gh:VAC.AfTepmH6qmkb9wypryJF4bj7dOkAOzjp"
    assert resp.absolute_copy_number.dict(by_alias=True) == expected


@pytest.mark.asyncio
async def test_genomic_del2_relative_cnv(test_handler, genomic_del2_rel_38,
                                         genomic_del2_rel_37):
    """Test that genomic deletion works correctly"""
    q = "NC_000003.12:g.10146595_10146613del"  # 38
    resp = await test_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="low-level gain", do_liftover=False)
    assert resp.relative_copy_number.dict(by_alias=True) == genomic_del2_rel_38

    q = "NC_000003.11:g.10188279_10188297del"  # 37
    resp = await test_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="low-level gain", do_liftover=False)
    assert resp.relative_copy_number.dict(by_alias=True) == genomic_del2_rel_37

    resp = await test_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="low-level gain", do_liftover=True)
    assert resp.relative_copy_number.dict(by_alias=True) == genomic_del2_rel_38


@pytest.mark.asyncio
async def test_genomic_del3_absolute_cnv(test_handler, genomic_del3_abs_38,
                                         genomic_del3_abs_37):
    """Test that genomic deletion works correctly"""
    q = "NC_000023.11:g.(31060227_31100351)_(33274278_33417151)del"  # 38
    resp = await test_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=3, do_liftover=False)
    assert resp.absolute_copy_number.dict(by_alias=True) == genomic_del3_abs_38

    q = "NC_000023.10:g.(31078344_31118468)_(33292395_33435268)del"  # 37
    resp = await test_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=3, do_liftover=False)
    assert resp.absolute_copy_number.dict(by_alias=True) == genomic_del3_abs_37

    resp = await test_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=3, do_liftover=True)
    assert resp.absolute_copy_number.dict(by_alias=True) == genomic_del3_abs_38

    resp = await test_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=2, do_liftover=True)
    expected = copy.deepcopy(genomic_del3_abs_38)
    expected["copies"] = {"value": 1, "type": "Number"}
    expected["_id"] = "ga4gh:VAC.ErMN6pFh9a0Qiirf-o3096XY3yjl73de"
    assert resp.absolute_copy_number.dict(by_alias=True) == expected


@pytest.mark.asyncio
async def test_genomic_del3_relative_cnv(test_handler, genomic_del3_rel_38,
                                         genomic_del3_rel_37):
    """Test that genomic deletion works correctly"""
    q = "NC_000023.11:g.(31060227_31100351)_(33274278_33417151)del"  # 38
    resp = await test_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="complete loss", do_liftover=False)
    assert resp.relative_copy_number.dict(by_alias=True) == genomic_del3_rel_38

    q = "NC_000023.10:g.(31078344_31118468)_(33292395_33435268)del"  # 37
    resp = await test_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="complete loss", do_liftover=False)
    assert resp.relative_copy_number.dict(by_alias=True) == genomic_del3_rel_37

    resp = await test_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="complete loss", do_liftover=True)
    assert resp.relative_copy_number.dict(by_alias=True) == genomic_del3_rel_38


@pytest.mark.asyncio
async def test_genomic_del4_absolute_cnv(test_handler, genomic_del4_abs_38,
                                         genomic_del4_abs_37):
    """Test that genomic deletion works correctly"""
    q = "NC_000023.11:g.(?_31120496)_(33339477_?)del"  # 38
    resp = await test_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=5, do_liftover=False)
    assert resp.absolute_copy_number.dict(by_alias=True) == genomic_del4_abs_38

    q = "NC_000023.10:g.(?_31138613)_(33357594_?)del"  # 37
    resp = await test_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=5, do_liftover=False)
    assert resp.absolute_copy_number.dict(by_alias=True) == genomic_del4_abs_37

    resp = await test_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=5, do_liftover=True)
    assert resp.absolute_copy_number.dict(by_alias=True) == genomic_del4_abs_38

    resp = await test_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=3, do_liftover=True)
    expected = copy.deepcopy(genomic_del4_abs_38)
    expected["copies"] = {"value": 2, "type": "Number"}
    expected["_id"] = "ga4gh:VAC.CK7qBFD64kIib5frYSFOxYXOqNbqfXVj"
    assert resp.absolute_copy_number.dict(by_alias=True) == expected


@pytest.mark.asyncio
async def test_genomic_del4_relative_cnv(test_handler, genomic_del4_rel_38,
                                         genomic_del4_rel_37):
    """Test that genomic deletion works correctly"""
    q = "NC_000023.11:g.(?_31120496)_(33339477_?)del"  # 38
    resp = await test_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="partial loss", do_liftover=False)
    assert resp.relative_copy_number.dict(by_alias=True) == genomic_del4_rel_38

    q = "NC_000023.10:g.(?_31138613)_(33357594_?)del"  # 37
    resp = await test_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="partial loss", do_liftover=False)
    assert resp.relative_copy_number.dict(by_alias=True) == genomic_del4_rel_37

    resp = await test_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="partial loss", do_liftover=True)
    assert resp.relative_copy_number.dict(by_alias=True) == genomic_del4_rel_38


@pytest.mark.asyncio
async def test_genomic_del5_absolute_cnv(test_handler, genomic_del5_abs_38,
                                         genomic_del5_abs_37):
    """Test that genomic deletion works correctly"""
    q = "NC_000023.11:g.(?_18575354)_18653629del"  # 38
    resp = await test_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=3, do_liftover=False)
    assert resp.absolute_copy_number.dict(by_alias=True) == genomic_del5_abs_38

    q = "NC_000023.10:g.(?_18593474)_18671749del"  # 37
    resp = await test_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=3, do_liftover=False)
    assert resp.absolute_copy_number.dict(by_alias=True) == genomic_del5_abs_37

    resp = await test_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=3, do_liftover=True)
    assert resp.absolute_copy_number.dict(by_alias=True) == genomic_del5_abs_38

    resp = await test_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=2, do_liftover=True)
    expected = copy.deepcopy(genomic_del5_abs_38)
    expected["copies"] = {"value": 1, "type": "Number"}
    expected["_id"] = "ga4gh:VAC.OAF6eGtER6jSRo2Y8srmitbxns5rMc7z"
    assert resp.absolute_copy_number.dict(by_alias=True) == expected


@pytest.mark.asyncio
async def test_genomic_del5_relative_cnv(test_handler, genomic_del5_rel_38,
                                         genomic_del5_rel_37):
    """Test that genomic deletion works correctly"""
    q = "NC_000023.11:g.(?_18575354)_18653629del"  # 38
    resp = await test_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="copy neutral", do_liftover=False)
    assert resp.relative_copy_number.dict(by_alias=True) == genomic_del5_rel_38

    q = "NC_000023.10:g.(?_18593474)_18671749del"  # 37
    resp = await test_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="copy neutral", do_liftover=False)
    assert resp.relative_copy_number.dict(by_alias=True) == genomic_del5_rel_37

    resp = await test_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="copy neutral", do_liftover=True)
    assert resp.relative_copy_number.dict(by_alias=True) == genomic_del5_rel_38


@pytest.mark.asyncio
async def test_genomic_del6_absolute_cnv(test_handler, genomic_del6_abs_38,
                                         genomic_del6_abs_37):
    """Test that genomic deletion works correctly"""
    q = "NC_000006.12:g.133462764_(133464858_?)del"  # 38
    resp = await test_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=2, do_liftover=False)
    assert resp.absolute_copy_number.dict(by_alias=True) == genomic_del6_abs_38

    q = "NC_000006.11:g.133783902_(133785996_?)del"  # 37
    resp = await test_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=2, do_liftover=False)
    assert resp.absolute_copy_number.dict(by_alias=True) == genomic_del6_abs_37

    resp = await test_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=2, do_liftover=True)
    assert resp.absolute_copy_number.dict(by_alias=True) == genomic_del6_abs_38

    resp = await test_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=3, do_liftover=True)
    expected = copy.deepcopy(genomic_del6_abs_38)
    expected["copies"]["value"] = 2
    expected["_id"] = "ga4gh:VAC.LZbbbtfCrpqAntJKaya0Qa9j0Lm5-u1R"
    assert resp.absolute_copy_number.dict(by_alias=True) == expected


@pytest.mark.asyncio
async def test_genomic_del6_relative_cnv(test_handler, genomic_del6_rel_38,
                                         genomic_del6_rel_37):
    """Test that genomic deletion works correctly"""
    q = "NC_000006.12:g.133462764_(133464858_?)del"  # 38
    resp = await test_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="low-level gain", do_liftover=False)
    assert resp.relative_copy_number.dict(by_alias=True) == genomic_del6_rel_38

    q = "NC_000006.11:g.133783902_(133785996_?)del"  # 37
    resp = await test_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="low-level gain", do_liftover=False)
    assert resp.relative_copy_number.dict(by_alias=True) == genomic_del6_rel_37

    resp = await test_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="low-level gain", do_liftover=True)
    assert resp.relative_copy_number.dict(by_alias=True) == genomic_del6_rel_38


@pytest.mark.asyncio
async def test_invalid_cnv_parameters(test_handler):
    """Check that invalid parameters return warnings"""
    q = "NC_000006.11:g.133783902_(133785996_?)del"
    resp, w = await test_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="low-level gains", do_liftover=True)
    assert resp is None
    assert w == ["low-level gains is not a valid relative copy class: ['complete loss', "  # noqa: E501
                 "'partial loss', 'copy neutral', 'low-level gain', 'high-level gain']"]


@pytest.mark.asyncio
async def test_invalid_cnv(test_handler):
    """Check that invalid input return warnings"""
    q = "DAG1 g.49568695dup"
    resp = await test_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="low-level gain", do_liftover=True,
        untranslatable_returns_text=True)
    assert set(resp.warnings) == {"Unable to translate DAG1 g.49568695dup to copy number variation",  # noqa: E501
                                  "DAG1 g.49568695dup is not a supported HGVS genomic duplication or deletion"}  # noqa: E501
    assert resp.relative_copy_number.type == "Text"

    q = "braf v600e"
    resp = await test_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="low-level gain", do_liftover=True)
    assert set(resp.warnings) == {"Unable to translate braf v600e to copy number variation",  # noqa: E501
                                  "braf v600e is not a supported HGVS genomic duplication or deletion"}  # noqa: E501
    assert resp.relative_copy_number is None
