"""Module for testing the hgvs to copy number count and copy number change endpoints"""
import copy

import pytest


@pytest.fixture(scope="module")
def genomic_dup1_cx_38(genomic_dup1_seq_loc_normalized):
    """Create test fixture copy number change variation"""
    return {
        "type": "CopyNumberChange",
        "_id": "ga4gh:CX.hGuvyiJmDtx4-MRjsLja0fb_DqOE2chN",
        "subject": genomic_dup1_seq_loc_normalized,
        "copy_change": "efo:0030069"
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
        "type": "SequenceLocation"
    }


@pytest.fixture(scope="module")
def genomic_dup1_cn_37(genomic_dup1_37_loc):
    """Create test fixture copy number count variation"""
    return {
        "type": "CopyNumberCount",
        "_id": "ga4gh:CN.pycd7qcmwmrqChKzm8VkIJUr6zYDoykP",
        "subject": genomic_dup1_37_loc,
        "copies": {"type": "Number", "value": 3}
    }


@pytest.fixture(scope="module")
def genomic_dup1_cx_37(genomic_dup1_37_loc):
    """Create test fixture copy number change variation"""
    return {
        "type": "CopyNumberChange",
        "_id": "ga4gh:CX.UiX_v-GZihnyYExWIno6WjN3ZYKABjTF",
        "subject": genomic_dup1_37_loc,
        "copy_change": "efo:0030069"
    }


@pytest.fixture(scope="module")
def genomic_dup2_cx_38(genomic_dup2_seq_loc_normalized):
    """Create test fixture copy number change variation"""
    return {
        "type": "CopyNumberChange",
        "_id": "ga4gh:CX.I2J6PFYfQMShQO72ga2QK9I52gbJTLoK",
        "subject": genomic_dup2_seq_loc_normalized,
        "copy_change": "efo:0030067"
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
def genomic_dup2_cn_37(genomic_dup2_37_loc):
    """Create test fixture copy number count variation"""
    return {
        "type": "CopyNumberCount",
        "_id": "ga4gh:CN.emp-0jz8WHt0ulMMUvon5iKRLxCKipLs",
        "subject": genomic_dup2_37_loc,
        "copies": {"type": "Number", "value": 3}
    }


@pytest.fixture(scope="module")
def genomic_dup2_cx_37(genomic_dup2_37_loc):
    """Create test fixture copy number change variation"""
    return {
        "type": "CopyNumberChange",
        "_id": "ga4gh:CX.6ld_A8Six_-TpnT_D9t694K4ZwNNzqIG",
        "subject": genomic_dup2_37_loc,
        "copy_change": "efo:0030067"
    }


@pytest.fixture(scope="module")
def genomic_dup3_cn_38(genomic_del3_dup3_loc_not_normalized):
    """Create test fixture copy number count variation"""
    return {
        "type": "CopyNumberCount",
        "_id": "ga4gh:CN.m0wq_fm3nMQDehJMtPve8OMyc880D0HE",
        "subject": genomic_del3_dup3_loc_not_normalized,
        "copies": {"type": "Number", "value": 2}
    }


@pytest.fixture(scope="module")
def genomic_dup3_cx_38(genomic_del3_dup3_loc_not_normalized):
    """Create test fixture copy number change variation"""
    return {
        "type": "CopyNumberChange",
        "_id": "ga4gh:CX.5uD6QzjG9zaNq-Ce-2elulzfiIJjXnjC",
        "subject": genomic_del3_dup3_loc_not_normalized,
        "copy_change": "efo:0030072"
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
def genomic_dup3_cn_37(genomic_dup3_37_loc):
    """Create test fixture copy number count variation"""
    return {
        "type": "CopyNumberCount",
        "_id": "ga4gh:CN.ZyEYeDs0BWbGGjOvkkSRxjthSdVbXi8q",
        "subject": genomic_dup3_37_loc,
        "copies": {"type": "Number", "value": 2}
    }


@pytest.fixture(scope="module")
def genomic_dup3_cx_37(genomic_dup3_37_loc):
    """Create test fixture copy number change variation"""
    return {
        "type": "CopyNumberChange",
        "_id": "ga4gh:CX.xe8q7BxvT0qPv8tr_G_FA_kOWd2FkED6",
        "subject": genomic_dup3_37_loc,
        "copy_change": "efo:0030072"
    }


@pytest.fixture(scope="module")
def genomic_dup4_cn_38(genomic_dup4_loc):
    """Create test fixture copy number count variation"""
    return {
        "type": "CopyNumberCount",
        "_id": "ga4gh:CN.ytM9tmgtdo6cU_qNw413JTXozc-QnZcp",
        "subject": genomic_dup4_loc,
        "copies": {"type": "Number", "value": 3}
    }


@pytest.fixture(scope="module")
def genomic_dup4_cx_38(genomic_dup4_loc):
    """Create test fixture copy number change variation"""
    return {
        "type": "CopyNumberChange",
        "_id": "ga4gh:CX.3daulOR3D7p_n8DsKmU-0YbS2U1iGDLI",
        "subject": genomic_dup4_loc,
        "copy_change": "efo:0030069"
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
        "type": "SequenceLocation"
    }


@pytest.fixture(scope="module")
def genomic_dup4_cn_37(genomic_dup4_37_loc):
    """Create test fixture copy number change variation"""
    return {
        "type": "CopyNumberCount",
        "_id": "ga4gh:CN.1nsSj3PtWTkgp4QvG9AalvexKLiIKwav",
        "subject": genomic_dup4_37_loc,
        "copies": {"type": "Number", "value": 3}
    }


@pytest.fixture(scope="module")
def genomic_dup4_cx_37(genomic_dup4_37_loc):
    """Create test fixture copy number change variation"""
    return {
        "type": "CopyNumberChange",
        "_id": "ga4gh:CX.EYHE5IrKx7jedSMFB0KVXA0aJgok2OnO",
        "subject": genomic_dup4_37_loc,
        "copy_change": "efo:0030069"
    }


@pytest.fixture(scope="module")
def genomic_dup5_cn_38(genomic_dup5_loc):
    """Create test fixture copy number count variation"""
    return {
        "type": "CopyNumberCount",
        "_id": "ga4gh:CN.yFFM4ooIWAr0cutUWwnU3yCh_wN64TJ6",
        "subject": genomic_dup5_loc,
        "copies": {"type": "Number", "value": 4}
    }


@pytest.fixture(scope="module")
def genomic_dup5_cx_38(genomic_dup5_loc):
    """Create test fixture copy number change variation"""
    return {
        "type": "CopyNumberChange",
        "_id": "ga4gh:CX.d1LMxDZ8tlxRBWtehc_ARpy_uFpcgcbW",
        "subject": genomic_dup5_loc,
        "copy_change": "efo:0030067"
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
        "type": "SequenceLocation"
    }


@pytest.fixture(scope="module")
def genomic_dup5_cn_37(genomic_dup5_37_loc):
    """Create test fixture copy number count variation"""
    return {
        "type": "CopyNumberCount",
        "_id": "ga4gh:CN.3C2rgy_TZSzAecP4vVGfWCVPEa8aD2Bf",
        "subject": genomic_dup5_37_loc,
        "copies": {"type": "Number", "value": 4}
    }


@pytest.fixture(scope="module")
def genomic_dup5_cx_37(genomic_dup5_37_loc):
    """Create test fixture copy number change variation"""
    return {
        "type": "CopyNumberChange",
        "_id": "ga4gh:CX.Gj3U_zuLnMnm-1-x-5lglYy7t-WecmSL",
        "subject": genomic_dup5_37_loc,
        "copy_change": "efo:0030067"
    }


@pytest.fixture(scope="module")
def genomic_dup6_cn_38(genomic_dup6_loc):
    """Create test fixture copy number count variation"""
    return {
        "type": "CopyNumberCount",
        "_id": "ga4gh:CN.tNxea8UWRp9ORzCDE2vtmJIqXEsUqp0j",
        "subject": genomic_dup6_loc,
        "copies": {"type": "Number", "value": 2}
    }


@pytest.fixture(scope="module")
def genomic_dup6_cx_38(genomic_dup6_loc):
    """Create test fixture copy number change variation"""
    return {
        "type": "CopyNumberChange",
        "_id": "ga4gh:CX.9ctgoxzJ15ehI6kE9TQIw2jnn4Tv-b8l",
        "subject": genomic_dup6_loc,
        "copy_change": "efo:0030064"
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
        "type": "SequenceLocation"
    }


@pytest.fixture(scope="module")
def genomic_dup6_cn_37(genomic_dup6_37_loc):
    """Create test fixture copy number count variation"""
    return {
        "type": "CopyNumberCount",
        "_id": "ga4gh:CN.7zbbnPO9ZXCA1pmqhbH2yHTdEtJ9W8dW",
        "subject": genomic_dup6_37_loc,
        "copies": {"type": "Number", "value": 2}
    }


@pytest.fixture(scope="module")
def genomic_dup6_cx_37(genomic_dup6_37_loc):
    """Create test fixture copy number change variation"""
    return {
        "type": "CopyNumberChange",
        "_id": "ga4gh:CX.zzqab9hKg2C0r7bAvPTvNiT0nO7kcE83",
        "subject": genomic_dup6_37_loc,
        "copy_change": "efo:0030064"
    }


@pytest.fixture(scope="module")
def genomic_del1_cx_38(genomic_del1_seq_loc):
    """Create test fixture copy number change variation"""
    return {
        "type": "CopyNumberChange",
        "_id": "ga4gh:CX.d00Awg55O8wqzMAO6lognOO0nZIu8Vvj",
        "subject": genomic_del1_seq_loc,
        "copy_change": "efo:0030064"
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
        "type": "SequenceLocation"
    }


@pytest.fixture(scope="module")
def genomic_del1_cn_37(genomic_del1_37_loc):
    """Create test fixture copy number count variation"""
    return {
        "type": "CopyNumberCount",
        "_id": "ga4gh:CN.JrxQrffVY_yGyYU6-fAT437Qe-MXOD1r",
        "subject": genomic_del1_37_loc,
        "copies": {"type": "Number", "value": 1}
    }


@pytest.fixture(scope="module")
def genomic_del1_cx_37(genomic_del1_37_loc):
    """Create test fixture copy number change variation"""
    return {
        "type": "CopyNumberChange",
        "_id": "ga4gh:CX.QYZkRsu-peJDFyVH8YrZzAReMmPepxcV",
        "subject": genomic_del1_37_loc,
        "copy_change": "efo:0030064"
    }


@pytest.fixture(scope="module")
def genomic_del2_cx_38(genomic_del2_seq_loc):
    """Create test fixture copy number change variation"""
    return {
        "type": "CopyNumberChange",
        "_id": "ga4gh:CX.AeXQgDz9pQF68hAJHvJWv5ZudDU59GLp",
        "subject": genomic_del2_seq_loc,
        "copy_change": "efo:0030071"
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
        "type": "SequenceLocation"
    }


@pytest.fixture(scope="module")
def genomic_del2_cn_37(genomic_del2_37_loc):
    """Create test fixture copy number count variation"""
    return {
        "type": "CopyNumberCount",
        "_id": "ga4gh:CN.iE0p_xeQu5Sp1JNU14seHF0zkqcH2TOP",
        "subject": genomic_del2_37_loc,
        "copies": {"type": "Number", "value": 1}
    }


@pytest.fixture(scope="module")
def genomic_del2_cx_37(genomic_del2_37_loc):
    """Create test fixture copy number change variation"""
    return {
        "type": "CopyNumberChange",
        "_id": "ga4gh:CX.kapN1JhyRh7De5rx0JtSSeTIfcBpGZW9",
        "subject": genomic_del2_37_loc,
        "copy_change": "efo:0030071"
    }


@pytest.fixture(scope="module")
def genomic_del3_cn_38(genomic_del3_dup3_loc_not_normalized):
    """Create test fixture copy number count variation"""
    return {
        "type": "CopyNumberCount",
        "_id": "ga4gh:CN.m0wq_fm3nMQDehJMtPve8OMyc880D0HE",
        "subject": genomic_del3_dup3_loc_not_normalized,
        "copies": {"type": "Number", "value": 2}
    }


@pytest.fixture(scope="module")
def genomic_del3_cx_38(genomic_del3_dup3_loc_not_normalized):
    """Create test fixture copy number change variation"""
    return {
        "type": "CopyNumberChange",
        "_id": "ga4gh:CX.LVmKnSkleRvcOVB-odZkNjDXsSIKme6a",
        "subject": genomic_del3_dup3_loc_not_normalized,
        "copy_change": "efo:0030069"
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
        "type": "SequenceLocation"
    }


@pytest.fixture(scope="module")
def genomic_del3_cn_37(genomic_del3_37_loc):
    """Create test fixture copy number count variation"""
    return {
        "type": "CopyNumberCount",
        "_id": "ga4gh:CN.ZyEYeDs0BWbGGjOvkkSRxjthSdVbXi8q",
        "subject": genomic_del3_37_loc,
        "copies": {"type": "Number", "value": 2}
    }


@pytest.fixture(scope="module")
def genomic_del3_cx_37(genomic_del3_37_loc):
    """Create test fixture copy number change variation"""
    return {
        "type": "CopyNumberChange",
        "_id": "ga4gh:CX.Uq6_nL4vwvXG7e5LPbbhv0Vp1GYF3Z33",
        "subject": genomic_del3_37_loc,
        "copy_change": "efo:0030069"
    }


@pytest.fixture(scope="module")
def genomic_del4_cn_38(genomic_del4_seq_loc):
    """Create test fixture copy number count variation"""
    return {
        "type": "CopyNumberCount",
        "_id": "ga4gh:CN.yYWqvpxkfj3VSGeVH5_c28omzar4XUBG",
        "subject": genomic_del4_seq_loc,
        "copies": {"type": "Number", "value": 4}
    }


@pytest.fixture(scope="module")
def genomic_del4_cx_38(genomic_del4_seq_loc):
    """Create test fixture copy number change variation"""
    return {
        "type": "CopyNumberChange",
        "_id": "ga4gh:CX.uHRB-mwjqAZlyOBE6Zi-T7QcKJVimIj-",
        "subject": genomic_del4_seq_loc,
        "copy_change": "efo:0030067"
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
def genomic_del4_cn_37(genomic_del4_37_loc):
    """Create test fixture copy number count variation"""
    return {
        "type": "CopyNumberCount",
        "_id": "ga4gh:CN.WTtNsuNX82V2xyuD6ff54C103A_E2Tf-",
        "subject": genomic_del4_37_loc,
        "copies": {"type": "Number", "value": 4}
    }


@pytest.fixture(scope="module")
def genomic_del4_cx_37(genomic_del4_37_loc):
    """Create test fixture copy number change variation"""
    return {
        "type": "CopyNumberChange",
        "_id": "ga4gh:CX.-RbBxVNVr1CpSRj2YT9ZNt8s4HnZ1Ntj",
        "subject": genomic_del4_37_loc,
        "copy_change": "efo:0030067"
    }


@pytest.fixture(scope="module")
def genomic_del5_cn_38(genomic_del5_seq_loc):
    """Create test fixture copy number count variation"""
    return {
        "type": "CopyNumberCount",
        "_id": "ga4gh:CN.4Md_wx8rtUUgr4fJfmIMxZqkKg9jkWIJ",
        "subject": genomic_del5_seq_loc,
        "copies": {"type": "Number", "value": 2}
    }


@pytest.fixture(scope="module")
def genomic_del5_cx_38(genomic_del5_seq_loc):
    """Create test fixture copy number change variation"""
    return {
        "type": "CopyNumberChange",
        "_id": "ga4gh:CX.pvZXOwOXLfKOQkw0dbZ_jshgtVvOBmHs",
        "subject": genomic_del5_seq_loc,
        "copy_change": "efo:0030064"
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
def genomic_del5_cn_37(genomic_del5_37_loc):
    """Create test fixture copy number count variation"""
    return {
        "type": "CopyNumberCount",
        "_id": "ga4gh:CN.plfngpTzJfSoPCWMCPIf24VnBIo1gPbT",
        "subject": genomic_del5_37_loc,
        "copies": {"type": "Number", "value": 2}
    }


@pytest.fixture(scope="module")
def genomic_del5_cx_37(genomic_del5_37_loc):
    """Create test fixture copy number change variation"""
    return {
        "type": "CopyNumberChange",
        "_id": "ga4gh:CX.WebbZzARTCpOxWL3jhbGqnutlRRi1rdN",
        "subject": genomic_del5_37_loc,
        "copy_change": "efo:0030064"
    }


@pytest.fixture(scope="module")
def genomic_del6_cn_38(genomic_del6_seq_loc):
    """Create test fixture copy number count variation"""
    return {
        "type": "CopyNumberCount",
        "_id": "ga4gh:CN.x-SGeFBpSWWI5qeMe7CIi5lTZxhnKvKJ",
        "subject": genomic_del6_seq_loc,
        "copies": {"type": "Number", "value": 1}
    }


@pytest.fixture(scope="module")
def genomic_del6_cx_38(genomic_del6_seq_loc):
    """Create test fixture copy number change variation"""
    return {
        "type": "CopyNumberChange",
        "_id": "ga4gh:CX.6egM2II947oerhaIGxXuX_lawZxS9p1e",
        "subject": genomic_del6_seq_loc,
        "copy_change": "efo:0030071"
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
def genomic_del6_cn_37(genomic_del6_37_loc):
    """Create test fixture copy number count variation"""
    return {
        "type": "CopyNumberCount",
        "_id": "ga4gh:CN.XQOdeq8E1UaSYHJ_rUlCng8Hfzowz9r8",
        "subject": genomic_del6_37_loc,
        "copies": {"type": "Number", "value": 1}
    }


@pytest.fixture(scope="module")
def genomic_del6_cx_37(genomic_del6_37_loc):
    """Create test fixture copy number change variation"""
    return {
        "type": "CopyNumberChange",
        "_id": "ga4gh:CX.4UHOph1nd5XlYxde6IdPxWhupLGWx5aL",
        "subject": genomic_del6_37_loc,
        "copy_change": "efo:0030071"
    }


@pytest.mark.asyncio
async def test_genomic_dup1_copy_number_count(test_cnv_handler, genomic_dup1_38_cn,
                                              genomic_dup1_cn_37):
    """Test that genomic duplication works correctly"""
    q = "NC_000003.12:g.49531262dup"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=2, do_liftover=False, )
    assert resp.copy_number_count.dict(by_alias=True) == genomic_dup1_38_cn

    q = "NC_000003.11:g.49568695dup"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=2, do_liftover=False)
    assert resp.copy_number_count.dict(by_alias=True) == genomic_dup1_cn_37

    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=2, do_liftover=True)
    assert resp.copy_number_count.dict(by_alias=True) == genomic_dup1_38_cn

    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=1, do_liftover=True)
    expected = copy.deepcopy(genomic_dup1_38_cn)
    expected["copies"]["value"] = 2
    expected["_id"] = "ga4gh:CN.lPsJh-HYw3fn3I0EF4geYUUxzCCf9nrU"
    assert resp.copy_number_count.dict(by_alias=True) == expected


@pytest.mark.asyncio
async def test_genomic_dup1_copy_number_change(test_cnv_handler, genomic_dup1_cx_38,
                                               genomic_dup1_cx_37):
    """Test that genomic duplication works correctly"""
    q = "NC_000003.12:g.49531262dup"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030069", do_liftover=False)
    assert resp.copy_number_change.dict(by_alias=True) == genomic_dup1_cx_38

    q = "NC_000003.11:g.49568695dup"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030069", do_liftover=False)
    assert resp.copy_number_change.dict(by_alias=True) == genomic_dup1_cx_37

    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030069", do_liftover=True)
    assert resp.copy_number_change.dict(by_alias=True) == genomic_dup1_cx_38


@pytest.mark.asyncio
async def test_genomic_dup2_copy_number_count(test_cnv_handler, genomic_dup2_38_cn,
                                              genomic_dup2_cn_37):
    """Test that genomic duplication works correctly"""
    q = "NC_000023.11:g.33211290_33211293dup"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=2, do_liftover=False)
    assert resp.copy_number_count.dict(by_alias=True) == genomic_dup2_38_cn

    q = "NC_000023.10:g.33229407_33229410dup"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=2, do_liftover=False)
    assert resp.copy_number_count.dict(by_alias=True) == genomic_dup2_cn_37

    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=2, do_liftover=True)
    assert resp.copy_number_count.dict(by_alias=True) == genomic_dup2_38_cn

    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=1, do_liftover=True)
    expected = copy.deepcopy(genomic_dup2_38_cn)
    expected["copies"]["value"] = 2
    expected["_id"] = "ga4gh:CN.CDu0H77EVAC_xR1STNJy0BwPzVP3vvCK"
    assert resp.copy_number_count.dict(by_alias=True) == expected


@pytest.mark.asyncio
async def test_genomic_dup2_copy_number_change(test_cnv_handler, genomic_dup2_cx_38,
                                               genomic_dup2_cx_37):
    """Test that genomic duplication works correctly"""
    q = "NC_000023.11:g.33211290_33211293dup"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030067", do_liftover=False)
    assert resp.copy_number_change.dict(by_alias=True) == genomic_dup2_cx_38

    q = "NC_000023.10:g.33229407_33229410dup"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030067", do_liftover=False)
    assert resp.copy_number_change.dict(by_alias=True) == genomic_dup2_cx_37

    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030067", do_liftover=True)
    assert resp.copy_number_change.dict(by_alias=True) == genomic_dup2_cx_38


@pytest.mark.asyncio
async def test_genomic_dup3_copy_number_count(test_cnv_handler, genomic_dup3_cn_38,
                                              genomic_dup3_cn_37):
    """Test that genomic duplication works correctly"""
    q = "NC_000023.11:g.(31060227_31100351)_(33274278_33417151)dup"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=1, do_liftover=False)
    assert resp.copy_number_count.dict(by_alias=True) == genomic_dup3_cn_38

    q = "NC_000023.10:g.(31078344_31118468)_(33292395_33435268)dup"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=1, do_liftover=False)
    assert resp.copy_number_count.dict(by_alias=True) == genomic_dup3_cn_37

    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=1, do_liftover=True)
    assert resp.copy_number_count.dict(by_alias=True) == genomic_dup3_cn_38

    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=2, do_liftover=True)
    expected = copy.deepcopy(genomic_dup3_cn_38)
    expected["copies"] = {"value": 3, "type": "Number"}
    expected["_id"] = "ga4gh:CN.nSFsagp2K8It3Pg-UbIMLAG3o7sAlwF8"
    assert resp.copy_number_count.dict(by_alias=True) == expected


@pytest.mark.asyncio
async def test_genomic_dup3_copy_number_change(test_cnv_handler, genomic_dup3_cx_38,
                                               genomic_dup3_cx_37):
    """Test that genomic duplication works correctly"""
    q = "NC_000023.11:g.(31060227_31100351)_(33274278_33417151)dup"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030072", do_liftover=False)
    assert resp.copy_number_change.dict(by_alias=True) == genomic_dup3_cx_38

    q = "NC_000023.10:g.(31078344_31118468)_(33292395_33435268)dup"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030072", do_liftover=False)
    assert resp.copy_number_change.dict(by_alias=True) == genomic_dup3_cx_37

    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030072", do_liftover=True)
    assert resp.copy_number_change.dict(by_alias=True) == genomic_dup3_cx_38


@pytest.mark.asyncio
async def test_genomic_dup4_copy_number_count(test_cnv_handler, genomic_dup4_cn_38,
                                              genomic_dup4_cn_37):
    """Test that genomic duplication works correctly"""
    q = "NC_000020.11:g.(?_30417576)_(31394018_?)dup"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=2, do_liftover=False)
    assert resp.copy_number_count.dict(by_alias=True) == genomic_dup4_cn_38

    q = "NC_000020.10:g.(?_29652252)_(29981821_?)dup"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=2, do_liftover=False)
    assert resp.copy_number_count.dict(by_alias=True) == genomic_dup4_cn_37

    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=2, do_liftover=True)
    assert resp.copy_number_count.dict(by_alias=True) == genomic_dup4_cn_38

    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=1, do_liftover=True)
    expected = copy.deepcopy(genomic_dup4_cn_38)
    expected["copies"]["value"] = 2
    expected["_id"] = "ga4gh:CN.V4cTNZe00EJzrR8VIRacCAQUhZc3wWuD"
    assert resp.copy_number_count.dict(by_alias=True) == expected


@pytest.mark.asyncio
async def test_genomic_dup4_copy_number_change(test_cnv_handler, genomic_dup4_cx_38,
                                               genomic_dup4_cx_37):
    """Test that genomic duplication works correctly"""
    q = "NC_000020.11:g.(?_30417576)_(31394018_?)dup"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030069", do_liftover=False)
    assert resp.copy_number_change.dict(by_alias=True) == genomic_dup4_cx_38

    q = "NC_000020.10:g.(?_29652252)_(29981821_?)dup"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030069", do_liftover=False)
    assert resp.copy_number_change.dict(by_alias=True) == genomic_dup4_cx_37

    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030069", do_liftover=True)
    assert resp.copy_number_change.dict(by_alias=True) == genomic_dup4_cx_38


@pytest.mark.asyncio
async def test_genomic_dup5_copy_number_count(test_cnv_handler, genomic_dup5_cn_38,
                                              genomic_dup5_cn_37):
    """Test that genomic duplication works correctly"""
    q = "NC_000023.11:g.(?_154021812)_154092209dup"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=3, do_liftover=False)
    assert resp.copy_number_count.dict(by_alias=True) == genomic_dup5_cn_38

    q = "NC_000023.10:g.(?_153287263)_153357667dup"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=3, do_liftover=False)
    assert resp.copy_number_count.dict(by_alias=True) == genomic_dup5_cn_37

    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=3, do_liftover=True)
    assert resp.copy_number_count.dict(by_alias=True) == genomic_dup5_cn_38

    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=4, do_liftover=True)
    expected = copy.deepcopy(genomic_dup5_cn_38)
    expected["copies"] = {"value": 5, "type": "Number"}
    expected["_id"] = "ga4gh:CN.5M3WmbZb-DJSHCX59VpEqoKjQXFws7Kz"
    assert resp.copy_number_count.dict(by_alias=True) == expected


@pytest.mark.asyncio
async def test_genomic_dup5_copy_number_change(test_cnv_handler, genomic_dup5_cx_38,
                                               genomic_dup5_cx_37):
    """Test that genomic duplication works correctly"""
    q = "NC_000023.11:g.(?_154021812)_154092209dup"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030067", do_liftover=False)
    assert resp.copy_number_change.dict(by_alias=True) == genomic_dup5_cx_38

    q = "NC_000023.10:g.(?_153287263)_153357667dup"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030067", do_liftover=False)
    assert resp.copy_number_change.dict(by_alias=True) == genomic_dup5_cx_37

    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030067", do_liftover=True)
    assert resp.copy_number_change.dict(by_alias=True) == genomic_dup5_cx_38


@pytest.mark.asyncio
async def test_genomic_dup6_copy_number_count(test_cnv_handler, genomic_dup6_cn_38,
                                              genomic_dup6_cn_37):
    """Test that genomic duplication works correctly"""
    q = "NC_000023.11:g.154021812_(154092209_?)dup"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=1, do_liftover=False)
    assert resp.copy_number_count.dict(by_alias=True) == genomic_dup6_cn_38

    q = "NC_000023.10:g.153287263_(153357667_?)dup"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=1, do_liftover=False)
    assert resp.copy_number_count.dict(by_alias=True) == genomic_dup6_cn_37

    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=1, do_liftover=True)
    assert resp.copy_number_count.dict(by_alias=True) == genomic_dup6_cn_38

    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=2, do_liftover=True)
    expected = copy.deepcopy(genomic_dup6_cn_38)
    expected["copies"] = {"value": 3, "type": "Number"}
    expected["_id"] = "ga4gh:CN.BfAbxCE8G4Ag3BDhwE0xU7z7l1eTOBtF"
    assert resp.copy_number_count.dict(by_alias=True) == expected


@pytest.mark.asyncio
async def test_genomic_dup6_copy_number_change(test_cnv_handler, genomic_dup6_cx_38,
                                               genomic_dup6_cx_37):
    """Test that genomic duplication works correctly"""
    q = "NC_000023.11:g.154021812_(154092209_?)dup"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030064", do_liftover=False)
    assert resp.copy_number_change.dict(by_alias=True) == genomic_dup6_cx_38

    q = "NC_000023.10:g.153287263_(153357667_?)dup"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030064", do_liftover=False)
    assert resp.copy_number_change.dict(by_alias=True) == genomic_dup6_cx_37

    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030064", do_liftover=True)
    assert resp.copy_number_change.dict(by_alias=True) == genomic_dup6_cx_38


@pytest.mark.asyncio
async def test_genomic_del1_copy_number_count(test_cnv_handler, genomic_del1_38_cn,
                                              genomic_del1_cn_37):
    """Test that genomic deletion works correctly"""
    q = "NC_000003.12:g.10149811del"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=2, do_liftover=False)
    assert resp.copy_number_count.dict(by_alias=True) == genomic_del1_38_cn

    q = "NC_000003.11:g.10191495del"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=2, do_liftover=False)
    assert resp.copy_number_count.dict(by_alias=True) == genomic_del1_cn_37

    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=2, do_liftover=True)
    assert resp.copy_number_count.dict(by_alias=True) == genomic_del1_38_cn

    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=3, do_liftover=True)
    expected = copy.deepcopy(genomic_del1_38_cn)
    expected["copies"]["value"] = 2
    expected["_id"] = "ga4gh:CN.V-JXhJ2ZmmZR2IxfXyvcVY3jbyJkv0Ku"
    assert resp.copy_number_count.dict(by_alias=True) == expected


@pytest.mark.asyncio
async def test_genomic_del1_copy_number_change(test_cnv_handler, genomic_del1_cx_38,
                                               genomic_del1_cx_37):
    """Test that genomic deletion works correctly"""
    q = "NC_000003.12:g.10149811del"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030064", do_liftover=False)
    assert resp.copy_number_change.dict(by_alias=True) == genomic_del1_cx_38

    q = "NC_000003.11:g.10191495del"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030064", do_liftover=False)
    assert resp.copy_number_change.dict(by_alias=True) == genomic_del1_cx_37

    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030064", do_liftover=True)
    assert resp.copy_number_change.dict(by_alias=True) == genomic_del1_cx_38


@pytest.mark.asyncio
async def test_genomic_del2_copy_number_count(test_cnv_handler, genomic_del2_38_cn,
                                              genomic_del2_cn_37):
    """Test that genomic deletion works correctly"""
    q = "NC_000003.12:g.10146595_10146613del"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=2, do_liftover=False)
    assert resp.copy_number_count.dict(by_alias=True) == genomic_del2_38_cn

    q = "NC_000003.11:g.10188279_10188297del"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=2, do_liftover=False)
    assert resp.copy_number_count.dict(by_alias=True) == genomic_del2_cn_37

    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=2, do_liftover=True)
    assert resp.copy_number_count.dict(by_alias=True) == genomic_del2_38_cn

    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=4, do_liftover=True)
    expected = copy.deepcopy(genomic_del2_38_cn)
    expected["copies"]["value"] = 3
    expected["_id"] = "ga4gh:CN.LnafQdwGW6J8wxSQUZ7VHB69mYvX5fnq"
    assert resp.copy_number_count.dict(by_alias=True) == expected


@pytest.mark.asyncio
async def test_genomic_del2_copy_number_change(test_cnv_handler, genomic_del2_cx_38,
                                               genomic_del2_cx_37):
    """Test that genomic deletion works correctly"""
    q = "NC_000003.12:g.10146595_10146613del"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030071", do_liftover=False)
    assert resp.copy_number_change.dict(by_alias=True) == genomic_del2_cx_38

    q = "NC_000003.11:g.10188279_10188297del"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030071", do_liftover=False)
    assert resp.copy_number_change.dict(by_alias=True) == genomic_del2_cx_37

    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030071", do_liftover=True)
    assert resp.copy_number_change.dict(by_alias=True) == genomic_del2_cx_38


@pytest.mark.asyncio
async def test_genomic_del3_copy_number_count(test_cnv_handler, genomic_del3_cn_38,
                                              genomic_del3_cn_37):
    """Test that genomic deletion works correctly"""
    q = "NC_000023.11:g.(31060227_31100351)_(33274278_33417151)del"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=3, do_liftover=False)
    assert resp.copy_number_count.dict(by_alias=True) == genomic_del3_cn_38

    q = "NC_000023.10:g.(31078344_31118468)_(33292395_33435268)del"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=3, do_liftover=False)
    assert resp.copy_number_count.dict(by_alias=True) == genomic_del3_cn_37

    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=3, do_liftover=True)
    assert resp.copy_number_count.dict(by_alias=True) == genomic_del3_cn_38

    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=2, do_liftover=True)
    expected = copy.deepcopy(genomic_del3_cn_38)
    expected["copies"] = {"value": 1, "type": "Number"}
    expected["_id"] = "ga4gh:CN.9ydcgFUfIXFnCbzke6Iu1Cg7VssMRawR"
    assert resp.copy_number_count.dict(by_alias=True) == expected


@pytest.mark.asyncio
async def test_genomic_del3_copy_number_change(test_cnv_handler, genomic_del3_cx_38,
                                               genomic_del3_cx_37):
    """Test that genomic deletion works correctly"""
    q = "NC_000023.11:g.(31060227_31100351)_(33274278_33417151)del"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030069", do_liftover=False)
    assert resp.copy_number_change.dict(by_alias=True) == genomic_del3_cx_38

    q = "NC_000023.10:g.(31078344_31118468)_(33292395_33435268)del"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030069", do_liftover=False)
    assert resp.copy_number_change.dict(by_alias=True) == genomic_del3_cx_37

    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030069", do_liftover=True)
    assert resp.copy_number_change.dict(by_alias=True) == genomic_del3_cx_38


@pytest.mark.asyncio
async def test_genomic_del4_copy_number_count(test_cnv_handler, genomic_del4_cn_38,
                                              genomic_del4_cn_37):
    """Test that genomic deletion works correctly"""
    q = "NC_000023.11:g.(?_31120496)_(33339477_?)del"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=5, do_liftover=False)
    assert resp.copy_number_count.dict(by_alias=True) == genomic_del4_cn_38

    q = "NC_000023.10:g.(?_31138613)_(33357594_?)del"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=5, do_liftover=False)
    assert resp.copy_number_count.dict(by_alias=True) == genomic_del4_cn_37

    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=5, do_liftover=True)
    assert resp.copy_number_count.dict(by_alias=True) == genomic_del4_cn_38

    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=3, do_liftover=True)
    expected = copy.deepcopy(genomic_del4_cn_38)
    expected["copies"] = {"value": 2, "type": "Number"}
    expected["_id"] = "ga4gh:CN.y2oNmJM4XzlVbAxcPWPq4JgCPuFQKGYQ"
    assert resp.copy_number_count.dict(by_alias=True) == expected


@pytest.mark.asyncio
async def test_genomic_del4_copy_number_change(test_cnv_handler, genomic_del4_cx_38,
                                               genomic_del4_cx_37):
    """Test that genomic deletion works correctly"""
    q = "NC_000023.11:g.(?_31120496)_(33339477_?)del"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030067", do_liftover=False)
    assert resp.copy_number_change.dict(by_alias=True) == genomic_del4_cx_38

    q = "NC_000023.10:g.(?_31138613)_(33357594_?)del"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030067", do_liftover=False)
    assert resp.copy_number_change.dict(by_alias=True) == genomic_del4_cx_37

    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030067", do_liftover=True)
    assert resp.copy_number_change.dict(by_alias=True) == genomic_del4_cx_38


@pytest.mark.asyncio
async def test_genomic_del5_copy_number_count(test_cnv_handler, genomic_del5_cn_38,
                                              genomic_del5_cn_37):
    """Test that genomic deletion works correctly"""
    q = "NC_000023.11:g.(?_18575354)_18653629del"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=3, do_liftover=False)
    assert resp.copy_number_count.dict(by_alias=True) == genomic_del5_cn_38

    q = "NC_000023.10:g.(?_18593474)_18671749del"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=3, do_liftover=False)
    assert resp.copy_number_count.dict(by_alias=True) == genomic_del5_cn_37

    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=3, do_liftover=True)
    assert resp.copy_number_count.dict(by_alias=True) == genomic_del5_cn_38

    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=2, do_liftover=True)
    expected = copy.deepcopy(genomic_del5_cn_38)
    expected["copies"] = {"value": 1, "type": "Number"}
    expected["_id"] = "ga4gh:CN.KuWPcxhww7j1W7nwDs5741_6d_TgJOlf"
    assert resp.copy_number_count.dict(by_alias=True) == expected


@pytest.mark.asyncio
async def test_genomic_del5_copy_number_change(test_cnv_handler, genomic_del5_cx_38,
                                               genomic_del5_cx_37):
    """Test that genomic deletion works correctly"""
    q = "NC_000023.11:g.(?_18575354)_18653629del"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030064", do_liftover=False)
    assert resp.copy_number_change.dict(by_alias=True) == genomic_del5_cx_38

    q = "NC_000023.10:g.(?_18593474)_18671749del"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030064", do_liftover=False)
    assert resp.copy_number_change.dict(by_alias=True) == genomic_del5_cx_37

    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030064", do_liftover=True)
    assert resp.copy_number_change.dict(by_alias=True) == genomic_del5_cx_38


@pytest.mark.asyncio
async def test_genomic_del6_copy_number_count(test_cnv_handler, genomic_del6_cn_38,
                                              genomic_del6_cn_37):
    """Test that genomic deletion works correctly"""
    q = "NC_000006.12:g.133462764_(133464858_?)del"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=2, do_liftover=False)
    assert resp.copy_number_count.dict(by_alias=True) == genomic_del6_cn_38

    q = "NC_000006.11:g.133783902_(133785996_?)del"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=2, do_liftover=False)
    assert resp.copy_number_count.dict(by_alias=True) == genomic_del6_cn_37

    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=2, do_liftover=True)
    assert resp.copy_number_count.dict(by_alias=True) == genomic_del6_cn_38

    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=3, do_liftover=True)
    expected = copy.deepcopy(genomic_del6_cn_38)
    expected["copies"]["value"] = 2
    expected["_id"] = "ga4gh:CN.vfPoDhlxDQsLJaxA9WiIbRtyUCgD3y1k"
    assert resp.copy_number_count.dict(by_alias=True) == expected


@pytest.mark.asyncio
async def test_genomic_del6_copy_number_change(test_cnv_handler, genomic_del6_cx_38,
                                               genomic_del6_cx_37):
    """Test that genomic deletion works correctly"""
    q = "NC_000006.12:g.133462764_(133464858_?)del"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030071", do_liftover=False)
    assert resp.copy_number_change.dict(by_alias=True) == genomic_del6_cx_38

    q = "NC_000006.11:g.133783902_(133785996_?)del"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030071", do_liftover=False)
    assert resp.copy_number_change.dict(by_alias=True) == genomic_del6_cx_37

    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030071", do_liftover=True)
    assert resp.copy_number_change.dict(by_alias=True) == genomic_del6_cx_38


@pytest.mark.asyncio
async def test_invalid_cnv_parameters(test_cnv_handler):
    """Check that invalid parameters return warnings"""
    q = "NC_000006.11:g.133783902_(133785996_?)del"
    resp, w = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="low-level gain", do_liftover=True)
    assert resp is None
    assert w == ["low-level gain is not a valid copy change: ['efo:0030069', "
                 "'efo:0020073', 'efo:0030068', 'efo:0030067', 'efo:0030064', "
                 "'efo:0030070', 'efo:0030071', 'efo:0030072']"]


@pytest.mark.asyncio
async def test_invalid_cnv(test_cnv_handler):
    """Check that invalid input return warnings"""
    q = "DAG1 g.49568695dup"
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030071", do_liftover=True,
        untranslatable_returns_text=True)
    assert set(resp.warnings) == {"Unable to translate DAG1 g.49568695dup to copy number variation",  # noqa: E501
                                  "DAG1 g.49568695dup is not a supported HGVS genomic duplication or deletion"}  # noqa: E501
    assert resp.copy_number_change.type == "Text"

    q = "braf V600E"
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030071", do_liftover=True)
    assert set(resp.warnings) == {"Unable to translate braf V600E to copy number variation",  # noqa: E501
                                  "braf V600E is not a supported HGVS genomic duplication or deletion"}  # noqa: E501
    assert resp.copy_number_change is None

    # Not yet supported
    for q in ["NC_000018.9:g.(48556994_48573289)_48573471dup",
              "NC_000018.9:g.48556994_(48573289_48573471)dup"]:
        resp = await test_cnv_handler.hgvs_to_copy_number_change(
            q, copy_change="efo:0030070"
        )
        assert resp.warnings == [f"Unable to find classification for: {q}"], q
        assert resp.copy_number_change is None, q
