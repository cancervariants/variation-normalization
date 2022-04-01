"""Module for testing the hgvs to absolute and relative copy number endpoints"""
import copy

import pytest


@pytest.fixture(scope="module")
def genomic_dup1_38_subject(genomic_dup1_seq_loc):
    """Create test fixture GRCh38 duplication subject"""
    return {
        "location": genomic_dup1_seq_loc,
        "reverse_complement": False,
        "type": "DerivedSequenceExpression"
    }


@pytest.fixture(scope="module")
def genomic_dup1_rel_38(genomic_dup1_38_subject):
    """Create test fixture relative copy number variation"""
    return {
        "type": "RelativeCopyNumber",
        "_id": "ga4gh:VRC.XiXamTGYJ43rc8xheleMKcjxEBOFp82l",
        "subject": genomic_dup1_38_subject,
        "relative_copy_class": "complete loss"
    }


@pytest.fixture(scope="module")
def genomic_dup1_37_subject():
    """Create test fixture GRCh37 duplication subject"""
    return {
        "location": {
            "_id": "ga4gh:VSL.CXcLL6RUPkro3dLXN0miGEzlzPYiqw2q",
            "sequence_id": "ga4gh:SQ.VNBualIltAyi2AI_uXcKU7M9XUOuA7MS",
            "interval": {
                "type": "SequenceInterval",
                "start": {"value": 49568693, "type": "Number"},
                "end": {"value": 49568695, "type": "Number"},
            },
            "type": "SequenceLocation",
        },
        "reverse_complement": False,
        "type": "DerivedSequenceExpression"
    }


@pytest.fixture(scope="module")
def genomic_dup1_abs_37(genomic_dup1_37_subject):
    """Create test fixture absolute copy number variation"""
    return {
        "type": "AbsoluteCopyNumber",
        "_id": "ga4gh:VAC.ReDrt_J0qRXOokz_32CWYdaWy7aATEHT",
        "subject": genomic_dup1_37_subject,
        "copies": {
            "type": "Number",
            "value": 3
        }
    }


@pytest.fixture(scope="module")
def genomic_dup1_rel_37(genomic_dup1_37_subject):
    """Create test fixture relative copy number variation"""
    return {
        "type": "RelativeCopyNumber",
        "_id": "ga4gh:VRC.Y7OyCA7q7mil6YFiHEh4pyZHPLkGaPDE",
        "subject": genomic_dup1_37_subject,
        "relative_copy_class": "complete loss"
    }


@pytest.fixture(scope="module")
def genomic_dup2_38_subject(genomic_dup2_seq_loc):
    """Create test fixture GRCh38 duplication subject"""
    return {
        "location": genomic_dup2_seq_loc,
        "reverse_complement": False,
        "type": "DerivedSequenceExpression"
    }


@pytest.fixture(scope="module")
def genomic_dup2_rel_38(genomic_dup2_38_subject):
    """Create test fixture relative copy number variation"""
    return {
        "type": "RelativeCopyNumber",
        "_id": "ga4gh:VRC.ZafHwmugzRCgqgy4jAS6W-l3y3gSEezQ",
        "subject": genomic_dup2_38_subject,
        "relative_copy_class": "partial loss"
    }


@pytest.fixture(scope="module")
def genomic_dup2_37_subject():
    """Create test fixture GRCh37 duplication subject"""
    return {
        "location": {
            "_id": "ga4gh:VSL.VY7qo3C2Nl-bpRtjOQ6MODJ2FfGUJRV1",
            "sequence_id": "ga4gh:SQ.W6wLoIFOn4G7cjopxPxYNk2lcEqhLQFb",
            "interval": {
                "type": "SequenceInterval",
                "start": {"value": 2137938, "type": "Number"},
                "end": {"value": 2137949, "type": "Number"},
            },
            "type": "SequenceLocation",
        },
        "reverse_complement": False,
        "type": "DerivedSequenceExpression"
    }


@pytest.fixture(scope="module")
def genomic_dup2_abs_37(genomic_dup2_37_subject):
    """Create test fixture absolute copy number variation"""
    return {
        "type": "AbsoluteCopyNumber",
        "_id": "ga4gh:VAC.8N1BeYiz3saze_AC26Gg4ocKF8IO7GqF",
        "subject": genomic_dup2_37_subject,
        "copies": {
            "type": "Number",
            "value": 3
        }
    }


@pytest.fixture(scope="module")
def genomic_dup2_rel_37(genomic_dup2_37_subject):
    """Create test fixture relative copy number variation"""
    return {
        "type": "RelativeCopyNumber",
        "_id": "ga4gh:VRC.1vUSqjdiWCEJrPYfGs_4PDEVLHohztdJ",
        "subject": genomic_dup2_37_subject,
        "relative_copy_class": "partial loss"
    }


@pytest.fixture(scope="module")
def genomic_dup3_38_subject(genomic_del3_dup3_loc):
    """Create test fixture GRCh38 duplication subject"""
    return {
        "location": genomic_del3_dup3_loc,
        "reverse_complement": False,
        "type": "DerivedSequenceExpression"
    }


@pytest.fixture(scope="module")
def genomic_dup3_abs_38(genomic_dup3_38_subject):
    """Create test fixture absolute copy number variation"""
    return {
        "type": "AbsoluteCopyNumber",
        "_id": "ga4gh:VAC.LT_jOeYNJg6UblsnItIhX0tppRff4fCh",
        "subject": genomic_dup3_38_subject,
        "copies": {
            "type": "DefiniteRange",
            "min": 2,
            "max": 3
        }
    }


@pytest.fixture(scope="module")
def genomic_dup3_rel_38(genomic_dup3_38_subject):
    """Create test fixture relative copy number variation"""
    return {
        "type": "RelativeCopyNumber",
        "_id": "ga4gh:VRC.-y-_z5rsOV0x-KsgPjYiDN6XDNV4MX7U",
        "subject": genomic_dup3_38_subject,
        "relative_copy_class": "high-level gain"
    }


@pytest.fixture(scope="module")
def genomic_dup3_37_subject():
    """Create test fixture GRCh37 duplication subject"""
    return {
        "location": {
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
        },
        "reverse_complement": False,
        "type": "DerivedSequenceExpression"
    }


@pytest.fixture(scope="module")
def genomic_dup3_abs_37(genomic_dup3_37_subject):
    """Create test fixture absolute copy number variation"""
    return {
        "type": "AbsoluteCopyNumber",
        "_id": "ga4gh:VAC.LAY9brjDTRS9v6EHIqorBJSWdPCdGehM",
        "subject": genomic_dup3_37_subject,
        "copies": {
            "type": "DefiniteRange",
            "min": 2,
            "max": 3
        }
    }


@pytest.fixture(scope="module")
def genomic_dup3_rel_37(genomic_dup3_37_subject):
    """Create test fixture relative copy number variation"""
    return {
        "type": "RelativeCopyNumber",
        "_id": "ga4gh:VRC.sPsCwxoV3zg4laO26t0iYm_FQoTThCBF",
        "subject": genomic_dup3_37_subject,
        "relative_copy_class": "high-level gain"
    }


@pytest.fixture(scope="module")
def genomic_dup4_38_subject(genoimc_dup4_loc):
    """Create test fixture GRCh38 duplication subject"""
    return {
        "location": genoimc_dup4_loc,
        "reverse_complement": False,
        "type": "DerivedSequenceExpression"
    }


@pytest.fixture(scope="module")
def genomic_dup4_abs_38(genomic_dup4_38_subject):
    """Create test fixture absolute copy number variation"""
    return {
        "type": "AbsoluteCopyNumber",
        "_id": "ga4gh:VAC.aHxuxlXZTVm-HroDd34Jh9BNPrgVHtML",
        "subject": genomic_dup4_38_subject,
        "copies": {
            "type": "Number",
            "value": 3
        }
    }


@pytest.fixture(scope="module")
def genomic_dup4_rel_38(genomic_dup4_38_subject):
    """Create test fixture relative copy number variation"""
    return {
        "type": "RelativeCopyNumber",
        "_id": "ga4gh:VRC.E2IxQfOBiCWvuXa7U-swwgbRoiV5aURe",
        "subject": genomic_dup4_38_subject,
        "relative_copy_class": "complete loss"
    }


@pytest.fixture(scope="module")
def genomic_dup4_37_subject():
    """Create test fixture GRCh37 duplication subject"""
    return {
        "location": {
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
        },
        "reverse_complement": False,
        "type": "DerivedSequenceExpression"
    }


@pytest.fixture(scope="module")
def genomic_dup4_abs_37(genomic_dup4_37_subject):
    """Create test fixture relative copy number variation"""
    return {
        "type": "AbsoluteCopyNumber",
        "_id": "ga4gh:VAC.2NR9ZYtmbFy31pIic_1cBZ9SEev0q-Sx",
        "subject": genomic_dup4_37_subject,
        "copies": {
            "type": "Number",
            "value": 3
        }
    }


@pytest.fixture(scope="module")
def genomic_dup4_rel_37(genomic_dup4_37_subject):
    """Create test fixture relative copy number variation"""
    return {
        "type": "RelativeCopyNumber",
        "_id": "ga4gh:VRC.JAk-UphGbw9-rFi0WoHJNW-CZKBIA81I",
        "subject": genomic_dup4_37_subject,
        "relative_copy_class": "complete loss"
    }


@pytest.fixture(scope="module")
def genomic_dup5_38_subject(genoimc_dup5_loc):
    """Create test fixture GRCh38 duplication subject"""
    return {
        "location": genoimc_dup5_loc,
        "reverse_complement": False,
        "type": "DerivedSequenceExpression"
    }


@pytest.fixture(scope="module")
def genomic_dup5_abs_38(genomic_dup5_38_subject):
    """Create test fixture absolute copy number variation"""
    return {
        "type": "AbsoluteCopyNumber",
        "_id": "ga4gh:VAC.aaz7BxN1V9TzvG934_kSQbbPNCT6vHU1",
        "subject": genomic_dup5_38_subject,
        "copies": {
            "type": "DefiniteRange",
            "min": 2,
            "max": 3
        }
    }


@pytest.fixture(scope="module")
def genomic_dup5_rel_38(genomic_dup5_38_subject):
    """Create test fixture relative copy number variation"""
    return {
        "type": "RelativeCopyNumber",
        "_id": "ga4gh:VRC.dfAwyWhv8FlJLXGG8fS86qp5kSmgK9tb",
        "subject": genomic_dup5_38_subject,
        "relative_copy_class": "partial loss"
    }


@pytest.fixture(scope="module")
def genomic_dup5_37_subject():
    """Create test fixture GRCh37 duplication subject"""
    return {
        "location": {
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
        },
        "reverse_complement": False,
        "type": "DerivedSequenceExpression"
    }


@pytest.fixture(scope="module")
def genomic_dup5_abs_37(genomic_dup5_37_subject):
    """Create test fixture absolute copy number variation"""
    return {
        "type": "AbsoluteCopyNumber",
        "_id": "ga4gh:VAC.99BliOHLTPSR7gqbPNKTCfRMXHFrwqgd",
        "subject": genomic_dup5_37_subject,
        "copies": {
            "type": "DefiniteRange",
            "min": 2,
            "max": 3
        }
    }


@pytest.fixture(scope="module")
def genomic_dup5_rel_37(genomic_dup5_37_subject):
    """Create test fixture relative copy number variation"""
    return {
        "type": "RelativeCopyNumber",
        "_id": "ga4gh:VRC.u7y0KCsAhXEy3bdcfYwOtGBBGNziLb3W",
        "subject": genomic_dup5_37_subject,
        "relative_copy_class": "partial loss"
    }


@pytest.fixture(scope="module")
def genomic_dup6_38_subject(genoimc_dup6_loc):
    """Create test fixture GRCh38 duplication subject"""
    return {
        "location": genoimc_dup6_loc,
        "reverse_complement": False,
        "type": "DerivedSequenceExpression"
    }


@pytest.fixture(scope="module")
def genomic_dup6_abs_38(genomic_dup6_38_subject):
    """Create test fixture absolute copy number variation"""
    return {
        "type": "AbsoluteCopyNumber",
        "_id": "ga4gh:VAC.v-dhCkHL7W67y3SHp-qfMg7F5DOaAbQn",
        "subject": genomic_dup6_38_subject,
        "copies": {
            "type": "DefiniteRange",
            "min": 2,
            "max": 3
        }
    }


@pytest.fixture(scope="module")
def genomic_dup6_rel_38(genomic_dup6_38_subject):
    """Create test fixture relative copy number variation"""
    return {
        "type": "RelativeCopyNumber",
        "_id": "ga4gh:VRC.dLmIPRSPDfaZsY6JkQZJPumqsZTHgvb8",
        "subject": genomic_dup6_38_subject,
        "relative_copy_class": "copy neutral"
    }


@pytest.fixture(scope="module")
def genomic_dup6_37_subject():
    """Create test fixture GRCh37 duplication subject"""
    return {
        "location": {
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
        },
        "reverse_complement": False,
        "type": "DerivedSequenceExpression"
    }


@pytest.fixture(scope="module")
def genomic_dup6_abs_37(genomic_dup6_37_subject):
    """Create test fixture absolute copy number variation"""
    return {
        "type": "AbsoluteCopyNumber",
        "_id": "ga4gh:VAC.aIHl5yt6YJPqrbh4FL6Fmlg7DTEqCh8y",
        "subject": genomic_dup6_37_subject,
        "copies": {
            "type": "DefiniteRange",
            "min": 2,
            "max": 3
        }
    }


@pytest.fixture(scope="module")
def genomic_dup6_rel_37(genomic_dup6_37_subject):
    """Create test fixture relative copy number variation"""
    return {
        "type": "RelativeCopyNumber",
        "_id": "ga4gh:VRC.LkA6PMk1quZ0K_bgGhI6XpmCylDocIR_",
        "subject": genomic_dup6_37_subject,
        "relative_copy_class": "copy neutral"
    }


@pytest.fixture(scope="module")
def genomic_del1_38_subject(genomic_del1_seq_loc):
    """Create test fixture GRCh38 deletion subject"""
    return {
        "location": genomic_del1_seq_loc,
        "reverse_complement": False,
        "type": "DerivedSequenceExpression"
    }


@pytest.fixture(scope="module")
def genomic_del1_rel_38(genomic_del1_38_subject):
    """Create test fixture relative copy number variation"""
    return {
        "type": "RelativeCopyNumber",
        "_id": "ga4gh:VRC.KQwe4TvYODBZeSUqC7K3qssIRQhxnevQ",
        "subject": genomic_del1_38_subject,
        "relative_copy_class": "copy neutral"
    }


@pytest.fixture(scope="module")
def genomic_del1_37_subject():
    """Create test fixture GRCh37 deletion subject"""
    return {
        "location": {
            "_id": "ga4gh:VSL.IpNE_OIaSJSeiiJTEwrxyiRsgpUay8_n",
            "sequence_id": "ga4gh:SQ.VNBualIltAyi2AI_uXcKU7M9XUOuA7MS",
            "interval": {
                "type": "SequenceInterval",
                "start": {"value": 10191494, "type": "Number"},
                "end": {"value": 10191495, "type": "Number"},
            },
            "type": "SequenceLocation",
        },
        "reverse_complement": False,
        "type": "DerivedSequenceExpression"
    }


@pytest.fixture(scope="module")
def genomic_del1_abs_37(genomic_del1_37_subject):
    """Create test fixture absolute copy number variation"""
    return {
        "type": "AbsoluteCopyNumber",
        "_id": "ga4gh:VAC.BJcNDL3rGhci1ftMCL0_exxDb45SWgxn",
        "subject": genomic_del1_37_subject,
        "copies": {
            "type": "Number",
            "value": 1
        }
    }


@pytest.fixture(scope="module")
def genomic_del1_rel_37(genomic_del1_37_subject):
    """Create test fixture relative copy number variation"""
    return {
        "type": "RelativeCopyNumber",
        "_id": "ga4gh:VRC.T-BAA72t1pfMsgxCaJGfS8vg_cdS-8dX",
        "subject": genomic_del1_37_subject,
        "relative_copy_class": "copy neutral"
    }


@pytest.fixture(scope="module")
def genomic_del2_38_subject(genomic_del2_seq_loc):
    """Create test fixture GRCh38 deletion subject"""
    return {
        "location": genomic_del2_seq_loc,
        "reverse_complement": False,
        "type": "DerivedSequenceExpression"
    }


@pytest.fixture(scope="module")
def genomic_del2_rel_38(genomic_del2_38_subject):
    """Create test fixture relative copy number variation"""
    return {
        "type": "RelativeCopyNumber",
        "_id": "ga4gh:VRC.DRcKvwgwmtxjpRUctWaFBgzK-3mu7p7L",
        "subject": genomic_del2_38_subject,
        "relative_copy_class": "low-level gain"
    }


@pytest.fixture(scope="module")
def genomic_del2_37_subject():
    """Create test fixture GRCh37 deletion subject"""
    return {
        "location": {
            "_id": "ga4gh:VSL.B-feFlnwavUnci2jFAqG8WOQrqFAZJ6N",
            "sequence_id": "ga4gh:SQ.VNBualIltAyi2AI_uXcKU7M9XUOuA7MS",
            "interval": {
                "type": "SequenceInterval",
                "start": {"value": 10188278, "type": "Number"},
                "end": {"value": 10188297, "type": "Number"},
            },
            "type": "SequenceLocation",
        },
        "reverse_complement": False,
        "type": "DerivedSequenceExpression"
    }


@pytest.fixture(scope="module")
def genomic_del2_abs_37(genomic_del2_37_subject):
    """Create test fixture absolute copy number variation"""
    return {
        "type": "AbsoluteCopyNumber",
        "_id": "ga4gh:VAC.PgSgggCYoZfyVIDT5LIYemgu6UQHAFOh",
        "subject": genomic_del2_37_subject,
        "copies": {
            "type": "Number",
            "value": 1
        }
    }


@pytest.fixture(scope="module")
def genomic_del2_rel_37(genomic_del2_37_subject):
    """Create test fixture relative copy number variation"""
    return {
        "type": "RelativeCopyNumber",
        "_id": "ga4gh:VRC.yV0_s897kIPggWedY4buNUzaokXTVzZG",
        "subject": genomic_del2_37_subject,
        "relative_copy_class": "low-level gain"
    }


@pytest.fixture(scope="module")
def genomic_del3_38_subject(genomic_del3_dup3_loc):
    """Create test fixture GRCh38 deletion subject"""
    return {
        "location": genomic_del3_dup3_loc,
        "reverse_complement": False,
        "type": "DerivedSequenceExpression"
    }


@pytest.fixture(scope="module")
def genomic_del3_abs_38(genomic_del3_38_subject):
    """Create test fixture absolute copy number variation"""
    return {
        "type": "AbsoluteCopyNumber",
        "_id": "ga4gh:VAC.o39DwfEb2I-LR9NFKR5cznWrRhI9Vx8b",
        "subject": genomic_del3_38_subject,
        "copies": {
            "type": "DefiniteRange",
            "min": 0,
            "max": 1
        }
    }


@pytest.fixture(scope="module")
def genomic_del3_rel_38(genomic_del3_38_subject):
    """Create test fixture relative copy number variation"""
    return {
        "type": "RelativeCopyNumber",
        "_id": "ga4gh:VRC.vpFW_FTFwvziXVlkOsxZTTFn0rvzHfnA",
        "subject": genomic_del3_38_subject,
        "relative_copy_class": "complete loss"
    }


@pytest.fixture(scope="module")
def genomic_del3_37_subject():
    """Create test fixture GRCh37 deletion subject"""
    return {
        "location": {
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
        },
        "reverse_complement": False,
        "type": "DerivedSequenceExpression"
    }


@pytest.fixture(scope="module")
def genomic_del3_abs_37(genomic_del3_37_subject):
    """Create test fixture absolute copy number variation"""
    return {
        "type": "AbsoluteCopyNumber",
        "_id": "ga4gh:VAC.0EC3Bi-mP2GtfIOhPMvYjhXVHIUhxYQP",
        "subject": genomic_del3_37_subject,
        "copies": {
            "type": "DefiniteRange",
            "min": 0,
            "max": 1
        }
    }


@pytest.fixture(scope="module")
def genomic_del3_rel_37(genomic_del3_37_subject):
    """Create test fixture relative copy number variation"""
    return {
        "type": "RelativeCopyNumber",
        "_id": "ga4gh:VRC.G1O15kvCHpTTTyAeuWNFG0rTj1SgvzFq",
        "subject": genomic_del3_37_subject,
        "relative_copy_class": "complete loss"
    }


@pytest.fixture(scope="module")
def genomic_del4_38_subject(genomic_del4_seq_loc):
    """Create test fixture GRCh38 deletion subject"""
    return {
        "location": genomic_del4_seq_loc,
        "reverse_complement": False,
        "type": "DerivedSequenceExpression"
    }


@pytest.fixture(scope="module")
def genomic_del4_abs_38(genomic_del4_38_subject):
    """Create test fixture absolute copy number variation"""
    return {
        "type": "AbsoluteCopyNumber",
        "_id": "ga4gh:VAC.oCrOJZcSJ-knvR1l6u-wpPNNMCLDCB0L",
        "subject": genomic_del4_38_subject,
        "copies": {
            "type": "DefiniteRange",
            "min": 0,
            "max": 1
        }
    }


@pytest.fixture(scope="module")
def genomic_del4_rel_38(genomic_del4_38_subject):
    """Create test fixture relative copy number variation"""
    return {
        "type": "RelativeCopyNumber",
        "_id": "ga4gh:VRC.u7we9nJeMlt-CAlQ47HRN7lWJ0UXay_W",
        "subject": genomic_del4_38_subject,
        "relative_copy_class": "partial loss"
    }


@pytest.fixture(scope="module")
def genomic_del4_37_subject():
    """Create test fixture GRCh37 deletion subject"""
    return {
        "location": {
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
            "type": "SequenceLocation",
        },
        "reverse_complement": False,
        "type": "DerivedSequenceExpression"
    }


@pytest.fixture(scope="module")
def genomic_del4_abs_37(genomic_del4_37_subject):
    """Create test fixture absolute copy number variation"""
    return {
        "type": "AbsoluteCopyNumber",
        "_id": "ga4gh:VAC.sBsyUfIqJcBr2Sl7rsscPcjlId-wX9_V",
        "subject": genomic_del4_37_subject,
        "copies": {
            "type": "DefiniteRange",
            "min": 0,
            "max": 1
        }
    }


@pytest.fixture(scope="module")
def genomic_del4_rel_37(genomic_del4_37_subject):
    """Create test fixture relative copy number variation"""
    return {
        "type": "RelativeCopyNumber",
        "_id": "ga4gh:VRC.E2X_1uJAfCD5LVHgv6kJKrq0bgUxZmEs",
        "subject": genomic_del4_37_subject,
        "relative_copy_class": "partial loss"
    }


@pytest.fixture(scope="module")
def genomic_del5_38_subject(genomic_del5_seq_loc):
    """Create test fixture GRCh38 deletion subject"""
    return {
        "location": genomic_del5_seq_loc,
        "reverse_complement": False,
        "type": "DerivedSequenceExpression"
    }


@pytest.fixture(scope="module")
def genomic_del5_abs_38(genomic_del5_38_subject):
    """Create test fixture absolute copy number variation"""
    return {
        "type": "AbsoluteCopyNumber",
        "_id": "ga4gh:VAC.WQkCH9rBzLPtB4f6LRb8Y-rhob68h-7M",
        "subject": genomic_del5_38_subject,
        "copies": {
            "type": "DefiniteRange",
            "min": 0,
            "max": 1
        }
    }


@pytest.fixture(scope="module")
def genomic_del5_rel_38(genomic_del5_38_subject):
    """Create test fixture relative copy number variation"""
    return {
        "type": "RelativeCopyNumber",
        "_id": "ga4gh:VRC.WNlJVKAmuN8BgZDvl584zG60asYKo7Iz",
        "subject": genomic_del5_38_subject,
        "relative_copy_class": "copy neutral"
    }


@pytest.fixture(scope="module")
def genomic_del5_37_subject():
    """Create test fixture GRCh37 deletion subject"""
    return {
        "location": {
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
            "type": "SequenceLocation",
        },
        "reverse_complement": False,
        "type": "DerivedSequenceExpression"
    }


@pytest.fixture(scope="module")
def genomic_del5_abs_37(genomic_del5_37_subject):
    """Create test fixture absolute copy number variation"""
    return {
        "type": "AbsoluteCopyNumber",
        "_id": "ga4gh:VAC.sA4s5CrkkknwUt64RUVNcy0AqY3uD3W4",
        "subject": genomic_del5_37_subject,
        "copies": {
            "type": "DefiniteRange",
            "min": 0,
            "max": 1
        }
    }


@pytest.fixture(scope="module")
def genomic_del5_rel_37(genomic_del5_37_subject):
    """Create test fixture relative copy number variation"""
    return {
        "type": "RelativeCopyNumber",
        "_id": "ga4gh:VRC.waFZlE_5ZXtD02BR4bWqYJGQ7LNf6I7a",
        "subject": genomic_del5_37_subject,
        "relative_copy_class": "copy neutral"
    }


@pytest.fixture(scope="module")
def genomic_del6_38_subject(genomic_del6_seq_loc):
    """Create test fixture GRCh38 deletion subject"""
    return {
        "location": genomic_del6_seq_loc,
        "reverse_complement": False,
        "type": "DerivedSequenceExpression"
    }


@pytest.fixture(scope="module")
def genomic_del6_abs_38(genomic_del6_38_subject):
    """Create test fixture absolute copy number variation"""
    return {
        "type": "AbsoluteCopyNumber",
        "_id": "ga4gh:VAC.YAyl7c-YMRkbqZymPFI2GKkYSQ9v7EeD",
        "subject": genomic_del6_38_subject,
        "copies": {
            "type": "Number",
            "value": 1
        }
    }


@pytest.fixture(scope="module")
def genomic_del6_rel_38(genomic_del6_38_subject):
    """Create test fixture relative copy number variation"""
    return {
        "type": "RelativeCopyNumber",
        "_id": "ga4gh:VRC.wS7t1-8Q0XfNKF9SbOEL_mD7X8joaslw",
        "subject": genomic_del6_38_subject,
        "relative_copy_class": "low-level gain"
    }


@pytest.fixture(scope="module")
def genomic_del6_37_subject():
    """Create test fixture GRCh37 deletion subject"""
    return {
        "location": {
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
            "type": "SequenceLocation",
        },
        "reverse_complement": False,
        "type": "DerivedSequenceExpression"
    }


@pytest.fixture(scope="module")
def genomic_del6_abs_37(genomic_del6_37_subject):
    """Create test fixture absolute copy number variation"""
    return {
        "type": "AbsoluteCopyNumber",
        "_id": "ga4gh:VAC.V9-mDixtWzx19X7dAdzTvXq4hXqQ9udx",
        "subject": genomic_del6_37_subject,
        "copies": {
            "type": "Number",
            "value": 1
        }
    }


@pytest.fixture(scope="module")
def genomic_del6_rel_37(genomic_del6_37_subject):
    """Create test fixture relative copy number variation"""
    return {
        "type": "RelativeCopyNumber",
        "_id": "ga4gh:VRC.ah-x-uIo3eTsIMXCvmXkEvhzm3d52X7Y",
        "subject": genomic_del6_37_subject,
        "relative_copy_class": "low-level gain"
    }


@pytest.mark.asyncio
async def test_genomic_dup1_absolute_cnv(test_query_handler, genomic_dup1_38_vac,
                                         genomic_dup1_abs_37):
    """Test that genomic duplication works correctly"""
    q = "NC_000003.12:g.49531262dup"  # 38
    resp, _ = await test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=None, do_liftover=False)
    assert resp.dict(by_alias=True) == genomic_dup1_38_vac

    q = "NC_000003.11:g.49568695dup"  # 37
    resp, _ = await test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=None, do_liftover=False)
    assert resp.dict(by_alias=True) == genomic_dup1_abs_37

    resp, _ = await test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=None, do_liftover=True)
    assert resp.dict(by_alias=True) == genomic_dup1_38_vac

    resp, _ = await test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=1, do_liftover=True)
    expected = copy.deepcopy(genomic_dup1_38_vac)
    expected["copies"]["value"] = 2
    expected["_id"] = "ga4gh:VAC.PDhQpQe82ssZiIMlghL8wDlf6xBIV_Ca"
    assert resp.dict(by_alias=True) == expected


@pytest.mark.asyncio
async def test_genomic_dup1_relative_cnv(test_query_handler, genomic_dup1_rel_38,
                                         genomic_dup1_rel_37):
    """Test that genomic duplication works correctly"""
    q = "NC_000003.12:g.49531262dup"  # 38
    resp, _ = await test_query_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="complete loss", do_liftover=False)
    assert resp.dict(by_alias=True) == genomic_dup1_rel_38

    q = "NC_000003.11:g.49568695dup"  # 37
    resp, _ = await test_query_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="complete loss", do_liftover=False)
    assert resp.dict(by_alias=True) == genomic_dup1_rel_37

    resp, _ = await test_query_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="complete loss", do_liftover=True)
    assert resp.dict(by_alias=True) == genomic_dup1_rel_38


@pytest.mark.asyncio
async def test_genomic_dup2_absolute_cnv(test_query_handler, genomic_dup2_38_vac,
                                         genomic_dup2_abs_37):
    """Test that genomic duplication works correctly"""
    q = "NC_000016.10:g.2087938_2087948dup"  # 38
    resp, _ = await test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=None, do_liftover=False)
    assert resp.dict(by_alias=True) == genomic_dup2_38_vac

    q = "NC_000016.9:g.2137939_2137949dup"  # 37
    resp, _ = await test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=None, do_liftover=False)
    assert resp.dict(by_alias=True) == genomic_dup2_abs_37

    resp, _ = await test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=None, do_liftover=True)
    assert resp.dict(by_alias=True) == genomic_dup2_38_vac

    resp, _ = await test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=1, do_liftover=True)
    expected = copy.deepcopy(genomic_dup2_38_vac)
    expected["copies"]["value"] = 2
    expected["_id"] = "ga4gh:VAC.Riu3BNPMzofH0ABV0Ulr1W1hxojyl4Qk"
    assert resp.dict(by_alias=True) == expected


@pytest.mark.asyncio
async def test_genomic_dup2_relative_cnv(test_query_handler, genomic_dup2_rel_38,
                                         genomic_dup2_rel_37):
    """Test that genomic duplication works correctly"""
    q = "NC_000016.10:g.2087938_2087948dup"  # 38
    resp, _ = await test_query_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="partial loss", do_liftover=False)
    assert resp.dict(by_alias=True) == genomic_dup2_rel_38

    q = "NC_000016.9:g.2137939_2137949dup"  # 37
    resp, _ = await test_query_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="partial loss", do_liftover=False)
    assert resp.dict(by_alias=True) == genomic_dup2_rel_37

    resp, _ = await test_query_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="partial loss", do_liftover=True)
    assert resp.dict(by_alias=True) == genomic_dup2_rel_38


@pytest.mark.asyncio
async def test_genomic_dup3_absolute_cnv(test_query_handler, genomic_dup3_abs_38,
                                         genomic_dup3_abs_37):
    """Test that genomic duplication works correctly"""
    q = "NC_000023.11:g.(31060227_31100351)_(33274278_33417151)dup"  # 38
    resp, _ = await test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=None, do_liftover=False)
    assert resp.dict(by_alias=True) == genomic_dup3_abs_38

    q = "NC_000023.10:g.(31078344_31118468)_(33292395_33435268)dup"  # 37
    resp, _ = await test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=None, do_liftover=False)
    assert resp.dict(by_alias=True) == genomic_dup3_abs_37

    resp, _ = await test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=None, do_liftover=True)
    assert resp.dict(by_alias=True) == genomic_dup3_abs_38

    resp, _ = await test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=2, do_liftover=True)
    expected = copy.deepcopy(genomic_dup3_abs_38)
    expected["copies"] = {"value": 3, "type": "Number"}
    expected["_id"] = "ga4gh:VAC.kuf5yPlF13yjhXKJMGc_gMiClDlvIFAY"
    assert resp.dict(by_alias=True) == expected


@pytest.mark.asyncio
async def test_genomic_dup3_relative_cnv(test_query_handler, genomic_dup3_rel_38,
                                         genomic_dup3_rel_37):
    """Test that genomic duplication works correctly"""
    q = "NC_000023.11:g.(31060227_31100351)_(33274278_33417151)dup"  # 38
    resp, _ = await test_query_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="high-level gain", do_liftover=False)
    assert resp.dict(by_alias=True) == genomic_dup3_rel_38

    q = "NC_000023.10:g.(31078344_31118468)_(33292395_33435268)dup"  # 37
    resp, _ = await test_query_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="high-level gain", do_liftover=False)
    assert resp.dict(by_alias=True) == genomic_dup3_rel_37

    resp, _ = await test_query_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="high-level gain", do_liftover=True)
    assert resp.dict(by_alias=True) == genomic_dup3_rel_38


@pytest.mark.asyncio
async def test_genomic_dup4_absolute_cnv(test_query_handler, genomic_dup4_abs_38,
                                         genomic_dup4_abs_37):
    """Test that genomic duplication works correctly"""
    q = "NC_000020.11:g.(?_30417576)_(31394018_?)dup"  # 38
    resp, _ = await test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=None, do_liftover=False)
    assert resp.dict(by_alias=True) == genomic_dup4_abs_38

    q = "NC_000020.10:g.(?_29652252)_(29981821_?)dup"  # 37
    resp, _ = await test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=None, do_liftover=False)
    assert resp.dict(by_alias=True) == genomic_dup4_abs_37

    resp, _ = await test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=None, do_liftover=True)
    assert resp.dict(by_alias=True) == genomic_dup4_abs_38

    resp, _ = await test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=1, do_liftover=True)
    expected = copy.deepcopy(genomic_dup4_abs_38)
    expected["copies"]["value"] = 2
    expected["_id"] = "ga4gh:VAC.qAzWMQEm93Q5FOSgsf8uWOVbiB9UZVc7"
    assert resp.dict(by_alias=True) == expected


@pytest.mark.asyncio
async def test_genomic_dup4_relative_cnv(test_query_handler, genomic_dup4_rel_38,
                                         genomic_dup4_rel_37):
    """Test that genomic duplication works correctly"""
    q = "NC_000020.11:g.(?_30417576)_(31394018_?)dup"  # 38
    resp, _ = await test_query_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="complete loss", do_liftover=False)
    assert resp.dict(by_alias=True) == genomic_dup4_rel_38

    q = "NC_000020.10:g.(?_29652252)_(29981821_?)dup"  # 37
    resp, _ = await test_query_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="complete loss", do_liftover=False)
    assert resp.dict(by_alias=True) == genomic_dup4_rel_37

    resp, _ = await test_query_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="complete loss", do_liftover=True)
    assert resp.dict(by_alias=True) == genomic_dup4_rel_38


@pytest.mark.asyncio
async def test_genomic_dup5_absolute_cnv(test_query_handler, genomic_dup5_abs_38,
                                         genomic_dup5_abs_37):
    """Test that genomic duplication works correctly"""
    q = "NC_000023.11:g.(?_154021812)_154092209dup"  # 38
    resp, _ = await test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=None, do_liftover=False)
    assert resp.dict(by_alias=True) == genomic_dup5_abs_38

    q = "NC_000023.10:g.(?_153287263)_153357667dup"  # 37
    resp, _ = await test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=None, do_liftover=False)
    assert resp.dict(by_alias=True) == genomic_dup5_abs_37

    resp, _ = await test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=None, do_liftover=True)
    assert resp.dict(by_alias=True) == genomic_dup5_abs_38

    resp, _ = await test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=4, do_liftover=True)
    expected = copy.deepcopy(genomic_dup5_abs_38)
    expected["copies"] = {"value": 5, "type": "Number"}
    expected["_id"] = "ga4gh:VAC.OCz5vEoeJ8rkAcJHiYiZt9yrscRYniVL"
    assert resp.dict(by_alias=True) == expected


@pytest.mark.asyncio
async def test_genomic_dup5_relative_cnv(test_query_handler, genomic_dup5_rel_38,
                                         genomic_dup5_rel_37):
    """Test that genomic duplication works correctly"""
    q = "NC_000023.11:g.(?_154021812)_154092209dup"  # 38
    resp, _ = await test_query_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="partial loss", do_liftover=False)
    assert resp.dict(by_alias=True) == genomic_dup5_rel_38

    q = "NC_000023.10:g.(?_153287263)_153357667dup"  # 37
    resp, _ = await test_query_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="partial loss", do_liftover=False)
    assert resp.dict(by_alias=True) == genomic_dup5_rel_37

    resp, _ = await test_query_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="partial loss", do_liftover=True)
    assert resp.dict(by_alias=True) == genomic_dup5_rel_38


@pytest.mark.asyncio
async def test_genomic_dup6_absolute_cnv(test_query_handler, genomic_dup6_abs_38,
                                         genomic_dup6_abs_37):
    """Test that genomic duplication works correctly"""
    q = "NC_000023.11:g.154021812_(154092209_?)dup"  # 38
    resp, _ = await test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=None, do_liftover=False)
    assert resp.dict(by_alias=True) == genomic_dup6_abs_38

    q = "NC_000023.10:g.153287263_(153357667_?)dup"  # 37
    resp, _ = await test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=None, do_liftover=False)
    assert resp.dict(by_alias=True) == genomic_dup6_abs_37

    resp, _ = await test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=None, do_liftover=True)
    assert resp.dict(by_alias=True) == genomic_dup6_abs_38

    resp, _ = await test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=2, do_liftover=True)
    expected = copy.deepcopy(genomic_dup6_abs_38)
    expected["copies"] = {"value": 3, "type": "Number"}
    expected["_id"] = "ga4gh:VAC.lfKTjaLLNWyis9O1lpn1aBfi_YSIFCz4"
    assert resp.dict(by_alias=True) == expected


@pytest.mark.asyncio
async def test_genomic_dup6_relative_cnv(test_query_handler, genomic_dup6_rel_38,
                                         genomic_dup6_rel_37):
    """Test that genomic duplication works correctly"""
    q = "NC_000023.11:g.154021812_(154092209_?)dup"  # 38
    resp, _ = await test_query_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="copy neutral", do_liftover=False)
    assert resp.dict(by_alias=True) == genomic_dup6_rel_38

    q = "NC_000023.10:g.153287263_(153357667_?)dup"  # 37
    resp, _ = await test_query_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="copy neutral", do_liftover=False)
    assert resp.dict(by_alias=True) == genomic_dup6_rel_37

    resp, _ = await test_query_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="copy neutral", do_liftover=True)
    assert resp.dict(by_alias=True) == genomic_dup6_rel_38


@pytest.mark.asyncio
async def test_genomic_del1_absolute_cnv(test_query_handler, genomic_del1_38_vac,
                                         genomic_del1_abs_37):
    """Test that genomic deletion works correctly"""
    q = "NC_000003.12:g.10149811del"  # 38
    resp, _ = await test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=None, do_liftover=False)
    assert resp.dict(by_alias=True) == genomic_del1_38_vac

    q = "NC_000003.11:g.10191495del"  # 37
    resp, _ = await test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=None, do_liftover=False)
    assert resp.dict(by_alias=True) == genomic_del1_abs_37

    resp, _ = await test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=None, do_liftover=True)
    assert resp.dict(by_alias=True) == genomic_del1_38_vac

    resp, _ = await test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=3, do_liftover=True)
    expected = copy.deepcopy(genomic_del1_38_vac)
    expected["copies"]["value"] = 2
    expected["_id"] = "ga4gh:VAC.YwwG6b9B7LDFsKlw_t1aRilDUu-hGajk"
    assert resp.dict(by_alias=True) == expected


@pytest.mark.asyncio
async def test_genomic_del1_relative_cnv(test_query_handler, genomic_del1_rel_38,
                                         genomic_del1_rel_37):
    """Test that genomic deletion works correctly"""
    q = "NC_000003.12:g.10149811del"  # 38
    resp, _ = await test_query_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="copy neutral", do_liftover=False)
    assert resp.dict(by_alias=True) == genomic_del1_rel_38

    q = "NC_000003.11:g.10191495del"  # 37
    resp, _ = await test_query_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="copy neutral", do_liftover=False)
    assert resp.dict(by_alias=True) == genomic_del1_rel_37

    resp, _ = await test_query_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="copy neutral", do_liftover=True)
    assert resp.dict(by_alias=True) == genomic_del1_rel_38


@pytest.mark.asyncio
async def test_genomic_del2_absolute_cnv(test_query_handler, genomic_del2_38_vac,
                                         genomic_del2_abs_37):
    """Test that genomic deletion works correctly"""
    q = "NC_000003.12:g.10146595_10146613del"  # 38
    resp, _ = await test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=None, do_liftover=False)
    assert resp.dict(by_alias=True) == genomic_del2_38_vac

    q = "NC_000003.11:g.10188279_10188297del"  # 37
    resp, _ = await test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=None, do_liftover=False)
    assert resp.dict(by_alias=True) == genomic_del2_abs_37

    resp, _ = await test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=None, do_liftover=True)
    assert resp.dict(by_alias=True) == genomic_del2_38_vac

    resp, _ = await test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=4, do_liftover=True)
    expected = copy.deepcopy(genomic_del2_38_vac)
    expected["copies"]["value"] = 3
    expected["_id"] = "ga4gh:VAC.lvGZeXqN0zTOD43dzLa9HmSjYz2-MBJf"
    assert resp.dict(by_alias=True) == expected


@pytest.mark.asyncio
async def test_genomic_del2_relative_cnv(test_query_handler, genomic_del2_rel_38,
                                         genomic_del2_rel_37):
    """Test that genomic deletion works correctly"""
    q = "NC_000003.12:g.10146595_10146613del"  # 38
    resp, _ = await test_query_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="low-level gain", do_liftover=False)
    assert resp.dict(by_alias=True) == genomic_del2_rel_38

    q = "NC_000003.11:g.10188279_10188297del"  # 37
    resp, _ = await test_query_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="low-level gain", do_liftover=False)
    assert resp.dict(by_alias=True) == genomic_del2_rel_37

    resp, _ = await test_query_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="low-level gain", do_liftover=True)
    assert resp.dict(by_alias=True) == genomic_del2_rel_38


@pytest.mark.asyncio
async def test_genomic_del3_absolute_cnv(test_query_handler, genomic_del3_abs_38,
                                         genomic_del3_abs_37):
    """Test that genomic deletion works correctly"""
    q = "NC_000023.11:g.(31060227_31100351)_(33274278_33417151)del"  # 38
    resp, _ = await test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=None, do_liftover=False)
    assert resp.dict(by_alias=True) == genomic_del3_abs_38

    q = "NC_000023.10:g.(31078344_31118468)_(33292395_33435268)del"  # 37
    resp, _ = await test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=None, do_liftover=False)
    assert resp.dict(by_alias=True) == genomic_del3_abs_37

    resp, _ = await test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=None, do_liftover=True)
    assert resp.dict(by_alias=True) == genomic_del3_abs_38

    resp, _ = await test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=2, do_liftover=True)
    expected = copy.deepcopy(genomic_del3_abs_38)
    expected["copies"] = {"value": 1, "type": "Number"}
    expected["_id"] = "ga4gh:VAC.ohgrbhR6eL_JXnXxXd-4mXTHot2HfPJQ"
    assert resp.dict(by_alias=True) == expected


@pytest.mark.asyncio
async def test_genomic_del3_relative_cnv(test_query_handler, genomic_del3_rel_38,
                                         genomic_del3_rel_37):
    """Test that genomic deletion works correctly"""
    q = "NC_000023.11:g.(31060227_31100351)_(33274278_33417151)del"  # 38
    resp, _ = await test_query_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="complete loss", do_liftover=False)
    assert resp.dict(by_alias=True) == genomic_del3_rel_38

    q = "NC_000023.10:g.(31078344_31118468)_(33292395_33435268)del"  # 37
    resp, _ = await test_query_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="complete loss", do_liftover=False)
    assert resp.dict(by_alias=True) == genomic_del3_rel_37

    resp, _ = await test_query_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="complete loss", do_liftover=True)
    assert resp.dict(by_alias=True) == genomic_del3_rel_38


@pytest.mark.asyncio
async def test_genomic_del4_absolute_cnv(test_query_handler, genomic_del4_abs_38,
                                         genomic_del4_abs_37):
    """Test that genomic deletion works correctly"""
    q = "NC_000023.11:g.(?_31120496)_(33339477_?)del"  # 38
    resp, _ = await test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=None, do_liftover=False)
    assert resp.dict(by_alias=True) == genomic_del4_abs_38

    q = "NC_000023.10:g.(?_31138613)_(33357594_?)del"  # 37
    resp, _ = await test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=None, do_liftover=False)
    assert resp.dict(by_alias=True) == genomic_del4_abs_37

    resp, _ = await test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=None, do_liftover=True)
    assert resp.dict(by_alias=True) == genomic_del4_abs_38

    resp, _ = await test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=3, do_liftover=True)
    expected = copy.deepcopy(genomic_del4_abs_38)
    expected["copies"] = {"value": 2, "type": "Number"}
    expected["_id"] = "ga4gh:VAC.QMZBxhPgywSZ-d1je0IKs_z6g0a2iN-X"
    assert resp.dict(by_alias=True) == expected


@pytest.mark.asyncio
async def test_genomic_del4_relative_cnv(test_query_handler, genomic_del4_rel_38,
                                         genomic_del4_rel_37):
    """Test that genomic deletion works correctly"""
    q = "NC_000023.11:g.(?_31120496)_(33339477_?)del"  # 38
    resp, _ = await test_query_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="partial loss", do_liftover=False)
    assert resp.dict(by_alias=True) == genomic_del4_rel_38

    q = "NC_000023.10:g.(?_31138613)_(33357594_?)del"  # 37
    resp, _ = await test_query_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="partial loss", do_liftover=False)
    assert resp.dict(by_alias=True) == genomic_del4_rel_37

    resp, _ = await test_query_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="partial loss", do_liftover=True)
    assert resp.dict(by_alias=True) == genomic_del4_rel_38


@pytest.mark.asyncio
async def test_genomic_del5_absolute_cnv(test_query_handler, genomic_del5_abs_38,
                                         genomic_del5_abs_37):
    """Test that genomic deletion works correctly"""
    q = "NC_000023.11:g.(?_18575354)_18653629del"  # 38
    resp, _ = await test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=None, do_liftover=False)
    assert resp.dict(by_alias=True) == genomic_del5_abs_38

    q = "NC_000023.10:g.(?_18593474)_18671749del"  # 37
    resp, _ = await test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=None, do_liftover=False)
    assert resp.dict(by_alias=True) == genomic_del5_abs_37

    resp, _ = await test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=None, do_liftover=True)
    assert resp.dict(by_alias=True) == genomic_del5_abs_38

    resp, _ = await test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=2, do_liftover=True)
    expected = copy.deepcopy(genomic_del5_abs_38)
    expected["copies"] = {"value": 1, "type": "Number"}
    expected["_id"] = "ga4gh:VAC.q4lg4INA2pSXs6cu3liKnJGacKofQDoB"
    assert resp.dict(by_alias=True) == expected


@pytest.mark.asyncio
async def test_genomic_del5_relative_cnv(test_query_handler, genomic_del5_rel_38,
                                         genomic_del5_rel_37):
    """Test that genomic deletion works correctly"""
    q = "NC_000023.11:g.(?_18575354)_18653629del"  # 38
    resp, _ = await test_query_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="copy neutral", do_liftover=False)
    assert resp.dict(by_alias=True) == genomic_del5_rel_38

    q = "NC_000023.10:g.(?_18593474)_18671749del"  # 37
    resp, _ = await test_query_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="copy neutral", do_liftover=False)
    assert resp.dict(by_alias=True) == genomic_del5_rel_37

    resp, _ = await test_query_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="copy neutral", do_liftover=True)
    assert resp.dict(by_alias=True) == genomic_del5_rel_38


@pytest.mark.asyncio
async def test_genomic_del6_absolute_cnv(test_query_handler, genomic_del6_abs_38,
                                         genomic_del6_abs_37):
    """Test that genomic deletion works correctly"""
    q = "NC_000006.12:g.133462764_(133464858_?)del"  # 38
    resp, _ = await test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=None, do_liftover=False)
    assert resp.dict(by_alias=True) == genomic_del6_abs_38

    q = "NC_000006.11:g.133783902_(133785996_?)del"  # 37
    resp, _ = await test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=None, do_liftover=False)
    assert resp.dict(by_alias=True) == genomic_del6_abs_37

    resp, _ = await test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=None, do_liftover=True)
    assert resp.dict(by_alias=True) == genomic_del6_abs_38

    resp, _ = await test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=3, do_liftover=True)
    expected = copy.deepcopy(genomic_del6_abs_38)
    expected["copies"]["value"] = 2
    expected["_id"] = "ga4gh:VAC.--AWt37YlF0oVaExdssLIrrLl5EEh_0k"
    assert resp.dict(by_alias=True) == expected


@pytest.mark.asyncio
async def test_genomic_del6_relative_cnv(test_query_handler, genomic_del6_rel_38,
                                         genomic_del6_rel_37):
    """Test that genomic deletion works correctly"""
    q = "NC_000006.12:g.133462764_(133464858_?)del"  # 38
    resp, _ = await test_query_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="low-level gain", do_liftover=False)
    assert resp.dict(by_alias=True) == genomic_del6_rel_38

    q = "NC_000006.11:g.133783902_(133785996_?)del"  # 37
    resp, _ = await test_query_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="low-level gain", do_liftover=False)
    assert resp.dict(by_alias=True) == genomic_del6_rel_37

    resp, _ = await test_query_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="low-level gain", do_liftover=True)
    assert resp.dict(by_alias=True) == genomic_del6_rel_38


@pytest.mark.asyncio
async def test_invalid_cnv_parameters(test_query_handler):
    """Check that invalid parameters return warnings"""
    q = "NC_000006.11:g.133783902_(133785996_?)del"
    resp, w = await test_query_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="low-level gains", do_liftover=True)
    assert resp is None
    assert w == ["low-level gains is not a valid relative copy class: ['complete loss', "  # noqa: E501
                 "'partial loss', 'copy neutral', 'low-level gain', 'high-level gain']"]


@pytest.mark.asyncio
async def test_invalid_cnv(test_query_handler):
    """Check that invalid input return warnings"""
    q = "DAG1 g.49568695dup"
    resp, w = await test_query_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="low-level gain", do_liftover=True)
    assert set(w) == {"Unable to translate DAG1 g.49568695dup to copy number variation",
                      "DAG1 g.49568695dup is not a supported HGVS genomic duplication or deletion"}  # noqa: E501
    assert resp.type == "Text"

    q = "braf v600e"
    resp, w = await test_query_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="low-level gain", do_liftover=True)
    assert set(w) == {"Unable to translate braf v600e to copy number variation",
                      "braf v600e is not a supported HGVS genomic duplication or deletion"}  # noqa: E501
    assert resp.type == "Text"
