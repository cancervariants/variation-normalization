"""Module for testing the hgvs to absolute and relative copy number endpoints"""
import pytest


@pytest.fixture(scope="module")
def genomic_dup1_38_subject():
    """Create test fixture GRCh38 duplication subject"""
    return {
        "location": {
            "_id": "ga4gh:VSL.G_J9WrfooiONRgjbmGPuCBYbBYFQnYOg",
            "sequence_id": "ga4gh:SQ.Zu7h9AggXxhTaGVsy7h_EZSChSZGcmgX",
            "interval": {
                "type": "SequenceInterval",
                "start": {"value": 49531260, "type": "Number"},
                "end": {"value": 49531262, "type": "Number"},
            },
            "type": "SequenceLocation",
        },
        "reverse_complement": False,
        "type": "DerivedSequenceExpression"
    }


@pytest.fixture(scope="module")
def genomic_dup1_abs_38(genomic_dup1_38_subject):
    """Create test fixture absolute copy number variation"""
    return {
        "type": "CopyNumber",
        "_id": "ga4gh:VCN.KdBguJLeiXM2yr3JaRQ2kxLxaAd4pPlq",
        "subject": genomic_dup1_38_subject,
        "copies": {
            "type": "Number",
            "value": 3
        }
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
        "type": "CopyNumber",
        "_id": "ga4gh:VCN.XIojCJVUI8SQFlV8jk7V-CewsQy8PQa1",
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
def genomic_dup2_38_subject():
    """Create test fixture GRCh38 duplication subject"""
    return {
        "location": {
            "_id": "ga4gh:VSL.4mH68huylkPmu6zyUwH4wiazIYr9cQUX",
            "sequence_id": "ga4gh:SQ.yC_0RBj3fgBlvgyAuycbzdubtLxq-rE0",
            "interval": {
                "type": "SequenceInterval",
                "start": {"value": 2087937, "type": "Number"},
                "end": {"value": 2087948, "type": "Number"},
            },
            "type": "SequenceLocation",
        },
        "reverse_complement": False,
        "type": "DerivedSequenceExpression"
    }


@pytest.fixture(scope="module")
def genomic_dup2_abs_38(genomic_dup2_38_subject):
    """Create test fixture absolute copy number variation"""
    return {
        "type": "CopyNumber",
        "_id": "ga4gh:VCN.rd1wobb8NXRxk6O__njJUQg_ekZUALGx",
        "subject": genomic_dup2_38_subject,
        "copies": {
            "type": "Number",
            "value": 3
        }
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
        "type": "CopyNumber",
        "_id": "ga4gh:VCN.Q4-OUpSNO6SXY4EyWEKehneRpHgjmpZo",
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
def genomic_dup3_38_subject():
    """Create test fixture GRCh38 duplication subject"""
    return {
        "location": {
            "_id": "ga4gh:VSL.DgEMxYt1AdPe-HZAQbT2AVz5OejICnOj",
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
    }


@pytest.fixture(scope="module")
def genomic_dup3_abs_38(genomic_dup3_38_subject):
    """Create test fixture absolute copy number variation"""
    return {
        "type": "CopyNumber",
        "_id": "ga4gh:VCN.IgQATuKrM_J5MDHm2VemKThFOkzz-7AZ",
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
        "type": "CopyNumber",
        "_id": "ga4gh:VCN.tGhEnpnYFCtcwq8oZRzbmp2ekCypl_PW",
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
def genomic_dup4_38_subject():
    """Create test fixture GRCh38 duplication subject"""
    return {
        "location": {
            "_id": "ga4gh:VSL.us51izImAQQWr-Hu6Q7HQm-vYvmb-jJo",
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
    }


@pytest.fixture(scope="module")
def genomic_dup4_abs_38(genomic_dup4_38_subject):
    """Create test fixture absolute copy number variation"""
    return {
        "type": "CopyNumber",
        "_id": "ga4gh:VCN.3rvfUmiIb4hSxVQhXKOonuOY6Q3xTkKx",
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
        "type": "CopyNumber",
        "_id": "ga4gh:VCN.bGinS0pRwqxl1OKRDdJ8itpGRMAL9O6O",
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
def genomic_dup5_38_subject():
    """Create test fixture GRCh38 duplication subject"""
    return {
        "location": {
            "_id": "ga4gh:VSL.k2FXLyqyS8pbtZxEHCpNd2SHD6iCtH9C",
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
    }


@pytest.fixture(scope="module")
def genomic_dup5_abs_38(genomic_dup5_38_subject):
    """Create test fixture absolute copy number variation"""
    return {
        "type": "CopyNumber",
        "_id": "ga4gh:VCN.eLAZZ-ht1h2dTtZqzhO9TVhBdFufv67-",
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
        "type": "CopyNumber",
        "_id": "ga4gh:VCN.AVN4KtF50NqzkOEX9On7z5TDEBo6t_CN",
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
def genomic_dup6_38_subject():
    """Create test fixture GRCh38 duplication subject"""
    return {
        "location": {
            "_id": "ga4gh:VSL.h0_xXu36uSnPEuLoxvVmTAFQCS1ZFuLN",
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
    }


@pytest.fixture(scope="module")
def genomic_dup6_abs_38(genomic_dup6_38_subject):
    """Create test fixture absolute copy number variation"""
    return {
        "type": "CopyNumber",
        "_id": "ga4gh:VCN.Rekk_MmUQ777V76S51x7nZGjh4U3LkLy",
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
        "type": "CopyNumber",
        "_id": "ga4gh:VCN.j1rT3jTVfMzblHZqGkQjMl3BhGzMa2ub",
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
def genomic_del1_38_subject():
    """Create test fixture GRCh38 deletion subject"""
    return {
        "location": {
            "_id": "ga4gh:VSL.Yg5B66zErDjK9Lqeaw-kuzAB9w5-uUaS",
            "sequence_id": "ga4gh:SQ.Zu7h9AggXxhTaGVsy7h_EZSChSZGcmgX",
            "interval": {
                "type": "SequenceInterval",
                "start": {"value": 10149810, "type": "Number"},
                "end": {"value": 10149811, "type": "Number"},
            },
            "type": "SequenceLocation",
        },
        "reverse_complement": False,
        "type": "DerivedSequenceExpression"
    }


@pytest.fixture(scope="module")
def genomic_del1_abs_38(genomic_del1_38_subject):
    """Create test fixture absolute copy number variation"""
    return {
        "type": "CopyNumber",
        "_id": "ga4gh:VCN._Iv1RBu8ctlHOaobb4emjxwbxPdkBIVF",
        "subject": genomic_del1_38_subject,
        "copies": {
            "type": "Number",
            "value": 1
        }
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
        "type": "CopyNumber",
        "_id": "ga4gh:VCN.ptB4eMIbnOpFqXRDbT5jQzcZVdyulFyi",
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
def genomic_del2_38_subject():
    """Create test fixture GRCh38 deletion subject"""
    return {
        "location": {
            "_id": "ga4gh:VSL.lksYAhEQvP8biy_nxoOJ_Zwu75a_kYtQ",
            "sequence_id": "ga4gh:SQ.Zu7h9AggXxhTaGVsy7h_EZSChSZGcmgX",
            "interval": {
                "type": "SequenceInterval",
                "start": {"value": 10146594, "type": "Number"},
                "end": {"value": 10146613, "type": "Number"},
            },
            "type": "SequenceLocation",
        },
        "reverse_complement": False,
        "type": "DerivedSequenceExpression"
    }


@pytest.fixture(scope="module")
def genomic_del2_abs_38(genomic_del2_38_subject):
    """Create test fixture absolute copy number variation"""
    return {
        "type": "CopyNumber",
        "_id": "ga4gh:VCN.gBHXvaw64pQg04DAhp_Gtzh8ADUf7HuI",
        "subject": genomic_del2_38_subject,
        "copies": {
            "type": "Number",
            "value": 1
        }
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
        "type": "CopyNumber",
        "_id": "ga4gh:VCN.XUiX4XhQy_VkW_7ShRbCDfdRRb685Ofm",
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
def genomic_del3_38_subject():
    """Create test fixture GRCh38 deletion subject"""
    return {
        "location": {
            "_id": "ga4gh:VSL.DgEMxYt1AdPe-HZAQbT2AVz5OejICnOj",
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
    }


@pytest.fixture(scope="module")
def genomic_del3_abs_38(genomic_del3_38_subject):
    """Create test fixture absolute copy number variation"""
    return {
        "type": "CopyNumber",
        "_id": "ga4gh:VCN.9h2LkajTwHBdXYMRyrD9HkYwU9d7fIBr",
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
        "type": "CopyNumber",
        "_id": "ga4gh:VCN.BIRFSIemPTP_HfRQl5o7DnRbuSGfivnz",
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
def genomic_del4_38_subject():
    """Create test fixture GRCh38 deletion subject"""
    return {
        "location": {
            "_id": "ga4gh:VSL.7OJ5EFgu_2C4zPFDUBgn-ziE6BZwsRcv",
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
    }


@pytest.fixture(scope="module")
def genomic_del4_abs_38(genomic_del4_38_subject):
    """Create test fixture absolute copy number variation"""
    return {
        "type": "CopyNumber",
        "_id": "ga4gh:VCN.yQJnQz12MXlZGWx6BuzccVGrCCic_tMk",
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
        "type": "CopyNumber",
        "_id": "ga4gh:VCN.gRshXhruFQw-QdKwU4xc2iKBNLIbFNzt",
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
def genomic_del5_38_subject():
    """Create test fixture GRCh38 deletion subject"""
    return {
        "location": {
            "_id": "ga4gh:VSL.jURzcCBf3kJVx19uuJJtwt78LuBbtfwD",
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
    }


@pytest.fixture(scope="module")
def genomic_del5_abs_38(genomic_del5_38_subject):
    """Create test fixture absolute copy number variation"""
    return {
        "type": "CopyNumber",
        "_id": "ga4gh:VCN._RIw5UC5bZeLeHnBLYAow7Ml-lv2nKJW",
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
        "type": "CopyNumber",
        "_id": "ga4gh:VCN.SwAiS0A3UJ6Up-ghFm1hKIeka2LOeQH3",
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
def genomic_del6_38_subject():
    """Create test fixture GRCh38 deletion subject"""
    return {
        "location": {
            "_id": "ga4gh:VSL.TPwsB5ymsNI7TynTlI8_8CI_NmNrBHUQ",
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
    }


@pytest.fixture(scope="module")
def genomic_del6_abs_38(genomic_del6_38_subject):
    """Create test fixture absolute copy number variation"""
    return {
        "type": "CopyNumber",
        "_id": "ga4gh:VCN.F3U6Rmov1WO2mhmRHWumJb-YALOMkeeI",
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
        "type": "CopyNumber",
        "_id": "ga4gh:VCN.tb5MqxC7Ljh7T1ZTg81ClfMyUnZuxKZl",
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


def test_genomic_dup1_absolute_cnv(test_query_handler, genomic_dup1_abs_38,
                                   genomic_dup1_abs_37):
    """Test that genomic duplication works correctly"""
    q = "NC_000003.12:g.49531262dup"  # 38
    resp, w = test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=None, do_liftover=False)
    assert resp == genomic_dup1_abs_38

    q = "NC_000003.11:g.49568695dup"  # 37
    resp, w = test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=None, do_liftover=False)
    assert resp == genomic_dup1_abs_37

    resp, w = test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=None, do_liftover=True)
    assert resp == genomic_dup1_abs_38


def test_genomic_dup1_relative_cnv(test_query_handler, genomic_dup1_rel_38,
                                   genomic_dup1_rel_37):
    """Test that genomic duplication works correctly"""
    q = "NC_000003.12:g.49531262dup"  # 38
    resp, w = test_query_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="complete loss", do_liftover=False)
    assert resp == genomic_dup1_rel_38

    q = "NC_000003.11:g.49568695dup"  # 37
    resp, w = test_query_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="complete loss", do_liftover=False)
    assert resp == genomic_dup1_rel_37

    resp, w = test_query_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="complete loss", do_liftover=True)
    assert resp == genomic_dup1_rel_38


def test_genomic_dup2_absolute_cnv(test_query_handler, genomic_dup2_abs_38,
                                   genomic_dup2_abs_37):
    """Test that genomic duplication works correctly"""
    q = "NC_000016.10:g.2087938_2087948dup"  # 38
    resp, w = test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=None, do_liftover=False)
    assert resp == genomic_dup2_abs_38

    q = "NC_000016.9:g.2137939_2137949dup"  # 37
    resp, w = test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=None, do_liftover=False)
    assert resp == genomic_dup2_abs_37

    resp, w = test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=None, do_liftover=True)
    assert resp == genomic_dup2_abs_38


def test_genomic_dup2_relative_cnv(test_query_handler, genomic_dup2_rel_38,
                                   genomic_dup2_rel_37):
    """Test that genomic duplication works correctly"""
    q = "NC_000016.10:g.2087938_2087948dup"  # 38
    resp, w = test_query_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="partial loss", do_liftover=False)
    assert resp == genomic_dup2_rel_38

    q = "NC_000016.9:g.2137939_2137949dup"  # 37
    resp, w = test_query_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="partial loss", do_liftover=False)
    assert resp == genomic_dup2_rel_37

    resp, w = test_query_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="partial loss", do_liftover=True)
    assert resp == genomic_dup2_rel_38


def test_genomic_dup3_absolute_cnv(test_query_handler, genomic_dup3_abs_38,
                                   genomic_dup3_abs_37):
    """Test that genomic duplication works correctly"""
    q = "NC_000023.11:g.(31060227_31100351)_(33274278_33417151)dup"  # 38
    resp, w = test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=None, do_liftover=False)
    assert resp == genomic_dup3_abs_38

    q = "NC_000023.10:g.(31078344_31118468)_(33292395_33435268)dup"  # 37
    resp, w = test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=None, do_liftover=False)
    assert resp == genomic_dup3_abs_37

    resp, w = test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=None, do_liftover=True)
    assert resp == genomic_dup3_abs_38


def test_genomic_dup3_relative_cnv(test_query_handler, genomic_dup3_rel_38,
                                   genomic_dup3_rel_37):
    """Test that genomic duplication works correctly"""
    q = "NC_000023.11:g.(31060227_31100351)_(33274278_33417151)dup"  # 38
    resp, w = test_query_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="high-level gain", do_liftover=False)
    assert resp == genomic_dup3_rel_38

    q = "NC_000023.10:g.(31078344_31118468)_(33292395_33435268)dup"  # 37
    resp, w = test_query_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="high-level gain", do_liftover=False)
    assert resp == genomic_dup3_rel_37

    resp, w = test_query_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="high-level gain", do_liftover=True)
    assert resp == genomic_dup3_rel_38


def test_genomic_dup4_absolute_cnv(test_query_handler, genomic_dup4_abs_38,
                                   genomic_dup4_abs_37):
    """Test that genomic duplication works correctly"""
    q = "NC_000020.11:g.(?_30417576)_(31394018_?)dup"  # 38
    resp, w = test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=None, do_liftover=False)
    assert resp == genomic_dup4_abs_38

    q = "NC_000020.10:g.(?_29652252)_(29981821_?)dup"  # 37
    resp, w = test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=None, do_liftover=False)
    assert resp == genomic_dup4_abs_37

    resp, w = test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=None, do_liftover=True)
    assert resp == genomic_dup4_abs_38


def test_genomic_dup4_relative_cnv(test_query_handler, genomic_dup4_rel_38,
                                   genomic_dup4_rel_37):
    """Test that genomic duplication works correctly"""
    q = "NC_000020.11:g.(?_30417576)_(31394018_?)dup"  # 38
    resp, w = test_query_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="complete loss", do_liftover=False)
    assert resp == genomic_dup4_rel_38

    q = "NC_000020.10:g.(?_29652252)_(29981821_?)dup"  # 37
    resp, w = test_query_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="complete loss", do_liftover=False)
    assert resp == genomic_dup4_rel_37

    resp, w = test_query_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="complete loss", do_liftover=True)
    assert resp == genomic_dup4_rel_38


def test_genomic_dup5_absolute_cnv(test_query_handler, genomic_dup5_abs_38,
                                   genomic_dup5_abs_37):
    """Test that genomic duplication works correctly"""
    q = "NC_000023.11:g.(?_154021812)_154092209dup"  # 38
    resp, w = test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=None, do_liftover=False)
    assert resp == genomic_dup5_abs_38

    q = "NC_000023.10:g.(?_153287263)_153357667dup"  # 37
    resp, w = test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=None, do_liftover=False)
    assert resp == genomic_dup5_abs_37

    resp, w = test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=None, do_liftover=True)
    assert resp == genomic_dup5_abs_38


def test_genomic_dup5_relative_cnv(test_query_handler, genomic_dup5_rel_38,
                                   genomic_dup5_rel_37):
    """Test that genomic duplication works correctly"""
    q = "NC_000023.11:g.(?_154021812)_154092209dup"  # 38
    resp, w = test_query_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="partial loss", do_liftover=False)
    assert resp == genomic_dup5_rel_38

    q = "NC_000023.10:g.(?_153287263)_153357667dup"  # 37
    resp, w = test_query_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="partial loss", do_liftover=False)
    assert resp == genomic_dup5_rel_37

    resp, w = test_query_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="partial loss", do_liftover=True)
    assert resp == genomic_dup5_rel_38


def test_genomic_dup6_absolute_cnv(test_query_handler, genomic_dup6_abs_38,
                                   genomic_dup6_abs_37):
    """Test that genomic duplication works correctly"""
    q = "NC_000023.11:g.154021812_(154092209_?)dup"  # 38
    resp, w = test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=None, do_liftover=False)
    assert resp == genomic_dup6_abs_38

    q = "NC_000023.10:g.153287263_(153357667_?)dup"  # 37
    resp, w = test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=None, do_liftover=False)
    assert resp == genomic_dup6_abs_37

    resp, w = test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=None, do_liftover=True)
    assert resp == genomic_dup6_abs_38


def test_genomic_dup6_relative_cnv(test_query_handler, genomic_dup6_rel_38,
                                   genomic_dup6_rel_37):
    """Test that genomic duplication works correctly"""
    q = "NC_000023.11:g.154021812_(154092209_?)dup"  # 38
    resp, w = test_query_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="copy neutral", do_liftover=False)
    assert resp == genomic_dup6_rel_38

    q = "NC_000023.10:g.153287263_(153357667_?)dup"  # 37
    resp, w = test_query_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="copy neutral", do_liftover=False)
    assert resp == genomic_dup6_rel_37

    resp, w = test_query_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="copy neutral", do_liftover=True)
    assert resp == genomic_dup6_rel_38


def test_genomic_del1_absolute_cnv(test_query_handler, genomic_del1_abs_38,
                                   genomic_del1_abs_37):
    """Test that genomic deletion works correctly"""
    q = "NC_000003.12:g.10149811del"  # 38
    resp, w = test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=None, do_liftover=False)
    assert resp == genomic_del1_abs_38

    q = "NC_000003.11:g.10191495del"  # 37
    resp, w = test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=None, do_liftover=False)
    assert resp == genomic_del1_abs_37

    resp, w = test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=None, do_liftover=True)
    assert resp == genomic_del1_abs_38


def test_genomic_del1_relative_cnv(test_query_handler, genomic_del1_rel_38,
                                   genomic_del1_rel_37):
    """Test that genomic deletion works correctly"""
    q = "NC_000003.12:g.10149811del"  # 38
    resp, w = test_query_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="copy neutral", do_liftover=False)
    assert resp == genomic_del1_rel_38

    q = "NC_000003.11:g.10191495del"  # 37
    resp, w = test_query_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="copy neutral", do_liftover=False)
    assert resp == genomic_del1_rel_37

    resp, w = test_query_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="copy neutral", do_liftover=True)
    assert resp == genomic_del1_rel_38


def test_genomic_del2_absolute_cnv(test_query_handler, genomic_del2_abs_38,
                                   genomic_del2_abs_37):
    """Test that genomic deletion works correctly"""
    q = "NC_000003.12:g.10146595_10146613del"  # 38
    resp, w = test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=None, do_liftover=False)
    assert resp == genomic_del2_abs_38

    q = "NC_000003.11:g.10188279_10188297del"  # 37
    resp, w = test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=None, do_liftover=False)
    assert resp == genomic_del2_abs_37

    resp, w = test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=None, do_liftover=True)
    assert resp == genomic_del2_abs_38


def test_genomic_del2_relative_cnv(test_query_handler, genomic_del2_rel_38,
                                   genomic_del2_rel_37):
    """Test that genomic deletion works correctly"""
    q = "NC_000003.12:g.10146595_10146613del"  # 38
    resp, w = test_query_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="low-level gain", do_liftover=False)
    assert resp == genomic_del2_rel_38

    q = "NC_000003.11:g.10188279_10188297del"  # 37
    resp, w = test_query_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="low-level gain", do_liftover=False)
    assert resp == genomic_del2_rel_37

    resp, w = test_query_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="low-level gain", do_liftover=True)
    assert resp == genomic_del2_rel_38


def test_genomic_del3_absolute_cnv(test_query_handler, genomic_del3_abs_38,
                                   genomic_del3_abs_37):
    """Test that genomic deletion works correctly"""
    q = "NC_000023.11:g.(31060227_31100351)_(33274278_33417151)del"  # 38
    resp, w = test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=None, do_liftover=False)
    assert resp == genomic_del3_abs_38

    q = "NC_000023.10:g.(31078344_31118468)_(33292395_33435268)del"  # 37
    resp, w = test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=None, do_liftover=False)
    assert resp == genomic_del3_abs_37

    resp, w = test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=None, do_liftover=True)
    assert resp == genomic_del3_abs_38


def test_genomic_del3_relative_cnv(test_query_handler, genomic_del3_rel_38,
                                   genomic_del3_rel_37):
    """Test that genomic deletion works correctly"""
    q = "NC_000023.11:g.(31060227_31100351)_(33274278_33417151)del"  # 38
    resp, w = test_query_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="complete loss", do_liftover=False)
    assert resp == genomic_del3_rel_38

    q = "NC_000023.10:g.(31078344_31118468)_(33292395_33435268)del"  # 37
    resp, w = test_query_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="complete loss", do_liftover=False)
    assert resp == genomic_del3_rel_37

    resp, w = test_query_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="complete loss", do_liftover=True)
    assert resp == genomic_del3_rel_38


def test_genomic_del4_absolute_cnv(test_query_handler, genomic_del4_abs_38,
                                   genomic_del4_abs_37):
    """Test that genomic deletion works correctly"""
    q = "NC_000023.11:g.(?_31120496)_(33339477_?)del"  # 38
    resp, w = test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=None, do_liftover=False)
    assert resp == genomic_del4_abs_38

    q = "NC_000023.10:g.(?_31138613)_(33357594_?)del"  # 37
    resp, w = test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=None, do_liftover=False)
    assert resp == genomic_del4_abs_37

    resp, w = test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=None, do_liftover=True)
    assert resp == genomic_del4_abs_38


def test_genomic_del4_relative_cnv(test_query_handler, genomic_del4_rel_38,
                                   genomic_del4_rel_37):
    """Test that genomic deletion works correctly"""
    q = "NC_000023.11:g.(?_31120496)_(33339477_?)del"  # 38
    resp, w = test_query_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="partial loss", do_liftover=False)
    assert resp == genomic_del4_rel_38

    q = "NC_000023.10:g.(?_31138613)_(33357594_?)del"  # 37
    resp, w = test_query_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="partial loss", do_liftover=False)
    assert resp == genomic_del4_rel_37

    resp, w = test_query_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="partial loss", do_liftover=True)
    assert resp == genomic_del4_rel_38


def test_genomic_del5_absolute_cnv(test_query_handler, genomic_del5_abs_38,
                                   genomic_del5_abs_37):
    """Test that genomic deletion works correctly"""
    q = "NC_000023.11:g.(?_18575354)_18653629del"  # 38
    resp, w = test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=None, do_liftover=False)
    assert resp == genomic_del5_abs_38

    q = "NC_000023.10:g.(?_18593474)_18671749del"  # 37
    resp, w = test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=None, do_liftover=False)
    assert resp == genomic_del5_abs_37

    resp, w = test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=None, do_liftover=True)
    assert resp == genomic_del5_abs_38


def test_genomic_del5_relative_cnv(test_query_handler, genomic_del5_rel_38,
                                   genomic_del5_rel_37):
    """Test that genomic deletion works correctly"""
    q = "NC_000023.11:g.(?_18575354)_18653629del"  # 38
    resp, w = test_query_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="copy neutral", do_liftover=False)
    assert resp == genomic_del5_rel_38

    q = "NC_000023.10:g.(?_18593474)_18671749del"  # 37
    resp, w = test_query_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="copy neutral", do_liftover=False)
    assert resp == genomic_del5_rel_37

    resp, w = test_query_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="copy neutral", do_liftover=True)
    assert resp == genomic_del5_rel_38


def test_genomic_del6_absolute_cnv(test_query_handler, genomic_del6_abs_38,
                                   genomic_del6_abs_37):
    """Test that genomic deletion works correctly"""
    q = "NC_000006.12:g.133462764_(133464858_?)del"  # 38
    resp, w = test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=None, do_liftover=False)
    assert resp == genomic_del6_abs_38

    q = "NC_000006.11:g.133783902_(133785996_?)del"  # 37
    resp, w = test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=None, do_liftover=False)
    assert resp == genomic_del6_abs_37

    resp, w = test_query_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=None, do_liftover=True)
    assert resp == genomic_del6_abs_38


def test_genomic_del6_relative_cnv(test_query_handler, genomic_del6_rel_38,
                                   genomic_del6_rel_37):
    """Test that genomic deletion works correctly"""
    q = "NC_000006.12:g.133462764_(133464858_?)del"  # 38
    resp, w = test_query_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="low-level gain", do_liftover=False)
    assert resp == genomic_del6_rel_38

    q = "NC_000006.11:g.133783902_(133785996_?)del"  # 37
    resp, w = test_query_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="low-level gain", do_liftover=False)
    assert resp == genomic_del6_rel_37

    resp, w = test_query_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="low-level gain", do_liftover=True)
    assert resp == genomic_del6_rel_38

# def test_absolute_cnv_parameters():
#     """Check that invalid parameters return warnings"""
#     pass

# def test_relative_cnv_parameters():
#     """Check that invalid parameters return warnings"""
#     pass

# def test_invalid_cnv(test_query_handler):
#     """Check that invalid input return warnings"""
#     pass
