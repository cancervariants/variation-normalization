"""Module for testing the hgvs to copy number count and copy number change endpoints"""
import copy

import pytest
from ga4gh.vrs import models

from tests.conftest import cnv_assertion_checks


@pytest.fixture(scope="module")
def genomic_dup1_cx_38(genomic_dup1_seq_loc_not_normalized):
    """Create test fixture copy number change variation"""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.t2ng-5-owqfItCFyz1uA_xNw3FSZTlYo",
        "subject": genomic_dup1_seq_loc_not_normalized,
        "copyChange": "efo:0030069",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_dup1_37_loc():
    """Create test fixture GRCh37 duplication subject"""
    return {
        "id": "ga4gh:SL.fvq_Ta_h7VUfmV_z45Fli7jrALAhqrk_",
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.VNBualIltAyi2AI_uXcKU7M9XUOuA7MS",
        },
        "start": 49568694,
        "end": 49568695,
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_dup1_cn_37(genomic_dup1_37_loc):
    """Create test fixture copy number count variation (not normalized)"""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN._aUTJcDsXX1uNJF1TpahgmjnndKBwsVm",
        "subject": genomic_dup1_37_loc,
        "copies": 3,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_dup1_cx_37(genomic_dup1_37_loc):
    """Create test fixture copy number change variation"""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.5MwKOwU8GqzAkIaVWeA9plU9mYviJHWW",
        "subject": genomic_dup1_37_loc,
        "copyChange": "efo:0030069",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_dup2_cx_38(genomic_dup2_seq_loc_normalized):
    """Create test fixture copy number change variation"""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.J92RlKs7pjcjRF5fzjsj43uNLl7YLm0O",
        "subject": genomic_dup2_seq_loc_normalized,
        "copyChange": "efo:0030067",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_dup2_37_loc():
    """Create test fixture GRCh37 duplication subject"""
    return {
        "id": "ga4gh:SL.dUhZPlPAMN4B1zt7Ww9E4lMuvHwJ-aBb",
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.v7noePfnNpK8ghYXEqZ9NukMXW7YeNsm",
        },
        "start": 33229406,
        "end": 33229410,
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_dup2_cn_37(genomic_dup2_37_loc):
    """Create test fixture copy number count variation"""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.N3Vm6lp21K3VqX8HlKNyKVaqHAHxoP01",
        "subject": genomic_dup2_37_loc,
        "copies": 3,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_dup2_cx_37(genomic_dup2_37_loc):
    """Create test fixture copy number change variation"""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.rbg77MAqTjtZ3qgpO6hq1E_1V-vufD_y",
        "subject": genomic_dup2_37_loc,
        "copyChange": "efo:0030067",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_dup3_cn_38(genomic_del3_dup3_loc_not_normalized):
    """Create test fixture copy number count variation"""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.BEAOoq3yyFS_3Mjy7u2d557vfS_EEtuM",
        "subject": genomic_del3_dup3_loc_not_normalized,
        "copies": 2,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_dup3_cx_38(genomic_del3_dup3_loc_not_normalized):
    """Create test fixture copy number change variation"""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.6GLzav8IHczludUusuNX-HYuf6Z8VTgv",
        "subject": genomic_del3_dup3_loc_not_normalized,
        "copyChange": "efo:0030072",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_dup3_37_loc():
    """Create test fixture GRCh37 duplication subject"""
    return {
        "id": "ga4gh:SL.VpH1hrh1I_xHHQHGbrexGa0VF33c7awQ",
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.v7noePfnNpK8ghYXEqZ9NukMXW7YeNsm",
        },
        "start": [31078343, 31118467],
        "end": [33292396, 33435269],
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_dup3_cn_37(genomic_dup3_37_loc):
    """Create test fixture copy number count variation"""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.VgkxQW-HOSeNKmz23B_g8RpCI4olkDp7",
        "subject": genomic_dup3_37_loc,
        "copies": 2,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_dup3_cx_37(genomic_dup3_37_loc):
    """Create test fixture copy number change variation"""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.A1I4UrYZEg0zX8R2alvBypgstReRDna3",
        "subject": genomic_dup3_37_loc,
        "copyChange": "efo:0030072",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_dup4_cn_38(genomic_dup4_loc):
    """Create test fixture copy number count variation"""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.IGupUKw7D4ErAXoMIL9gcuE-aIhqGHCn",
        "subject": genomic_dup4_loc,
        "copies": 3,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_dup4_cx_38(genomic_dup4_loc):
    """Create test fixture copy number change variation"""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.-o60HFLhclY_B-J_hxd-uWbOIAgvjzNZ",
        "subject": genomic_dup4_loc,
        "copyChange": "efo:0030069",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_dup4_37_loc():
    """Create test fixture GRCh37 duplication subject"""
    return {
        "id": "ga4gh:SL.3npNZqZCncn_X1IYWr7Ch0pBE2iAFvfU",
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.iy_UbUrvECxFRX5LPTH_KPojdlT7BKsf",
        },
        "start": [None, 29652251],
        "end": [29981821, None],
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_dup4_cn_37(genomic_dup4_37_loc):
    """Create test fixture copy number change variation"""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.COBsV4umQseKdzN3QEc431hpVytFJAaM",
        "subject": genomic_dup4_37_loc,
        "copies": 3,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_dup4_cx_37(genomic_dup4_37_loc):
    """Create test fixture copy number change variation"""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.XkcOlZOMrgPoo75reJxKqiQm1bWH6RKG",
        "subject": genomic_dup4_37_loc,
        "copyChange": "efo:0030069",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_dup5_cn_38(genomic_dup5_loc):
    """Create test fixture copy number count variation"""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.t8BziFJ42omY_ojbY9NphXxORWAkx_oC",
        "subject": genomic_dup5_loc,
        "copies": 4,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_dup5_cx_38(genomic_dup5_loc):
    """Create test fixture copy number change variation"""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.VK_C9npe_G00xzeIbLVSWsLjuHmkIkP8",
        "subject": genomic_dup5_loc,
        "copyChange": "efo:0030067",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_dup5_37_loc():
    """Create test fixture GRCh37 duplication subject"""
    return {
        "id": "ga4gh:SL.HGR3fSQPnhcTUnWc27LqqxXCiNCAtR4K",
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.v7noePfnNpK8ghYXEqZ9NukMXW7YeNsm",
        },
        "start": [None, 153287262],
        "end": 153357667,
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_dup5_cn_37(genomic_dup5_37_loc):
    """Create test fixture copy number count variation"""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.AzqR_6YYj7a6H5ZbyR8FLwPVPzf4NpuU",
        "subject": genomic_dup5_37_loc,
        "copies": 4,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_dup5_cx_37(genomic_dup5_37_loc):
    """Create test fixture copy number change variation"""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.9ps5hp71hwVJ_FGx4Mq2fvIh9CBI-jDy",
        "subject": genomic_dup5_37_loc,
        "copyChange": "efo:0030067",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_dup6_cn_38(genomic_dup6_loc):
    """Create test fixture copy number count variation"""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.6tOzN9-2xi301VUy0Rf_X7UzhS5rx5VS",
        "subject": genomic_dup6_loc,
        "copies": 2,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_dup6_cx_38(genomic_dup6_loc):
    """Create test fixture copy number change variation"""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.LzLOEMjzFoS2KrSIjp9kKWvDFDMkom6F",
        "subject": genomic_dup6_loc,
        "copyChange": "efo:0030064",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_dup6_37_loc():
    """Create test fixture GRCh37 duplication subject"""
    return {
        "id": "ga4gh:SL.0y_Tb1ed8YpOjD2NeiB__HGXV71bx6lv",
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.v7noePfnNpK8ghYXEqZ9NukMXW7YeNsm",
        },
        "start": 153287262,
        "end": [153357667, None],
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_dup6_cn_37(genomic_dup6_37_loc):
    """Create test fixture copy number count variation"""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.H8iMhIPfP4K3yrog6UG_QRbSld2-kC90",
        "subject": genomic_dup6_37_loc,
        "copies": 2,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_dup6_cx_37(genomic_dup6_37_loc):
    """Create test fixture copy number change variation"""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.3lAcmnj47emOukkuWnzYrlqcwxOcxiG9",
        "subject": genomic_dup6_37_loc,
        "copyChange": "efo:0030064",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_del1_cx_38(genomic_del1_seq_loc):
    """Create test fixture copy number change variation"""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.yG6nbQytOHuPMn3Jx8Av2bJBK3qvjhNF",
        "subject": genomic_del1_seq_loc,
        "copyChange": "efo:0030064",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_del1_37_loc():
    """Create test fixture GRCh37 deletion subject"""
    return {
        "id": "ga4gh:SL.jFRhLvxA6d-UrFd8Z-gBQDkjAfoyAxwY",
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.VNBualIltAyi2AI_uXcKU7M9XUOuA7MS",
        },
        "start": 10191494,
        "end": 10191495,
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_del1_cn_37(genomic_del1_37_loc):
    """Create test fixture copy number count variation"""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.yCFnwvrFAOYDf3q3v6oiTx_MPChQ2mZQ",
        "subject": genomic_del1_37_loc,
        "copies": 1,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_del1_cx_37(genomic_del1_37_loc):
    """Create test fixture copy number change variation"""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.gfkDdqRnbJQNFx1oiueBro3hKTIHboyY",
        "subject": genomic_del1_37_loc,
        "copyChange": "efo:0030064",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_del2_cx_38(genomic_del2_seq_loc):
    """Create test fixture copy number change variation"""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.uPLKIvgFxL0VeDpquZ0zm8_CQFcYwlK7",
        "subject": genomic_del2_seq_loc,
        "copyChange": "efo:0030071",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_del2_37_loc():
    """Create test fixture GRCh37 deletion subject"""
    return {
        "id": "ga4gh:SL.m1eIfrgB4ZfrHuM5N7dAGWuFAFU39dso",
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.VNBualIltAyi2AI_uXcKU7M9XUOuA7MS",
        },
        "start": 10188278,
        "end": 10188297,
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_del2_cn_37(genomic_del2_37_loc):
    """Create test fixture copy number count variation"""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.XOcI_zOSqtlVEDc4nMYqTL2T6mZHVNdC",
        "subject": genomic_del2_37_loc,
        "copies": 1,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_del2_cx_37(genomic_del2_37_loc):
    """Create test fixture copy number change variation"""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.JTwe-amwER_ZCEPgxkurRWoTALEgokI-",
        "subject": genomic_del2_37_loc,
        "copyChange": "efo:0030071",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_del3_cn_38(genomic_del3_dup3_loc_not_normalized):
    """Create test fixture copy number count variation"""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.BEAOoq3yyFS_3Mjy7u2d557vfS_EEtuM",
        "subject": genomic_del3_dup3_loc_not_normalized,
        "copies": 2,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_del3_cx_38(genomic_del3_dup3_loc_not_normalized):
    """Create test fixture copy number change variation"""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.d_9kg95iLP6vfqD64ZpwdKFQ66lQaSGs",
        "subject": genomic_del3_dup3_loc_not_normalized,
        "copyChange": "efo:0030069",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_del3_37_loc():
    """Create test fixture GRCh37 deletion subject"""
    return {
        "id": "ga4gh:SL.VpH1hrh1I_xHHQHGbrexGa0VF33c7awQ",
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.v7noePfnNpK8ghYXEqZ9NukMXW7YeNsm",
        },
        "start": [31078343, 31118467],
        "end": [33292396, 33435269],
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_del3_cn_37(genomic_del3_37_loc):
    """Create test fixture copy number count variation"""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.VgkxQW-HOSeNKmz23B_g8RpCI4olkDp7",
        "subject": genomic_del3_37_loc,
        "copies": 2,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_del3_cx_37(genomic_del3_37_loc):
    """Create test fixture copy number change variation"""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.ubn9DGt10rIRsdDwotP-IjeKxeS5bjyW",
        "subject": genomic_del3_37_loc,
        "copyChange": "efo:0030069",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_del4_cn_38(genomic_del4_seq_loc):
    """Create test fixture copy number count variation"""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.iV7Cjhe_DiokCGT_p8xbQIKHQ7lzEuZn",
        "subject": genomic_del4_seq_loc,
        "copies": 4,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_del4_cx_38(genomic_del4_seq_loc):
    """Create test fixture copy number change variation"""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.4gr9Ar1_evzZ6Q1f51oP7DkGM0LfxqGS",
        "subject": genomic_del4_seq_loc,
        "copyChange": "efo:0030067",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_del4_37_loc():
    """Create test fixture GRCh37 deletion subject"""
    return {
        "id": "ga4gh:SL.JJ04d-LYT8qrReNiZDdvj8J_uoBVwL0n",
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.v7noePfnNpK8ghYXEqZ9NukMXW7YeNsm",
        },
        "start": [None, 31138612],
        "end": [33357594, None],
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_del4_cn_37(genomic_del4_37_loc):
    """Create test fixture copy number count variation"""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.rvy2iu1XQ7zf2I1nKHOjGA_sRZGqNqk4",
        "subject": genomic_del4_37_loc,
        "copies": 4,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_del4_cx_37(genomic_del4_37_loc):
    """Create test fixture copy number change variation"""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.q_iLtX3wYfg8TeNI9JTR4pdoiJMzgtNN",
        "subject": genomic_del4_37_loc,
        "copyChange": "efo:0030067",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_del5_cn_38(genomic_del5_seq_loc):
    """Create test fixture copy number count variation"""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.ezuTVsAnNGPlwCRzcAZpcAMffIcKKZaA",
        "subject": genomic_del5_seq_loc,
        "copies": 2,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_del5_cx_38(genomic_del5_seq_loc):
    """Create test fixture copy number change variation"""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.UJpKxY2DFpV_CD_aI_df6s7YwRfuaBUc",
        "subject": genomic_del5_seq_loc,
        "copyChange": "efo:0030064",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_del5_37_loc():
    """Create test fixture GRCh37 deletion subject"""
    return {
        "id": "ga4gh:SL.Xhto4K6EdNUd-rYX1LYF-KXz2Lv2ZJfK",
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.v7noePfnNpK8ghYXEqZ9NukMXW7YeNsm",
        },
        "start": [None, 18593473],
        "end": 18671749,
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_del5_cn_37(genomic_del5_37_loc):
    """Create test fixture copy number count variation"""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.A9bYHmbYHHiEmWWu3cIZRAF7ZRL0Bykf",
        "subject": genomic_del5_37_loc,
        "copies": 2,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_del5_cx_37(genomic_del5_37_loc):
    """Create test fixture copy number change variation"""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.w2EqBSFkdKg1NREkHxdRYOO7UyZC9B9I",
        "subject": genomic_del5_37_loc,
        "copyChange": "efo:0030064",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_del6_cn_38(genomic_del6_seq_loc):
    """Create test fixture copy number count variation"""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.UoYMn-iOOJRz5OF1IdS1fJDhlHhAH-oR",
        "subject": genomic_del6_seq_loc,
        "copies": 1,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_del6_cx_38(genomic_del6_seq_loc):
    """Create test fixture copy number change variation"""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.5dfKF4nVUiKKHdrDPI0mBf2HzJGwNWhN",
        "subject": genomic_del6_seq_loc,
        "copyChange": "efo:0030071",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_del6_37_loc():
    """Create test fixture GRCh37 deletion subject"""
    return {
        "id": "ga4gh:SL.cKF6MGu-iUp74rlPikMYUyZ9A7ex1RyL",
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.KqaUhJMW3CDjhoVtBetdEKT1n6hM-7Ek",
        },
        "start": 133783901,
        "end": [133785996, None],
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_del6_cn_37(genomic_del6_37_loc):
    """Create test fixture copy number count variation"""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.SD_tu46e05joK4QvurSYa6f-UFRDMfKO",
        "subject": genomic_del6_37_loc,
        "copies": 1,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_del6_cx_37(genomic_del6_37_loc):
    """Create test fixture copy number change variation"""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.SIy1Piu11pRW5g75uFOEs0FmZNWz_4RZ",
        "subject": genomic_del6_37_loc,
        "copyChange": "efo:0030071",
    }
    return models.CopyNumberChange(**params)


@pytest.mark.asyncio
async def test_genomic_dup1_copy_number_count(
    test_cnv_handler, genomic_dup1_38_cn, genomic_dup1_cn_37
):
    """Test that genomic duplication works correctly"""
    q = "NC_000003.12:g.49531262dup"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q,
        baseline_copies=2,
        do_liftover=False,
    )
    cnv_assertion_checks(resp, genomic_dup1_38_cn)

    q = "NC_000003.11:g.49568695dup"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=2, do_liftover=False
    )
    cnv_assertion_checks(resp, genomic_dup1_cn_37)

    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=2, do_liftover=True
    )
    cnv_assertion_checks(resp, genomic_dup1_38_cn)

    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=1, do_liftover=True
    )
    expected = copy.deepcopy(genomic_dup1_38_cn)
    expected.copies = 2
    expected.id = "ga4gh:CN.pzlDkCIQevs-uQwEb8kcj0lEqctO1yQb"
    cnv_assertion_checks(resp, expected)


@pytest.mark.asyncio
async def test_genomic_dup1_copy_number_change(
    test_cnv_handler, genomic_dup1_cx_38, genomic_dup1_cx_37
):
    """Test that genomic duplication works correctly"""
    q = "NC_000003.12:g.49531262dup"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030069", do_liftover=False
    )
    cnv_assertion_checks(resp, genomic_dup1_cx_38)

    q = "NC_000003.11:g.49568695dup"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030069", do_liftover=False
    )
    cnv_assertion_checks(resp, genomic_dup1_cx_37)

    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030069", do_liftover=True
    )
    cnv_assertion_checks(resp, genomic_dup1_cx_38)


@pytest.mark.asyncio
async def test_genomic_dup2_copy_number_count(
    test_cnv_handler, genomic_dup2_38_cn, genomic_dup2_cn_37
):
    """Test that genomic duplication works correctly"""
    q = "NC_000023.11:g.33211290_33211293dup"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=2, do_liftover=False
    )
    cnv_assertion_checks(resp, genomic_dup2_38_cn)

    q = "NC_000023.10:g.33229407_33229410dup"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=2, do_liftover=False
    )
    cnv_assertion_checks(resp, genomic_dup2_cn_37)

    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=2, do_liftover=True
    )
    cnv_assertion_checks(resp, genomic_dup2_38_cn)

    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=1, do_liftover=True
    )
    expected = copy.deepcopy(genomic_dup2_38_cn)
    expected.copies = 2
    expected.id = "ga4gh:CN.6aH6YUxdUYU6Rbjrt-n74K8STpv4-dO3"
    cnv_assertion_checks(resp, expected)


@pytest.mark.asyncio
async def test_genomic_dup2_copy_number_change(
    test_cnv_handler, genomic_dup2_cx_38, genomic_dup2_cx_37
):
    """Test that genomic duplication works correctly"""
    q = "NC_000023.11:g.33211290_33211293dup"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030067", do_liftover=False
    )
    cnv_assertion_checks(resp, genomic_dup2_cx_38)

    q = "NC_000023.10:g.33229407_33229410dup"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030067", do_liftover=False
    )
    cnv_assertion_checks(resp, genomic_dup2_cx_37)

    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030067", do_liftover=True
    )
    cnv_assertion_checks(resp, genomic_dup2_cx_38)


@pytest.mark.asyncio
async def test_genomic_dup3_copy_number_count(
    test_cnv_handler, genomic_dup3_cn_38, genomic_dup3_cn_37
):
    """Test that genomic duplication works correctly"""
    q = "NC_000023.11:g.(31060227_31100351)_(33274278_33417151)dup"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=1, do_liftover=False
    )
    cnv_assertion_checks(resp, genomic_dup3_cn_38)

    q = "NC_000023.10:g.(31078344_31118468)_(33292395_33435268)dup"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=1, do_liftover=False
    )
    cnv_assertion_checks(resp, genomic_dup3_cn_37)

    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=1, do_liftover=True
    )
    cnv_assertion_checks(resp, genomic_dup3_cn_38)

    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=2, do_liftover=True
    )
    expected = copy.deepcopy(genomic_dup3_cn_38)
    expected.copies = 3
    expected.id = "ga4gh:CN.K-4wE5WJjZxv4vAegv8IeHenrBYfjJ46"
    cnv_assertion_checks(resp, expected)


@pytest.mark.asyncio
async def test_genomic_dup3_copy_number_change(
    test_cnv_handler, genomic_dup3_cx_38, genomic_dup3_cx_37
):
    """Test that genomic duplication works correctly"""
    q = "NC_000023.11:g.(31060227_31100351)_(33274278_33417151)dup"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030072", do_liftover=False
    )
    cnv_assertion_checks(resp, genomic_dup3_cx_38)

    q = "NC_000023.10:g.(31078344_31118468)_(33292395_33435268)dup"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030072", do_liftover=False
    )
    cnv_assertion_checks(resp, genomic_dup3_cx_37)

    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030072", do_liftover=True
    )
    cnv_assertion_checks(resp, genomic_dup3_cx_38)


@pytest.mark.asyncio
async def test_genomic_dup4_copy_number_count(
    test_cnv_handler, genomic_dup4_cn_38, genomic_dup4_cn_37
):
    """Test that genomic duplication works correctly"""
    q = "NC_000020.11:g.(?_30417576)_(31394018_?)dup"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=2, do_liftover=False
    )
    cnv_assertion_checks(resp, genomic_dup4_cn_38)

    q = "NC_000020.10:g.(?_29652252)_(29981821_?)dup"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=2, do_liftover=False
    )
    cnv_assertion_checks(resp, genomic_dup4_cn_37)

    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=2, do_liftover=True
    )
    cnv_assertion_checks(resp, genomic_dup4_cn_38)

    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=1, do_liftover=True
    )
    expected = copy.deepcopy(genomic_dup4_cn_38)
    expected.copies = 2
    expected.id = "ga4gh:CN.wmR2Ux1h5CqHsHPePxaUd4vVoEb6taDr"
    cnv_assertion_checks(resp, expected)


@pytest.mark.asyncio
async def test_genomic_dup4_copy_number_change(
    test_cnv_handler, genomic_dup4_cx_38, genomic_dup4_cx_37
):
    """Test that genomic duplication works correctly"""
    q = "NC_000020.11:g.(?_30417576)_(31394018_?)dup"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030069", do_liftover=False
    )
    cnv_assertion_checks(resp, genomic_dup4_cx_38)

    q = "NC_000020.10:g.(?_29652252)_(29981821_?)dup"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030069", do_liftover=False
    )
    cnv_assertion_checks(resp, genomic_dup4_cx_37)

    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030069", do_liftover=True
    )
    cnv_assertion_checks(resp, genomic_dup4_cx_38)


@pytest.mark.asyncio
async def test_genomic_dup5_copy_number_count(
    test_cnv_handler, genomic_dup5_cn_38, genomic_dup5_cn_37
):
    """Test that genomic duplication works correctly"""
    q = "NC_000023.11:g.(?_154021812)_154092209dup"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=3, do_liftover=False
    )
    cnv_assertion_checks(resp, genomic_dup5_cn_38)

    q = "NC_000023.10:g.(?_153287263)_153357667dup"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=3, do_liftover=False
    )
    cnv_assertion_checks(resp, genomic_dup5_cn_37)

    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=3, do_liftover=True
    )
    cnv_assertion_checks(resp, genomic_dup5_cn_38)

    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=4, do_liftover=True
    )
    expected = copy.deepcopy(genomic_dup5_cn_38)
    expected.copies = 5
    expected.id = "ga4gh:CN.5D9zEmtdVcOWlUCThRX30e5z811qHSQV"
    cnv_assertion_checks(resp, expected)


@pytest.mark.asyncio
async def test_genomic_dup5_copy_number_change(
    test_cnv_handler, genomic_dup5_cx_38, genomic_dup5_cx_37
):
    """Test that genomic duplication works correctly"""
    q = "NC_000023.11:g.(?_154021812)_154092209dup"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030067", do_liftover=False
    )
    cnv_assertion_checks(resp, genomic_dup5_cx_38)

    q = "NC_000023.10:g.(?_153287263)_153357667dup"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030067", do_liftover=False
    )
    cnv_assertion_checks(resp, genomic_dup5_cx_37)

    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030067", do_liftover=True
    )
    cnv_assertion_checks(resp, genomic_dup5_cx_38)


@pytest.mark.asyncio
async def test_genomic_dup6_copy_number_count(
    test_cnv_handler, genomic_dup6_cn_38, genomic_dup6_cn_37
):
    """Test that genomic duplication works correctly"""
    q = "NC_000023.11:g.154021812_(154092209_?)dup"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=1, do_liftover=False
    )
    cnv_assertion_checks(resp, genomic_dup6_cn_38)

    q = "NC_000023.10:g.153287263_(153357667_?)dup"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=1, do_liftover=False
    )
    cnv_assertion_checks(resp, genomic_dup6_cn_37)

    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=1, do_liftover=True
    )
    cnv_assertion_checks(resp, genomic_dup6_cn_38)

    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=2, do_liftover=True
    )
    expected = copy.deepcopy(genomic_dup6_cn_38)
    expected.copies = 3
    expected.id = "ga4gh:CN.R1Uqv_hy2X7ybuaLcqA0j-3-kPL9rK4a"
    cnv_assertion_checks(resp, expected)


@pytest.mark.asyncio
async def test_genomic_dup6_copy_number_change(
    test_cnv_handler, genomic_dup6_cx_38, genomic_dup6_cx_37
):
    """Test that genomic duplication works correctly"""
    q = "NC_000023.11:g.154021812_(154092209_?)dup"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030064", do_liftover=False
    )
    cnv_assertion_checks(resp, genomic_dup6_cx_38)

    q = "NC_000023.10:g.153287263_(153357667_?)dup"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030064", do_liftover=False
    )
    cnv_assertion_checks(resp, genomic_dup6_cx_37)

    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030064", do_liftover=True
    )
    cnv_assertion_checks(resp, genomic_dup6_cx_38)


@pytest.mark.asyncio
async def test_genomic_del1_copy_number_count(
    test_cnv_handler, genomic_del1_38_cn, genomic_del1_cn_37
):
    """Test that genomic deletion works correctly"""
    q = "NC_000003.12:g.10149811del"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=2, do_liftover=False
    )
    cnv_assertion_checks(resp, genomic_del1_38_cn)

    q = "NC_000003.11:g.10191495del"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=2, do_liftover=False
    )
    cnv_assertion_checks(resp, genomic_del1_cn_37)

    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=2, do_liftover=True
    )
    cnv_assertion_checks(resp, genomic_del1_38_cn)

    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=3, do_liftover=True
    )
    expected = copy.deepcopy(genomic_del1_38_cn)
    expected.copies = 2
    expected.id = "ga4gh:CN.XhOq5QVWHUmGEt3zZ77CGUIButVptkMI"
    cnv_assertion_checks(resp, expected)


@pytest.mark.asyncio
async def test_genomic_del1_copy_number_change(
    test_cnv_handler, genomic_del1_cx_38, genomic_del1_cx_37
):
    """Test that genomic deletion works correctly"""
    q = "NC_000003.12:g.10149811del"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030064", do_liftover=False
    )
    cnv_assertion_checks(resp, genomic_del1_cx_38)

    q = "NC_000003.11:g.10191495del"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030064", do_liftover=False
    )
    cnv_assertion_checks(resp, genomic_del1_cx_37)

    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030064", do_liftover=True
    )
    cnv_assertion_checks(resp, genomic_del1_cx_38)


@pytest.mark.asyncio
async def test_genomic_del2_copy_number_count(
    test_cnv_handler, genomic_del2_38_cn, genomic_del2_cn_37
):
    """Test that genomic deletion works correctly"""
    q = "NC_000003.12:g.10146595_10146613del"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=2, do_liftover=False
    )
    cnv_assertion_checks(resp, genomic_del2_38_cn)

    q = "NC_000003.11:g.10188279_10188297del"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=2, do_liftover=False
    )
    cnv_assertion_checks(resp, genomic_del2_cn_37)

    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=2, do_liftover=True
    )
    cnv_assertion_checks(resp, genomic_del2_38_cn)

    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=4, do_liftover=True
    )
    expected = copy.deepcopy(genomic_del2_38_cn)
    expected.copies = 3
    expected.id = "ga4gh:CN.8M5so67MgpZom6UdaIoDUB115poKMOap"
    cnv_assertion_checks(resp, expected)


@pytest.mark.asyncio
async def test_genomic_del2_copy_number_change(
    test_cnv_handler, genomic_del2_cx_38, genomic_del2_cx_37
):
    """Test that genomic deletion works correctly"""
    q = "NC_000003.12:g.10146595_10146613del"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030071", do_liftover=False
    )
    cnv_assertion_checks(resp, genomic_del2_cx_38)

    q = "NC_000003.11:g.10188279_10188297del"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030071", do_liftover=False
    )
    cnv_assertion_checks(resp, genomic_del2_cx_37)

    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030071", do_liftover=True
    )
    cnv_assertion_checks(resp, genomic_del2_cx_38)


@pytest.mark.asyncio
async def test_genomic_del3_copy_number_count(
    test_cnv_handler, genomic_del3_cn_38, genomic_del3_cn_37
):
    """Test that genomic deletion works correctly"""
    q = "NC_000023.11:g.(31060227_31100351)_(33274278_33417151)del"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=3, do_liftover=False
    )
    cnv_assertion_checks(resp, genomic_del3_cn_38)

    q = "NC_000023.10:g.(31078344_31118468)_(33292395_33435268)del"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=3, do_liftover=False
    )
    cnv_assertion_checks(resp, genomic_del3_cn_37)

    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=3, do_liftover=True
    )
    cnv_assertion_checks(resp, genomic_del3_cn_38)

    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=2, do_liftover=True
    )
    expected = copy.deepcopy(genomic_del3_cn_38)
    expected.copies = 1
    expected.id = "ga4gh:CN.mYLtO_1CdyfkKfby__CF1UYbX9KNb5Iv"
    cnv_assertion_checks(resp, expected)


@pytest.mark.asyncio
async def test_genomic_del3_copy_number_change(
    test_cnv_handler, genomic_del3_cx_38, genomic_del3_cx_37
):
    """Test that genomic deletion works correctly"""
    q = "NC_000023.11:g.(31060227_31100351)_(33274278_33417151)del"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030069", do_liftover=False
    )
    cnv_assertion_checks(resp, genomic_del3_cx_38)

    q = "NC_000023.10:g.(31078344_31118468)_(33292395_33435268)del"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030069", do_liftover=False
    )
    cnv_assertion_checks(resp, genomic_del3_cx_37)

    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030069", do_liftover=True
    )
    cnv_assertion_checks(resp, genomic_del3_cx_38)


@pytest.mark.asyncio
async def test_genomic_del4_copy_number_count(
    test_cnv_handler, genomic_del4_cn_38, genomic_del4_cn_37
):
    """Test that genomic deletion works correctly"""
    q = "NC_000023.11:g.(?_31120496)_(33339477_?)del"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=5, do_liftover=False
    )
    cnv_assertion_checks(resp, genomic_del4_cn_38)

    q = "NC_000023.10:g.(?_31138613)_(33357594_?)del"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=5, do_liftover=False
    )
    cnv_assertion_checks(resp, genomic_del4_cn_37)

    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=5, do_liftover=True
    )
    cnv_assertion_checks(resp, genomic_del4_cn_38)

    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=3, do_liftover=True
    )
    expected = copy.deepcopy(genomic_del4_cn_38)
    expected.copies = 2
    expected.id = "ga4gh:CN.NvFf_XtkTqqxU3LpMoLbo41UAFAINHH8"
    cnv_assertion_checks(resp, expected)


@pytest.mark.asyncio
async def test_genomic_del4_copy_number_change(
    test_cnv_handler, genomic_del4_cx_38, genomic_del4_cx_37
):
    """Test that genomic deletion works correctly"""
    q = "NC_000023.11:g.(?_31120496)_(33339477_?)del"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030067", do_liftover=False
    )
    cnv_assertion_checks(resp, genomic_del4_cx_38)

    q = "NC_000023.10:g.(?_31138613)_(33357594_?)del"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030067", do_liftover=False
    )
    cnv_assertion_checks(resp, genomic_del4_cx_37)

    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030067", do_liftover=True
    )
    cnv_assertion_checks(resp, genomic_del4_cx_38)


@pytest.mark.asyncio
async def test_genomic_del5_copy_number_count(
    test_cnv_handler, genomic_del5_cn_38, genomic_del5_cn_37
):
    """Test that genomic deletion works correctly"""
    q = "NC_000023.11:g.(?_18575354)_18653629del"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=3, do_liftover=False
    )
    cnv_assertion_checks(resp, genomic_del5_cn_38)

    q = "NC_000023.10:g.(?_18593474)_18671749del"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=3, do_liftover=False
    )
    cnv_assertion_checks(resp, genomic_del5_cn_37)

    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=3, do_liftover=True
    )
    cnv_assertion_checks(resp, genomic_del5_cn_38)

    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=2, do_liftover=True
    )
    expected = copy.deepcopy(genomic_del5_cn_38)
    expected.copies = 1
    expected.id = "ga4gh:CN.704pxEQyFZ466GOl8IxEKAJFtP5xLb8J"
    cnv_assertion_checks(resp, expected)


@pytest.mark.asyncio
async def test_genomic_del5_copy_number_change(
    test_cnv_handler, genomic_del5_cx_38, genomic_del5_cx_37
):
    """Test that genomic deletion works correctly"""
    q = "NC_000023.11:g.(?_18575354)_18653629del"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030064", do_liftover=False
    )
    cnv_assertion_checks(resp, genomic_del5_cx_38)

    q = "NC_000023.10:g.(?_18593474)_18671749del"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030064", do_liftover=False
    )
    cnv_assertion_checks(resp, genomic_del5_cx_37)

    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030064", do_liftover=True
    )
    cnv_assertion_checks(resp, genomic_del5_cx_38)


@pytest.mark.asyncio
async def test_genomic_del6_copy_number_count(
    test_cnv_handler, genomic_del6_cn_38, genomic_del6_cn_37
):
    """Test that genomic deletion works correctly"""
    q = "NC_000006.12:g.133462764_(133464858_?)del"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=2, do_liftover=False
    )
    cnv_assertion_checks(resp, genomic_del6_cn_38)

    q = "NC_000006.11:g.133783902_(133785996_?)del"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=2, do_liftover=False
    )
    cnv_assertion_checks(resp, genomic_del6_cn_37)

    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=2, do_liftover=True
    )
    cnv_assertion_checks(resp, genomic_del6_cn_38)

    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=3, do_liftover=True
    )
    expected = copy.deepcopy(genomic_del6_cn_38)
    expected.copies = 2
    expected.id = "ga4gh:CN._ZaD3a4RsGWvUcw3i6e75-TIO2nbyHYv"
    cnv_assertion_checks(resp, expected)


@pytest.mark.asyncio
async def test_genomic_del6_copy_number_change(
    test_cnv_handler, genomic_del6_cx_38, genomic_del6_cx_37
):
    """Test that genomic deletion works correctly"""
    q = "NC_000006.12:g.133462764_(133464858_?)del"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030071", do_liftover=False
    )
    cnv_assertion_checks(resp, genomic_del6_cx_38)

    q = "NC_000006.11:g.133783902_(133785996_?)del"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030071", do_liftover=False
    )
    cnv_assertion_checks(resp, genomic_del6_cx_37)

    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030071", do_liftover=True
    )
    cnv_assertion_checks(resp, genomic_del6_cx_38)


@pytest.mark.asyncio
async def test_invalid_cnv(test_cnv_handler):
    """Check that invalid input return warnings"""
    q = "DAG1 g.49568695dup"
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q,
        copy_change="efo:0030071",
        do_liftover=True,
    )
    assert set(resp.warnings) == {
        "DAG1 g.49568695dup is not a supported HGVS genomic duplication or deletion"
    }
    assert resp.copy_number_change is None

    q = "braf V600E"
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030071", do_liftover=True
    )
    assert set(resp.warnings) == {
        "braf V600E is not a supported HGVS genomic duplication or deletion"
    }
    assert resp.copy_number_change is None

    # Not yet supported
    for q in [
        "NC_000018.9:g.(48556994_48573289)_48573471dup",
        "NC_000018.9:g.48556994_(48573289_48573471)dup",
    ]:
        resp = await test_cnv_handler.hgvs_to_copy_number_change(
            q, copy_change="efo:0030070"
        )
        assert resp.warnings == [f"Unable to find classification for: {q}"], q
        assert resp.copy_number_change is None, q
