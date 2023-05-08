"""Module for testing the hgvs to copy number count and copy number change endpoints"""
import copy

import pytest


@pytest.fixture(scope="module")
def genomic_dup1_cx_38(genomic_dup1_seq_loc):
    """Create test fixture copy number change variation"""
    return {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.ENQD2_J-4FE964fFrO9cjBDBee09ORWH",
        "subject": genomic_dup1_seq_loc,
        "copy_change": "efo:0030069"
    }


@pytest.fixture(scope="module")
def genomic_dup1_37_loc():
    """Create test fixture GRCh37 duplication location"""
    return {
        "id": "ga4gh:SL.PoXklQujGEI2KSuYEhM_boWkdZNe6XVp",
        "sequence_id": "ga4gh:SQ.VNBualIltAyi2AI_uXcKU7M9XUOuA7MS",
        "start": {"value": 49568693, "type": "Number"},
        "end": {"value": 49568695, "type": "Number"},
        "type": "SequenceLocation"
    }


@pytest.fixture(scope="module")
def genomic_dup1_cn_37(genomic_dup1_37_loc):
    """Create test fixture copy number count variation"""
    return {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.FmZdc1ZTOJoawPld_zSYtWDmJH_N0-5N",
        "subject": genomic_dup1_37_loc,
        "copies": {"type": "Number", "value": 3}
    }


@pytest.fixture(scope="module")
def genomic_dup1_cx_37(genomic_dup1_37_loc):
    """Create test fixture copy number change variation"""
    return {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.hrzxucrCYF5kLU3_1WwmnN5F6erRDnQe",
        "subject": genomic_dup1_37_loc,
        "copy_change": "efo:0030069"
    }


@pytest.fixture(scope="module")
def genomic_dup2_cx_38(genomic_dup2_seq_loc):
    """Create test fixture copy number change variation"""
    return {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.dCxABnpPBpZs05zbf4tBA6TmYaGReCF5",
        "subject": genomic_dup2_seq_loc,
        "copy_change": "efo:0030068"
    }


@pytest.fixture(scope="module")
def genomic_dup2_37_loc():
    """Create test fixture GRCh37 duplication location"""
    return {
        "id": "ga4gh:SL.rfaH4tbIPn7o3ZnmMqjxtZPq_57vxjc0",
        "sequence_id": "ga4gh:SQ.W6wLoIFOn4G7cjopxPxYNk2lcEqhLQFb",
        "start": {"value": 2137938, "type": "Number"},
        "end": {"value": 2137949, "type": "Number"},
        "type": "SequenceLocation"
    }


@pytest.fixture(scope="module")
def genomic_dup2_cn_37(genomic_dup2_37_loc):
    """Create test fixture copy number count variation"""
    return {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.H1o0rrnS2XekKUvVyMWIHyenk4vvgFLv",
        "subject": genomic_dup2_37_loc,
        "copies": {"type": "Number", "value": 3}
    }


@pytest.fixture(scope="module")
def genomic_dup2_cx_37(genomic_dup2_37_loc):
    """Create test fixture copy number change variation"""
    return {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.10_4sZHdFwDA9ekxHcPMpmNbJRgvUIEl",
        "subject": genomic_dup2_37_loc,
        "copy_change": "efo:0030068"
    }


@pytest.fixture(scope="module")
def genomic_dup3_cn_38(genomic_del3_dup3_loc):
    """Create test fixture copy number count variation"""
    return {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.a2-YlA01Q-y2zAOgn4r-cP9hiIEGaaqq",
        "subject": genomic_del3_dup3_loc,
        "copies": {"type": "Number", "value": 2}
    }


@pytest.fixture(scope="module")
def genomic_dup3_cx_38(genomic_del3_dup3_loc):
    """Create test fixture copy number change variation"""
    return {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.MGhVuAEVskM_cAhATCYT_tRcOsqDbHNd",
        "subject": genomic_del3_dup3_loc,
        "copy_change": "efo:0030072"
    }


@pytest.fixture(scope="module")
def genomic_dup3_37_loc():
    """Create test fixture GRCh37 duplication location"""
    return {
        "id": "ga4gh:SL.MmBA0qlC0J6yxWTyqhpsd4H_oErccpUD",
        "sequence_id": "ga4gh:SQ.v7noePfnNpK8ghYXEqZ9NukMXW7YeNsm",
        "start": {"min": 31078343, "max": 31118467, "type": "DefiniteRange"},
        "end": {"min": 33292396, "max": 33435269, "type": "DefiniteRange"},
        "type": "SequenceLocation"
    }


@pytest.fixture(scope="module")
def genomic_dup3_cn_37(genomic_dup3_37_loc):
    """Create test fixture copy number count variation"""
    return {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.GOhhUsvbhFqaXgq3hNH3nuJldMIwVMpO",
        "subject": genomic_dup3_37_loc,
        "copies": {"type": "Number", "value": 2}
    }


@pytest.fixture(scope="module")
def genomic_dup3_cx_37(genomic_dup3_37_loc):
    """Create test fixture copy number change variation"""
    return {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.fpx5MzUfd0Lll7kOwHxQ_znt5k-KJBau",
        "subject": genomic_dup3_37_loc,
        "copy_change": "efo:0030072"
    }


@pytest.fixture(scope="module")
def genomic_dup4_cn_38(genomic_dup4_loc):
    """Create test fixture copy number count variation"""
    return {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.Fvt2VSCSLNROTjZd20vKZUzdMYxz-cfZ",
        "subject": genomic_dup4_loc,
        "copies": {"type": "Number", "value": 3}
    }


@pytest.fixture(scope="module")
def genomic_dup4_cx_38(genomic_dup4_loc):
    """Create test fixture copy number change variation"""
    return {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.WlcvrOPHx67R2Qyw8O_s3MyHTJLX8uJD",
        "subject": genomic_dup4_loc,
        "copy_change": "efo:0030069"
    }


@pytest.fixture(scope="module")
def genomic_dup4_37_loc():
    """Create test fixture GRCh37 duplication location"""
    return {
        "id": "ga4gh:SL.RUcHRUFPNBJhkwXz-WzUZLbT7XOlddAv",
        "sequence_id": "ga4gh:SQ.iy_UbUrvECxFRX5LPTH_KPojdlT7BKsf",
        "start": {"value": 29652251, "comparator": "<=", "type": "IndefiniteRange"},
        "end": {"value": 29981821, "comparator": ">=", "type": "IndefiniteRange"},
        "type": "SequenceLocation"
    }


@pytest.fixture(scope="module")
def genomic_dup4_cn_37(genomic_dup4_37_loc):
    """Create test fixture copy number change variation"""
    return {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN._KLpF-U9nyddrZvZott7sQpJmdMRAh9E",
        "subject": genomic_dup4_37_loc,
        "copies": {"type": "Number", "value": 3}
    }


@pytest.fixture(scope="module")
def genomic_dup4_cx_37(genomic_dup4_37_loc):
    """Create test fixture copy number change variation"""
    return {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.q0kxeY8Nm8I9hLkrlOldryNy8Qv7prje",
        "subject": genomic_dup4_37_loc,
        "copy_change": "efo:0030069"
    }


@pytest.fixture(scope="module")
def genomic_dup5_cn_38(genomic_dup5_loc):
    """Create test fixture copy number count variation"""
    return {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.cegGbGJ5faMNTnlr6Oc4r93BOgklfftc",
        "subject": genomic_dup5_loc,
        "copies": {"type": "Number", "value": 4}
    }


@pytest.fixture(scope="module")
def genomic_dup5_cx_38(genomic_dup5_loc):
    """Create test fixture copy number change variation"""
    return {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.VWyaC_x9kGtwvzkNYL17Bf0lpf1RdHyg",
        "subject": genomic_dup5_loc,
        "copy_change": "efo:0030068"
    }


@pytest.fixture(scope="module")
def genomic_dup5_37_loc():
    """Create test fixture GRCh37 duplication location"""
    return {
        "id": "ga4gh:SL.1ELvCCFb5kq7_k47XzpatrJMgtvW74Co",
        "sequence_id": "ga4gh:SQ.v7noePfnNpK8ghYXEqZ9NukMXW7YeNsm",
        "start": {"value": 153287262, "comparator": "<=", "type": "IndefiniteRange"},
        "end": {"value": 153357667, "type": "Number"},
        "type": "SequenceLocation"
    }


@pytest.fixture(scope="module")
def genomic_dup5_cn_37(genomic_dup5_37_loc):
    """Create test fixture copy number count variation"""
    return {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.wWg3b2U9qwdvS7ASQaG-6Gh6rhBIS5AA",
        "subject": genomic_dup5_37_loc,
        "copies": {"type": "Number", "value": 4}
    }


@pytest.fixture(scope="module")
def genomic_dup5_cx_37(genomic_dup5_37_loc):
    """Create test fixture copy number change variation"""
    return {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.ggLmm0rU4lgCL9_HjkYCOMu-mCsiWFBy",
        "subject": genomic_dup5_37_loc,
        "copy_change": "efo:0030068"
    }


@pytest.fixture(scope="module")
def genomic_dup6_cn_38(genomic_dup6_loc):
    """Create test fixture copy number count variation"""
    return {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.Tuvj2C1xYVUF0UL3CSB1GZjCRC6K_eSt",
        "subject": genomic_dup6_loc,
        "copies": {"type": "Number", "value": 2}
    }


@pytest.fixture(scope="module")
def genomic_dup6_cx_38(genomic_dup6_loc):
    """Create test fixture copy number change variation"""
    return {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.RKJvR7nY_nCn0rb70DzTDPHAdcwq1gaG",
        "subject": genomic_dup6_loc,
        "copy_change": "efo:0030067"
    }


@pytest.fixture(scope="module")
def genomic_dup6_37_loc():
    """Create test fixture GRCh37 duplication location"""
    return {
        "id": "ga4gh:SL.GdeHSPT6RDTaMycAXoD_rqhMjumKGd-4",
        "sequence_id": "ga4gh:SQ.v7noePfnNpK8ghYXEqZ9NukMXW7YeNsm",
        "start": {"value": 153287262, "type": "Number"},
        "end": {"value": 153357667, "comparator": ">=", "type": "IndefiniteRange"},
        "type": "SequenceLocation"
    }


@pytest.fixture(scope="module")
def genomic_dup6_cn_37(genomic_dup6_37_loc):
    """Create test fixture copy number count variation"""
    return {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.JvO8QAavKtd38hFB-wpzsm_9JCrFCM7i",
        "subject": genomic_dup6_37_loc,
        "copies": {"type": "Number", "value": 2}
    }


@pytest.fixture(scope="module")
def genomic_dup6_cx_37(genomic_dup6_37_loc):
    """Create test fixture copy number change variation"""
    return {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.d27oYKn13QK6oLQGcO_dM-KyFYgVD6Os",
        "subject": genomic_dup6_37_loc,
        "copy_change": "efo:0030067"
    }


@pytest.fixture(scope="module")
def genomic_del1_cx_38(genomic_del1_seq_loc):
    """Create test fixture copy number change variation"""
    return {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.LIMo4t9dxh_lxtppith_0xULcbIoGYi0",
        "subject": genomic_del1_seq_loc,
        "copy_change": "efo:0030067"
    }


@pytest.fixture(scope="module")
def genomic_del1_37_loc():
    """Create test fixture GRCh37 deletion location"""
    return {
        "id": "ga4gh:SL.QtY9ACIT6TJAae_KDSAZKJfBLerBbFfq",
        "sequence_id": "ga4gh:SQ.VNBualIltAyi2AI_uXcKU7M9XUOuA7MS",
        "start": {"value": 10191494, "type": "Number"},
        "end": {"value": 10191495, "type": "Number"},
        "type": "SequenceLocation"
    }


@pytest.fixture(scope="module")
def genomic_del1_cn_37(genomic_del1_37_loc):
    """Create test fixture copy number count variation"""
    return {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.rkgfXxRiPAPhf8b8M_7JZ_Sd_v--xTSp",
        "subject": genomic_del1_37_loc,
        "copies": {"type": "Number", "value": 1}
    }


@pytest.fixture(scope="module")
def genomic_del1_cx_37(genomic_del1_37_loc):
    """Create test fixture copy number change variation"""
    return {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.TUVvFFsCWkRQW2wmkEFE6WguYBY8OVyb",
        "subject": genomic_del1_37_loc,
        "copy_change": "efo:0030067"
    }


@pytest.fixture(scope="module")
def genomic_del2_cx_38(genomic_del2_seq_loc):
    """Create test fixture copy number change variation"""
    return {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.Y_94z-QmAZgbQfSRRFxBbaBe1GZ4TBCG",
        "subject": genomic_del2_seq_loc,
        "copy_change": "efo:0030071"
    }


@pytest.fixture(scope="module")
def genomic_del2_37_loc():
    """Create test fixture GRCh37 deletion location"""
    return {
        "id": "ga4gh:SL.XXW-NlnipMp1hTv0xGdTY-RsjDUUZWht",
        "sequence_id": "ga4gh:SQ.VNBualIltAyi2AI_uXcKU7M9XUOuA7MS",
        "start": {"value": 10188278, "type": "Number"},
        "end": {"value": 10188297, "type": "Number"},
        "type": "SequenceLocation"
    }


@pytest.fixture(scope="module")
def genomic_del2_cn_37(genomic_del2_37_loc):
    """Create test fixture copy number count variation"""
    return {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.pn37YpWXjDcIPLhuTUtAb9-pwFKrnM4M",
        "subject": genomic_del2_37_loc,
        "copies": {"type": "Number", "value": 1}
    }


@pytest.fixture(scope="module")
def genomic_del2_cx_37(genomic_del2_37_loc):
    """Create test fixture copy number change variation"""
    return {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.pWZIXKKaRGj8LwUJBcdEPnJKtMzLAUJ_",
        "subject": genomic_del2_37_loc,
        "copy_change": "efo:0030071"
    }


@pytest.fixture(scope="module")
def genomic_del3_cn_38(genomic_del3_dup3_loc):
    """Create test fixture copy number count variation"""
    return {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.a2-YlA01Q-y2zAOgn4r-cP9hiIEGaaqq",
        "subject": genomic_del3_dup3_loc,
        "copies": {"type": "Number", "value": 2}
    }


@pytest.fixture(scope="module")
def genomic_del3_cx_38(genomic_del3_dup3_loc):
    """Create test fixture copy number change variation"""
    return {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.Xffhy7i-ZkZtQSoMz6eytI6QQlShmfPo",
        "subject": genomic_del3_dup3_loc,
        "copy_change": "efo:0030069"
    }


@pytest.fixture(scope="module")
def genomic_del3_37_loc():
    """Create test fixture GRCh37 deletion location"""
    return {
        "id": "ga4gh:SL.MmBA0qlC0J6yxWTyqhpsd4H_oErccpUD",
        "sequence_id": "ga4gh:SQ.v7noePfnNpK8ghYXEqZ9NukMXW7YeNsm",
        "start": {"min": 31078343, "max": 31118467, "type": "DefiniteRange"},
        "end": {"min": 33292396, "max": 33435269, "type": "DefiniteRange"},
        "type": "SequenceLocation"
    }


@pytest.fixture(scope="module")
def genomic_del3_cn_37(genomic_del3_37_loc):
    """Create test fixture copy number count variation"""
    return {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.GOhhUsvbhFqaXgq3hNH3nuJldMIwVMpO",
        "subject": genomic_del3_37_loc,
        "copies": {"type": "Number", "value": 2}
    }


@pytest.fixture(scope="module")
def genomic_del3_cx_37(genomic_del3_37_loc):
    """Create test fixture copy number change variation"""
    return {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.5xJqmmaJTg0gv_BVa5gFHW_krgAGXviI",
        "subject": genomic_del3_37_loc,
        "copy_change": "efo:0030069"
    }


@pytest.fixture(scope="module")
def genomic_del4_cn_38(genomic_del4_seq_loc):
    """Create test fixture copy number count variation"""
    return {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.R0U1Bg1lhPYRVyJx_4LVE3JmR_5V7e8e",
        "subject": genomic_del4_seq_loc,
        "copies": {"type": "Number", "value": 4}
    }


@pytest.fixture(scope="module")
def genomic_del4_cx_38(genomic_del4_seq_loc):
    """Create test fixture copy number change variation"""
    return {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.bHHDJw8nQXVLnR8HM8t2jsCxY1exj-69",
        "subject": genomic_del4_seq_loc,
        "copy_change": "efo:0030068"
    }


@pytest.fixture(scope="module")
def genomic_del4_37_loc():
    """Create test fixture GRCh37 deletion location"""
    return {
        "id": "ga4gh:SL.ghgq_DzgZ93L1gVXbnHKEQwBzwnNdfSJ",
        "sequence_id": "ga4gh:SQ.v7noePfnNpK8ghYXEqZ9NukMXW7YeNsm",
        "start": {"value": 31138612, "comparator": "<=", "type": "IndefiniteRange"},
        "end": {"value": 33357594, "comparator": ">=", "type": "IndefiniteRange"},
        "type": "SequenceLocation"
    }


@pytest.fixture(scope="module")
def genomic_del4_cn_37(genomic_del4_37_loc):
    """Create test fixture copy number count variation"""
    return {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.lI6I5Tf7rRjN4P_T1euZgAn0YucxrMj4",
        "subject": genomic_del4_37_loc,
        "copies": {"type": "Number", "value": 4}
    }


@pytest.fixture(scope="module")
def genomic_del4_cx_37(genomic_del4_37_loc):
    """Create test fixture copy number change variation"""
    return {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.9JHPpA1EUtPvK7ta6S5oSL7IIfu_Aun2",
        "subject": genomic_del4_37_loc,
        "copy_change": "efo:0030068"
    }


@pytest.fixture(scope="module")
def genomic_del5_cn_38(genomic_del5_seq_loc):
    """Create test fixture copy number count variation"""
    return {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.UhXP8ENfA3T-FhNvKSqu-9UwrL5GZFvF",
        "subject": genomic_del5_seq_loc,
        "copies": {"type": "Number", "value": 2}
    }


@pytest.fixture(scope="module")
def genomic_del5_cx_38(genomic_del5_seq_loc):
    """Create test fixture copy number change variation"""
    return {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.0TqknNhM5e5i-5tAKxGXv4yDJlaQgDRR",
        "subject": genomic_del5_seq_loc,
        "copy_change": "efo:0030067"
    }


@pytest.fixture(scope="module")
def genomic_del5_37_loc():
    """Create test fixture GRCh37 deletion location"""
    return {
        "id": "ga4gh:SL.kJXgk4nQIwS53UN5gxgIFIPcNSt53n8O",
        "sequence_id": "ga4gh:SQ.v7noePfnNpK8ghYXEqZ9NukMXW7YeNsm",
        "start": {"value": 18593473, "comparator": "<=", "type": "IndefiniteRange"},
        "end": {"value": 18671749, "type": "Number"},
        "type": "SequenceLocation"
    }


@pytest.fixture(scope="module")
def genomic_del5_cn_37(genomic_del5_37_loc):
    """Create test fixture copy number count variation"""
    return {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.3Yu1OOpMjnzt6yrt0OLjKWnhD2YD9pkg",
        "subject": genomic_del5_37_loc,
        "copies": {"type": "Number", "value": 2}
    }


@pytest.fixture(scope="module")
def genomic_del5_cx_37(genomic_del5_37_loc):
    """Create test fixture copy number change variation"""
    return {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.njnotyulNBRZMFuJ-rjACLwZS1MH6Ft2",
        "subject": genomic_del5_37_loc,
        "copy_change": "efo:0030067"
    }


@pytest.fixture(scope="module")
def genomic_del6_cn_38(genomic_del6_seq_loc):
    """Create test fixture copy number count variation"""
    return {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.lfI0nj_cOjRwul_R_4tohZK_vqiZ-LSz",
        "subject": genomic_del6_seq_loc,
        "copies": {"type": "Number", "value": 1}
    }


@pytest.fixture(scope="module")
def genomic_del6_cx_38(genomic_del6_seq_loc):
    """Create test fixture copy number change variation"""
    return {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.HbonS2GEi_jBQA_2gvaiwONy_KHqN_G0",
        "subject": genomic_del6_seq_loc,
        "copy_change": "efo:0030071"
    }


@pytest.fixture(scope="module")
def genomic_del6_37_loc():
    """Create test fixture GRCh37 deletion location"""
    return {
        "id": "ga4gh:SL.LU0WrrgeZ-HIG7f0LF_eS0yGvhDtiHuT",
        "sequence_id": "ga4gh:SQ.KqaUhJMW3CDjhoVtBetdEKT1n6hM-7Ek",
        "start": {"value": 133783901, "type": "Number"},
        "end": {"value": 133785996, "comparator": ">=", "type": "IndefiniteRange"},
        "type": "SequenceLocation"
    }


@pytest.fixture(scope="module")
def genomic_del6_cn_37(genomic_del6_37_loc):
    """Create test fixture copy number count variation"""
    return {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.UclTlSHjepPFKjkBdZjTPquoStwlr24H",
        "subject": genomic_del6_37_loc,
        "copies": {"type": "Number", "value": 1}
    }


@pytest.fixture(scope="module")
def genomic_del6_cx_37(genomic_del6_37_loc):
    """Create test fixture copy number change variation"""
    return {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.j6BFcKCBMSfz3Vh-ekevUx8GqLCGJCU6",
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
    expected["id"] = "ga4gh:CN.XbKkhTdqvXkyIT8rJwotqYW0v2Ka1zOU"
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
    q = "NC_000016.10:g.2087938_2087948dup"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=2, do_liftover=False)
    assert resp.copy_number_count.dict(by_alias=True) == genomic_dup2_38_cn

    q = "NC_000016.9:g.2137939_2137949dup"  # 37
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
    expected["id"] = "ga4gh:CN.Nt2omwS2KB2JfskoeX4v3eBBL2l7T-s-"
    assert resp.copy_number_count.dict(by_alias=True) == expected


@pytest.mark.asyncio
async def test_genomic_dup2_copy_number_change(test_cnv_handler, genomic_dup2_cx_38,
                                               genomic_dup2_cx_37):
    """Test that genomic duplication works correctly"""
    q = "NC_000016.10:g.2087938_2087948dup"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030068", do_liftover=False)
    assert resp.copy_number_change.dict(by_alias=True) == genomic_dup2_cx_38

    q = "NC_000016.9:g.2137939_2137949dup"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030068", do_liftover=False)
    assert resp.copy_number_change.dict(by_alias=True) == genomic_dup2_cx_37

    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030068", do_liftover=True)
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
    expected["id"] = "ga4gh:CN.-4AcyJ5imRkXzgeULrzWVTraqAe4YRYc"
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
    expected["id"] = "ga4gh:CN.ObqoAUmFaGiMknzr8GurBYuoBNeE2HtV"
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
    expected["id"] = "ga4gh:CN.KxU8n8a8UXCQnnwHDX9daUXd0zdxtFtB"
    assert resp.copy_number_count.dict(by_alias=True) == expected


@pytest.mark.asyncio
async def test_genomic_dup5_copy_number_change(test_cnv_handler, genomic_dup5_cx_38,
                                               genomic_dup5_cx_37):
    """Test that genomic duplication works correctly"""
    q = "NC_000023.11:g.(?_154021812)_154092209dup"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030068", do_liftover=False)
    assert resp.copy_number_change.dict(by_alias=True) == genomic_dup5_cx_38

    q = "NC_000023.10:g.(?_153287263)_153357667dup"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030068", do_liftover=False)
    assert resp.copy_number_change.dict(by_alias=True) == genomic_dup5_cx_37

    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030068", do_liftover=True)
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
    expected["id"] = "ga4gh:CN.TQ36aHJOnfKG3a8hpVMR1POBTOVBgYcS"
    assert resp.copy_number_count.dict(by_alias=True) == expected


@pytest.mark.asyncio
async def test_genomic_dup6_copy_number_change(test_cnv_handler, genomic_dup6_cx_38,
                                               genomic_dup6_cx_37):
    """Test that genomic duplication works correctly"""
    q = "NC_000023.11:g.154021812_(154092209_?)dup"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030067", do_liftover=False)
    assert resp.copy_number_change.dict(by_alias=True) == genomic_dup6_cx_38

    q = "NC_000023.10:g.153287263_(153357667_?)dup"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030067", do_liftover=False)
    assert resp.copy_number_change.dict(by_alias=True) == genomic_dup6_cx_37

    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030067", do_liftover=True)
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
    expected["id"] = "ga4gh:CN.WwgB7O1Ztcn6QqN2hjfjHucKOaWXic8n"
    assert resp.copy_number_count.dict(by_alias=True) == expected


@pytest.mark.asyncio
async def test_genomic_del1_copy_number_change(test_cnv_handler, genomic_del1_cx_38,
                                               genomic_del1_cx_37):
    """Test that genomic deletion works correctly"""
    q = "NC_000003.12:g.10149811del"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030067", do_liftover=False)
    assert resp.copy_number_change.dict(by_alias=True) == genomic_del1_cx_38

    q = "NC_000003.11:g.10191495del"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030067", do_liftover=False)
    assert resp.copy_number_change.dict(by_alias=True) == genomic_del1_cx_37

    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030067", do_liftover=True)
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
    expected["id"] = "ga4gh:CN.GiZEIUkpM1Yl_O9v0Xz3X9Uza1BzdYwn"
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
    expected["id"] = "ga4gh:CN.a7prdJg9CK7KAFug-0VMjZC7xDxiSByV"
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
    expected["id"] = "ga4gh:CN.PlqY0WzSBD2KV2mD0L9QkyOGf1z2zZLu"
    assert resp.copy_number_count.dict(by_alias=True) == expected


@pytest.mark.asyncio
async def test_genomic_del4_copy_number_change(test_cnv_handler, genomic_del4_cx_38,
                                               genomic_del4_cx_37):
    """Test that genomic deletion works correctly"""
    q = "NC_000023.11:g.(?_31120496)_(33339477_?)del"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030068", do_liftover=False)
    assert resp.copy_number_change.dict(by_alias=True) == genomic_del4_cx_38

    q = "NC_000023.10:g.(?_31138613)_(33357594_?)del"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030068", do_liftover=False)
    assert resp.copy_number_change.dict(by_alias=True) == genomic_del4_cx_37

    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030068", do_liftover=True)
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
    expected["id"] = "ga4gh:CN.HzpGrUmI6mfMtPGC0H8ElVcNc0-V1ldA"
    assert resp.copy_number_count.dict(by_alias=True) == expected


@pytest.mark.asyncio
async def test_genomic_del5_copy_number_change(test_cnv_handler, genomic_del5_cx_38,
                                               genomic_del5_cx_37):
    """Test that genomic deletion works correctly"""
    q = "NC_000023.11:g.(?_18575354)_18653629del"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030067", do_liftover=False)
    assert resp.copy_number_change.dict(by_alias=True) == genomic_del5_cx_38

    q = "NC_000023.10:g.(?_18593474)_18671749del"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030067", do_liftover=False)
    assert resp.copy_number_change.dict(by_alias=True) == genomic_del5_cx_37

    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030067", do_liftover=True)
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
    expected["id"] = "ga4gh:CN.NWXsPHBTkOaY0LgfZK6bWFoJ3lvf3fOH"
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
        q, copy_change="low-level gains", do_liftover=True)
    assert resp is None
    assert w == ["low-level gains is not a valid copy change: ['efo:0030069', "
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

    q = "braf v600e"
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030071", do_liftover=True)
    assert set(resp.warnings) == {"Unable to translate braf v600e to copy number variation",  # noqa: E501
                                  "braf v600e is not a supported HGVS genomic duplication or deletion"}  # noqa: E501
    assert resp.copy_number_change is None

    # Not yet supported
    for q in ["NC_000018.9:g.(48556994_48573289)_48573471dup",
              "NC_000018.9:g.48556994_(48573289_48573471)dup"]:
        resp = await test_cnv_handler.hgvs_to_copy_number_change(
            q, copy_change="efo:0030070"
        )
        assert resp.warnings == ["Unable to find valid result for classifications: "
                                 "{'genomic duplication'}"], q
        assert resp.copy_number_change is None, q
