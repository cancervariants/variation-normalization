"""Module for testing the hgvs to absolute and relative copy number endpoints"""
import copy

import pytest


@pytest.fixture(scope="module")
def genomic_dup1_rel_38(genomic_dup1_seq_loc):
    """Create test fixture relative copy number variation"""
    return {
        "type": "RelativeCopyNumber",
        "id": "ga4gh:RCN._1Nz0yj2g9Q6cRO8j6oRi5peUJYjTAga",
        "location": genomic_dup1_seq_loc,
        "relative_copy_class": "complete loss"
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
def genomic_dup1_abs_37(genomic_dup1_37_loc):
    """Create test fixture absolute copy number variation"""
    return {
        "type": "AbsoluteCopyNumber",
        "id": "ga4gh:ACN.JWc1vKWUXuxG6l5CdftVGD6axJ91geGj",
        "location": genomic_dup1_37_loc,
        "copies": {"type": "Number", "value": 3}
    }


@pytest.fixture(scope="module")
def genomic_dup1_rel_37(genomic_dup1_37_loc):
    """Create test fixture relative copy number variation"""
    return {
        "type": "RelativeCopyNumber",
        "id": "ga4gh:RCN.dfqRkwlXaJqc8ZtG5mORZU4Cdsp3DTcz",
        "location": genomic_dup1_37_loc,
        "relative_copy_class": "complete loss"
    }


@pytest.fixture(scope="module")
def genomic_dup2_rel_38(genomic_dup2_seq_loc):
    """Create test fixture relative copy number variation"""
    return {
        "type": "RelativeCopyNumber",
        "id": "ga4gh:RCN.U9sPvq7Ggxf3jzcJlTD_53dAaesWZ6-o",
        "location": genomic_dup2_seq_loc,
        "relative_copy_class": "partial loss"
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
def genomic_dup2_abs_37(genomic_dup2_37_loc):
    """Create test fixture absolute copy number variation"""
    return {
        "type": "AbsoluteCopyNumber",
        "id": "ga4gh:ACN.ElLSWlxkk5iSfG0E3IONQRXJUzIcBhED",
        "location": genomic_dup2_37_loc,
        "copies": {"type": "Number", "value": 3}
    }


@pytest.fixture(scope="module")
def genomic_dup2_rel_37(genomic_dup2_37_loc):
    """Create test fixture relative copy number variation"""
    return {
        "type": "RelativeCopyNumber",
        "id": "ga4gh:RCN.BmS2zuzMCgnrbgU5zZgqwaPEBCGV6Wxo",
        "location": genomic_dup2_37_loc,
        "relative_copy_class": "partial loss"
    }


@pytest.fixture(scope="module")
def genomic_dup3_abs_38(genomic_del3_dup3_loc):
    """Create test fixture absolute copy number variation"""
    return {
        "type": "AbsoluteCopyNumber",
        "id": "ga4gh:ACN.2ZUQcccwvtoGZ5LZZRUoDZp6218Y6sQK",
        "location": genomic_del3_dup3_loc,
        "copies": {"type": "Number", "value": 2}
    }


@pytest.fixture(scope="module")
def genomic_dup3_rel_38(genomic_del3_dup3_loc):
    """Create test fixture relative copy number variation"""
    return {
        "type": "RelativeCopyNumber",
        "id": "ga4gh:RCN.FoK9bxEWUcnG6yb4MQwuEOKnyvmJQHWQ",
        "location": genomic_del3_dup3_loc,
        "relative_copy_class": "high-level gain"
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
def genomic_dup3_abs_37(genomic_dup3_37_loc):
    """Create test fixture absolute copy number variation"""
    return {
        "type": "AbsoluteCopyNumber",
        "id": "ga4gh:ACN.3ZeHMw-dKq2hTNEL9zzihfCPP2ZQV6zl",
        "location": genomic_dup3_37_loc,
        "copies": {"type": "Number", "value": 2}
    }


@pytest.fixture(scope="module")
def genomic_dup3_rel_37(genomic_dup3_37_loc):
    """Create test fixture relative copy number variation"""
    return {
        "type": "RelativeCopyNumber",
        "id": "ga4gh:RCN.oBHqFRckDEI_Z5urrQ_oexexObWv39FS",
        "location": genomic_dup3_37_loc,
        "relative_copy_class": "high-level gain"
    }


@pytest.fixture(scope="module")
def genomic_dup4_abs_38(genoimc_dup4_loc):
    """Create test fixture absolute copy number variation"""
    return {
        "type": "AbsoluteCopyNumber",
        "id": "ga4gh:ACN.OBSQ6SY9waWJ7YOyYA_qK0te_mTfOT4A",
        "location": genoimc_dup4_loc,
        "copies": {"type": "Number", "value": 3}
    }


@pytest.fixture(scope="module")
def genomic_dup4_rel_38(genoimc_dup4_loc):
    """Create test fixture relative copy number variation"""
    return {
        "type": "RelativeCopyNumber",
        "id": "ga4gh:RCN.4aCUMyIGHAaqGLGnWdGF3pU81nMPiRMf",
        "location": genoimc_dup4_loc,
        "relative_copy_class": "complete loss"
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
def genomic_dup4_abs_37(genomic_dup4_37_loc):
    """Create test fixture relative copy number variation"""
    return {
        "type": "AbsoluteCopyNumber",
        "id": "ga4gh:ACN.rloHun4f85F8QJ__SxoKMJFxA5IiNTMB",
        "location": genomic_dup4_37_loc,
        "copies": {"type": "Number", "value": 3}
    }


@pytest.fixture(scope="module")
def genomic_dup4_rel_37(genomic_dup4_37_loc):
    """Create test fixture relative copy number variation"""
    return {
        "type": "RelativeCopyNumber",
        "id": "ga4gh:RCN.U5eCTbzeX2JOyWFfQ9xOYacjbBKcz4lM",
        "location": genomic_dup4_37_loc,
        "relative_copy_class": "complete loss"
    }


@pytest.fixture(scope="module")
def genomic_dup5_abs_38(genomic_dup5_loc):
    """Create test fixture absolute copy number variation"""
    return {
        "type": "AbsoluteCopyNumber",
        "id": "ga4gh:ACN.oW400VdyoT8TChxzOLVl4oQLonMStzkK",
        "location": genomic_dup5_loc,
        "copies": {"type": "Number", "value": 4}
    }


@pytest.fixture(scope="module")
def genomic_dup5_rel_38(genomic_dup5_loc):
    """Create test fixture relative copy number variation"""
    return {
        "type": "RelativeCopyNumber",
        "id": "ga4gh:RCN.y4ia336ms3p6s51U36K5-kkXeW9PFXMz",
        "location": genomic_dup5_loc,
        "relative_copy_class": "partial loss"
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
def genomic_dup5_abs_37(genomic_dup5_37_loc):
    """Create test fixture absolute copy number variation"""
    return {
        "type": "AbsoluteCopyNumber",
        "id": "ga4gh:ACN.1k68ebEgFbNe9-FonMK620m_wbBkFIyZ",
        "location": genomic_dup5_37_loc,
        "copies": {"type": "Number", "value": 4}
    }


@pytest.fixture(scope="module")
def genomic_dup5_rel_37(genomic_dup5_37_loc):
    """Create test fixture relative copy number variation"""
    return {
        "type": "RelativeCopyNumber",
        "id": "ga4gh:RCN.AQqRwFtHFlvZG9Lcl28flKTITe_9wR3Z",
        "location": genomic_dup5_37_loc,
        "relative_copy_class": "partial loss"
    }


@pytest.fixture(scope="module")
def genomic_dup6_abs_38(genoimc_dup6_loc):
    """Create test fixture absolute copy number variation"""
    return {
        "type": "AbsoluteCopyNumber",
        "id": "ga4gh:ACN.6q8D5d_ie9MO0HdNEJmRJmGMg8C5LdAM",
        "location": genoimc_dup6_loc,
        "copies": {"type": "Number", "value": 2}
    }


@pytest.fixture(scope="module")
def genomic_dup6_rel_38(genoimc_dup6_loc):
    """Create test fixture relative copy number variation"""
    return {
        "type": "RelativeCopyNumber",
        "id": "ga4gh:RCN.popti3moTOvPqsqPopElu7-TqgINIq6I",
        "location": genoimc_dup6_loc,
        "relative_copy_class": "copy neutral"
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
def genomic_dup6_abs_37(genomic_dup6_37_loc):
    """Create test fixture absolute copy number variation"""
    return {
        "type": "AbsoluteCopyNumber",
        "id": "ga4gh:ACN.3RqSmp3UqxFy_1V7M8bYWEVi3vmIvM_C",
        "location": genomic_dup6_37_loc,
        "copies": {"type": "Number", "value": 2}
    }


@pytest.fixture(scope="module")
def genomic_dup6_rel_37(genomic_dup6_37_loc):
    """Create test fixture relative copy number variation"""
    return {
        "type": "RelativeCopyNumber",
        "id": "ga4gh:RCN.Lvhuup_Md0tsrN6eqO5kV3sUO4vRtMsW",
        "location": genomic_dup6_37_loc,
        "relative_copy_class": "copy neutral"
    }


@pytest.fixture(scope="module")
def genomic_del1_rel_38(genomic_del1_seq_loc):
    """Create test fixture relative copy number variation"""
    return {
        "type": "RelativeCopyNumber",
        "id": "ga4gh:RCN.z7MU8QUSR_aeWG7MP161H4jwPGoyo1No",
        "location": genomic_del1_seq_loc,
        "relative_copy_class": "copy neutral"
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
def genomic_del1_abs_37(genomic_del1_37_loc):
    """Create test fixture absolute copy number variation"""
    return {
        "type": "AbsoluteCopyNumber",
        "id": "ga4gh:ACN.oZ1GlVH9riWL7autRckjfTdef7MqVHex",
        "location": genomic_del1_37_loc,
        "copies": {"type": "Number", "value": 1}
    }


@pytest.fixture(scope="module")
def genomic_del1_rel_37(genomic_del1_37_loc):
    """Create test fixture relative copy number variation"""
    return {
        "type": "RelativeCopyNumber",
        "id": "ga4gh:RCN.bT-NgstX8QjwhOR3FMbQmfXozLQ9wdHB",
        "location": genomic_del1_37_loc,
        "relative_copy_class": "copy neutral"
    }


@pytest.fixture(scope="module")
def genomic_del2_rel_38(genomic_del2_seq_loc):
    """Create test fixture relative copy number variation"""
    return {
        "type": "RelativeCopyNumber",
        "id": "ga4gh:RCN.xXa8l2bXplY26DKNBDTttUUc0aN1Pwpo",
        "location": genomic_del2_seq_loc,
        "relative_copy_class": "low-level gain"
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
def genomic_del2_abs_37(genomic_del2_37_loc):
    """Create test fixture absolute copy number variation"""
    return {
        "type": "AbsoluteCopyNumber",
        "id": "ga4gh:ACN.doD8TUapoh2ekJQalFu2pqhU6RgMg2VI",
        "location": genomic_del2_37_loc,
        "copies": {"type": "Number", "value": 1}
    }


@pytest.fixture(scope="module")
def genomic_del2_rel_37(genomic_del2_37_loc):
    """Create test fixture relative copy number variation"""
    return {
        "type": "RelativeCopyNumber",
        "id": "ga4gh:RCN.Kw5y22t_9O3R_DBfP43ZgXkuLtrjbK2e",
        "location": genomic_del2_37_loc,
        "relative_copy_class": "low-level gain"
    }


@pytest.fixture(scope="module")
def genomic_del3_abs_38(genomic_del3_dup3_loc):
    """Create test fixture absolute copy number variation"""
    return {
        "type": "AbsoluteCopyNumber",
        "id": "ga4gh:ACN.2ZUQcccwvtoGZ5LZZRUoDZp6218Y6sQK",
        "location": genomic_del3_dup3_loc,
        "copies": {"type": "Number", "value": 2}
    }


@pytest.fixture(scope="module")
def genomic_del3_rel_38(genomic_del3_dup3_loc):
    """Create test fixture relative copy number variation"""
    return {
        "type": "RelativeCopyNumber",
        "id": "ga4gh:RCN.RoVV-V_2zSMFA8IUJSs3Ah8Y-3jh6ktV",
        "location": genomic_del3_dup3_loc,
        "relative_copy_class": "complete loss"
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
def genomic_del3_abs_37(genomic_del3_37_loc):
    """Create test fixture absolute copy number variation"""
    return {
        "type": "AbsoluteCopyNumber",
        "id": "ga4gh:ACN.3ZeHMw-dKq2hTNEL9zzihfCPP2ZQV6zl",
        "location": genomic_del3_37_loc,
        "copies": {"type": "Number", "value": 2}
    }


@pytest.fixture(scope="module")
def genomic_del3_rel_37(genomic_del3_37_loc):
    """Create test fixture relative copy number variation"""
    return {
        "type": "RelativeCopyNumber",
        "id": "ga4gh:RCN.rKupZ9FRkAXrYtZmaMd_7RF6vU4fOxLS",
        "location": genomic_del3_37_loc,
        "relative_copy_class": "complete loss"
    }


@pytest.fixture(scope="module")
def genomic_del4_abs_38(genomic_del4_seq_loc):
    """Create test fixture absolute copy number variation"""
    return {
        "type": "AbsoluteCopyNumber",
        "id": "ga4gh:ACN.RTb0O1AM8twFPfcArsuu_Hiffs98GcXm",
        "location": genomic_del4_seq_loc,
        "copies": {"type": "Number", "value": 4}
    }


@pytest.fixture(scope="module")
def genomic_del4_rel_38(genomic_del4_seq_loc):
    """Create test fixture relative copy number variation"""
    return {
        "type": "RelativeCopyNumber",
        "id": "ga4gh:RCN.BwZOFAfo5u8TcwbR3DMi8qbIImv96VQU",
        "location": genomic_del4_seq_loc,
        "relative_copy_class": "partial loss"
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
def genomic_del4_abs_37(genomic_del4_37_loc):
    """Create test fixture absolute copy number variation"""
    return {
        "type": "AbsoluteCopyNumber",
        "id": "ga4gh:ACN.4Lg4kT2VZ4yCjntMcVXy7wyrbAGflrTp",
        "location": genomic_del4_37_loc,
        "copies": {"type": "Number", "value": 4}
    }


@pytest.fixture(scope="module")
def genomic_del4_rel_37(genomic_del4_37_loc):
    """Create test fixture relative copy number variation"""
    return {
        "type": "RelativeCopyNumber",
        "id": "ga4gh:RCN.xeNugTAGZ5HPU5hyoOa6Jk_lzQKFj-2S",
        "location": genomic_del4_37_loc,
        "relative_copy_class": "partial loss"
    }


@pytest.fixture(scope="module")
def genomic_del5_abs_38(genomic_del5_seq_loc):
    """Create test fixture absolute copy number variation"""
    return {
        "type": "AbsoluteCopyNumber",
        "id": "ga4gh:ACN.2su0y7vCqZmbTBU1cqKTOAR7m9vkNqBW",
        "location": genomic_del5_seq_loc,
        "copies": {"type": "Number", "value": 2}
    }


@pytest.fixture(scope="module")
def genomic_del5_rel_38(genomic_del5_seq_loc):
    """Create test fixture relative copy number variation"""
    return {
        "type": "RelativeCopyNumber",
        "id": "ga4gh:RCN.y5LOigojbNT1BqtAv9BsKQ7-2i-iL8jA",
        "location": genomic_del5_seq_loc,
        "relative_copy_class": "copy neutral"
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
def genomic_del5_abs_37(genomic_del5_37_loc):
    """Create test fixture absolute copy number variation"""
    return {
        "type": "AbsoluteCopyNumber",
        "id": "ga4gh:ACN.GxN5Fwl4w9NpVZSOXdxOMvrpPLmGWEd_",
        "location": genomic_del5_37_loc,
        "copies": {"type": "Number", "value": 2}
    }


@pytest.fixture(scope="module")
def genomic_del5_rel_37(genomic_del5_37_loc):
    """Create test fixture relative copy number variation"""
    return {
        "type": "RelativeCopyNumber",
        "id": "ga4gh:RCN.jExtRl5UIhAlKPw68TtcA05-Gfi6cXQZ",
        "location": genomic_del5_37_loc,
        "relative_copy_class": "copy neutral"
    }


@pytest.fixture(scope="module")
def genomic_del6_abs_38(genomic_del6_seq_loc):
    """Create test fixture absolute copy number variation"""
    return {
        "type": "AbsoluteCopyNumber",
        "id": "ga4gh:ACN.ZnnJNutwCrHNzFQamAWXMbLC7PfILmqA",
        "location": genomic_del6_seq_loc,
        "copies": {"type": "Number", "value": 1}
    }


@pytest.fixture(scope="module")
def genomic_del6_rel_38(genomic_del6_seq_loc):
    """Create test fixture relative copy number variation"""
    return {
        "type": "RelativeCopyNumber",
        "id": "ga4gh:RCN.nmDpVRjFPJZF4K-HxmIgcAewOKtCOfGC",
        "location": genomic_del6_seq_loc,
        "relative_copy_class": "low-level gain"
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
def genomic_del6_abs_37(genomic_del6_37_loc):
    """Create test fixture absolute copy number variation"""
    return {
        "type": "AbsoluteCopyNumber",
        "id": "ga4gh:ACN.aOQqe7zijdBB8VKaMSW97V2YVPeBnYMv",
        "location": genomic_del6_37_loc,
        "copies": {"type": "Number", "value": 1}
    }


@pytest.fixture(scope="module")
def genomic_del6_rel_37(genomic_del6_37_loc):
    """Create test fixture relative copy number variation"""
    return {
        "type": "RelativeCopyNumber",
        "id": "ga4gh:RCN.C597L06IPEjWaoP-ktRkPxbbayBWBg12",
        "location": genomic_del6_37_loc,
        "relative_copy_class": "low-level gain"
    }


@pytest.mark.asyncio
async def test_genomic_dup1_absolute_cnv(test_cnv_handler, genomic_dup1_38_vac,
                                         genomic_dup1_abs_37):
    """Test that genomic duplication works correctly"""
    q = "NC_000003.12:g.49531262dup"  # 38
    resp = await test_cnv_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=2, do_liftover=False, )
    assert resp.absolute_copy_number.dict(by_alias=True) == genomic_dup1_38_vac

    q = "NC_000003.11:g.49568695dup"  # 37
    resp = await test_cnv_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=2, do_liftover=False)
    assert resp.absolute_copy_number.dict(by_alias=True) == genomic_dup1_abs_37

    resp = await test_cnv_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=2, do_liftover=True)
    assert resp.absolute_copy_number.dict(by_alias=True) == genomic_dup1_38_vac

    resp = await test_cnv_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=1, do_liftover=True)
    expected = copy.deepcopy(genomic_dup1_38_vac)
    expected["copies"]["value"] = 2
    expected["id"] = "ga4gh:ACN.1IRAq_PwqIELC3n5kjZsr5X-xLmzGxrL"
    assert resp.absolute_copy_number.dict(by_alias=True) == expected


@pytest.mark.asyncio
async def test_genomic_dup1_relative_cnv(test_cnv_handler, genomic_dup1_rel_38,
                                         genomic_dup1_rel_37):
    """Test that genomic duplication works correctly"""
    q = "NC_000003.12:g.49531262dup"  # 38
    resp = await test_cnv_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="complete loss", do_liftover=False)
    assert resp.relative_copy_number.dict(by_alias=True) == genomic_dup1_rel_38

    q = "NC_000003.11:g.49568695dup"  # 37
    resp = await test_cnv_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="complete loss", do_liftover=False)
    assert resp.relative_copy_number.dict(by_alias=True) == genomic_dup1_rel_37

    resp = await test_cnv_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="complete loss", do_liftover=True)
    assert resp.relative_copy_number.dict(by_alias=True) == genomic_dup1_rel_38


@pytest.mark.asyncio
async def test_genomic_dup2_absolute_cnv(test_cnv_handler, genomic_dup2_38_vac,
                                         genomic_dup2_abs_37):
    """Test that genomic duplication works correctly"""
    q = "NC_000016.10:g.2087938_2087948dup"  # 38
    resp = await test_cnv_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=2, do_liftover=False)
    assert resp.absolute_copy_number.dict(by_alias=True) == genomic_dup2_38_vac

    q = "NC_000016.9:g.2137939_2137949dup"  # 37
    resp = await test_cnv_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=2, do_liftover=False)
    assert resp.absolute_copy_number.dict(by_alias=True) == genomic_dup2_abs_37

    resp = await test_cnv_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=2, do_liftover=True)
    assert resp.absolute_copy_number.dict(by_alias=True) == genomic_dup2_38_vac

    resp = await test_cnv_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=1, do_liftover=True)
    expected = copy.deepcopy(genomic_dup2_38_vac)
    expected["copies"]["value"] = 2
    expected["id"] = "ga4gh:ACN.rDCO3xAHcwYBhdNPP8qqzy6t62jrNFYs"
    assert resp.absolute_copy_number.dict(by_alias=True) == expected


@pytest.mark.asyncio
async def test_genomic_dup2_relative_cnv(test_cnv_handler, genomic_dup2_rel_38,
                                         genomic_dup2_rel_37):
    """Test that genomic duplication works correctly"""
    q = "NC_000016.10:g.2087938_2087948dup"  # 38
    resp = await test_cnv_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="partial loss", do_liftover=False)
    assert resp.relative_copy_number.dict(by_alias=True) == genomic_dup2_rel_38

    q = "NC_000016.9:g.2137939_2137949dup"  # 37
    resp = await test_cnv_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="partial loss", do_liftover=False)
    assert resp.relative_copy_number.dict(by_alias=True) == genomic_dup2_rel_37

    resp = await test_cnv_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="partial loss", do_liftover=True)
    assert resp.relative_copy_number.dict(by_alias=True) == genomic_dup2_rel_38


@pytest.mark.asyncio
async def test_genomic_dup3_absolute_cnv(test_cnv_handler, genomic_dup3_abs_38,
                                         genomic_dup3_abs_37):
    """Test that genomic duplication works correctly"""
    q = "NC_000023.11:g.(31060227_31100351)_(33274278_33417151)dup"  # 38
    resp = await test_cnv_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=1, do_liftover=False)
    assert resp.absolute_copy_number.dict(by_alias=True) == genomic_dup3_abs_38

    q = "NC_000023.10:g.(31078344_31118468)_(33292395_33435268)dup"  # 37
    resp = await test_cnv_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=1, do_liftover=False)
    assert resp.absolute_copy_number.dict(by_alias=True) == genomic_dup3_abs_37

    resp = await test_cnv_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=1, do_liftover=True)
    assert resp.absolute_copy_number.dict(by_alias=True) == genomic_dup3_abs_38

    resp = await test_cnv_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=2, do_liftover=True)
    expected = copy.deepcopy(genomic_dup3_abs_38)
    expected["copies"] = {"value": 3, "type": "Number"}
    expected["id"] = "ga4gh:ACN.I1zTu1vvN_OIrSWh9mbPLFyZQnyl1jsu"
    assert resp.absolute_copy_number.dict(by_alias=True) == expected


@pytest.mark.asyncio
async def test_genomic_dup3_relative_cnv(test_cnv_handler, genomic_dup3_rel_38,
                                         genomic_dup3_rel_37):
    """Test that genomic duplication works correctly"""
    q = "NC_000023.11:g.(31060227_31100351)_(33274278_33417151)dup"  # 38
    resp = await test_cnv_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="high-level gain", do_liftover=False)
    assert resp.relative_copy_number.dict(by_alias=True) == genomic_dup3_rel_38

    q = "NC_000023.10:g.(31078344_31118468)_(33292395_33435268)dup"  # 37
    resp = await test_cnv_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="high-level gain", do_liftover=False)
    assert resp.relative_copy_number.dict(by_alias=True) == genomic_dup3_rel_37

    resp = await test_cnv_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="high-level gain", do_liftover=True)
    assert resp.relative_copy_number.dict(by_alias=True) == genomic_dup3_rel_38


@pytest.mark.asyncio
async def test_genomic_dup4_absolute_cnv(test_cnv_handler, genomic_dup4_abs_38,
                                         genomic_dup4_abs_37):
    """Test that genomic duplication works correctly"""
    q = "NC_000020.11:g.(?_30417576)_(31394018_?)dup"  # 38
    resp = await test_cnv_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=2, do_liftover=False)
    assert resp.absolute_copy_number.dict(by_alias=True) == genomic_dup4_abs_38

    q = "NC_000020.10:g.(?_29652252)_(29981821_?)dup"  # 37
    resp = await test_cnv_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=2, do_liftover=False)
    assert resp.absolute_copy_number.dict(by_alias=True) == genomic_dup4_abs_37

    resp = await test_cnv_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=2, do_liftover=True)
    assert resp.absolute_copy_number.dict(by_alias=True) == genomic_dup4_abs_38

    resp = await test_cnv_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=1, do_liftover=True)
    expected = copy.deepcopy(genomic_dup4_abs_38)
    expected["copies"]["value"] = 2
    expected["id"] = "ga4gh:ACN.jLffc6gpShAdrPSDMbn4r9-a3KOBIk0h"
    assert resp.absolute_copy_number.dict(by_alias=True) == expected


@pytest.mark.asyncio
async def test_genomic_dup4_relative_cnv(test_cnv_handler, genomic_dup4_rel_38,
                                         genomic_dup4_rel_37):
    """Test that genomic duplication works correctly"""
    q = "NC_000020.11:g.(?_30417576)_(31394018_?)dup"  # 38
    resp = await test_cnv_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="complete loss", do_liftover=False)
    assert resp.relative_copy_number.dict(by_alias=True) == genomic_dup4_rel_38

    q = "NC_000020.10:g.(?_29652252)_(29981821_?)dup"  # 37
    resp = await test_cnv_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="complete loss", do_liftover=False)
    assert resp.relative_copy_number.dict(by_alias=True) == genomic_dup4_rel_37

    resp = await test_cnv_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="complete loss", do_liftover=True)
    assert resp.relative_copy_number.dict(by_alias=True) == genomic_dup4_rel_38


@pytest.mark.asyncio
async def test_genomic_dup5_absolute_cnv(test_cnv_handler, genomic_dup5_abs_38,
                                         genomic_dup5_abs_37):
    """Test that genomic duplication works correctly"""
    q = "NC_000023.11:g.(?_154021812)_154092209dup"  # 38
    resp = await test_cnv_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=3, do_liftover=False)
    assert resp.absolute_copy_number.dict(by_alias=True) == genomic_dup5_abs_38

    q = "NC_000023.10:g.(?_153287263)_153357667dup"  # 37
    resp = await test_cnv_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=3, do_liftover=False)
    assert resp.absolute_copy_number.dict(by_alias=True) == genomic_dup5_abs_37

    resp = await test_cnv_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=3, do_liftover=True)
    assert resp.absolute_copy_number.dict(by_alias=True) == genomic_dup5_abs_38

    resp = await test_cnv_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=4, do_liftover=True)
    expected = copy.deepcopy(genomic_dup5_abs_38)
    expected["copies"] = {"value": 5, "type": "Number"}
    expected["id"] = "ga4gh:ACN.74oGYuYCD0rAx1tjbTCG6EeZkfH2hLx9"
    assert resp.absolute_copy_number.dict(by_alias=True) == expected


@pytest.mark.asyncio
async def test_genomic_dup5_relative_cnv(test_cnv_handler, genomic_dup5_rel_38,
                                         genomic_dup5_rel_37):
    """Test that genomic duplication works correctly"""
    q = "NC_000023.11:g.(?_154021812)_154092209dup"  # 38
    resp = await test_cnv_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="partial loss", do_liftover=False)
    assert resp.relative_copy_number.dict(by_alias=True) == genomic_dup5_rel_38

    q = "NC_000023.10:g.(?_153287263)_153357667dup"  # 37
    resp = await test_cnv_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="partial loss", do_liftover=False)
    assert resp.relative_copy_number.dict(by_alias=True) == genomic_dup5_rel_37

    resp = await test_cnv_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="partial loss", do_liftover=True)
    assert resp.relative_copy_number.dict(by_alias=True) == genomic_dup5_rel_38


@pytest.mark.asyncio
async def test_genomic_dup6_absolute_cnv(test_cnv_handler, genomic_dup6_abs_38,
                                         genomic_dup6_abs_37):
    """Test that genomic duplication works correctly"""
    q = "NC_000023.11:g.154021812_(154092209_?)dup"  # 38
    resp = await test_cnv_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=1, do_liftover=False)
    assert resp.absolute_copy_number.dict(by_alias=True) == genomic_dup6_abs_38

    q = "NC_000023.10:g.153287263_(153357667_?)dup"  # 37
    resp = await test_cnv_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=1, do_liftover=False)
    assert resp.absolute_copy_number.dict(by_alias=True) == genomic_dup6_abs_37

    resp = await test_cnv_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=1, do_liftover=True)
    assert resp.absolute_copy_number.dict(by_alias=True) == genomic_dup6_abs_38

    resp = await test_cnv_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=2, do_liftover=True)
    expected = copy.deepcopy(genomic_dup6_abs_38)
    expected["copies"] = {"value": 3, "type": "Number"}
    expected["id"] = "ga4gh:ACN.uxA-iMw_FRIvuh9EBS8E-Loqr-Kzjd0c"
    assert resp.absolute_copy_number.dict(by_alias=True) == expected


@pytest.mark.asyncio
async def test_genomic_dup6_relative_cnv(test_cnv_handler, genomic_dup6_rel_38,
                                         genomic_dup6_rel_37):
    """Test that genomic duplication works correctly"""
    q = "NC_000023.11:g.154021812_(154092209_?)dup"  # 38
    resp = await test_cnv_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="copy neutral", do_liftover=False)
    assert resp.relative_copy_number.dict(by_alias=True) == genomic_dup6_rel_38

    q = "NC_000023.10:g.153287263_(153357667_?)dup"  # 37
    resp = await test_cnv_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="copy neutral", do_liftover=False)
    assert resp.relative_copy_number.dict(by_alias=True) == genomic_dup6_rel_37

    resp = await test_cnv_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="copy neutral", do_liftover=True)
    assert resp.relative_copy_number.dict(by_alias=True) == genomic_dup6_rel_38


@pytest.mark.asyncio
async def test_genomic_del1_absolute_cnv(test_cnv_handler, genomic_del1_38_vac,
                                         genomic_del1_abs_37):
    """Test that genomic deletion works correctly"""
    q = "NC_000003.12:g.10149811del"  # 38
    resp = await test_cnv_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=2, do_liftover=False)
    assert resp.absolute_copy_number.dict(by_alias=True) == genomic_del1_38_vac

    q = "NC_000003.11:g.10191495del"  # 37
    resp = await test_cnv_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=2, do_liftover=False)
    assert resp.absolute_copy_number.dict(by_alias=True) == genomic_del1_abs_37

    resp = await test_cnv_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=2, do_liftover=True)
    assert resp.absolute_copy_number.dict(by_alias=True) == genomic_del1_38_vac

    resp = await test_cnv_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=3, do_liftover=True)
    expected = copy.deepcopy(genomic_del1_38_vac)
    expected["copies"]["value"] = 2
    expected["id"] = "ga4gh:ACN.6fvLDf6cBMYpIAu9xbjSWn8YL4ru0eM-"
    assert resp.absolute_copy_number.dict(by_alias=True) == expected


@pytest.mark.asyncio
async def test_genomic_del1_relative_cnv(test_cnv_handler, genomic_del1_rel_38,
                                         genomic_del1_rel_37):
    """Test that genomic deletion works correctly"""
    q = "NC_000003.12:g.10149811del"  # 38
    resp = await test_cnv_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="copy neutral", do_liftover=False)
    assert resp.relative_copy_number.dict(by_alias=True) == genomic_del1_rel_38

    q = "NC_000003.11:g.10191495del"  # 37
    resp = await test_cnv_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="copy neutral", do_liftover=False)
    assert resp.relative_copy_number.dict(by_alias=True) == genomic_del1_rel_37

    resp = await test_cnv_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="copy neutral", do_liftover=True)
    assert resp.relative_copy_number.dict(by_alias=True) == genomic_del1_rel_38


@pytest.mark.asyncio
async def test_genomic_del2_absolute_cnv(test_cnv_handler, genomic_del2_38_vac,
                                         genomic_del2_abs_37):
    """Test that genomic deletion works correctly"""
    q = "NC_000003.12:g.10146595_10146613del"  # 38
    resp = await test_cnv_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=2, do_liftover=False)
    assert resp.absolute_copy_number.dict(by_alias=True) == genomic_del2_38_vac

    q = "NC_000003.11:g.10188279_10188297del"  # 37
    resp = await test_cnv_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=2, do_liftover=False)
    assert resp.absolute_copy_number.dict(by_alias=True) == genomic_del2_abs_37

    resp = await test_cnv_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=2, do_liftover=True)
    assert resp.absolute_copy_number.dict(by_alias=True) == genomic_del2_38_vac

    resp = await test_cnv_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=4, do_liftover=True)
    expected = copy.deepcopy(genomic_del2_38_vac)
    expected["copies"]["value"] = 3
    expected["id"] = "ga4gh:ACN.1mf7xNFTZu8JsC1twohPFCEF4IY_OPhq"
    assert resp.absolute_copy_number.dict(by_alias=True) == expected


@pytest.mark.asyncio
async def test_genomic_del2_relative_cnv(test_cnv_handler, genomic_del2_rel_38,
                                         genomic_del2_rel_37):
    """Test that genomic deletion works correctly"""
    q = "NC_000003.12:g.10146595_10146613del"  # 38
    resp = await test_cnv_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="low-level gain", do_liftover=False)
    assert resp.relative_copy_number.dict(by_alias=True) == genomic_del2_rel_38

    q = "NC_000003.11:g.10188279_10188297del"  # 37
    resp = await test_cnv_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="low-level gain", do_liftover=False)
    assert resp.relative_copy_number.dict(by_alias=True) == genomic_del2_rel_37

    resp = await test_cnv_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="low-level gain", do_liftover=True)
    assert resp.relative_copy_number.dict(by_alias=True) == genomic_del2_rel_38


@pytest.mark.asyncio
async def test_genomic_del3_absolute_cnv(test_cnv_handler, genomic_del3_abs_38,
                                         genomic_del3_abs_37):
    """Test that genomic deletion works correctly"""
    q = "NC_000023.11:g.(31060227_31100351)_(33274278_33417151)del"  # 38
    resp = await test_cnv_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=3, do_liftover=False)
    assert resp.absolute_copy_number.dict(by_alias=True) == genomic_del3_abs_38

    q = "NC_000023.10:g.(31078344_31118468)_(33292395_33435268)del"  # 37
    resp = await test_cnv_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=3, do_liftover=False)
    assert resp.absolute_copy_number.dict(by_alias=True) == genomic_del3_abs_37

    resp = await test_cnv_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=3, do_liftover=True)
    assert resp.absolute_copy_number.dict(by_alias=True) == genomic_del3_abs_38

    resp = await test_cnv_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=2, do_liftover=True)
    expected = copy.deepcopy(genomic_del3_abs_38)
    expected["copies"] = {"value": 1, "type": "Number"}
    expected["id"] = "ga4gh:ACN.Q41L3gg0DJBgQj4BSVIdCEt41Hkh0U_S"
    assert resp.absolute_copy_number.dict(by_alias=True) == expected


@pytest.mark.asyncio
async def test_genomic_del3_relative_cnv(test_cnv_handler, genomic_del3_rel_38,
                                         genomic_del3_rel_37):
    """Test that genomic deletion works correctly"""
    q = "NC_000023.11:g.(31060227_31100351)_(33274278_33417151)del"  # 38
    resp = await test_cnv_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="complete loss", do_liftover=False)
    assert resp.relative_copy_number.dict(by_alias=True) == genomic_del3_rel_38

    q = "NC_000023.10:g.(31078344_31118468)_(33292395_33435268)del"  # 37
    resp = await test_cnv_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="complete loss", do_liftover=False)
    assert resp.relative_copy_number.dict(by_alias=True) == genomic_del3_rel_37

    resp = await test_cnv_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="complete loss", do_liftover=True)
    assert resp.relative_copy_number.dict(by_alias=True) == genomic_del3_rel_38


@pytest.mark.asyncio
async def test_genomic_del4_absolute_cnv(test_cnv_handler, genomic_del4_abs_38,
                                         genomic_del4_abs_37):
    """Test that genomic deletion works correctly"""
    q = "NC_000023.11:g.(?_31120496)_(33339477_?)del"  # 38
    resp = await test_cnv_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=5, do_liftover=False)
    assert resp.absolute_copy_number.dict(by_alias=True) == genomic_del4_abs_38

    q = "NC_000023.10:g.(?_31138613)_(33357594_?)del"  # 37
    resp = await test_cnv_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=5, do_liftover=False)
    assert resp.absolute_copy_number.dict(by_alias=True) == genomic_del4_abs_37

    resp = await test_cnv_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=5, do_liftover=True)
    assert resp.absolute_copy_number.dict(by_alias=True) == genomic_del4_abs_38

    resp = await test_cnv_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=3, do_liftover=True)
    expected = copy.deepcopy(genomic_del4_abs_38)
    expected["copies"] = {"value": 2, "type": "Number"}
    expected["id"] = "ga4gh:ACN.YwVdXUJiRe_yLUbJT1WZ_4vKXOPhUk_-"
    assert resp.absolute_copy_number.dict(by_alias=True) == expected


@pytest.mark.asyncio
async def test_genomic_del4_relative_cnv(test_cnv_handler, genomic_del4_rel_38,
                                         genomic_del4_rel_37):
    """Test that genomic deletion works correctly"""
    q = "NC_000023.11:g.(?_31120496)_(33339477_?)del"  # 38
    resp = await test_cnv_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="partial loss", do_liftover=False)
    assert resp.relative_copy_number.dict(by_alias=True) == genomic_del4_rel_38

    q = "NC_000023.10:g.(?_31138613)_(33357594_?)del"  # 37
    resp = await test_cnv_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="partial loss", do_liftover=False)
    assert resp.relative_copy_number.dict(by_alias=True) == genomic_del4_rel_37

    resp = await test_cnv_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="partial loss", do_liftover=True)
    assert resp.relative_copy_number.dict(by_alias=True) == genomic_del4_rel_38


@pytest.mark.asyncio
async def test_genomic_del5_absolute_cnv(test_cnv_handler, genomic_del5_abs_38,
                                         genomic_del5_abs_37):
    """Test that genomic deletion works correctly"""
    q = "NC_000023.11:g.(?_18575354)_18653629del"  # 38
    resp = await test_cnv_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=3, do_liftover=False)
    assert resp.absolute_copy_number.dict(by_alias=True) == genomic_del5_abs_38

    q = "NC_000023.10:g.(?_18593474)_18671749del"  # 37
    resp = await test_cnv_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=3, do_liftover=False)
    assert resp.absolute_copy_number.dict(by_alias=True) == genomic_del5_abs_37

    resp = await test_cnv_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=3, do_liftover=True)
    assert resp.absolute_copy_number.dict(by_alias=True) == genomic_del5_abs_38

    resp = await test_cnv_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=2, do_liftover=True)
    expected = copy.deepcopy(genomic_del5_abs_38)
    expected["copies"] = {"value": 1, "type": "Number"}
    expected["id"] = "ga4gh:ACN.MPb2oPn0bAT-Cqz3UNjtf9EQ4JT2sIyR"
    assert resp.absolute_copy_number.dict(by_alias=True) == expected


@pytest.mark.asyncio
async def test_genomic_del5_relative_cnv(test_cnv_handler, genomic_del5_rel_38,
                                         genomic_del5_rel_37):
    """Test that genomic deletion works correctly"""
    q = "NC_000023.11:g.(?_18575354)_18653629del"  # 38
    resp = await test_cnv_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="copy neutral", do_liftover=False)
    assert resp.relative_copy_number.dict(by_alias=True) == genomic_del5_rel_38

    q = "NC_000023.10:g.(?_18593474)_18671749del"  # 37
    resp = await test_cnv_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="copy neutral", do_liftover=False)
    assert resp.relative_copy_number.dict(by_alias=True) == genomic_del5_rel_37

    resp = await test_cnv_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="copy neutral", do_liftover=True)
    assert resp.relative_copy_number.dict(by_alias=True) == genomic_del5_rel_38


@pytest.mark.asyncio
async def test_genomic_del6_absolute_cnv(test_cnv_handler, genomic_del6_abs_38,
                                         genomic_del6_abs_37):
    """Test that genomic deletion works correctly"""
    q = "NC_000006.12:g.133462764_(133464858_?)del"  # 38
    resp = await test_cnv_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=2, do_liftover=False)
    assert resp.absolute_copy_number.dict(by_alias=True) == genomic_del6_abs_38

    q = "NC_000006.11:g.133783902_(133785996_?)del"  # 37
    resp = await test_cnv_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=2, do_liftover=False)
    assert resp.absolute_copy_number.dict(by_alias=True) == genomic_del6_abs_37

    resp = await test_cnv_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=2, do_liftover=True)
    assert resp.absolute_copy_number.dict(by_alias=True) == genomic_del6_abs_38

    resp = await test_cnv_handler.hgvs_to_absolute_copy_number(
        q, baseline_copies=3, do_liftover=True)
    expected = copy.deepcopy(genomic_del6_abs_38)
    expected["copies"]["value"] = 2
    expected["id"] = "ga4gh:ACN.FIfMZ5asuZULdupfFV5FmUxLeeszG4Vl"
    assert resp.absolute_copy_number.dict(by_alias=True) == expected


@pytest.mark.asyncio
async def test_genomic_del6_relative_cnv(test_cnv_handler, genomic_del6_rel_38,
                                         genomic_del6_rel_37):
    """Test that genomic deletion works correctly"""
    q = "NC_000006.12:g.133462764_(133464858_?)del"  # 38
    resp = await test_cnv_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="low-level gain", do_liftover=False)
    assert resp.relative_copy_number.dict(by_alias=True) == genomic_del6_rel_38

    q = "NC_000006.11:g.133783902_(133785996_?)del"  # 37
    resp = await test_cnv_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="low-level gain", do_liftover=False)
    assert resp.relative_copy_number.dict(by_alias=True) == genomic_del6_rel_37

    resp = await test_cnv_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="low-level gain", do_liftover=True)
    assert resp.relative_copy_number.dict(by_alias=True) == genomic_del6_rel_38


@pytest.mark.asyncio
async def test_invalid_cnv_parameters(test_cnv_handler):
    """Check that invalid parameters return warnings"""
    q = "NC_000006.11:g.133783902_(133785996_?)del"
    resp, w = await test_cnv_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="low-level gains", do_liftover=True)
    assert resp is None
    assert w == ["low-level gains is not a valid relative copy class: ['complete loss', "  # noqa: E501
                 "'partial loss', 'copy neutral', 'low-level gain', 'high-level gain']"]


@pytest.mark.asyncio
async def test_invalid_cnv(test_cnv_handler):
    """Check that invalid input return warnings"""
    q = "DAG1 g.49568695dup"
    resp = await test_cnv_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="low-level gain", do_liftover=True,
        untranslatable_returns_text=True)
    assert set(resp.warnings) == {"Unable to translate DAG1 g.49568695dup to copy number variation",  # noqa: E501
                                  "DAG1 g.49568695dup is not a supported HGVS genomic duplication or deletion"}  # noqa: E501
    assert resp.relative_copy_number.type == "Text"

    q = "braf v600e"
    resp = await test_cnv_handler.hgvs_to_relative_copy_number(
        q, relative_copy_class="low-level gain", do_liftover=True)
    assert set(resp.warnings) == {"Unable to translate braf v600e to copy number variation",  # noqa: E501
                                  "braf v600e is not a supported HGVS genomic duplication or deletion"}  # noqa: E501
    assert resp.relative_copy_number is None
