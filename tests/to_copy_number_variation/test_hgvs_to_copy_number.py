"""Module for testing the hgvs to copy number count and copy number change endpoints"""
import copy

import pytest


@pytest.fixture(scope="module")
def genomic_dup1_cx_38(genomic_dup1_seq_loc_not_normalized):
    """Create test fixture copy number change variation"""
    return {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.ZZY--6eqtK9AVqxFpjGSW6VIrDADfqJB",
        "subject": genomic_dup1_seq_loc_not_normalized,
        "copy_change": "efo:0030069",
    }


@pytest.fixture(scope="module")
def genomic_dup1_37_loc():
    """Create test fixture GRCh37 duplication subject"""
    return {
        "id": "ga4gh:SL.ppeF8R8M_I21kT6EvmHYMFW99ieBPZhG",
        "sequence_id": "ga4gh:SQ.VNBualIltAyi2AI_uXcKU7M9XUOuA7MS",
        "start": {"value": 49568694, "type": "Number"},
        "end": {"value": 49568695, "type": "Number"},
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_dup1_cn_37(genomic_dup1_37_loc):
    """Create test fixture copy number count variation (not normalized)"""
    return {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.eOYBWrgITvZ0Rsf8q_fpZnsSOkvOBB1v",
        "subject": genomic_dup1_37_loc,
        "copies": {"type": "Number", "value": 3},
    }


@pytest.fixture(scope="module")
def genomic_dup1_cx_37(genomic_dup1_37_loc):
    """Create test fixture copy number change variation"""
    return {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.NKYonqhqFZEqVxXJYoQQXf98FM618OP_",
        "subject": genomic_dup1_37_loc,
        "copy_change": "efo:0030069",
    }


@pytest.fixture(scope="module")
def genomic_dup2_cx_38(genomic_dup2_seq_loc_normalized):
    """Create test fixture copy number change variation"""
    return {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.dY9cDcR2xhOC1N9lEKuggXxBny42teDF",
        "subject": genomic_dup2_seq_loc_normalized,
        "copy_change": "efo:0030067",
    }


@pytest.fixture(scope="module")
def genomic_dup2_37_loc():
    """Create test fixture GRCh37 duplication subject"""
    return {
        "id": "ga4gh:SL.m2vBWKNkl-LU5mIJVCS7VoEeRyj5yxpf",
        "sequence_id": "ga4gh:SQ.v7noePfnNpK8ghYXEqZ9NukMXW7YeNsm",
        "start": {"value": 33229406, "type": "Number"},
        "end": {"value": 33229410, "type": "Number"},
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_dup2_cn_37(genomic_dup2_37_loc):
    """Create test fixture copy number count variation"""
    return {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.8SMv30enQyCLcbey8b5P4Zc_xolCle93",
        "subject": genomic_dup2_37_loc,
        "copies": {"type": "Number", "value": 3},
    }


@pytest.fixture(scope="module")
def genomic_dup2_cx_37(genomic_dup2_37_loc):
    """Create test fixture copy number change variation"""
    return {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.lKn0mkjm3OPx2IqbSYGSBxDiIXtthC1h",
        "subject": genomic_dup2_37_loc,
        "copy_change": "efo:0030067",
    }


@pytest.fixture(scope="module")
def genomic_dup3_cn_38(genomic_del3_dup3_loc_not_normalized):
    """Create test fixture copy number count variation"""
    return {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.a2-YlA01Q-y2zAOgn4r-cP9hiIEGaaqq",
        "subject": genomic_del3_dup3_loc_not_normalized,
        "copies": {"type": "Number", "value": 2},
    }


@pytest.fixture(scope="module")
def genomic_dup3_cx_38(genomic_del3_dup3_loc_not_normalized):
    """Create test fixture copy number change variation"""
    return {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.MGhVuAEVskM_cAhATCYT_tRcOsqDbHNd",
        "subject": genomic_del3_dup3_loc_not_normalized,
        "copy_change": "efo:0030072",
    }


@pytest.fixture(scope="module")
def genomic_dup3_37_loc():
    """Create test fixture GRCh37 duplication subject"""
    return {
        "id": "ga4gh:SL.MmBA0qlC0J6yxWTyqhpsd4H_oErccpUD",
        "sequence_id": "ga4gh:SQ.v7noePfnNpK8ghYXEqZ9NukMXW7YeNsm",
        "start": {"min": 31078343, "max": 31118467, "type": "DefiniteRange"},
        "end": {"min": 33292396, "max": 33435269, "type": "DefiniteRange"},
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_dup3_cn_37(genomic_dup3_37_loc):
    """Create test fixture copy number count variation"""
    return {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.GOhhUsvbhFqaXgq3hNH3nuJldMIwVMpO",
        "subject": genomic_dup3_37_loc,
        "copies": {"type": "Number", "value": 2},
    }


@pytest.fixture(scope="module")
def genomic_dup3_cx_37(genomic_dup3_37_loc):
    """Create test fixture copy number change variation"""
    return {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.fpx5MzUfd0Lll7kOwHxQ_znt5k-KJBau",
        "subject": genomic_dup3_37_loc,
        "copy_change": "efo:0030072",
    }


@pytest.fixture(scope="module")
def genomic_dup4_cn_38(genomic_dup4_loc):
    """Create test fixture copy number count variation"""
    return {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.Fvt2VSCSLNROTjZd20vKZUzdMYxz-cfZ",
        "subject": genomic_dup4_loc,
        "copies": {"type": "Number", "value": 3},
    }


@pytest.fixture(scope="module")
def genomic_dup4_cx_38(genomic_dup4_loc):
    """Create test fixture copy number change variation"""
    return {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.WlcvrOPHx67R2Qyw8O_s3MyHTJLX8uJD",
        "subject": genomic_dup4_loc,
        "copy_change": "efo:0030069",
    }


@pytest.fixture(scope="module")
def genomic_dup4_37_loc():
    """Create test fixture GRCh37 duplication subject"""
    return {
        "id": "ga4gh:SL.RUcHRUFPNBJhkwXz-WzUZLbT7XOlddAv",
        "sequence_id": "ga4gh:SQ.iy_UbUrvECxFRX5LPTH_KPojdlT7BKsf",
        "start": {"value": 29652251, "comparator": "<=", "type": "IndefiniteRange"},
        "end": {"value": 29981821, "comparator": ">=", "type": "IndefiniteRange"},
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_dup4_cn_37(genomic_dup4_37_loc):
    """Create test fixture copy number change variation"""
    return {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN._KLpF-U9nyddrZvZott7sQpJmdMRAh9E",
        "subject": genomic_dup4_37_loc,
        "copies": {"type": "Number", "value": 3},
    }


@pytest.fixture(scope="module")
def genomic_dup4_cx_37(genomic_dup4_37_loc):
    """Create test fixture copy number change variation"""
    return {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.q0kxeY8Nm8I9hLkrlOldryNy8Qv7prje",
        "subject": genomic_dup4_37_loc,
        "copy_change": "efo:0030069",
    }


@pytest.fixture(scope="module")
def genomic_dup5_cn_38(genomic_dup5_loc):
    """Create test fixture copy number count variation"""
    return {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.cegGbGJ5faMNTnlr6Oc4r93BOgklfftc",
        "subject": genomic_dup5_loc,
        "copies": {"type": "Number", "value": 4},
    }


@pytest.fixture(scope="module")
def genomic_dup5_cx_38(genomic_dup5_loc):
    """Create test fixture copy number change variation"""
    return {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.OetWeFHjI2_WQ_2FIdmxE0KWdkw22TYL",
        "subject": genomic_dup5_loc,
        "copy_change": "efo:0030067",
    }


@pytest.fixture(scope="module")
def genomic_dup5_37_loc():
    """Create test fixture GRCh37 duplication subject"""
    return {
        "id": "ga4gh:SL.1ELvCCFb5kq7_k47XzpatrJMgtvW74Co",
        "sequence_id": "ga4gh:SQ.v7noePfnNpK8ghYXEqZ9NukMXW7YeNsm",
        "start": {"value": 153287262, "comparator": "<=", "type": "IndefiniteRange"},
        "end": {"value": 153357667, "type": "Number"},
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_dup5_cn_37(genomic_dup5_37_loc):
    """Create test fixture copy number count variation"""
    return {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.wWg3b2U9qwdvS7ASQaG-6Gh6rhBIS5AA",
        "subject": genomic_dup5_37_loc,
        "copies": {"type": "Number", "value": 4},
    }


@pytest.fixture(scope="module")
def genomic_dup5_cx_37(genomic_dup5_37_loc):
    """Create test fixture copy number change variation"""
    return {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.scvCnG-ZfmUdGn3TBMMFysuv8xRq0VG1",
        "subject": genomic_dup5_37_loc,
        "copy_change": "efo:0030067",
    }


@pytest.fixture(scope="module")
def genomic_dup6_cn_38(genomic_dup6_loc):
    """Create test fixture copy number count variation"""
    return {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.Tuvj2C1xYVUF0UL3CSB1GZjCRC6K_eSt",
        "subject": genomic_dup6_loc,
        "copies": {"type": "Number", "value": 2},
    }


@pytest.fixture(scope="module")
def genomic_dup6_cx_38(genomic_dup6_loc):
    """Create test fixture copy number change variation"""
    return {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.zjYLwGA9ICJ1quFCO7rVvDHcwCCX9uNF",
        "subject": genomic_dup6_loc,
        "copy_change": "efo:0030064",
    }


@pytest.fixture(scope="module")
def genomic_dup6_37_loc():
    """Create test fixture GRCh37 duplication subject"""
    return {
        "id": "ga4gh:SL.GdeHSPT6RDTaMycAXoD_rqhMjumKGd-4",
        "sequence_id": "ga4gh:SQ.v7noePfnNpK8ghYXEqZ9NukMXW7YeNsm",
        "start": {"value": 153287262, "type": "Number"},
        "end": {"value": 153357667, "comparator": ">=", "type": "IndefiniteRange"},
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_dup6_cn_37(genomic_dup6_37_loc):
    """Create test fixture copy number count variation"""
    return {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.JvO8QAavKtd38hFB-wpzsm_9JCrFCM7i",
        "subject": genomic_dup6_37_loc,
        "copies": {"type": "Number", "value": 2},
    }


@pytest.fixture(scope="module")
def genomic_dup6_cx_37(genomic_dup6_37_loc):
    """Create test fixture copy number change variation"""
    return {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.2mfbKzavCrZA50piqh5tauTNWLkuZSCq",
        "subject": genomic_dup6_37_loc,
        "copy_change": "efo:0030064",
    }


@pytest.fixture(scope="module")
def genomic_del1_cx_38(genomic_del1_seq_loc):
    """Create test fixture copy number change variation"""
    return {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.p8k6YurGbJALwJRjvKfq156CRUetenB6",
        "subject": genomic_del1_seq_loc,
        "copy_change": "efo:0030064",
    }


@pytest.fixture(scope="module")
def genomic_del1_37_loc():
    """Create test fixture GRCh37 deletion subject"""
    return {
        "id": "ga4gh:SL.QtY9ACIT6TJAae_KDSAZKJfBLerBbFfq",
        "sequence_id": "ga4gh:SQ.VNBualIltAyi2AI_uXcKU7M9XUOuA7MS",
        "start": {"value": 10191494, "type": "Number"},
        "end": {"value": 10191495, "type": "Number"},
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_del1_cn_37(genomic_del1_37_loc):
    """Create test fixture copy number count variation"""
    return {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.rkgfXxRiPAPhf8b8M_7JZ_Sd_v--xTSp",
        "subject": genomic_del1_37_loc,
        "copies": {"type": "Number", "value": 1},
    }


@pytest.fixture(scope="module")
def genomic_del1_cx_37(genomic_del1_37_loc):
    """Create test fixture copy number change variation"""
    return {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.ES3st8kvfHmoaGjX3pCsNqnVunz8-prn",
        "subject": genomic_del1_37_loc,
        "copy_change": "efo:0030064",
    }


@pytest.fixture(scope="module")
def genomic_del2_cx_38(genomic_del2_seq_loc):
    """Create test fixture copy number change variation"""
    return {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.Y_94z-QmAZgbQfSRRFxBbaBe1GZ4TBCG",
        "subject": genomic_del2_seq_loc,
        "copy_change": "efo:0030071",
    }


@pytest.fixture(scope="module")
def genomic_del2_37_loc():
    """Create test fixture GRCh37 deletion subject"""
    return {
        "id": "ga4gh:SL.XXW-NlnipMp1hTv0xGdTY-RsjDUUZWht",
        "sequence_id": "ga4gh:SQ.VNBualIltAyi2AI_uXcKU7M9XUOuA7MS",
        "start": {"value": 10188278, "type": "Number"},
        "end": {"value": 10188297, "type": "Number"},
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_del2_cn_37(genomic_del2_37_loc):
    """Create test fixture copy number count variation"""
    return {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.pn37YpWXjDcIPLhuTUtAb9-pwFKrnM4M",
        "subject": genomic_del2_37_loc,
        "copies": {"type": "Number", "value": 1},
    }


@pytest.fixture(scope="module")
def genomic_del2_cx_37(genomic_del2_37_loc):
    """Create test fixture copy number change variation"""
    return {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.pWZIXKKaRGj8LwUJBcdEPnJKtMzLAUJ_",
        "subject": genomic_del2_37_loc,
        "copy_change": "efo:0030071",
    }


@pytest.fixture(scope="module")
def genomic_del3_cn_38(genomic_del3_dup3_loc_not_normalized):
    """Create test fixture copy number count variation"""
    return {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.a2-YlA01Q-y2zAOgn4r-cP9hiIEGaaqq",
        "subject": genomic_del3_dup3_loc_not_normalized,
        "copies": {"type": "Number", "value": 2},
    }


@pytest.fixture(scope="module")
def genomic_del3_cx_38(genomic_del3_dup3_loc_not_normalized):
    """Create test fixture copy number change variation"""
    return {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.Xffhy7i-ZkZtQSoMz6eytI6QQlShmfPo",
        "subject": genomic_del3_dup3_loc_not_normalized,
        "copy_change": "efo:0030069",
    }


@pytest.fixture(scope="module")
def genomic_del3_37_loc():
    """Create test fixture GRCh37 deletion subject"""
    return {
        "id": "ga4gh:SL.MmBA0qlC0J6yxWTyqhpsd4H_oErccpUD",
        "sequence_id": "ga4gh:SQ.v7noePfnNpK8ghYXEqZ9NukMXW7YeNsm",
        "start": {"min": 31078343, "max": 31118467, "type": "DefiniteRange"},
        "end": {"min": 33292396, "max": 33435269, "type": "DefiniteRange"},
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_del3_cn_37(genomic_del3_37_loc):
    """Create test fixture copy number count variation"""
    return {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.GOhhUsvbhFqaXgq3hNH3nuJldMIwVMpO",
        "subject": genomic_del3_37_loc,
        "copies": {"type": "Number", "value": 2},
    }


@pytest.fixture(scope="module")
def genomic_del3_cx_37(genomic_del3_37_loc):
    """Create test fixture copy number change variation"""
    return {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.5xJqmmaJTg0gv_BVa5gFHW_krgAGXviI",
        "subject": genomic_del3_37_loc,
        "copy_change": "efo:0030069",
    }


@pytest.fixture(scope="module")
def genomic_del4_cn_38(genomic_del4_seq_loc):
    """Create test fixture copy number count variation"""
    return {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.R0U1Bg1lhPYRVyJx_4LVE3JmR_5V7e8e",
        "subject": genomic_del4_seq_loc,
        "copies": {"type": "Number", "value": 4},
    }


@pytest.fixture(scope="module")
def genomic_del4_cx_38(genomic_del4_seq_loc):
    """Create test fixture copy number change variation"""
    return {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.fYi0dG6Q8kACkyY-ICBzzvslv-ONWrPF",
        "subject": genomic_del4_seq_loc,
        "copy_change": "efo:0030067",
    }


@pytest.fixture(scope="module")
def genomic_del4_37_loc():
    """Create test fixture GRCh37 deletion subject"""
    return {
        "id": "ga4gh:SL.ghgq_DzgZ93L1gVXbnHKEQwBzwnNdfSJ",
        "sequence_id": "ga4gh:SQ.v7noePfnNpK8ghYXEqZ9NukMXW7YeNsm",
        "start": {"value": 31138612, "comparator": "<=", "type": "IndefiniteRange"},
        "end": {"value": 33357594, "comparator": ">=", "type": "IndefiniteRange"},
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_del4_cn_37(genomic_del4_37_loc):
    """Create test fixture copy number count variation"""
    return {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.lI6I5Tf7rRjN4P_T1euZgAn0YucxrMj4",
        "subject": genomic_del4_37_loc,
        "copies": {"type": "Number", "value": 4},
    }


@pytest.fixture(scope="module")
def genomic_del4_cx_37(genomic_del4_37_loc):
    """Create test fixture copy number change variation"""
    return {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.PcA3_b1KEy4vXkhfUHr5lM5WeWBY4UQY",
        "subject": genomic_del4_37_loc,
        "copy_change": "efo:0030067",
    }


@pytest.fixture(scope="module")
def genomic_del5_cn_38(genomic_del5_seq_loc):
    """Create test fixture copy number count variation"""
    return {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.UhXP8ENfA3T-FhNvKSqu-9UwrL5GZFvF",
        "subject": genomic_del5_seq_loc,
        "copies": {"type": "Number", "value": 2},
    }


@pytest.fixture(scope="module")
def genomic_del5_cx_38(genomic_del5_seq_loc):
    """Create test fixture copy number change variation"""
    return {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.BpxaMFwIjonSNZaS8hYdgeGngz4bOlkV",
        "subject": genomic_del5_seq_loc,
        "copy_change": "efo:0030064",
    }


@pytest.fixture(scope="module")
def genomic_del5_37_loc():
    """Create test fixture GRCh37 deletion subject"""
    return {
        "id": "ga4gh:SL.kJXgk4nQIwS53UN5gxgIFIPcNSt53n8O",
        "sequence_id": "ga4gh:SQ.v7noePfnNpK8ghYXEqZ9NukMXW7YeNsm",
        "start": {"value": 18593473, "comparator": "<=", "type": "IndefiniteRange"},
        "end": {"value": 18671749, "type": "Number"},
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_del5_cn_37(genomic_del5_37_loc):
    """Create test fixture copy number count variation"""
    return {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.3Yu1OOpMjnzt6yrt0OLjKWnhD2YD9pkg",
        "subject": genomic_del5_37_loc,
        "copies": {"type": "Number", "value": 2},
    }


@pytest.fixture(scope="module")
def genomic_del5_cx_37(genomic_del5_37_loc):
    """Create test fixture copy number change variation"""
    return {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.gEe40vGaYb0D0rdWwtmX2ODJepkSJv15",
        "subject": genomic_del5_37_loc,
        "copy_change": "efo:0030064",
    }


@pytest.fixture(scope="module")
def genomic_del6_cn_38(genomic_del6_seq_loc):
    """Create test fixture copy number count variation"""
    return {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.lfI0nj_cOjRwul_R_4tohZK_vqiZ-LSz",
        "subject": genomic_del6_seq_loc,
        "copies": {"type": "Number", "value": 1},
    }


@pytest.fixture(scope="module")
def genomic_del6_cx_38(genomic_del6_seq_loc):
    """Create test fixture copy number change variation"""
    return {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.HbonS2GEi_jBQA_2gvaiwONy_KHqN_G0",
        "subject": genomic_del6_seq_loc,
        "copy_change": "efo:0030071",
    }


@pytest.fixture(scope="module")
def genomic_del6_37_loc():
    """Create test fixture GRCh37 deletion subject"""
    return {
        "id": "ga4gh:SL.LU0WrrgeZ-HIG7f0LF_eS0yGvhDtiHuT",
        "sequence_id": "ga4gh:SQ.KqaUhJMW3CDjhoVtBetdEKT1n6hM-7Ek",
        "start": {"value": 133783901, "type": "Number"},
        "end": {"value": 133785996, "comparator": ">=", "type": "IndefiniteRange"},
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_del6_cn_37(genomic_del6_37_loc):
    """Create test fixture copy number count variation"""
    return {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.UclTlSHjepPFKjkBdZjTPquoStwlr24H",
        "subject": genomic_del6_37_loc,
        "copies": {"type": "Number", "value": 1},
    }


@pytest.fixture(scope="module")
def genomic_del6_cx_37(genomic_del6_37_loc):
    """Create test fixture copy number change variation"""
    return {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.j6BFcKCBMSfz3Vh-ekevUx8GqLCGJCU6",
        "subject": genomic_del6_37_loc,
        "copy_change": "efo:0030071",
    }


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
    assert resp.copy_number_count.dict(by_alias=True) == genomic_dup1_38_cn

    q = "NC_000003.11:g.49568695dup"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=2, do_liftover=False
    )
    assert resp.copy_number_count.dict(by_alias=True) == genomic_dup1_cn_37

    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=2, do_liftover=True
    )
    assert resp.copy_number_count.dict(by_alias=True) == genomic_dup1_38_cn

    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=1, do_liftover=True
    )
    expected = copy.deepcopy(genomic_dup1_38_cn)
    expected["copies"]["value"] = 2
    expected["id"] = "ga4gh:CN.AjWJMk7pirTUVxS0yiSnhnfLg_ql79X1"
    assert resp.copy_number_count.dict(by_alias=True) == expected


@pytest.mark.asyncio
async def test_genomic_dup1_copy_number_change(
    test_cnv_handler, genomic_dup1_cx_38, genomic_dup1_cx_37
):
    """Test that genomic duplication works correctly"""
    q = "NC_000003.12:g.49531262dup"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030069", do_liftover=False
    )
    assert resp.copy_number_change.dict(by_alias=True) == genomic_dup1_cx_38

    q = "NC_000003.11:g.49568695dup"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030069", do_liftover=False
    )
    assert resp.copy_number_change.dict(by_alias=True) == genomic_dup1_cx_37

    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030069", do_liftover=True
    )
    assert resp.copy_number_change.dict(by_alias=True) == genomic_dup1_cx_38


@pytest.mark.asyncio
async def test_genomic_dup2_copy_number_count(
    test_cnv_handler, genomic_dup2_38_cn, genomic_dup2_cn_37
):
    """Test that genomic duplication works correctly"""
    q = "NC_000023.11:g.33211290_33211293dup"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=2, do_liftover=False
    )
    assert resp.copy_number_count.dict(by_alias=True) == genomic_dup2_38_cn

    q = "NC_000023.10:g.33229407_33229410dup"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=2, do_liftover=False
    )
    assert resp.copy_number_count.dict(by_alias=True) == genomic_dup2_cn_37

    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=2, do_liftover=True
    )
    assert resp.copy_number_count.dict(by_alias=True) == genomic_dup2_38_cn

    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=1, do_liftover=True
    )
    expected = copy.deepcopy(genomic_dup2_38_cn)
    expected["copies"]["value"] = 2
    expected["id"] = "ga4gh:CN.-MMzWFQXvQJg1_b3wwCyIre2KhN2Nbso"
    assert resp.copy_number_count.dict(by_alias=True) == expected


@pytest.mark.asyncio
async def test_genomic_dup2_copy_number_change(
    test_cnv_handler, genomic_dup2_cx_38, genomic_dup2_cx_37
):
    """Test that genomic duplication works correctly"""
    q = "NC_000023.11:g.33211290_33211293dup"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030067", do_liftover=False
    )
    assert resp.copy_number_change.dict(by_alias=True) == genomic_dup2_cx_38

    q = "NC_000023.10:g.33229407_33229410dup"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030067", do_liftover=False
    )
    assert resp.copy_number_change.dict(by_alias=True) == genomic_dup2_cx_37

    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030067", do_liftover=True
    )
    assert resp.copy_number_change.dict(by_alias=True) == genomic_dup2_cx_38


@pytest.mark.asyncio
async def test_genomic_dup3_copy_number_count(
    test_cnv_handler, genomic_dup3_cn_38, genomic_dup3_cn_37
):
    """Test that genomic duplication works correctly"""
    q = "NC_000023.11:g.(31060227_31100351)_(33274278_33417151)dup"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=1, do_liftover=False
    )
    assert resp.copy_number_count.dict(by_alias=True) == genomic_dup3_cn_38

    q = "NC_000023.10:g.(31078344_31118468)_(33292395_33435268)dup"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=1, do_liftover=False
    )
    assert resp.copy_number_count.dict(by_alias=True) == genomic_dup3_cn_37

    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=1, do_liftover=True
    )
    assert resp.copy_number_count.dict(by_alias=True) == genomic_dup3_cn_38

    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=2, do_liftover=True
    )
    expected = copy.deepcopy(genomic_dup3_cn_38)
    expected["copies"] = {"value": 3, "type": "Number"}
    expected["id"] = "ga4gh:CN.-4AcyJ5imRkXzgeULrzWVTraqAe4YRYc"
    assert resp.copy_number_count.dict(by_alias=True) == expected


@pytest.mark.asyncio
async def test_genomic_dup3_copy_number_change(
    test_cnv_handler, genomic_dup3_cx_38, genomic_dup3_cx_37
):
    """Test that genomic duplication works correctly"""
    q = "NC_000023.11:g.(31060227_31100351)_(33274278_33417151)dup"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030072", do_liftover=False
    )
    assert resp.copy_number_change.dict(by_alias=True) == genomic_dup3_cx_38

    q = "NC_000023.10:g.(31078344_31118468)_(33292395_33435268)dup"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030072", do_liftover=False
    )
    assert resp.copy_number_change.dict(by_alias=True) == genomic_dup3_cx_37

    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030072", do_liftover=True
    )
    assert resp.copy_number_change.dict(by_alias=True) == genomic_dup3_cx_38


@pytest.mark.asyncio
async def test_genomic_dup4_copy_number_count(
    test_cnv_handler, genomic_dup4_cn_38, genomic_dup4_cn_37
):
    """Test that genomic duplication works correctly"""
    q = "NC_000020.11:g.(?_30417576)_(31394018_?)dup"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=2, do_liftover=False
    )
    assert resp.copy_number_count.dict(by_alias=True) == genomic_dup4_cn_38

    q = "NC_000020.10:g.(?_29652252)_(29981821_?)dup"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=2, do_liftover=False
    )
    assert resp.copy_number_count.dict(by_alias=True) == genomic_dup4_cn_37

    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=2, do_liftover=True
    )
    assert resp.copy_number_count.dict(by_alias=True) == genomic_dup4_cn_38

    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=1, do_liftover=True
    )
    expected = copy.deepcopy(genomic_dup4_cn_38)
    expected["copies"]["value"] = 2
    expected["id"] = "ga4gh:CN.ObqoAUmFaGiMknzr8GurBYuoBNeE2HtV"
    assert resp.copy_number_count.dict(by_alias=True) == expected


@pytest.mark.asyncio
async def test_genomic_dup4_copy_number_change(
    test_cnv_handler, genomic_dup4_cx_38, genomic_dup4_cx_37
):
    """Test that genomic duplication works correctly"""
    q = "NC_000020.11:g.(?_30417576)_(31394018_?)dup"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030069", do_liftover=False
    )
    assert resp.copy_number_change.dict(by_alias=True) == genomic_dup4_cx_38

    q = "NC_000020.10:g.(?_29652252)_(29981821_?)dup"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030069", do_liftover=False
    )
    assert resp.copy_number_change.dict(by_alias=True) == genomic_dup4_cx_37

    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030069", do_liftover=True
    )
    assert resp.copy_number_change.dict(by_alias=True) == genomic_dup4_cx_38


@pytest.mark.asyncio
async def test_genomic_dup5_copy_number_count(
    test_cnv_handler, genomic_dup5_cn_38, genomic_dup5_cn_37
):
    """Test that genomic duplication works correctly"""
    q = "NC_000023.11:g.(?_154021812)_154092209dup"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=3, do_liftover=False
    )
    assert resp.copy_number_count.dict(by_alias=True) == genomic_dup5_cn_38

    q = "NC_000023.10:g.(?_153287263)_153357667dup"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=3, do_liftover=False
    )
    assert resp.copy_number_count.dict(by_alias=True) == genomic_dup5_cn_37

    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=3, do_liftover=True
    )
    assert resp.copy_number_count.dict(by_alias=True) == genomic_dup5_cn_38

    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=4, do_liftover=True
    )
    expected = copy.deepcopy(genomic_dup5_cn_38)
    expected["copies"] = {"value": 5, "type": "Number"}
    expected["id"] = "ga4gh:CN.KxU8n8a8UXCQnnwHDX9daUXd0zdxtFtB"
    assert resp.copy_number_count.dict(by_alias=True) == expected


@pytest.mark.asyncio
async def test_genomic_dup5_copy_number_change(
    test_cnv_handler, genomic_dup5_cx_38, genomic_dup5_cx_37
):
    """Test that genomic duplication works correctly"""
    q = "NC_000023.11:g.(?_154021812)_154092209dup"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030067", do_liftover=False
    )
    assert resp.copy_number_change.dict(by_alias=True) == genomic_dup5_cx_38

    q = "NC_000023.10:g.(?_153287263)_153357667dup"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030067", do_liftover=False
    )
    assert resp.copy_number_change.dict(by_alias=True) == genomic_dup5_cx_37

    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030067", do_liftover=True
    )
    assert resp.copy_number_change.dict(by_alias=True) == genomic_dup5_cx_38


@pytest.mark.asyncio
async def test_genomic_dup6_copy_number_count(
    test_cnv_handler, genomic_dup6_cn_38, genomic_dup6_cn_37
):
    """Test that genomic duplication works correctly"""
    q = "NC_000023.11:g.154021812_(154092209_?)dup"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=1, do_liftover=False
    )
    assert resp.copy_number_count.dict(by_alias=True) == genomic_dup6_cn_38

    q = "NC_000023.10:g.153287263_(153357667_?)dup"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=1, do_liftover=False
    )
    assert resp.copy_number_count.dict(by_alias=True) == genomic_dup6_cn_37

    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=1, do_liftover=True
    )
    assert resp.copy_number_count.dict(by_alias=True) == genomic_dup6_cn_38

    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=2, do_liftover=True
    )
    expected = copy.deepcopy(genomic_dup6_cn_38)
    expected["copies"] = {"value": 3, "type": "Number"}
    expected["id"] = "ga4gh:CN.TQ36aHJOnfKG3a8hpVMR1POBTOVBgYcS"
    assert resp.copy_number_count.dict(by_alias=True) == expected


@pytest.mark.asyncio
async def test_genomic_dup6_copy_number_change(
    test_cnv_handler, genomic_dup6_cx_38, genomic_dup6_cx_37
):
    """Test that genomic duplication works correctly"""
    q = "NC_000023.11:g.154021812_(154092209_?)dup"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030064", do_liftover=False
    )
    assert resp.copy_number_change.dict(by_alias=True) == genomic_dup6_cx_38

    q = "NC_000023.10:g.153287263_(153357667_?)dup"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030064", do_liftover=False
    )
    assert resp.copy_number_change.dict(by_alias=True) == genomic_dup6_cx_37

    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030064", do_liftover=True
    )
    assert resp.copy_number_change.dict(by_alias=True) == genomic_dup6_cx_38


@pytest.mark.asyncio
async def test_genomic_del1_copy_number_count(
    test_cnv_handler, genomic_del1_38_cn, genomic_del1_cn_37
):
    """Test that genomic deletion works correctly"""
    q = "NC_000003.12:g.10149811del"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=2, do_liftover=False
    )
    assert resp.copy_number_count.dict(by_alias=True) == genomic_del1_38_cn

    q = "NC_000003.11:g.10191495del"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=2, do_liftover=False
    )
    assert resp.copy_number_count.dict(by_alias=True) == genomic_del1_cn_37

    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=2, do_liftover=True
    )
    assert resp.copy_number_count.dict(by_alias=True) == genomic_del1_38_cn

    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=3, do_liftover=True
    )
    expected = copy.deepcopy(genomic_del1_38_cn)
    expected["copies"]["value"] = 2
    expected["id"] = "ga4gh:CN.WwgB7O1Ztcn6QqN2hjfjHucKOaWXic8n"
    assert resp.copy_number_count.dict(by_alias=True) == expected


@pytest.mark.asyncio
async def test_genomic_del1_copy_number_change(
    test_cnv_handler, genomic_del1_cx_38, genomic_del1_cx_37
):
    """Test that genomic deletion works correctly"""
    q = "NC_000003.12:g.10149811del"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030064", do_liftover=False
    )
    assert resp.copy_number_change.dict(by_alias=True) == genomic_del1_cx_38

    q = "NC_000003.11:g.10191495del"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030064", do_liftover=False
    )
    assert resp.copy_number_change.dict(by_alias=True) == genomic_del1_cx_37

    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030064", do_liftover=True
    )
    assert resp.copy_number_change.dict(by_alias=True) == genomic_del1_cx_38


@pytest.mark.asyncio
async def test_genomic_del2_copy_number_count(
    test_cnv_handler, genomic_del2_38_cn, genomic_del2_cn_37
):
    """Test that genomic deletion works correctly"""
    q = "NC_000003.12:g.10146595_10146613del"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=2, do_liftover=False
    )
    assert resp.copy_number_count.dict(by_alias=True) == genomic_del2_38_cn

    q = "NC_000003.11:g.10188279_10188297del"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=2, do_liftover=False
    )
    assert resp.copy_number_count.dict(by_alias=True) == genomic_del2_cn_37

    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=2, do_liftover=True
    )
    assert resp.copy_number_count.dict(by_alias=True) == genomic_del2_38_cn

    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=4, do_liftover=True
    )
    expected = copy.deepcopy(genomic_del2_38_cn)
    expected["copies"]["value"] = 3
    expected["id"] = "ga4gh:CN.GiZEIUkpM1Yl_O9v0Xz3X9Uza1BzdYwn"
    assert resp.copy_number_count.dict(by_alias=True) == expected


@pytest.mark.asyncio
async def test_genomic_del2_copy_number_change(
    test_cnv_handler, genomic_del2_cx_38, genomic_del2_cx_37
):
    """Test that genomic deletion works correctly"""
    q = "NC_000003.12:g.10146595_10146613del"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030071", do_liftover=False
    )
    assert resp.copy_number_change.dict(by_alias=True) == genomic_del2_cx_38

    q = "NC_000003.11:g.10188279_10188297del"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030071", do_liftover=False
    )
    assert resp.copy_number_change.dict(by_alias=True) == genomic_del2_cx_37

    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030071", do_liftover=True
    )
    assert resp.copy_number_change.dict(by_alias=True) == genomic_del2_cx_38


@pytest.mark.asyncio
async def test_genomic_del3_copy_number_count(
    test_cnv_handler, genomic_del3_cn_38, genomic_del3_cn_37
):
    """Test that genomic deletion works correctly"""
    q = "NC_000023.11:g.(31060227_31100351)_(33274278_33417151)del"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=3, do_liftover=False
    )
    assert resp.copy_number_count.dict(by_alias=True) == genomic_del3_cn_38

    q = "NC_000023.10:g.(31078344_31118468)_(33292395_33435268)del"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=3, do_liftover=False
    )
    assert resp.copy_number_count.dict(by_alias=True) == genomic_del3_cn_37

    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=3, do_liftover=True
    )
    assert resp.copy_number_count.dict(by_alias=True) == genomic_del3_cn_38

    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=2, do_liftover=True
    )
    expected = copy.deepcopy(genomic_del3_cn_38)
    expected["copies"] = {"value": 1, "type": "Number"}
    expected["id"] = "ga4gh:CN.a7prdJg9CK7KAFug-0VMjZC7xDxiSByV"
    assert resp.copy_number_count.dict(by_alias=True) == expected


@pytest.mark.asyncio
async def test_genomic_del3_copy_number_change(
    test_cnv_handler, genomic_del3_cx_38, genomic_del3_cx_37
):
    """Test that genomic deletion works correctly"""
    q = "NC_000023.11:g.(31060227_31100351)_(33274278_33417151)del"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030069", do_liftover=False
    )
    assert resp.copy_number_change.dict(by_alias=True) == genomic_del3_cx_38

    q = "NC_000023.10:g.(31078344_31118468)_(33292395_33435268)del"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030069", do_liftover=False
    )
    assert resp.copy_number_change.dict(by_alias=True) == genomic_del3_cx_37

    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030069", do_liftover=True
    )
    assert resp.copy_number_change.dict(by_alias=True) == genomic_del3_cx_38


@pytest.mark.asyncio
async def test_genomic_del4_copy_number_count(
    test_cnv_handler, genomic_del4_cn_38, genomic_del4_cn_37
):
    """Test that genomic deletion works correctly"""
    q = "NC_000023.11:g.(?_31120496)_(33339477_?)del"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=5, do_liftover=False
    )
    assert resp.copy_number_count.dict(by_alias=True) == genomic_del4_cn_38

    q = "NC_000023.10:g.(?_31138613)_(33357594_?)del"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=5, do_liftover=False
    )
    assert resp.copy_number_count.dict(by_alias=True) == genomic_del4_cn_37

    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=5, do_liftover=True
    )
    assert resp.copy_number_count.dict(by_alias=True) == genomic_del4_cn_38

    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=3, do_liftover=True
    )
    expected = copy.deepcopy(genomic_del4_cn_38)
    expected["copies"] = {"value": 2, "type": "Number"}
    expected["id"] = "ga4gh:CN.PlqY0WzSBD2KV2mD0L9QkyOGf1z2zZLu"
    assert resp.copy_number_count.dict(by_alias=True) == expected


@pytest.mark.asyncio
async def test_genomic_del4_copy_number_change(
    test_cnv_handler, genomic_del4_cx_38, genomic_del4_cx_37
):
    """Test that genomic deletion works correctly"""
    q = "NC_000023.11:g.(?_31120496)_(33339477_?)del"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030067", do_liftover=False
    )
    assert resp.copy_number_change.dict(by_alias=True) == genomic_del4_cx_38

    q = "NC_000023.10:g.(?_31138613)_(33357594_?)del"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030067", do_liftover=False
    )
    assert resp.copy_number_change.dict(by_alias=True) == genomic_del4_cx_37

    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030067", do_liftover=True
    )
    assert resp.copy_number_change.dict(by_alias=True) == genomic_del4_cx_38


@pytest.mark.asyncio
async def test_genomic_del5_copy_number_count(
    test_cnv_handler, genomic_del5_cn_38, genomic_del5_cn_37
):
    """Test that genomic deletion works correctly"""
    q = "NC_000023.11:g.(?_18575354)_18653629del"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=3, do_liftover=False
    )
    assert resp.copy_number_count.dict(by_alias=True) == genomic_del5_cn_38

    q = "NC_000023.10:g.(?_18593474)_18671749del"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=3, do_liftover=False
    )
    assert resp.copy_number_count.dict(by_alias=True) == genomic_del5_cn_37

    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=3, do_liftover=True
    )
    assert resp.copy_number_count.dict(by_alias=True) == genomic_del5_cn_38

    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=2, do_liftover=True
    )
    expected = copy.deepcopy(genomic_del5_cn_38)
    expected["copies"] = {"value": 1, "type": "Number"}
    expected["id"] = "ga4gh:CN.HzpGrUmI6mfMtPGC0H8ElVcNc0-V1ldA"
    assert resp.copy_number_count.dict(by_alias=True) == expected


@pytest.mark.asyncio
async def test_genomic_del5_copy_number_change(
    test_cnv_handler, genomic_del5_cx_38, genomic_del5_cx_37
):
    """Test that genomic deletion works correctly"""
    q = "NC_000023.11:g.(?_18575354)_18653629del"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030064", do_liftover=False
    )
    assert resp.copy_number_change.dict(by_alias=True) == genomic_del5_cx_38

    q = "NC_000023.10:g.(?_18593474)_18671749del"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030064", do_liftover=False
    )
    assert resp.copy_number_change.dict(by_alias=True) == genomic_del5_cx_37

    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030064", do_liftover=True
    )
    assert resp.copy_number_change.dict(by_alias=True) == genomic_del5_cx_38


@pytest.mark.asyncio
async def test_genomic_del6_copy_number_count(
    test_cnv_handler, genomic_del6_cn_38, genomic_del6_cn_37
):
    """Test that genomic deletion works correctly"""
    q = "NC_000006.12:g.133462764_(133464858_?)del"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=2, do_liftover=False
    )
    assert resp.copy_number_count.dict(by_alias=True) == genomic_del6_cn_38

    q = "NC_000006.11:g.133783902_(133785996_?)del"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=2, do_liftover=False
    )
    assert resp.copy_number_count.dict(by_alias=True) == genomic_del6_cn_37

    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=2, do_liftover=True
    )
    assert resp.copy_number_count.dict(by_alias=True) == genomic_del6_cn_38

    resp = await test_cnv_handler.hgvs_to_copy_number_count(
        q, baseline_copies=3, do_liftover=True
    )
    expected = copy.deepcopy(genomic_del6_cn_38)
    expected["copies"]["value"] = 2
    expected["id"] = "ga4gh:CN.NWXsPHBTkOaY0LgfZK6bWFoJ3lvf3fOH"
    assert resp.copy_number_count.dict(by_alias=True) == expected


@pytest.mark.asyncio
async def test_genomic_del6_copy_number_change(
    test_cnv_handler, genomic_del6_cx_38, genomic_del6_cx_37
):
    """Test that genomic deletion works correctly"""
    q = "NC_000006.12:g.133462764_(133464858_?)del"  # 38
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030071", do_liftover=False
    )
    assert resp.copy_number_change.dict(by_alias=True) == genomic_del6_cx_38

    q = "NC_000006.11:g.133783902_(133785996_?)del"  # 37
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030071", do_liftover=False
    )
    assert resp.copy_number_change.dict(by_alias=True) == genomic_del6_cx_37

    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030071", do_liftover=True
    )
    assert resp.copy_number_change.dict(by_alias=True) == genomic_del6_cx_38


@pytest.mark.asyncio
async def test_invalid_cnv(test_cnv_handler):
    """Check that invalid input return warnings"""
    q = "DAG1 g.49568695dup"
    resp = await test_cnv_handler.hgvs_to_copy_number_change(
        q, copy_change="efo:0030071", do_liftover=True, untranslatable_returns_text=True
    )
    assert set(resp.warnings) == {
        "DAG1 g.49568695dup is not a supported HGVS genomic duplication or deletion"
    }
    assert resp.copy_number_change.type == "Text"

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
