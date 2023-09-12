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
        "id": "ga4gh:CX.d8BWSLNKN0K4n8ySG0jWPCr4cJIqEf5g",
        "location": genomic_dup1_seq_loc_not_normalized,
        "copyChange": "efo:0030069",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_dup1_37_loc():
    """Create test fixture GRCh37 duplication location"""
    return {
        "id": "ga4gh:SL.ntKfbY4eZVFNOAMuZPb4RBRhINxvOmM9",
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
        "id": "ga4gh:CN.-TsNBiFhHFcWLb7pz9GCWcJunlBmb_B4",
        "location": genomic_dup1_37_loc,
        "copies": 3,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_dup1_cx_37(genomic_dup1_37_loc):
    """Create test fixture copy number change variation"""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.inPVJ3ANsN-Xb22HakFf_BmMg73gZiFo",
        "location": genomic_dup1_37_loc,
        "copyChange": "efo:0030069",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_dup2_cx_38(genomic_dup2_seq_loc_normalized):
    """Create test fixture copy number change variation"""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.uchcAvP6DXLirT7zASWHcWwmGWPLc8ye",
        "location": genomic_dup2_seq_loc_normalized,
        "copyChange": "efo:0030067",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_dup2_37_loc():
    """Create test fixture GRCh37 duplication location"""
    return {
        "id": "ga4gh:SL.g5_YcFgvTQSCuirJLMviwlue4NTb9EJ-",
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
        "id": "ga4gh:CN.jQdYqHKs7Y7BPQq5lBIddKd208iJrskw",
        "location": genomic_dup2_37_loc,
        "copies": 3,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_dup2_cx_37(genomic_dup2_37_loc):
    """Create test fixture copy number change variation"""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.6Nzma3vnsYb7ZqJ0XSiqv2F8-XFDKMyP",
        "location": genomic_dup2_37_loc,
        "copyChange": "efo:0030067",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_dup3_cn_38(genomic_del3_dup3_loc_not_normalized):
    """Create test fixture copy number count variation"""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.Ka-Wsibx4GHmHgurHCdk0W4deqZt26y4",
        "location": genomic_del3_dup3_loc_not_normalized,
        "copies": 2,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_dup3_cx_38(genomic_del3_dup3_loc_not_normalized):
    """Create test fixture copy number change variation"""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.VDYux3JxzXmKQU5Ll4Knt6Y4RW5UMUbJ",
        "location": genomic_del3_dup3_loc_not_normalized,
        "copyChange": "efo:0030072",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_dup3_37_loc():
    """Create test fixture GRCh37 duplication location"""
    return {
        "id": "ga4gh:SL.8aHcfpqhpWQ0kigsdn9mY61Alod63tDP",
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
        "id": "ga4gh:CN.gOq-IVcswisgIfsx9PpZ-BLxzcMVeJeR",
        "location": genomic_dup3_37_loc,
        "copies": 2,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_dup3_cx_37(genomic_dup3_37_loc):
    """Create test fixture copy number change variation"""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.8UZ5UZWGhsRU293nzuMfIpqG_6NBzHbs",
        "location": genomic_dup3_37_loc,
        "copyChange": "efo:0030072",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_dup4_cn_38(genomic_dup4_loc):
    """Create test fixture copy number count variation"""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.KqbQewUgZYfmottbgn1xYq58DiPVU5SZ",
        "location": genomic_dup4_loc,
        "copies": 3,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_dup4_cx_38(genomic_dup4_loc):
    """Create test fixture copy number change variation"""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.X0nbb6qzkjqisoT8Ls_7-OS9KkYfTCCu",
        "location": genomic_dup4_loc,
        "copyChange": "efo:0030069",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_dup4_37_loc():
    """Create test fixture GRCh37 duplication location"""
    return {
        "id": "ga4gh:SL.WVXqGHNVaD96semkKQVfTuEo4TN-yMGE",
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
        "id": "ga4gh:CN.m19rHAtbxPm8ojn1pZm7Pq7e00WJkE75",
        "location": genomic_dup4_37_loc,
        "copies": 3,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_dup4_cx_37(genomic_dup4_37_loc):
    """Create test fixture copy number change variation"""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.5cVAesP5EKfh3uOTXe-BWyHdsFEvW7a0",
        "location": genomic_dup4_37_loc,
        "copyChange": "efo:0030069",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_dup5_cn_38(genomic_dup5_loc):
    """Create test fixture copy number count variation"""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.hTq6Z10Mg7rS_T_uoODhF24HmvRdEkLk",
        "location": genomic_dup5_loc,
        "copies": 4,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_dup5_cx_38(genomic_dup5_loc):
    """Create test fixture copy number change variation"""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.bR_i_ztx0rLrEbML1B_xoqJ50W39qlKN",
        "location": genomic_dup5_loc,
        "copyChange": "efo:0030067",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_dup5_37_loc():
    """Create test fixture GRCh37 duplication location"""
    return {
        "id": "ga4gh:SL.hS8rIdNliH9F4YuLlHhxLIdrHVVrCEXz",
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
        "id": "ga4gh:CN.QdhteI3bI3qRFEZOk9pmICPsKFI6MCXd",
        "location": genomic_dup5_37_loc,
        "copies": 4,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_dup5_cx_37(genomic_dup5_37_loc):
    """Create test fixture copy number change variation"""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.PJWAoRYtFYcewmPtM5W2eG6TSef7gAyJ",
        "location": genomic_dup5_37_loc,
        "copyChange": "efo:0030067",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_dup6_cn_38(genomic_dup6_loc):
    """Create test fixture copy number count variation"""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.KSFn5KQIPuPVJ6FjWaF0vzl7eRwwHbX9",
        "location": genomic_dup6_loc,
        "copies": 2,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_dup6_cx_38(genomic_dup6_loc):
    """Create test fixture copy number change variation"""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.5irbnyk3aqsHCEbs2WAgiNpChn3asZgQ",
        "location": genomic_dup6_loc,
        "copyChange": "efo:0030064",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_dup6_37_loc():
    """Create test fixture GRCh37 duplication location"""
    return {
        "id": "ga4gh:SL.xoQ6AO6YVSUdlvcO5WtpXGJykEbMktY6",
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
        "id": "ga4gh:CN.SijB_fXpFimGmPnWe2YHNcVv4NyEp9Uo",
        "location": genomic_dup6_37_loc,
        "copies": 2,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_dup6_cx_37(genomic_dup6_37_loc):
    """Create test fixture copy number change variation"""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.IkUmxq_kjLmo2mPv4eg9fgcd3B12tT4G",
        "location": genomic_dup6_37_loc,
        "copyChange": "efo:0030064",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_del1_cx_38(genomic_del1_seq_loc):
    """Create test fixture copy number change variation"""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.LWRBNtBgcETMXEKezrr7WUPjO9WoOaqL",
        "location": genomic_del1_seq_loc,
        "copyChange": "efo:0030064",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_del1_37_loc():
    """Create test fixture GRCh37 deletion location"""
    return {
        "id": "ga4gh:SL.zz9AUA6ANv3OPftr3dI-7GxGeaZKeADW",
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
        "id": "ga4gh:CN._n2CTBQdDMSObUf6DPdj5t7XLPAQ-Ojt",
        "location": genomic_del1_37_loc,
        "copies": 1,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_del1_cx_37(genomic_del1_37_loc):
    """Create test fixture copy number change variation"""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.KqiBf1gspLI7WQu-wVOFYrG1HhLyqFlV",
        "location": genomic_del1_37_loc,
        "copyChange": "efo:0030064",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_del2_cx_38(genomic_del2_seq_loc):
    """Create test fixture copy number change variation"""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.ehFxhcx6VUozuQxTd4vskjpaUIptU0Qg",
        "location": genomic_del2_seq_loc,
        "copyChange": "efo:0030071",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_del2_37_loc():
    """Create test fixture GRCh37 deletion location"""
    return {
        "id": "ga4gh:SL.BwPpPS55EYrek4Gs7xdd8sNjYJvcZRi7",
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
        "id": "ga4gh:CN.ic8DuXVLz-cFut46rneGHbkACsEQg832",
        "location": genomic_del2_37_loc,
        "copies": 1,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_del2_cx_37(genomic_del2_37_loc):
    """Create test fixture copy number change variation"""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.0i-lR3bpBiuYAGC4D7UErkrrgrawJahz",
        "location": genomic_del2_37_loc,
        "copyChange": "efo:0030071",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_del3_cn_38(genomic_del3_dup3_loc_not_normalized):
    """Create test fixture copy number count variation"""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.Ka-Wsibx4GHmHgurHCdk0W4deqZt26y4",
        "location": genomic_del3_dup3_loc_not_normalized,
        "copies": 2,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_del3_cx_38(genomic_del3_dup3_loc_not_normalized):
    """Create test fixture copy number change variation"""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.JQvddf8xWhkCjdux8NLWNs4MdEMqDiLX",
        "location": genomic_del3_dup3_loc_not_normalized,
        "copyChange": "efo:0030069",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_del3_37_loc():
    """Create test fixture GRCh37 deletion location"""
    return {
        "id": "ga4gh:SL.8aHcfpqhpWQ0kigsdn9mY61Alod63tDP",
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
        "id": "ga4gh:CN.gOq-IVcswisgIfsx9PpZ-BLxzcMVeJeR",
        "location": genomic_del3_37_loc,
        "copies": 2,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_del3_cx_37(genomic_del3_37_loc):
    """Create test fixture copy number change variation"""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.-s0MQI48uGQ97FAGOVDeFgGAAL5J4wak",
        "location": genomic_del3_37_loc,
        "copyChange": "efo:0030069",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_del4_cn_38(genomic_del4_seq_loc):
    """Create test fixture copy number count variation"""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.2k3RnvQtBIIIqSqVlTT7Uh0KeuD8LIpK",
        "location": genomic_del4_seq_loc,
        "copies": 4,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_del4_cx_38(genomic_del4_seq_loc):
    """Create test fixture copy number change variation"""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.1DiUzraiKZLJb8oF8ynARS816fthsJpV",
        "location": genomic_del4_seq_loc,
        "copyChange": "efo:0030067",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_del4_37_loc():
    """Create test fixture GRCh37 deletion location"""
    return {
        "id": "ga4gh:SL.2ORImLGRcezhWXDYI9sZvLwFGPDH9WKS",
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
        "id": "ga4gh:CN.Mv42RTfmPVQV4nSoytbDbK4cfsUO0wyD",
        "location": genomic_del4_37_loc,
        "copies": 4,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_del4_cx_37(genomic_del4_37_loc):
    """Create test fixture copy number change variation"""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.CJdXMZXSEE4hCIwjGxke4EWY7lMENYPj",
        "location": genomic_del4_37_loc,
        "copyChange": "efo:0030067",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_del5_cn_38(genomic_del5_seq_loc):
    """Create test fixture copy number count variation"""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.VIYxNpNj8c5ulSfr8mFkxTBCGMunoUIN",
        "location": genomic_del5_seq_loc,
        "copies": 2,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_del5_cx_38(genomic_del5_seq_loc):
    """Create test fixture copy number change variation"""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.w-H5jcMK0AijJomgeu_LhN5IQ_0Z0GIi",
        "location": genomic_del5_seq_loc,
        "copyChange": "efo:0030064",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_del5_37_loc():
    """Create test fixture GRCh37 deletion location"""
    return {
        "id": "ga4gh:SL.xYPrcJoKhU-xWyeBB4-DzYmFYOEynMBy",
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
        "id": "ga4gh:CN.9qwSqDhzXCHPSWycZ-qePwwEXrCpQ02A",
        "location": genomic_del5_37_loc,
        "copies": 2,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_del5_cx_37(genomic_del5_37_loc):
    """Create test fixture copy number change variation"""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.ICA0LcW-0gaOgJziPZJcGdXQ_MR3OiDu",
        "location": genomic_del5_37_loc,
        "copyChange": "efo:0030064",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_del6_cn_38(genomic_del6_seq_loc):
    """Create test fixture copy number count variation"""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.CZEc44pX7Dh9yJARvvz6EW9oQvgkbwYf",
        "location": genomic_del6_seq_loc,
        "copies": 1,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_del6_cx_38(genomic_del6_seq_loc):
    """Create test fixture copy number change variation"""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.02f8o7Gz9y11bMOWl7Vacc3V5J5M82bH",
        "location": genomic_del6_seq_loc,
        "copyChange": "efo:0030071",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_del6_37_loc():
    """Create test fixture GRCh37 deletion location"""
    return {
        "id": "ga4gh:SL.nENEk628TbRKyB9H9n12ssAQZnpwJDUo",
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
        "id": "ga4gh:CN.dPP5zGqiwSPNnuj49POIaul6_3msLwgF",
        "location": genomic_del6_37_loc,
        "copies": 1,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_del6_cx_37(genomic_del6_37_loc):
    """Create test fixture copy number change variation"""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.2YOJ93--D5aER5HsmnM4JP8GLWCcKKN-",
        "location": genomic_del6_37_loc,
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
    expected.id = "ga4gh:CN.BvXjVTvckz3zAp8yATcQjLp-1sWyEO9d"
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
    expected.id = "ga4gh:CN.40RjBTowl-97BT5vsPUgqdLJKNvL583c"
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
    expected.id = "ga4gh:CN.SiA0poCKmElUkXA-vTsNnuUiEBkzG6v_"
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
    expected.id = "ga4gh:CN.214AwcxGQiu1rY8UYQpud23sQI5DJbm1"
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
    expected.id = "ga4gh:CN.XvDElgE55k4blDLkzlDhwh5xcGyZrucn"
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
    expected.id = "ga4gh:CN.px3rDGLGlOmJLGKBoojl1UrFKu6Rhb1P"
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
    expected.id = "ga4gh:CN.S0viz7yEPYtrHIBxBLlCKnIc9x8FspQp"
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
    expected.id = "ga4gh:CN.d8thjAuvDC3acnDwFnz2Mg6PEMipVNJk"
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
    expected.id = "ga4gh:CN.JM71BLpaKrASgTbtoKaSXZ6GZ4FpeDsZ"
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
    expected.id = "ga4gh:CN.9qA8nLBcResblQhR0xfz16vPHP-tYXIA"
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
    expected.id = "ga4gh:CN.7tln8gxFt8FfLrJH7XlcX7PVSQCQqnNZ"
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
    expected.id = "ga4gh:CN.ydRovIJjiYNRUGv0w3iJ-pgK08MuYqGR"
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
