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
        "id": "ga4gh:CX.8GScA3atWKldpyaCkglqmMW5j4DkJO2a",
        "subject": genomic_dup1_seq_loc_not_normalized,
        "copyChange": "efo:0030069",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_dup1_37_loc():
    """Create test fixture GRCh37 duplication subject"""
    return {
        "id": "ga4gh:SL.7YabErC1ILGrve59Y8N6cOqhF7b2er4v",
        "sequence": "ga4gh:SQ.VNBualIltAyi2AI_uXcKU7M9XUOuA7MS",
        "start": 49568694,
        "end": 49568695,
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_dup1_cn_37(genomic_dup1_37_loc):
    """Create test fixture copy number count variation (not normalized)"""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.g7A0uaj_-VEV90tUW2g4VAnW690YNXeA",
        "subject": genomic_dup1_37_loc,
        "copies": 3,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_dup1_cx_37(genomic_dup1_37_loc):
    """Create test fixture copy number change variation"""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.QYFu49s3GBXjvWgkRrTYnSXRx3izj9dY",
        "subject": genomic_dup1_37_loc,
        "copyChange": "efo:0030069",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_dup2_cx_38(genomic_dup2_seq_loc_normalized):
    """Create test fixture copy number change variation"""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.bswUIc9CeAntdH9-d7jJv92iFpKy9o7d",
        "subject": genomic_dup2_seq_loc_normalized,
        "copyChange": "efo:0030067",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_dup2_37_loc():
    """Create test fixture GRCh37 duplication subject"""
    return {
        "id": "ga4gh:SL.32CVFHZY7OtujZPT08RU3F5LyqmBjsGn",
        "sequence": "ga4gh:SQ.v7noePfnNpK8ghYXEqZ9NukMXW7YeNsm",
        "start": 33229406,
        "end": 33229410,
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_dup2_cn_37(genomic_dup2_37_loc):
    """Create test fixture copy number count variation"""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.ky8FZkw6mvpTbL7FOJyBhY9h-A99-eZt",
        "subject": genomic_dup2_37_loc,
        "copies": 3,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_dup2_cx_37(genomic_dup2_37_loc):
    """Create test fixture copy number change variation"""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.evtjNtc7abS7fGQNLdrwelVUoAb4lxJK",
        "subject": genomic_dup2_37_loc,
        "copyChange": "efo:0030067",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_dup3_cn_38(genomic_del3_dup3_loc_not_normalized):
    """Create test fixture copy number count variation"""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.s0IG7djsMJhfP-RfuUMP9sArG6yABGTL",
        "subject": genomic_del3_dup3_loc_not_normalized,
        "copies": 2,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_dup3_cx_38(genomic_del3_dup3_loc_not_normalized):
    """Create test fixture copy number change variation"""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.muSgNlXbctZ4tUfjHO4Z_R6oRvQmiv6B",
        "subject": genomic_del3_dup3_loc_not_normalized,
        "copyChange": "efo:0030072",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_dup3_37_loc():
    """Create test fixture GRCh37 duplication subject"""
    return {
        "id": "ga4gh:SL.412An18pXfftZhNIsf3uzIO8yGwUiU2v",
        "sequence": "ga4gh:SQ.v7noePfnNpK8ghYXEqZ9NukMXW7YeNsm",
        "start": [31078343, 31118467],
        "end": [33292396, 33435269],
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_dup3_cn_37(genomic_dup3_37_loc):
    """Create test fixture copy number count variation"""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.pSMP4c6tFLc2GTpfuAs6VPnQmT1N0tXv",
        "subject": genomic_dup3_37_loc,
        "copies": 2,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_dup3_cx_37(genomic_dup3_37_loc):
    """Create test fixture copy number change variation"""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.fWhKH9pjazIuHwPB-LGQQLyUuocRUy9G",
        "subject": genomic_dup3_37_loc,
        "copyChange": "efo:0030072",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_dup4_cn_38(genomic_dup4_loc):
    """Create test fixture copy number count variation"""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.J1YKi-YzLjS1-BSGNrvXmk1ny8Xa-R5A",
        "subject": genomic_dup4_loc,
        "copies": 3,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_dup4_cx_38(genomic_dup4_loc):
    """Create test fixture copy number change variation"""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.3Fga1qCQXHP2_kmLcDliH2CyZV4X7P9f",
        "subject": genomic_dup4_loc,
        "copyChange": "efo:0030069",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_dup4_37_loc():
    """Create test fixture GRCh37 duplication subject"""
    return {
        "id": "ga4gh:SL.X_DGlcQ7MdEGL6qp8_K-V_jxLr_5obLl",
        "sequence": "ga4gh:SQ.iy_UbUrvECxFRX5LPTH_KPojdlT7BKsf",
        "start": [None, 29652251],
        "end": [29981821, None],
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_dup4_cn_37(genomic_dup4_37_loc):
    """Create test fixture copy number change variation"""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.a-P2c5-6SXjwHkuQosgSvK-6mC1-7KpR",
        "subject": genomic_dup4_37_loc,
        "copies": 3,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_dup4_cx_37(genomic_dup4_37_loc):
    """Create test fixture copy number change variation"""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.5D9lKosOsJcyra3ll4a_x3X4O1UMBTB5",
        "subject": genomic_dup4_37_loc,
        "copyChange": "efo:0030069",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_dup5_cn_38(genomic_dup5_loc):
    """Create test fixture copy number count variation"""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.i9rCqRSvXSFNsUf65g54Zh1G-cp-WqeS",
        "subject": genomic_dup5_loc,
        "copies": 4,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_dup5_cx_38(genomic_dup5_loc):
    """Create test fixture copy number change variation"""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.ivbv2JfqPLnf99e0xdfHxBfRDMav8wUS",
        "subject": genomic_dup5_loc,
        "copyChange": "efo:0030067",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_dup5_37_loc():
    """Create test fixture GRCh37 duplication subject"""
    return {
        "id": "ga4gh:SL.syVEBRrdR_BEM70UTDRglclDKsnkxcB3",
        "sequence": "ga4gh:SQ.v7noePfnNpK8ghYXEqZ9NukMXW7YeNsm",
        "start": [None, 153287262],
        "end": 153357667,
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_dup5_cn_37(genomic_dup5_37_loc):
    """Create test fixture copy number count variation"""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.YJPXMzBfc1dFt_nF6w76pclkZ4SwgbGu",
        "subject": genomic_dup5_37_loc,
        "copies": 4,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_dup5_cx_37(genomic_dup5_37_loc):
    """Create test fixture copy number change variation"""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.pnWcJ6hCG463Sdm7_S09buc7OYBLQtLV",
        "subject": genomic_dup5_37_loc,
        "copyChange": "efo:0030067",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_dup6_cn_38(genomic_dup6_loc):
    """Create test fixture copy number count variation"""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.5RLW5eq2NW-oLvszAvJmg2_2QnqGYvKb",
        "subject": genomic_dup6_loc,
        "copies": 2,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_dup6_cx_38(genomic_dup6_loc):
    """Create test fixture copy number change variation"""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.L2UM699rAKlzdTKGkUx4WTU9CXnXnyqt",
        "subject": genomic_dup6_loc,
        "copyChange": "efo:0030064",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_dup6_37_loc():
    """Create test fixture GRCh37 duplication subject"""
    return {
        "id": "ga4gh:SL.wBR6YIurIYwlEi5n2MzyZDC5auO9JkIO",
        "sequence": "ga4gh:SQ.v7noePfnNpK8ghYXEqZ9NukMXW7YeNsm",
        "start": 153287262,
        "end": [153357667, None],
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_dup6_cn_37(genomic_dup6_37_loc):
    """Create test fixture copy number count variation"""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.T1q-ngv4nvfItg2Cj1MUmNSzsCTTlw-B",
        "subject": genomic_dup6_37_loc,
        "copies": 2,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_dup6_cx_37(genomic_dup6_37_loc):
    """Create test fixture copy number change variation"""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.LDMRe4O_Mf7muelwWn8VYR0VF2kmK-sT",
        "subject": genomic_dup6_37_loc,
        "copyChange": "efo:0030064",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_del1_cx_38(genomic_del1_seq_loc):
    """Create test fixture copy number change variation"""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.WkRsG6w5UfrBXNOS2yHk6VlI10Y-VeSE",
        "subject": genomic_del1_seq_loc,
        "copyChange": "efo:0030064",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_del1_37_loc():
    """Create test fixture GRCh37 deletion subject"""
    return {
        "id": "ga4gh:SL.pIzhJvbahc_sELvQlgd8h4NKXBgStfyx",
        "sequence": "ga4gh:SQ.VNBualIltAyi2AI_uXcKU7M9XUOuA7MS",
        "start": 10191494,
        "end": 10191495,
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_del1_cn_37(genomic_del1_37_loc):
    """Create test fixture copy number count variation"""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.ELYZsSswa0uRwNkO6jxB-iM-d7bmivaI",
        "subject": genomic_del1_37_loc,
        "copies": 1,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_del1_cx_37(genomic_del1_37_loc):
    """Create test fixture copy number change variation"""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.krxGdPCUUdl3QY3l7okgYRDULGeemnIA",
        "subject": genomic_del1_37_loc,
        "copyChange": "efo:0030064",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_del2_cx_38(genomic_del2_seq_loc):
    """Create test fixture copy number change variation"""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.3kx1MS2vjs5lhy50AbxW5gfG20Xt5iDE",
        "subject": genomic_del2_seq_loc,
        "copyChange": "efo:0030071",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_del2_37_loc():
    """Create test fixture GRCh37 deletion subject"""
    return {
        "id": "ga4gh:SL.5soMcT3HXUuNrOUlh5aDZ9JMjEgL7twY",
        "sequence": "ga4gh:SQ.VNBualIltAyi2AI_uXcKU7M9XUOuA7MS",
        "start": 10188278,
        "end": 10188297,
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_del2_cn_37(genomic_del2_37_loc):
    """Create test fixture copy number count variation"""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.7LhqNGXQVrPpuhdkucNovfmTRixp07L7",
        "subject": genomic_del2_37_loc,
        "copies": 1,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_del2_cx_37(genomic_del2_37_loc):
    """Create test fixture copy number change variation"""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.KgVgL2txdV5GFYVvR0MezDSTBIfKk2-n",
        "subject": genomic_del2_37_loc,
        "copyChange": "efo:0030071",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_del3_cn_38(genomic_del3_dup3_loc_not_normalized):
    """Create test fixture copy number count variation"""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.s0IG7djsMJhfP-RfuUMP9sArG6yABGTL",
        "subject": genomic_del3_dup3_loc_not_normalized,
        "copies": 2,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_del3_cx_38(genomic_del3_dup3_loc_not_normalized):
    """Create test fixture copy number change variation"""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.p4r77sGjW2uZUsWreZu9_QbldmmZCQ_Y",
        "subject": genomic_del3_dup3_loc_not_normalized,
        "copyChange": "efo:0030069",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_del3_37_loc():
    """Create test fixture GRCh37 deletion subject"""
    return {
        "id": "ga4gh:SL.412An18pXfftZhNIsf3uzIO8yGwUiU2v",
        "sequence": "ga4gh:SQ.v7noePfnNpK8ghYXEqZ9NukMXW7YeNsm",
        "start": [31078343, 31118467],
        "end": [33292396, 33435269],
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_del3_cn_37(genomic_del3_37_loc):
    """Create test fixture copy number count variation"""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.pSMP4c6tFLc2GTpfuAs6VPnQmT1N0tXv",
        "subject": genomic_del3_37_loc,
        "copies": 2,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_del3_cx_37(genomic_del3_37_loc):
    """Create test fixture copy number change variation"""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.WcaUFEEEFwveQnTSX62zdPvnYhPAeuA6",
        "subject": genomic_del3_37_loc,
        "copyChange": "efo:0030069",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_del4_cn_38(genomic_del4_seq_loc):
    """Create test fixture copy number count variation"""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.QeImfnMeJu1XlRY5yDlUUMdoCzeH3tux",
        "subject": genomic_del4_seq_loc,
        "copies": 4,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_del4_cx_38(genomic_del4_seq_loc):
    """Create test fixture copy number change variation"""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.FalxH2l6ld0FAcF69ylzOguaDPdp_mKJ",
        "subject": genomic_del4_seq_loc,
        "copyChange": "efo:0030067",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_del4_37_loc():
    """Create test fixture GRCh37 deletion subject"""
    return {
        "id": "ga4gh:SL.uTvw_gJLMS3I9j6IjhzsQqY0Yx6HUZ9g",
        "sequence": "ga4gh:SQ.v7noePfnNpK8ghYXEqZ9NukMXW7YeNsm",
        "start": [None, 31138612],
        "end": [33357594, None],
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_del4_cn_37(genomic_del4_37_loc):
    """Create test fixture copy number count variation"""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.dNcALkFlv6tm0MJTqfc7xREUvvB-KTW7",
        "subject": genomic_del4_37_loc,
        "copies": 4,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_del4_cx_37(genomic_del4_37_loc):
    """Create test fixture copy number change variation"""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.BzLCu4la3PUs6uvxQhBbEEAb3GXl_Dw_",
        "subject": genomic_del4_37_loc,
        "copyChange": "efo:0030067",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_del5_cn_38(genomic_del5_seq_loc):
    """Create test fixture copy number count variation"""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.Jv2-gjHlHkHw-UG-FXFG24KBeStB4oz7",
        "subject": genomic_del5_seq_loc,
        "copies": 2,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_del5_cx_38(genomic_del5_seq_loc):
    """Create test fixture copy number change variation"""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.hJI3MKVVfQ-IZuOnQ7BUbTi1rZ7cUEAz",
        "subject": genomic_del5_seq_loc,
        "copyChange": "efo:0030064",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_del5_37_loc():
    """Create test fixture GRCh37 deletion subject"""
    return {
        "id": "ga4gh:SL.mLDP5XN1HEr13tpKGmjMEnb1joojz-_d",
        "sequence": "ga4gh:SQ.v7noePfnNpK8ghYXEqZ9NukMXW7YeNsm",
        "start": [None, 18593473],
        "end": 18671749,
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_del5_cn_37(genomic_del5_37_loc):
    """Create test fixture copy number count variation"""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.X8rLxAHB7V5eL4aa0rgczlR1tViay_dK",
        "subject": genomic_del5_37_loc,
        "copies": 2,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_del5_cx_37(genomic_del5_37_loc):
    """Create test fixture copy number change variation"""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.Y0e06LvtrQVte3TB_kgZrUA1QVNG_qt8",
        "subject": genomic_del5_37_loc,
        "copyChange": "efo:0030064",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_del6_cn_38(genomic_del6_seq_loc):
    """Create test fixture copy number count variation"""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.z-cD_rqwLYx1dEmNtpM7t-27cEx7DPep",
        "subject": genomic_del6_seq_loc,
        "copies": 1,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_del6_cx_38(genomic_del6_seq_loc):
    """Create test fixture copy number change variation"""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.iSNMsibFyqI4kSIB-C-OQFjIb41haGZt",
        "subject": genomic_del6_seq_loc,
        "copyChange": "efo:0030071",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_del6_37_loc():
    """Create test fixture GRCh37 deletion subject"""
    return {
        "id": "ga4gh:SL.DQ0Y1rdwcUevAxueB1_qH029AcwLwnS-",
        "sequence": "ga4gh:SQ.KqaUhJMW3CDjhoVtBetdEKT1n6hM-7Ek",
        "start": 133783901,
        "end": [133785996, None],
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_del6_cn_37(genomic_del6_37_loc):
    """Create test fixture copy number count variation"""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.nLZ_2pSNzw8aM8kfc3ga4RaU6uO_ZfYe",
        "subject": genomic_del6_37_loc,
        "copies": 1,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_del6_cx_37(genomic_del6_37_loc):
    """Create test fixture copy number change variation"""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.GP51Wy7zQmwgRCfB8E_Ifd6vBIk9P1XT",
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
    expected.id = "ga4gh:CN.gZmQQZBNf4SwopDTEx87ggIYkUhu06Rb"
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
    expected.id = "ga4gh:CN.bP3qlci2dIJNmI8L1cqfxRE97U3oVvmc"
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
    expected.id = "ga4gh:CN.IddtKtPorUXk4dIRenAn7Yj-zJv_wEyX"
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
    expected.id = "ga4gh:CN.ItTxiV7Ctua6Gi17JEClouBC3hbDOy5P"
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
    expected.id = "ga4gh:CN.nnLG0bgtgrcL-E2jhm_bZOl589vvC76G"
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
    expected.id = "ga4gh:CN.BWEXvQ-rjeKaylhUNXYvAEGrMw97ZgDx"
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
    expected.id = "ga4gh:CN.eQv1lrBI0S1WLw95eRieuvSHRffWXNmS"
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
    expected.id = "ga4gh:CN.C_wesdHpdhAB4aNRoCI87tHoUBTDwmDx"
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
    expected.id = "ga4gh:CN.K4DEqpFr_80NtRAbjSnTTXvQqsYJJxIv"
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
    expected.id = "ga4gh:CN.bQkERTvaMlQ635IamrQXjNSK3qamIYK5"
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
    expected.id = "ga4gh:CN.jvwmKulUMO6WKNyxL1tB9IEaZsKM8Lu3"
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
    expected.id = "ga4gh:CN.03lwEhXRvdNopMxHg0OsQFyizzO1Dabe"
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
