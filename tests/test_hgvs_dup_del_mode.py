"""Module for testing HGVS Dup Del mode."""
import pytest
from ga4gh.vrs import models

from tests.conftest import assertion_checks
from variation.schemas.normalize_response_schema import HGVSDupDelModeOption


@pytest.fixture(scope="module")
def test_handler(test_query_handler):
    """Create test fixture for normalize handler"""
    return test_query_handler.normalize_handler


@pytest.fixture(scope="module")
def genomic_dup1_lse(genomic_dup1_seq_loc_normalized):
    """Create a test fixture for genomic dup LSE."""
    params = {
        "type": "Allele",
        "id": "ga4gh:VA.CHNQRjx52keAGF5WcbvKORtfLiitZKE4",
        "location": genomic_dup1_seq_loc_normalized,
        "state": {
            "type": "ReferenceLengthExpression",
            "repeatSubunitLength": 1,
            "length": 3,
            "sequence": "GGG",
        },
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def genomic_dup1_cx(genomic_dup1_seq_loc_not_normalized):
    """Create a test fixture for genomic dup copy number change."""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.7WKEz2E_jwZZdyRc2Gw-_LIbHDJyRXwr",
        "location": genomic_dup1_seq_loc_not_normalized,
        "copyChange": "efo:0030072",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_dup1_free_text_seq_loc_normalized():
    """Create genomic dup1 free text sequence location"""
    return {
        "id": "ga4gh:SL.iyddzpD5lYY2Ayv87Np462l6P8QH7rH9",
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.tpvbnWsfEGqip8gJQZnWJAF8-bWDUDKd",
        },
        "start": 1032,
        "end": 1034,
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_dup1_free_text_seq_loc_not_normalized():
    """Create genomic dup1 free text sequence location"""
    return {
        "id": "ga4gh:SL.L89XFOyAxF-wdQHXUV8OAAkx80Mltokc",
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.tpvbnWsfEGqip8gJQZnWJAF8-bWDUDKd",
        },
        "start": 1033,
        "end": 1034,
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_dup1_free_text_lse(genomic_dup1_free_text_seq_loc_normalized):
    """Create a test fixture for genomic dup LSE."""
    params = {
        "type": "Allele",
        "id": "ga4gh:VA.muSdvI3Q126oFLKx3DrVkbzGfQ40kFhx",
        "location": genomic_dup1_free_text_seq_loc_normalized,
        "state": {
            "type": "ReferenceLengthExpression",
            "repeatSubunitLength": 1,
            "length": 3,
            "sequence": "GGG",
        },
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def genomic_dup1_free_text_cn(genomic_dup1_free_text_seq_loc_not_normalized):
    """Create a test fixture for genomic dup copy number count."""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.-97yN7Nq98gCVA7s0VslwuNdDVFLW6Af",
        "location": genomic_dup1_free_text_seq_loc_not_normalized,
        "copies": 3,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_dup2_lse(genomic_dup2_seq_loc_normalized):
    """Create a test fixture for genomic dup LSE."""
    params = {
        "type": "Allele",
        "id": "ga4gh:VA.u4ffOvroo0SV1X13zWMA41EOdu1QSO9B",
        "location": genomic_dup2_seq_loc_normalized,
        "state": {
            "type": "ReferenceLengthExpression",
            "repeatSubunitLength": 4,
            "length": 8,
            "sequence": "TCTATCTA",
        },
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def genomic_dup2_cx(genomic_dup2_seq_loc_normalized):
    """Create a test fixture for genomic dup copy number change."""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.g4l6d1hb3Rd1slsYWSe4Z4x3ocKdCB3w",
        "location": genomic_dup2_seq_loc_normalized,
        "copyChange": "efo:0030070",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def seq_loc_gt_100_bp():
    """Create seq loc for positions 33211290, 33211490 on NC_000023.11"""
    return {
        "id": "ga4gh:SL.HYv7UB8dh8paRuy_Sb3g4sHQaTqJ3m8Q",
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
        },
        "start": 33211289,
        "end": 33211490,
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_dup2_rle2(seq_loc_gt_100_bp):
    """Create a test fixture for genomic dup RSE where bp > 100."""
    params = {
        "type": "Allele",
        "id": "ga4gh:VA.gOA4L7Juk4KcUZnq4CBOk32-gkuz5keM",
        "location": seq_loc_gt_100_bp,
        "state": {
            "type": "ReferenceLengthExpression",
            "repeatSubunitLength": 201,
            "length": 402,
        },
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def genomic_dup2_free_text_seq_loc():
    """Create genomic dup2 free text sequence location"""
    return {
        "id": "ga4gh:SL.D4MxySRp4-wlbC3whkRZIhcfON2pKKgx",
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.1DeZLYHMnd-smp3GDlpRxETb9_0AokO7",
        },
        "start": 256,
        "end": 260,
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_dup2_free_text_default(genomic_dup2_free_text_seq_loc):
    """Create a test fixture for genomic dup default and LSE."""
    params = {
        "type": "Allele",
        "id": "ga4gh:VA.rEgZI-6A3SdNKhqqNapVYlcF_mzUeUGg",
        "location": genomic_dup2_free_text_seq_loc,
        "state": {
            "type": "ReferenceLengthExpression",
            "repeatSubunitLength": 4,
            "length": 8,
            "sequence": "TAGATAGA",
        },
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def genomic_dup2_free_text_cn(genomic_dup2_free_text_seq_loc):
    """Create a test fixture for genomic dup copy number count."""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.Dg3jMd1lsFKpXiJAPfzrh_50PQM2g1C3",
        "location": genomic_dup2_free_text_seq_loc,
        "copies": 3,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_dup3_cn(genomic_del3_dup3_loc_not_normalized):
    """Create a test fixture for genomic dup copy number count."""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.Ka-Wsibx4GHmHgurHCdk0W4deqZt26y4",
        "location": genomic_del3_dup3_loc_not_normalized,
        "copies": 2,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_dup3_cx(genomic_del3_dup3_loc_not_normalized):
    """Create a test fixture for genomic dup copy number change."""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.Gq5d-uH01bgf7m703dOSfSp_29wFqpsb",
        "location": genomic_del3_dup3_loc_not_normalized,
        "copyChange": "efo:0030070",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_dup3_free_text_subject():
    """Create test fixture for genomic dup3 free text location"""
    return {
        "id": "ga4gh:SL.h-0akgFon48yzCHKVdRwU5ImXRux4huN",
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
        },
        "start": [31147273, 31147277],
        "end": [31182738, 31182740],
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_dup3_free_text_cx(genomic_dup3_free_text_subject):
    """Create a test fixture for genomic dup copy number change."""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.aJYe9oIrxrwI0VTkEnvC0fKIT5x-vIne",
        "location": genomic_dup3_free_text_subject,
        "copyChange": "efo:0030070",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_dup3_free_text_cn(genomic_dup3_free_text_subject):
    """Create a test fixture for genomic dup copy number count."""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.KXh3T716GG9uGwq9ejGmflT_67fNN8n3",
        "location": genomic_dup3_free_text_subject,
        "copies": 4,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_dup4_cx(genomic_dup4_loc):
    """Create a test fixture for genomic dup copy number change."""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.4JLs2ICAAvj5JgG0xHJk1voSKLb8gNQ9",
        "location": genomic_dup4_loc,
        "copyChange": "efo:0030070",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_dup4_cn(genomic_dup4_loc):
    """Create a test fixture for genomic dup copy number count."""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.KqbQewUgZYfmottbgn1xYq58DiPVU5SZ",
        "location": genomic_dup4_loc,
        "copies": 3,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_dup4_free_text_subject():
    """Create test fixture for genomic dup4 free text location"""
    return {
        "id": "ga4gh:SL.SIeDb2iPT5pM-1SDKM9ew8NjzZAgF8nb",
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.dLZ15tNO1Ur0IcGjwc3Sdi_0A6Yf4zm7",
        },
        "start": [None, 1674441],
        "end": [1684571, None],
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_dup4_free_text_cx(genomic_dup4_free_text_subject):
    """Create a test fixture for genomic dup copy number change."""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.A7iAltzEFjlPJPCBfdMwjsis-vt51o3L",
        "location": genomic_dup4_free_text_subject,
        "copyChange": "efo:0030070",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_dup4_free_text_cn(genomic_dup4_free_text_subject):
    """Create a test fixture for genomic dup copy number count."""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.hnSfHkympkuTJQlfJHjHUqvUMU-EM2_Z",
        "location": genomic_dup4_free_text_subject,
        "copies": 3,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_dup5_cx(genomic_dup5_loc):
    """Create a test fixture for genomic dup5 copy number change."""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.6KvwSUu1Vp3FcC2VzbZxqpLAouOMCPi9",
        "location": genomic_dup5_loc,
        "copyChange": "efo:0030070",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_dup5_cn(genomic_dup5_loc):
    """Create a test fixture for genomic dup5 copy number count."""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.nlDhmSyOYeLZ8Fv2_F0niIraPbHUvpOU",
        "location": genomic_dup5_loc,
        "copies": 3,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_dup6_cx(genomic_dup6_loc):
    """Create a test fixture for genomic dup copy number change."""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.yMUFkF1QBwq3mA2tUS8wTLH3--dEHJJD",
        "location": genomic_dup6_loc,
        "copyChange": "efo:0030070",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_dup6_cn(genomic_dup6_loc):
    """Create a test fixture for genomic dup copy number count."""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.KSFn5KQIPuPVJ6FjWaF0vzl7eRwwHbX9",
        "location": genomic_dup6_loc,
        "copies": 2,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_del1_lse(genomic_del1_seq_loc):
    """Create a test fixture for genomic del LSE."""
    params = {
        "type": "Allele",
        "id": "ga4gh:VA.gztc0BFS6p5V1_QVnEYIJ6DwzZQeDCd2",
        "location": genomic_del1_seq_loc,
        "state": {
            "type": "ReferenceLengthExpression",
            "repeatSubunitLength": 1,
            "length": 0,
            "sequence": "",
        },
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def genomic_del1_cx(genomic_del1_seq_loc):
    """Create a test fixture for genomic del copy number change."""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.LWRBNtBgcETMXEKezrr7WUPjO9WoOaqL",
        "location": genomic_del1_seq_loc,
        "copyChange": "efo:0030064",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_del1_rle(genomic_del1_seq_loc):
    """Create a test fixture for genomic del RSE."""
    params = {
        "type": "Allele",
        "id": "ga4gh:VA.Kg0FrJBjKRtIDsIKO0LxAwOPiXIOowoc",
        "location": genomic_del1_seq_loc,
        "state": {
            "type": "ReferenceLengthExpression",
            "repeatSubunitLength": 1,
            "length": 2,
        },
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def genomic_del1_free_text_seq_loc():
    """Create genomic del1 free text sequence location"""
    return {
        "id": "ga4gh:SL.072FoTQ7ZWLfOOOdyTI3Vj5pc2qwDii6",
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.xBKOKptLLDr-k4hTyCetvARn16pDS_rW",
        },
        "start": 557,
        "end": 558,
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_del1_free_text_lse(genomic_del1_free_text_seq_loc):
    """Create a test fixture for genomic del LSE."""
    params = {
        "type": "Allele",
        "id": "ga4gh:VA.8xbobLhnVeBLQ6ANUur7BcPNdXrLsSja",
        "location": genomic_del1_free_text_seq_loc,
        "state": {
            "type": "ReferenceLengthExpression",
            "repeatSubunitLength": 1,
            "length": 0,
            "sequence": "",
        },
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def genomic_del1_free_text_cn(genomic_del1_free_text_seq_loc):
    """Create a test fixture for genomic del copy number count."""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.OakJW5ITpO4m1ffP4tGoBK72_IIqEBM6",
        "location": genomic_del1_free_text_seq_loc,
        "copies": 1,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_del1_free_text_rle(genomic_del1_free_text_seq_loc):
    """Create a test fixture for genomic del RSE."""
    params = {
        "type": "Allele",
        "id": "ga4gh:VA.XZMMF_xhn76bLMxN5RnewNgrXkYuK-ni",
        "location": genomic_del1_free_text_seq_loc,
        "state": {
            "type": "ReferenceLengthExpression",
            "repeatSubunitLength": 1,
            "length": 0,
        },
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def genomic_del2_lse(genomic_del2_seq_loc):
    """Create a test fixture for genomic del LSE."""
    params = {
        "type": "Allele",
        "id": "ga4gh:VA.9NmH0sRYerurt-CE6WlF9UaxZiujByIE",
        "location": genomic_del2_seq_loc,
        "state": {
            "type": "ReferenceLengthExpression",
            "repeatSubunitLength": 19,
            "length": 0,
            "sequence": "",
        },
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def genomic_del2_lse2(seq_loc_gt_100_bp):
    """Create a test fixture for genomic del LSE where bp > 100."""
    params = {
        "type": "Allele",
        "id": "ga4gh:VA.1_cveYe6e74MEUt8EdTQmEtW5t6nA5bU",
        "location": seq_loc_gt_100_bp,
        "state": {
            "type": "ReferenceLengthExpression",
            "repeatSubunitLength": 201,
            "length": 0,
            "sequence": "",
        },
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def genomic_del2_cx(genomic_del2_seq_loc):
    """Create a test fixture for genomic del copy number change."""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.xOZeCpcgWTj-xTYJdIeXbRy8h48qfbQ5",
        "location": genomic_del2_seq_loc,
        "copyChange": "efo:0030069",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_del2_rle(genomic_del2_seq_loc):
    """Create a test fixture for genomic del RSE."""
    params = {
        "type": "Allele",
        "id": "ga4gh:VA.J6hHiLw-qq27H8CZ8aQRdJwBGHqd3BvB",
        "location": genomic_del2_seq_loc,
        "state": {
            "type": "ReferenceLengthExpression",
            "repeatSubunitLength": 19,
            "length": 0,
        },
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def genomic_del2_free_text_seq_loc():
    """Create genomic del2 free text sequence location"""
    return {
        "id": "ga4gh:SL.b06yJ2UPwSSo-4bmYE8ZqHkDfo6_KZuu",
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.xBKOKptLLDr-k4hTyCetvARn16pDS_rW",
        },
        "start": 491,
        "end": 510,
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_del2_free_text_default(genomic_del2_free_text_seq_loc):
    """Create a test fixture for genomic del default and LSE."""
    params = {
        "type": "Allele",
        "id": "ga4gh:VA.ZmmZ_3it-b0nl8pdxaIG5ROYwTYhhRfk",
        "location": genomic_del2_free_text_seq_loc,
        "state": {
            "type": "ReferenceLengthExpression",
            "repeatSubunitLength": 19,
            "length": 0,
            "sequence": "",
        },
    }
    return models.Allele(**params)


@pytest.fixture(scope="module")
def genomic_del2_free_text_cnv(genomic_del2_free_text_seq_loc):
    """Create a test fixture for genomic del CNV."""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.DdDmUshRGGSxugHHqGI8agdffFmvwjFm",
        "location": genomic_del2_free_text_seq_loc,
        "copies": 1,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_del3_cx(genomic_del3_dup3_loc_not_normalized):
    """Create a test fixture for genomic del copy number change."""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.h9krodUtVK--XmizkspBdOrptRNqrDHm",
        "location": genomic_del3_dup3_loc_not_normalized,
        "copyChange": "efo:0030067",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_del3_cn(genomic_del3_dup3_loc_not_normalized):
    """Create a test fixture for genomic del copy number count."""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.Ka-Wsibx4GHmHgurHCdk0W4deqZt26y4",
        "location": genomic_del3_dup3_loc_not_normalized,
        "copies": 2,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_del3_free_text_subject():
    """Create test fixture for genomic del3 free text location"""
    return {
        "id": "ga4gh:SL.HMjpquCLV9iYib972N0_3tn9TvnevIga",
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
        },
        "start": [68839264, 68839267],
        "end": [68841121, 68841126],
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_del3_free_text_cx(genomic_del3_free_text_subject):
    """Create a test fixture for genomic del copy number change."""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.qHYL7nZaKXl-h9D7d4JvEmXDZtfV9M2A",
        "location": genomic_del3_free_text_subject,
        "copyChange": "efo:0030067",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_del3_free_text_cn(genomic_del3_free_text_subject):
    """Create a test fixture for genomic del copy number count."""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.ESfNw8Y3Vil7UpRiV2UY-zfMoZNRchA3",
        "location": genomic_del3_free_text_subject,
        "copies": 2,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_del4_cx(genomic_del4_seq_loc):
    """Create a test fixture for genomic del copy number change."""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.1DiUzraiKZLJb8oF8ynARS816fthsJpV",
        "location": genomic_del4_seq_loc,
        "copyChange": "efo:0030067",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_del4_cn(genomic_del4_seq_loc):
    """Create a test fixture for genomic del copy number count."""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.6RKML7P4zTx1U8EpJ1q7L23OXDEKFihS",
        "location": genomic_del4_seq_loc,
        "copies": 1,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_del4_free_text_subject():
    """Create test fixture for genomic del4 free text location"""
    return {
        "id": "ga4gh:SL.ebOW5blAtyPPVH512rIYi6cGsyKI2990",
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": "SQ.pnAqCRBrTsUoBghSD1yp_jXWSmlbdh4g",
        },
        "start": [None, 227022027],
        "end": [227025830, None],
        "type": "SequenceLocation",
    }


@pytest.fixture(scope="module")
def genomic_del4_free_text_cx(genomic_del4_free_text_subject):
    """Create a test fixture for genomic del copy number change."""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.iI0rwBfPht6XCBhYPsfWUUbTwDIRycFi",
        "location": genomic_del4_free_text_subject,
        "copyChange": "efo:0030067",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_del4_free_text_cn(genomic_del4_free_text_subject):
    """Create a test fixture for genomic del copy number count."""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.ps-GD4HSeJZjxnS1dhQrn4ntJFaA97a3",
        "location": genomic_del4_free_text_subject,
        "copies": 1,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_uncertain_del_2():
    """Create a genomic uncertain deletion on chr 2 test fixture."""
    params = {
        "id": "ga4gh:CX.q3OPPp2fWM5uM60RNHY_jDThCyxV3URW",
        "location": {
            "id": "ga4gh:SL.aRKiRW6-lS9CCLfcPJpQIGihZqoIOCZ_",
            "sequenceReference": {
                "type": "SequenceReference",
                "refgetAccession": "SQ.pnAqCRBrTsUoBghSD1yp_jXWSmlbdh4g",
            },
            "start": [None, 110104899],
            "end": [110207160, None],
            "type": "SequenceLocation",
        },
        "copyChange": "efo:0030067",
        "type": "CopyNumberChange",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_uncertain_del_y():
    """Create a genomic uncertain deletion on chr Y test fixture."""
    params = {
        "id": "ga4gh:CX.vR12PHS1zCnoYUi9CSX3ZwhGG38xa-RA",
        "location": {
            "id": "ga4gh:SL.N44ez-5301ZoNdLoiblcUvm__BS4-4Jv",
            "sequenceReference": {
                "type": "SequenceReference",
                "refgetAccession": "SQ.8_liLu1aycC0tPQPFmUaGXJLDs5SbPZ5",
            },
            "start": [None, 14076801],
            "end": [57165209, None],
            "type": "SequenceLocation",
        },
        "copyChange": "efo:0030067",
        "type": "CopyNumberChange",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_del5_cn_var(genomic_del5_seq_loc):
    """Create genomic del5 copy number count"""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.B2qGS47pqvyvBjFRQEj3MjdsqfXpnhhC",
        "location": genomic_del5_seq_loc,
        "copies": 3,
    }
    return models.CopyNumberCount(**params)


@pytest.fixture(scope="module")
def genomic_del5_cx_var(genomic_del5_seq_loc):
    """Create genomic del5 copy number change"""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.oUhToHxDGpH5NkuFaQmKTmbijF9z_Esb",
        "location": genomic_del5_seq_loc,
        "copyChange": "efo:0030067",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_del6_cx_var(genomic_del6_seq_loc):
    """Create genomic del6 copy number change"""
    params = {
        "type": "CopyNumberChange",
        "id": "ga4gh:CX.gmwszbknrMmklvVuu2yOqu5nOKV_fp72",
        "location": genomic_del6_seq_loc,
        "copyChange": "efo:0030067",
    }
    return models.CopyNumberChange(**params)


@pytest.fixture(scope="module")
def genomic_del6_cn_var(genomic_del6_seq_loc):
    """Create genomic del6 copy number count"""
    params = {
        "type": "CopyNumberCount",
        "id": "ga4gh:CN.CZEc44pX7Dh9yJARvvz6EW9oQvgkbwYf",
        "location": genomic_del6_seq_loc,
        "copies": 1,
    }
    return models.CopyNumberCount(**params)


def no_variation_check(resp, q):
    """Check that variation is None in normalize response"""
    assert resp.variation is None, q


@pytest.mark.asyncio
async def invalid_query_list_checks(query_list, test_handler):
    """Check that invalid queries in query list do not normalize"""
    for q in query_list:
        resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
        no_variation_check(resp, q)


@pytest.mark.asyncio
async def test_genomic_dup1(
    test_handler,
    genomic_dup1_lse,
    genomic_dup1_38_cn,
    genomic_dup1_cx,
    genomic_dup1_free_text_lse,
    genomic_dup1_free_text_cn,
):
    """Test that genomic duplication works correctly."""
    # https://reg.clinicalgenome.org/redmine/projects/registry/genboree_registry/allele?hgvsOrDescriptor=NC_000003.12%3Ag.49531262dup
    q = "NC_000003.12:g.49531262dup"  # 38
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp, genomic_dup1_lse)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp, genomic_dup1_lse)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
    )
    assertion_checks(resp, genomic_dup1_38_cn)

    resp = await test_handler.normalize(
        q,
        HGVSDupDelModeOption.COPY_NUMBER_CHANGE,
        copy_change=models.CopyChange.EFO_0030072,
    )
    assertion_checks(resp, genomic_dup1_cx)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.ALLELE)
    assertion_checks(resp, genomic_dup1_lse)

    q = "NC_000003.11:g.49568695dup"  # 37
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp, genomic_dup1_lse)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
    )
    assertion_checks(resp, genomic_dup1_38_cn)

    resp = await test_handler.normalize(
        q,
        HGVSDupDelModeOption.COPY_NUMBER_CHANGE,
        copy_change=models.CopyChange.EFO_0030072,
    )
    assertion_checks(resp, genomic_dup1_cx)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.ALLELE)
    assertion_checks(resp, genomic_dup1_lse)

    # Free Text
    for q in ["DAG1 g.49568695dup", "DAG1 g.49531262dup"]:  # 37  # 38
        resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
        assertion_checks(resp, genomic_dup1_free_text_lse)

        resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
        assertion_checks(resp, genomic_dup1_free_text_lse)

        resp = await test_handler.normalize(
            q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
        )
        assertion_checks(resp, genomic_dup1_free_text_cn)

        resp = await test_handler.normalize(q, HGVSDupDelModeOption.ALLELE)
        assertion_checks(resp, genomic_dup1_free_text_lse)

    # Invalid
    invalid_queries = [
        "NC_000007.13:g.159138670dup",
        "NC_000007.14:g.159345976dup",
        "BRAF g.140219337dup",
        "BRAF g.141024929dup",
    ]
    await invalid_query_list_checks(invalid_queries, test_handler)


@pytest.mark.asyncio
async def test_genomic_dup2(
    test_handler,
    genomic_dup2_lse,
    genomic_dup2_38_cn,
    genomic_dup2_cx,
    genomic_dup2_free_text_default,
    genomic_dup2_free_text_cn,
    genomic_dup2_rle2,
):
    """Test that genomic duplication works correctly."""
    # https://reg.clinicalgenome.org/redmine/projects/registry/genboree_registry/allele?hgvsOrDescriptor=NM_004006.2%3Ac.20_23dup
    q = "NC_000023.11:g.33211290_33211293dup"  # 38
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp, genomic_dup2_lse)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
    )
    assertion_checks(resp, genomic_dup2_38_cn)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE)
    assertion_checks(resp, genomic_dup2_cx)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.ALLELE)
    assertion_checks(resp, genomic_dup2_lse)

    q = "NC_000023.10:g.33229407_33229410dup"  # 37
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp, genomic_dup2_lse)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
    )
    assertion_checks(resp, genomic_dup2_38_cn)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE)
    assertion_checks(resp, genomic_dup2_cx)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.ALLELE)
    assertion_checks(resp, genomic_dup2_lse)

    # Free text
    for q in ["DMD g.33211290_33211293dup", "DMD g.33229407_33229410dup"]:  # 37  # 38
        resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
        assertion_checks(resp, genomic_dup2_free_text_default)

        resp = await test_handler.normalize(
            q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
        )
        assertion_checks(resp, genomic_dup2_free_text_cn)

        resp = await test_handler.normalize(q, HGVSDupDelModeOption.ALLELE)
        assertion_checks(resp, genomic_dup2_free_text_default)

    # Greater than 100 bps -> rse
    q = "NC_000023.11:g.33211290_33211490dup"
    resp = await test_handler.normalize(q)
    assertion_checks(resp, genomic_dup2_rle2)

    # Invalid
    invalid_queries = [
        "NC_000007.13:g.140413127_159138670dup",
        "NC_000007.14:g.140413127_159345976dup",
        "BRAF g.140219337_140924929dup",
        "BRAF g.140719326_141024929dup",
    ]
    await invalid_query_list_checks(invalid_queries, test_handler)


@pytest.mark.asyncio
async def test_genomic_dup3(
    test_handler,
    genomic_dup3_cx,
    genomic_dup3_cn,
    genomic_dup3_free_text_cn,
    genomic_dup3_free_text_cx,
):
    """Test that genomic duplication works correctly."""
    q = "NC_000023.11:g.(31060227_31100351)_(33274278_33417151)dup"  # 38
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp, genomic_dup3_cx)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=1
    )
    assertion_checks(resp, genomic_dup3_cn)

    resp = await test_handler.normalize(
        q,
        HGVSDupDelModeOption.COPY_NUMBER_CHANGE,
        copy_change=models.CopyChange.EFO_0030070,
    )
    assertion_checks(resp, genomic_dup3_cx)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.ALLELE)
    no_variation_check(resp, q)

    q = "NC_000023.10:g.(31078344_31118468)_(33292395_33435268)dup"  # 37
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp, genomic_dup3_cx)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=1
    )
    assertion_checks(resp, genomic_dup3_cn)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE)
    assertion_checks(resp, genomic_dup3_cx)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.ALLELE)
    no_variation_check(resp, q)

    # Free Text
    for q in ["DMD g.(31147274_31147278)_(31182737_31182739)dup"]:  # 38
        resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
        assertion_checks(resp, genomic_dup3_free_text_cx)

        resp = await test_handler.normalize(
            q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=3
        )
        assertion_checks(resp, genomic_dup3_free_text_cn)

        resp = await test_handler.normalize(q, HGVSDupDelModeOption.ALLELE)
        no_variation_check(resp, q)

    # Invalid
    invalid_queries = [
        "NC_000023.10:g.(31119221_31119227)_(31119300_155270562)dup",
        "NC_000023.11:g.(31119221_31119227)_(31119300_156040899)dup",
        "DMD g.(31060227_31100351)_(33274278_33417151)dup",
    ]
    await invalid_query_list_checks(invalid_queries, test_handler)


@pytest.mark.asyncio
async def test_genomic_dup4(
    test_handler,
    genomic_dup4_cn,
    genomic_dup4_cx,
    genomic_dup4_free_text_cn,
    genomic_dup4_free_text_cx,
):
    """Test that genomic duplication works correctly."""
    q = "NC_000020.11:g.(?_30417576)_(31394018_?)dup"  # 38
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp, genomic_dup4_cx)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
    )
    assertion_checks(resp, genomic_dup4_cn)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE)
    assertion_checks(resp, genomic_dup4_cx)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.ALLELE)
    no_variation_check(resp, q)

    q = "NC_000020.10:g.(?_29652252)_(29981821_?)dup"  # 37
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp, genomic_dup4_cx)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
    )
    assertion_checks(resp, genomic_dup4_cn)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE)
    assertion_checks(resp, genomic_dup4_cx)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.ALLELE)
    no_variation_check(resp, q)

    # Free Text
    for q in [
        "PRPF8 g.(?_1577736)_(1587865_?)dup",  # 37
        "PRPF8 g.(?_1674442)_(1684571_?)dup",  # 38
    ]:
        resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
        assertion_checks(resp, genomic_dup4_free_text_cx)

        resp = await test_handler.normalize(
            q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
        )
        assertion_checks(resp, genomic_dup4_free_text_cn)

        resp = await test_handler.normalize(q, HGVSDupDelModeOption.ALLELE)
        no_variation_check(resp, q)

    # Invalid
    invalid_queries = [
        "NC_000020.10:g.(?_29652252)_(63025530_?)dup",
        "NC_000020.11:g.(?_29652252)_(64444169_?)dup",
        "PRPF8 g.(?_1650628)_(1684571_?)dup",
    ]
    await invalid_query_list_checks(invalid_queries, test_handler)


@pytest.mark.asyncio
async def test_genomic_dup5(
    test_handler,
    genomic_dup5_cn,
    genomic_dup5_cx,
):
    """Test that genomic duplication works correctly."""
    q = "NC_000023.11:g.(?_154021812)_154092209dup"  # 38
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp, genomic_dup5_cx)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
    )
    assertion_checks(resp, genomic_dup5_cn)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE)
    assertion_checks(resp, genomic_dup5_cx)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.ALLELE)
    no_variation_check(resp, q)

    q = "NC_000023.10:g.(?_153287263)_153357667dup"  # 37
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp, genomic_dup5_cx)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
    )
    assertion_checks(resp, genomic_dup5_cn)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE)
    assertion_checks(resp, genomic_dup5_cx)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.ALLELE)
    no_variation_check(resp, q)

    # Free Text
    for q in [
        "MECP2 g.(?_153287263)_153357667dup",  # 37
        "MECP2 g.(?_154021812)_154092209dup",  # 38
    ]:
        resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
        assertion_checks(resp, genomic_dup5_cx)

        resp = await test_handler.normalize(
            q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
        )
        assertion_checks(resp, genomic_dup5_cn)

        resp = await test_handler.normalize(q, HGVSDupDelModeOption.ALLELE)
        no_variation_check(resp, q)

    # Invalid
    for q in [
        "NC_000023.10:g.(?_153287263)_155270561dup",
        "NC_000023.11:g.(?_154021812)_156040896dup",
        "MECP2 g.(?_154021812)_154097733dup"  # 37
        "MECP2 g.(?_154021572)_154092209dup",  # 38
    ]:
        resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
        assert resp.variation is None, q


@pytest.mark.asyncio
async def test_genomic_dup6(
    test_handler,
    genomic_dup6_cn,
    genomic_dup6_cx,
):
    """Test that genomic duplication works correctly."""
    q = "NC_000023.11:g.154021812_(154092209_?)dup"  # 38
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp, genomic_dup6_cx)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=1
    )
    assertion_checks(resp, genomic_dup6_cn)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE)
    assertion_checks(resp, genomic_dup6_cx)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.ALLELE)
    no_variation_check(resp, q)

    q = "NC_000023.10:g.153287263_(153357667_?)dup"  # 37
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp, genomic_dup6_cx)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=1
    )
    assertion_checks(resp, genomic_dup6_cn)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE)
    assertion_checks(resp, genomic_dup6_cx)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.ALLELE)
    no_variation_check(resp, q)

    # Free Text
    for q in [
        "MECP2 g.153287263_(153357667_?)dup",  # 37
        "MECP2 g.154021812_(154092209_?)dup",  # 38
    ]:
        resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
        assertion_checks(resp, genomic_dup6_cx)

        resp = await test_handler.normalize(
            q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=1
        )
        assertion_checks(resp, genomic_dup6_cn)

        resp = await test_handler.normalize(q, HGVSDupDelModeOption.ALLELE)
        no_variation_check(resp, q)

    # Invalid
    for q in [
        "NC_000023.10:g.153287263_(155270561_?)dup",
        "NC_000023.11:g.154021812_(156040896_?)dup",
        "MECP2 g.154021812_(154097733_?)dup"  # 37
        "MECP2 g.154021572_(154092209_?)dup",  # 38
    ]:
        resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
        assert resp.variation is None, q


@pytest.mark.asyncio
async def test_genomic_del1(
    test_handler,
    genomic_del1_lse,
    genomic_del1_38_cn,
    genomic_del1_cx,
    genomic_del1_free_text_lse,
    genomic_del1_free_text_cn,
):
    """Test that genomic deletion works correctly."""
    q = "NC_000003.12:g.10149811del"  # 38
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp, genomic_del1_lse)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
    )
    assertion_checks(resp, genomic_del1_38_cn)

    resp = await test_handler.normalize(
        q,
        HGVSDupDelModeOption.COPY_NUMBER_CHANGE,
        copy_change=models.CopyChange.EFO_0030064,
    )
    assertion_checks(resp, genomic_del1_cx)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.ALLELE)
    assertion_checks(resp, genomic_del1_lse)

    q = "NC_000003.11:g.10191495del"  # 37
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp, genomic_del1_lse)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
    )
    assertion_checks(resp, genomic_del1_38_cn)

    resp = await test_handler.normalize(
        q,
        HGVSDupDelModeOption.COPY_NUMBER_CHANGE,
        copy_change=models.CopyChange.EFO_0030064,
    )
    assertion_checks(resp, genomic_del1_cx)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.ALLELE)
    assertion_checks(resp, genomic_del1_lse)

    # Free text
    for q in ["VHL g.10191495del", "VHL g.10149811del"]:  # 37  # 38
        resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
        assertion_checks(resp, genomic_del1_free_text_lse)

        resp = await test_handler.normalize(
            q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
        )
        assertion_checks(resp, genomic_del1_free_text_cn)

        resp = await test_handler.normalize(q, HGVSDupDelModeOption.ALLELE)
        assertion_checks(resp, genomic_del1_free_text_lse)

    # Invalid
    invalid_queries = [
        "NC_000003.11:g.198022431del",
        "NC_000003.12:g.198295567del",
        "BRAF g.140413127del",
        "BRAF g.141024929del",
    ]
    await invalid_query_list_checks(invalid_queries, test_handler)


@pytest.mark.asyncio
async def test_genomic_del2(
    test_handler,
    genomic_del2_lse,
    genomic_del2_38_cn,
    genomic_del2_cx,
    genomic_del2_free_text_default,
    genomic_del2_free_text_cnv,
    genomic_del2_lse2,
):
    """Test that genomic deletion works correctly."""
    q = "NC_000003.12:g.10146595_10146613del"  # 38
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp, genomic_del2_lse)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
    )
    assertion_checks(resp, genomic_del2_38_cn)

    resp = await test_handler.normalize(
        q,
        HGVSDupDelModeOption.COPY_NUMBER_CHANGE,
        copy_change=models.CopyChange.EFO_0030069,
    )
    assertion_checks(resp, genomic_del2_cx)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.ALLELE)
    assertion_checks(resp, genomic_del2_lse)

    q = "NC_000003.11:g.10188279_10188297del"  # 37
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp, genomic_del2_lse)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
    )
    assertion_checks(resp, genomic_del2_38_cn)

    resp = await test_handler.normalize(
        q,
        HGVSDupDelModeOption.COPY_NUMBER_CHANGE,
        copy_change=models.CopyChange.EFO_0030069,
    )
    assertion_checks(resp, genomic_del2_cx)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.ALLELE)
    assertion_checks(resp, genomic_del2_lse)

    # Free text
    for q in ["VHL g.10188279_10188297del", "VHL g.10146595_10146613del"]:  # 37  # 38
        resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
        assertion_checks(resp, genomic_del2_free_text_default)

        resp = await test_handler.normalize(
            q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
        )
        assertion_checks(resp, genomic_del2_free_text_cnv)

        resp = await test_handler.normalize(q, HGVSDupDelModeOption.ALLELE)
        assertion_checks(resp, genomic_del2_free_text_default)

    # Check that del > 100 bps returns LSE
    q = "NC_000023.11:g.33211290_33211490del"
    resp = await test_handler.normalize(q)
    assertion_checks(resp, genomic_del2_lse2)

    # gnomad vcf
    q = "3-10146594-AATGTTGACGGACAGCCTAT-A"
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp, genomic_del2_lse)

    q = "3-10188278-AATGTTGACGGACAGCCTAT-A"
    resp = await test_handler.normalize(q)
    assertion_checks(resp, genomic_del2_lse)

    # Invalid
    invalid_queries = [
        "NC_000003.12:g.10146595_198295580del",
        "NC_000003.11:g.198022435_198022437del",
        "BRAF g.140413127_140419136del",
        "BRAF g.140719326_141024929del",
    ]
    await invalid_query_list_checks(invalid_queries, test_handler)


@pytest.mark.asyncio
async def test_genomic_del3(
    test_handler,
    genomic_del3_cn,
    genomic_del3_cx,
    genomic_del3_free_text_cn,
    genomic_del3_free_text_cx,
):
    """Test that genomic deletion works correctly."""
    q = "NC_000023.11:g.(31060227_31100351)_(33274278_33417151)del"  # 38
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp, genomic_del3_cx)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=3
    )
    assertion_checks(resp, genomic_del3_cn)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE)
    assertion_checks(resp, genomic_del3_cx)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.ALLELE)
    no_variation_check(resp, q)

    q = "NC_000023.10:g.(31078344_31118468)_(33292395_33435268)del"  # 37
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp, genomic_del3_cx)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=3
    )
    assertion_checks(resp, genomic_del3_cn)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE)
    assertion_checks(resp, genomic_del3_cx)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.ALLELE)
    no_variation_check(resp, q)

    # Free Text
    for q in [
        "EFNB1 g.(68059108_68059111)_(68060963_68060968)del",  # 37
        "EFNB1 g.(68839265_68839268)_(68841120_68841125)del",  # 38
    ]:
        resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
        assertion_checks(resp, genomic_del3_free_text_cx)

        resp = await test_handler.normalize(
            q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=3
        )
        assertion_checks(resp, genomic_del3_free_text_cn)

        resp = await test_handler.normalize(q, HGVSDupDelModeOption.ALLELE)
        no_variation_check(resp, q)

    # Invalid
    invalid_queries = [
        "NC_000023.11:g.(156040880_156040883)_(156040896_156040899)del",
        "NC_000023.10:g.(155270550_155270555)_(155270560_155270562)del",
        "EFNB1 g.(68048863_68048870)_(68842150_68842152)del",  # 37
        "EFNB1 g.(68829022_68829030)_(68842150_68842161)del",  # 38
    ]
    await invalid_query_list_checks(invalid_queries, test_handler)


@pytest.mark.asyncio
async def test_genomic_del4(
    test_handler,
    genomic_del4_cn,
    genomic_del4_cx,
    genomic_uncertain_del_2,
    genomic_uncertain_del_y,
    genomic_del4_free_text_cn,
    genomic_del4_free_text_cx,
):
    """Test that genomic deletion works correctly."""
    q = "NC_000023.11:g.(?_31120496)_(33339477_?)del"  # 38
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp, genomic_del4_cx)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
    )
    assertion_checks(resp, genomic_del4_cn)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE)
    assertion_checks(resp, genomic_del4_cx)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.ALLELE)
    no_variation_check(resp, q)

    q = "NC_000023.10:g.(?_31138613)_(33357594_?)del"  # 37
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp, genomic_del4_cx)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
    )
    assertion_checks(resp, genomic_del4_cn)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE)
    assertion_checks(resp, genomic_del4_cx)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.ALLELE)
    no_variation_check(resp, q)

    q = "NC_000002.12:g.(?_110104900)_(110207160_?)del"
    resp = await test_handler.normalize(q)
    assertion_checks(resp, genomic_uncertain_del_2)

    q = "NC_000024.10:g.(?_14076802)_(57165209_?)del"
    resp = await test_handler.normalize(q)
    assertion_checks(resp, genomic_uncertain_del_y)

    # Free Text
    for q in ["COL4A4 g.(?_227022028)_(227025830_?)del"]:  # 38
        resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
        assertion_checks(resp, genomic_del4_free_text_cx)

        resp = await test_handler.normalize(
            q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
        )
        assertion_checks(resp, genomic_del4_free_text_cn)

        resp = await test_handler.normalize(q, HGVSDupDelModeOption.ALLELE)
        no_variation_check(resp, q)

    # Invalid
    invalid_queries = [
        "NC_000023.11:g.(?_156040899)_(156040900_?)del",
        "NC_000024.10:g.(?_155270565)_(155270568_?)del",
        "COL4A4 g.(?_227002710)_(227003710_?)del",
        "COL4A4 g.(?_227867430)_(228029276_?)del",
    ]
    await invalid_query_list_checks(invalid_queries, test_handler)


@pytest.mark.asyncio
async def test_genomic_del5(
    test_handler,
    genomic_del5_cn_var,
    genomic_del5_cx_var,
):
    """Test that genomic deletion works correctly."""
    q = "NC_000023.11:g.(?_18575354)_18653629del"  # 38
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp, genomic_del5_cx_var)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=4
    )
    assertion_checks(resp, genomic_del5_cn_var)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE)
    assertion_checks(resp, genomic_del5_cx_var)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.ALLELE)
    no_variation_check(resp, q)

    q = "NC_000023.10:g.(?_18593474)_18671749del"  # 37
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp, genomic_del5_cx_var)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=4
    )
    assertion_checks(resp, genomic_del5_cn_var)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE)
    assertion_checks(resp, genomic_del5_cx_var)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.ALLELE)
    no_variation_check(resp, q)

    # Free text
    for q in ["CDKL5 g.(?_18575354)_18653629del"]:
        resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
        assertion_checks(resp, genomic_del5_cx_var)

        resp = await test_handler.normalize(
            q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=4
        )
        assertion_checks(resp, genomic_del5_cn_var)

        resp = await test_handler.normalize(q, HGVSDupDelModeOption.ALLELE)
        no_variation_check(resp, q)

    # Invalid
    invalid_queries = [
        "NC_000023.10:g.(?_155270550)_155270570del",
        "NC_000023.11:g.(?_18593474)_18671749del"
        "CDKL5  g.(?_18443702)_18671700del",  # 37
        "CDKL5  g.(?_18425585)_18653631del",  # 38
        "CDKL5  g.(?_18425582)_18653500del",  # 38
    ]
    await invalid_query_list_checks(invalid_queries, test_handler)


@pytest.mark.asyncio
async def test_genomic_del6(
    test_handler,
    genomic_del6_cn_var,
    genomic_del6_cx_var,
):
    """Test that genomic deletion works correctly."""
    q = "NC_000006.12:g.133462764_(133464858_?)del"  # 38
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp, genomic_del6_cx_var)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
    )
    assertion_checks(resp, genomic_del6_cn_var)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE)
    assertion_checks(resp, genomic_del6_cx_var)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.ALLELE)
    no_variation_check(resp, q)

    q = "NC_000006.11:g.133783902_(133785996_?)del"  # 37
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
    assertion_checks(resp, genomic_del6_cx_var)

    resp = await test_handler.normalize(
        q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
    )
    assertion_checks(resp, genomic_del6_cn_var)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.COPY_NUMBER_CHANGE)
    assertion_checks(resp, genomic_del6_cx_var)

    resp = await test_handler.normalize(q, HGVSDupDelModeOption.ALLELE)
    no_variation_check(resp, q)

    # Free text
    for q in ["EYA4 g.133462764_(133464858_?)del"]:  # 38
        resp = await test_handler.normalize(q, HGVSDupDelModeOption.DEFAULT)
        assertion_checks(resp, genomic_del6_cx_var)

        resp = await test_handler.normalize(
            q, HGVSDupDelModeOption.COPY_NUMBER_COUNT, baseline_copies=2
        )
        assertion_checks(resp, genomic_del6_cn_var)

        resp = await test_handler.normalize(q, HGVSDupDelModeOption.ALLELE)
        no_variation_check(resp, q)

    # Invalid
    invalid_queries = [
        "NC_000006.11:g.171115069_(171115080_?)del",
        "NC_000006.12:g.170805981_(170805989_?)del"
        "EYA4 g.133561700_(133853270_?)del",  # 37
        "EYA4 g.133561651_(133561708_?)del",  # 37
        "EYA4 g.133240513_(133240600_?)del",  # 38
        "EYA4 g.133240515_(133532130_?)del",  # 38
    ]
    await invalid_query_list_checks(invalid_queries, test_handler)


@pytest.mark.asyncio
async def test_parameters(test_handler):
    """Check that valid and invalid parameters work as intended."""
    resp = await test_handler.normalize("7-140453136-A-T")
    assert resp.variation
    assert resp.warnings == []

    q = "NC_000003.12:g.49531262dup"
    resp = await test_handler.normalize(q)
    assert resp.variation
    assert resp.warnings == []

    resp = await test_handler.normalize(q, hgvs_dup_del_mode=None)
    assert resp.variation
    assert resp.warnings == []

    resp = await test_handler.normalize(
        q, hgvs_dup_del_mode=HGVSDupDelModeOption.COPY_NUMBER_COUNT
    )
    assert resp.variation is None
    assert resp.warnings == ["copy_number_count mode requires `baseline_copies`"]
