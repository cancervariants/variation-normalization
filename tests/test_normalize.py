"""Module for testing the normalize endpoint."""
import copy
from datetime import datetime

import pytest
from ga4gh.vrsatile.pydantic.vrsatile_models import VariationDescriptor

from tests.conftest import assertion_checks
from variation.main import normalize as normalize_get_response
from variation.main import to_vrs as to_vrs_get_response
from variation.schemas.normalize_response_schema import HGVSDupDelModeOption


@pytest.fixture(scope="module")
def test_handler(test_query_handler):
    """Create test fixture for normalize handler"""
    return test_query_handler.normalize_handler


@pytest.fixture(scope="module")
def dis3_p63a():
    """Create DIS3 P63A test fixture."""
    params = {
        "id": "normalize.variation:DIS3%20P63A",
        "type": "VariationDescriptor",
        "variation": {
            "id": "ga4gh:VA.a6x8Rwd4ZkbfpC7ze6sJLPIrTWh-idrJ",
            "location": {
                "id": "ga4gh:SL.ehf9oZcXnqRu2ZK4-yBYfCJ6HUgfozFL",
                "end": {"value": 63, "type": "Number"},
                "start": {"value": 62, "type": "Number"},
                "sequence_id": "ga4gh:SQ.mlWsxfPKINN3o300stAI8oqN5U7P6kEu",
                "type": "SequenceLocation",
            },
            "state": {"sequence": "A", "type": "LiteralSequenceExpression"},
            "type": "Allele",
        },
        "molecule_context": "protein",
        "gene_context": "hgnc:20604",
        "vrs_ref_allele_seq": "P",
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def vhl():
    """Create VHL Tyr185Ter fixture."""
    params = {
        "id": "normalize.variation:NP_000542.1%3Ap.Tyr185Ter",
        "type": "VariationDescriptor",
        "variation": {
            "id": "ga4gh:VA.1QkysTg-l3UDx4vjlf0rPaYfJWPKl40T",
            "location": {
                "id": "ga4gh:SL.MdNTceNgzChb2uKAcXnJwKBW_QCbjStq",
                "end": {"value": 185, "type": "Number"},
                "start": {"value": 184, "type": "Number"},
                "sequence_id": "ga4gh:SQ.z-Oa0pZkJ6GHJHOYM7h5mY_umc0SJzTu",
                "type": "SequenceLocation",
            },
            "state": {"sequence": "*", "type": "LiteralSequenceExpression"},
            "type": "Allele",
        },
        "molecule_context": "protein",
        "vrs_ref_allele_seq": "Y",
        "gene_context": "hgnc:12687",
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def braf_v600e_nucleotide(braf_nuc_value):
    """Create a test fixture for BRAF V600E MANE select nucleotide hgvs."""
    variation = copy.deepcopy(braf_nuc_value)
    variation["id"] = "ga4gh:VA.yZv6UrLSFN_tRphhvO-6sty1cQDXvrYD"
    params = {
        "id": "normalize.variation:NM_004333.4%3Ac.1799T%3EA",
        "type": "VariationDescriptor",
        "variation": variation,
        "molecule_context": "transcript",
        "vrs_ref_allele_seq": "T",
        "gene_context": "hgnc:1097",
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def nm_004448_cdna_delins():
    """Create test fixture for NM_004448.4:c.2326_2327delinsCT."""
    params = {
        "id": "normalize.variation:NM_004448.4%3Ac.2326_2327delinsCT",
        "type": "VariationDescriptor",
        "variation": {
            "id": "ga4gh:VA.8pRTvdPPwD2BiMd0hrZ9jeJrGJjKODCl",
            "location": {
                "id": "ga4gh:SL.leD_McoqLU_NHEFSDAAjotnnRchhwjEV",
                "end": {"value": 2502, "type": "Number"},
                "start": {"value": 2500, "type": "Number"},
                "sequence_id": "ga4gh:SQ.y9b4LVMiCXpZxOg9Xt1NwRtssA03MwWM",
                "type": "SequenceLocation",
            },
            "state": {"sequence": "CT", "type": "LiteralSequenceExpression"},
            "type": "Allele",
        },
        "molecule_context": "transcript",
        "vrs_ref_allele_seq": "GG",
        "gene_context": "hgnc:3430",
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def nm_000551():
    """Create test fixture for NM_000551.4:c.615delinsAA."""
    params = {
        "id": "normalize.variation:temp",
        "type": "VariationDescriptor",
        "variation": {
            "id": "ga4gh:VA.1dJA9Ok8-NQYd34NFH3FEffcupe3D7Af",
            "location": {
                "id": "ga4gh:SL.SE972_512zN0a2MGyou9a3wOUaFvcUvE",
                "end": {"value": 685, "type": "Number"},
                "start": {"value": 684, "type": "Number"},
                "sequence_id": "ga4gh:SQ.xBKOKptLLDr-k4hTyCetvARn16pDS_rW",
                "type": "SequenceLocation",
            },
            "state": {"sequence": "AA", "type": "LiteralSequenceExpression"},
            "type": "Allele",
        },
        "molecule_context": "transcript",
        "vrs_ref_allele_seq": "C",
        "gene_context": "hgnc:12687",
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def braf_nuc_value():
    """Create test fixture for BRAF V600E value on c. coordinate."""
    return {
        "location": {
            "id": "ga4gh:SL.ikSRwwSnAg_AysHDQHwp7ISteYo4gXyS",
            "end": {"value": 2025, "type": "Number"},
            "start": {"value": 2024, "type": "Number"},
            "sequence_id": "ga4gh:SQ.aKMPEJgmlZXt_F6gRY5cUG3THH2n-GUa",
            "type": "SequenceLocation",
        },
        "state": {"sequence": "A", "type": "LiteralSequenceExpression"},
        "type": "Allele",
    }


@pytest.fixture(scope="module")
def cdna_reference_agree(braf_nuc_value):
    """Create test fixture for NM_004333.4:c.1799=."""
    value = copy.deepcopy(braf_nuc_value)
    value["state"]["sequence"] = "T"
    value["id"] = "ga4gh:VA.X1IJdIcz3Ftnlj8GSnjXr-zoJjuHQv_m"
    params = {
        "id": "normalize.variation:NM_004333.4%3Ac.1799%3D",
        "type": "VariationDescriptor",
        "variation": value,
        "molecule_context": "transcript",
        "vrs_ref_allele_seq": "T",
        "gene_context": "hgnc:1097",
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def nc_000007_reference_agree(braf_nuc_value):
    """Create test fixture for NC_000007.13:g.140453136=."""
    value = copy.deepcopy(braf_nuc_value)
    value["state"]["sequence"] = "T"
    value["id"] = "ga4gh:VA.X1IJdIcz3Ftnlj8GSnjXr-zoJjuHQv_m"
    params = {
        "id": "normalize.variation:NC_000007.13%3Ag.140453136%3D",
        "type": "VariationDescriptor",
        "variation": value,
        "molecule_context": "transcript",
        "vrs_ref_allele_seq": "T",
        "gene_context": "hgnc:1097",
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def protein_delins():
    """Create test fixture for protein delins."""
    params = {
        "id": "normalize.variation:NP_001333827.1%3Ap.Leu747_Thr751delinsPro",
        "type": "VariationDescriptor",
        "variation": {
            "id": "ga4gh:VA.l8IoVM1nO9s8D0o5UBsJS43IjUepddtR",
            "location": {
                "id": "ga4gh:SL.4utThcVn7rgC3ZhNudavmKw42cFnwVym",
                "end": {"value": 751, "type": "Number"},
                "start": {"value": 746, "type": "Number"},
                "sequence_id": "ga4gh:SQ.vyo55F6mA6n2LgN4cagcdRzOuh38V4mE",
                "type": "SequenceLocation",
            },
            "state": {"sequence": "P", "type": "LiteralSequenceExpression"},
            "type": "Allele",
        },
        "molecule_context": "protein",
        "vrs_ref_allele_seq": "LREAT",
        "gene_context": "hgnc:3236",
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def cdna_deletion():
    """Create test fixture for cdna deletion range with deleted
    sequence.
    """
    params = {
        "id": "normalize.variation:NM_004448.3%3Ac.2264_2278delTGAGGGAAAACACAT",
        "type": "VariationDescriptor",
        "variation": {
            "id": "ga4gh:VA.tQXi6GzpBF2RxR7qUHZvUoVibSU4Wb5c",
            "location": {
                "id": "ga4gh:SL._eX9k4gxNFIcH50Gm7QxkdbWjhIlLsXT",
                "end": {"value": 2453, "type": "Number"},
                "start": {"value": 2437, "type": "Number"},
                "sequence_id": "ga4gh:SQ.y9b4LVMiCXpZxOg9Xt1NwRtssA03MwWM",
                "type": "SequenceLocation",
            },
            "state": {"sequence": "T", "type": "LiteralSequenceExpression"},
            "type": "Allele",
        },
        "molecule_context": "transcript",
        "vrs_ref_allele_seq": "TTGAGGGAAAACACAT",
        "gene_context": "hgnc:3430",
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_deletion():
    """Create test fixture for genomic deletion range with deleted sequence.
    (CA915940709)
    """
    params = {
        "id": "normalize.variation:NC_000003.12%3Ag.10146527_10146528del",
        "type": "VariationDescriptor",
        "variation": {
            "id": "ga4gh:VA.AoFRR6KWkw6_YfrGAxGkT9SJsXieSk93",
            "location": {
                "id": "ga4gh:SL.6D7Wbq7XOvga3y-057BKIra4g9RgAFy9",
                "end": {"value": 10146528, "type": "Number"},
                "start": {"value": 10146524, "type": "Number"},
                "sequence_id": "ga4gh:SQ.Zu7h9AggXxhTaGVsy7h_EZSChSZGcmgX",
                "type": "SequenceLocation",
            },
            "state": {"sequence": "CT", "type": "LiteralSequenceExpression"},
            "type": "Allele",
        },
        "molecule_context": "genomic",
        "vrs_ref_allele_seq": "CTCT",
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def cdna_insertion():
    """Create test fixture for coding DNA insertion."""
    params = {
        "id": "normalize.variation:ENST00000331728.9%3Ac.2049_2050insA",
        "type": "VariationDescriptor",
        "variation": {
            "id": "ga4gh:VA.NuvjUrfrwALQWlI0qkWEiv4IngxNZO8f",
            "location": {
                "id": "ga4gh:SL.4WN2rqLPYke5PGXtesOkolvS9bYvN6KY",
                "end": {"value": 2160, "type": "Number"},
                "start": {"value": 2160, "type": "Number"},
                "sequence_id": "ga4gh:SQ.7_mlQyDN-uWH0RlxTQFvFEv6ykd2D-xF",
                "type": "SequenceLocation",
            },
            "state": {"sequence": "A", "type": "LiteralSequenceExpression"},
            "type": "Allele",
        },
        "molecule_context": "transcript",
        "gene_context": "hgnc:6614",
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_insertion():
    """Create a gene insertion test fixture."""
    params = {
        "id": "normalize.variation:NC_000017.10%3Ag.37880993_37880994insGCTTACGTGATG",
        "type": "VariationDescriptor",
        "variation": {
            "id": "ga4gh:VA.mIIHczyLBr_Z6Vr8_ja3VawRrDgOr3R9",
            "location": {
                "id": "ga4gh:SL.wWXXeIyW78C5WtbR-lWMeo3DXs5eKdk7",
                "end": {"value": 2500, "type": "Number"},
                "start": {"value": 2488, "type": "Number"},
                "sequence_id": "ga4gh:SQ.y9b4LVMiCXpZxOg9Xt1NwRtssA03MwWM",
                "type": "SequenceLocation",
            },
            "state": {
                "sequence": "TACGTGATGGCTTACGTGATGGCT",
                "type": "LiteralSequenceExpression",
            },
            "type": "Allele",
        },
        "molecule_context": "transcript",
        "vrs_ref_allele_seq": "TACGTGATGGCT",
        "gene_context": "hgnc:3430",
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_substitution():
    """Create a gene insertion test fixture."""
    params = {
        "id": "normalize.variation:NC_000007.13%3Ag.55249071C%3ET",
        "type": "VariationDescriptor",
        "variation": {
            "id": "ga4gh:VA.zKeMa2Pee-Z9NFDbDNhe2-LWfaH553YQ",
            "location": {
                "id": "ga4gh:SL.1yYNWks0-ne4DiCHOqaMqOQ_ZBjNzl4Y",
                "end": {"value": 2630, "type": "Number"},
                "start": {"value": 2629, "type": "Number"},
                "sequence_id": "ga4gh:SQ.d_QsP29RWJi6bac7GOC9cJ9AO7s_HUMN",
                "type": "SequenceLocation",
            },
            "state": {"sequence": "T", "type": "LiteralSequenceExpression"},
            "type": "Allele",
        },
        "molecule_context": "transcript",
        "vrs_ref_allele_seq": "C",
        "gene_context": "hgnc:3236",
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def gnomad_vcf_genomic_sub_mnv():
    """Create a genomic substitution mnv test fixture for 5-112175770-GGAA-AGAA."""
    params = {
        "id": "normalize.variation:5-112175770-GGAA-AGAA",
        "type": "VariationDescriptor",
        "variation": {
            "id": "ga4gh:VA.vq6AaaH72SnS7JK_z88F1W8zsAh_2FUh",
            "location": {
                "id": "ga4gh:SL.8_XkoE5ZOq89JNqR2m1PhnOWGP1oy-pW",
                "end": {"value": 112840073, "type": "Number"},
                "start": {"value": 112840072, "type": "Number"},
                "sequence_id": "ga4gh:SQ.aUiQCzCPZ2d0csHbMSbh2NzInhonSXwI",
                "type": "SequenceLocation",
            },
            "state": {"sequence": "A", "type": "LiteralSequenceExpression"},
            "type": "Allele",
        },
        "molecule_context": "genomic",
        "vrs_ref_allele_seq": "G",
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_sub_grch38():
    """Create a genomic substitution GRCh38 test fixture."""
    params = {
        "id": "normalize.variation:NC_000007.13%3Ag.55249071C%3ET",
        "type": "VariationDescriptor",
        "variation": {
            "id": "ga4gh:VA.AlXJ7wD4w6XiHdZeGyF2eLTuBj7PdChR",
            "location": {
                "id": "ga4gh:SL.1i5o_fcj4le345xy_o4zu8de3loank56",
                "end": {"value": 55181378, "type": "Number"},
                "start": {"value": 55181377, "type": "Number"},
                "sequence_id": "ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
                "type": "SequenceLocation",
            },
            "state": {"sequence": "T", "type": "LiteralSequenceExpression"},
            "type": "Allele",
        },
        "molecule_context": "genomic",
        "vrs_ref_allele_seq": "C",
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def grch38_braf_genom_sub(braf_v600e_genomic_sub):
    """Create a genomic substitution GRCh38 test fixture for BRAF."""
    params = {
        "id": "normalize.variation:NC_000007.13%3Ag.140453136A%3ET",
        "type": "VariationDescriptor",
        "variation": braf_v600e_genomic_sub,
        "molecule_context": "genomic",
        "vrs_ref_allele_seq": "A",
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def grch38_braf_genom_reference_agree():
    """Create a genomic reference agree GRCh38 test fixture for BRAF."""
    params = {
        "id": "normalize.variation:NC_000007.13%3Ag.140453136%3D",
        "type": "VariationDescriptor",
        "variation": {
            "id": "ga4gh:VA.XsQzkHmPY63Zsvn1NgrgDtA51ESlA_s9",
            "location": {
                "id": "ga4gh:SL.WBLxdkoypnRME6b8tJtlOWqZKU1ruqY1",
                "end": {"value": 140753336, "type": "Number"},
                "start": {"value": 140753335, "type": "Number"},
                "sequence_id": "ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
                "type": "SequenceLocation",
            },
            "state": {"sequence": "A", "type": "LiteralSequenceExpression"},
            "type": "Allele",
        },
        "molecule_context": "genomic",
        "vrs_ref_allele_seq": "A",
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def grch38_genomic_delins1():
    """Create a test fixture for NC_000007.13:g.140453135_140453136delinsAT."""
    params = {
        "id": "normalize.variation:NC_000007.13%3Ag.140453135_140453136delinsAT",
        "type": "VariationDescriptor",
        "variation": {
            "id": "ga4gh:VA.yeG4TlAgUCuYgV2azqbwflAThd-1SqJG",
            "location": {
                "id": "ga4gh:SL.j0odygfrzBUPJPeS2PTeGpuDq8sHLrkO",
                "end": {"value": 140753336, "type": "Number"},
                "start": {"value": 140753334, "type": "Number"},
                "sequence_id": "ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
                "type": "SequenceLocation",
            },
            "state": {"sequence": "AT", "type": "LiteralSequenceExpression"},
            "type": "Allele",
        },
        "molecule_context": "genomic",
        "vrs_ref_allele_seq": "CA",
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def grch38_genomic_delins2():
    """Create a test fixture for NC_000003.12:g.10149938delinsAA."""
    params = {
        "id": "normalize.variation:NC_000003.12%3Ag.10149938delinsAA",
        "type": "VariationDescriptor",
        "variation": {
            "id": "ga4gh:VA.NdPljo96MLPteFfhvZV5LnqrFpdDa1W5",
            "location": {
                "id": "ga4gh:SL.Hr9iLa3Zykj0J4unZC9GL9Cbxerbr0sa",
                "start": {"value": 10149937, "type": "Number"},
                "end": {"value": 10149938, "type": "Number"},
                "sequence_id": "ga4gh:SQ.Zu7h9AggXxhTaGVsy7h_EZSChSZGcmgX",
                "type": "SequenceLocation",
            },
            "state": {"sequence": "AA", "type": "LiteralSequenceExpression"},
            "type": "Allele",
        },
        "molecule_context": "genomic",
        "vrs_ref_allele_seq": "C",
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def genomic_delins_gene():
    """Create a test fixture for BRAF g.140453135_140453136delinsAT (CA16602419)."""
    params = {
        "id": "normalize.variation:BRAF%20g.140453135_140453136delinsAT",
        "type": "VariationDescriptor",
        "variation": {
            "id": "ga4gh:VA.8JSTEalx0qAj07DlwD1QkrhiZ0PqVRVH",
            "location": {
                "id": "ga4gh:SL.XA1l43TqxRVMdiZ7NhBfBj4Y7Vn4x-LH",
                "start": {"value": 2024, "type": "Number"},
                "end": {"value": 2026, "type": "Number"},
                "sequence_id": "ga4gh:SQ.aKMPEJgmlZXt_F6gRY5cUG3THH2n-GUa",
                "type": "SequenceLocation",
            },
            "state": {"sequence": "AT", "type": "LiteralSequenceExpression"},
            "type": "Allele",
        },
        "molecule_context": "transcript",
        "vrs_ref_allele_seq": "TG",
        "gene_context": "hgnc:1097",
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def gnomad_vcf_genomic_delins1():
    """Create a test fixture for 3-37050340-AAAAGCTTTA-GAGGCTTT.

    https://reg.clinicalgenome.org/redmine/projects/registry/genboree_registry/
    allele?hgvsOrDescriptor=NM_000249.3%3Ac.489_498delinsGAGGCTTT
    """
    params = {
        "id": "normalize.variation:3-37050340-AAAAGCTTTA-GAGGCTTT",
        "type": "VariationDescriptor",
        "variation": {
            "id": "ga4gh:VA.LbnRwxx_4Q_-UFCWm6jmIsglwYGQAESg",
            "location": {
                "id": "ga4gh:SL.KIIypD7gQdLFme54xc7N--zKUdD8YqrM",
                "start": {"value": 37008848, "type": "Number"},
                "end": {"value": 37008858, "type": "Number"},
                "sequence_id": "ga4gh:SQ.Zu7h9AggXxhTaGVsy7h_EZSChSZGcmgX",
                "type": "SequenceLocation",
            },
            "state": {"sequence": "GAGGCTTT", "type": "LiteralSequenceExpression"},
            "type": "Allele",
        },
        "molecule_context": "genomic",
        "vrs_ref_allele_seq": "AAAAGCTTTA",
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def gnomad_vcf_genomic_delins2():
    """Create a test fixture for 16-68846036-AG-TGAGTTT (CA396459910)"""
    params = {
        "id": "normalize.variation:16-68846036-AG-TGAGTTT",
        "type": "VariationDescriptor",
        "variation": {
            "id": "ga4gh:VA.VJRr7mheDEV9bhLg-t7_-k2Q6lJSa5-o",
            "location": {
                "id": "ga4gh:SL.Vf1g-9vg6g8K4_AQhzmvwXFTUn2cvdzm",
                "start": {"value": 68812132, "type": "Number"},
                "end": {"value": 68812134, "type": "Number"},
                "sequence_id": "ga4gh:SQ.yC_0RBj3fgBlvgyAuycbzdubtLxq-rE0",
                "type": "SequenceLocation",
            },
            "state": {"sequence": "TGAGTTT", "type": "LiteralSequenceExpression"},
            "type": "Allele",
        },
        "molecule_context": "genomic",
        "vrs_ref_allele_seq": "AG",
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def grch38_genomic_insertion(grch38_genomic_insertion_variation):
    """Create a test fixture for NC_000017.10:g.37880993_37880994insGCTTACGTGATG."""
    params = {
        "id": "normalize.variation:NC_000017.10%3Ag.37880993_37880994insGCTTACGTGATG",
        "type": "VariationDescriptor",
        "variation": grch38_genomic_insertion_variation,
        "molecule_context": "genomic",
        "vrs_ref_allele_seq": "TACGTGATGGCT",
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def gnomad_vcf_genomic_delins3():
    """Create a test fixture for X-70350063-AG-AGGCAGCGCATAAAGCGCATTCTCCG"""
    params = {
        "id": "normalize.variation:X-70350063-AG-AGGCAGCGCATAAAGCGCATTCTCCG",
        "type": "VariationDescriptor",
        "variation": {
            "id": "ga4gh:VA.DqeBlMxo9Z5MkJwBsboqzxcWqQjJfj8U",
            "location": {
                "id": "ga4gh:SL.JooNMvtPGU7DxZCWB4cat1GQQh3CaOUJ",
                "start": {"value": 71130213, "type": "Number"},
                "end": {"value": 71130215, "type": "Number"},
                "sequence_id": "ga4gh:SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
                "type": "SequenceLocation",
            },
            "state": {
                "sequence": "GGCAGCGCATAAAGCGCATTCTCCGG",
                "type": "LiteralSequenceExpression",
            },
            "type": "Allele",
        },
        "molecule_context": "genomic",
        "vrs_ref_allele_seq": "GG",
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def gnomad_vcf_genomic_delins4():
    """Create a test fixture for 1-55509715-AC-A"""
    params = {
        "id": "normalize.variation:1-55509715-AC-A",
        "type": "VariationDescriptor",
        "variation": {
            "id": "ga4gh:VA.FV9ABAdtGeNsHVqRth9ggLAma-KgwqIL",
            "location": {
                "id": "ga4gh:SL.EnaCJXjNspm_Pgx2hlHe2mrCE-9EXwFr",
                "end": {"value": 55044045, "type": "Number"},
                "start": {"value": 55044042, "type": "Number"},
                "sequence_id": "ga4gh:SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
                "type": "SequenceLocation",
            },
            "state": {"sequence": "CC", "type": "LiteralSequenceExpression"},
            "type": "Allele",
        },
        "molecule_context": "genomic",
        "vrs_ref_allele_seq": "CCC",
    }
    return VariationDescriptor(**params)


@pytest.fixture(scope="module")
def gnomad_vcf_genomic_delins5():
    """Create test fixture for 17-7578455-CGCGG-CGCG (CA497925643)"""
    params = {
        "id": "normalize.variation:17-7578455-CGCGG-CGCG",
        "type": "VariationDescriptor",
        "variation": {
            "id": "ga4gh:VA.V2BRiqZMnzwnIp9CaqXxoiMHuvSuUiJi",
            "type": "Allele",
            "location": {
                "id": "ga4gh:SL.3L904-SBWK3Hj-N-QZY6qFnbsap48ICP",
                "type": "SequenceLocation",
                "sequence_id": "ga4gh:SQ.dLZ15tNO1Ur0IcGjwc3Sdi_0A6Yf4zm7",
                "start": {"type": "Number", "value": 7675139},
                "end": {"type": "Number", "value": 7675141},
            },
            "state": {"type": "LiteralSequenceExpression", "sequence": "G"},
        },
        "molecule_context": "genomic",
        "vrs_ref_allele_seq": "GG",
    }
    return VariationDescriptor(**params)


@pytest.mark.asyncio
async def test_protein_substitution(test_handler, braf_v600e, dis3_p63a):
    """Test that protein substitutions normalize correctly."""
    resp = await test_handler.normalize("     BRAF      V600E    ")
    assertion_checks(resp.variation_descriptor, braf_v600e, "BRAF      V600E")

    braf_id = "normalize.variation:BRAF%20V600E"

    resp = await test_handler.normalize("NP_004324.2:p.Val600Glu")
    assert (
        resp.variation_descriptor.id == "normalize.variation:NP_004324.2%3Ap.Val600Glu"
    )
    resp.variation_descriptor.id = braf_id
    assertion_checks(resp.variation_descriptor, braf_v600e, "NP_004324.2:p.Val600Glu")

    resp = await test_handler.normalize("braf V512E")
    assert resp.variation_descriptor.id == "normalize.variation:braf%20V512E"
    resp.variation_descriptor.id = braf_id
    assertion_checks(resp.variation_descriptor, braf_v600e, "braf V512E")

    resp = await test_handler.normalize(" NP_001365404.1:p.Val512Glu  ")
    assert (
        resp.variation_descriptor.id
        == "normalize.variation:NP_001365404.1%3Ap.Val512Glu"
    )
    resp.variation_descriptor.id = braf_id
    assertion_checks(
        resp.variation_descriptor, braf_v600e, "NP_001365404.1:p.Val512Glu"
    )

    resp = await test_handler.normalize("DIS3 P63A")
    assertion_checks(resp.variation_descriptor, dis3_p63a, "DIS3 P63A")


@pytest.mark.asyncio
async def test_polypeptide_truncation(test_handler, vhl):
    """Test that polypeptide truncations normalize correctly."""
    resp = await test_handler.normalize("NP_000542.1:p.Tyr185Ter")
    assertion_checks(resp.variation_descriptor, vhl, "NP_000542.1:p.Tyr185Ter")


@pytest.mark.asyncio
async def test_reference_agree(test_handler, vhl_reference_agree):
    """Test that reference agrees normalize correctly."""
    resp = await test_handler.normalize("NP_000542.1:p.Pro61=")
    assertion_checks(
        resp.variation_descriptor, vhl_reference_agree, "NP_000542.1:p.Pro61="
    )


@pytest.mark.asyncio
async def test_cdna_and_genomic_substitution(
    test_handler,
    braf_v600e_nucleotide,
    genomic_substitution,
    genomic_sub_grch38,
    grch38_braf_genom_sub,
    gnomad_vcf_genomic_sub_mnv,
):
    """Test that cdna and genomic substitutions normalize correctly."""
    resp = await test_handler.normalize("NM_004333.4:c.1799T>A")
    assertion_checks(
        resp.variation_descriptor, braf_v600e_nucleotide, "NM_004333.4:c.1799T>A"
    )

    # MANE transcript
    refseq_id = "normalize.variation:NM_004333.4%3Ac.1799T%3EA"

    resp = await test_handler.normalize("ENST00000288602.10:c.1799T>A")
    assert (
        resp.variation_descriptor.id
        == "normalize.variation:ENST00000288602.10%3Ac.1799T%3EA"
    )
    resp.variation_descriptor.id = refseq_id
    assertion_checks(
        resp.variation_descriptor, braf_v600e_nucleotide, "ENST00000288602.10:c.1799T>A"
    )

    resp = await test_handler.normalize("BRAF V600E c.1799T>A")
    assert (
        resp.variation_descriptor.id == "normalize.variation:BRAF%20V600E%20c.1799T%3EA"
    )
    resp.variation_descriptor.id = refseq_id
    assertion_checks(
        resp.variation_descriptor, braf_v600e_nucleotide, "BRAF V600E c.1799T>A"
    )

    resp = await test_handler.normalize("BRAF V600E (c.1799T>A)")
    assert (
        resp.variation_descriptor.id
        == "normalize.variation:BRAF%20V600E%20%28c.1799T%3EA%29"
    )
    resp.variation_descriptor.id = refseq_id
    assertion_checks(
        resp.variation_descriptor, braf_v600e_nucleotide, "BRAF V600E (c.1799T>A)"
    )

    resp = await test_handler.normalize("BRAF c.1799T>A")
    assert resp.variation_descriptor.id == "normalize.variation:BRAF%20c.1799T%3EA"
    resp.variation_descriptor.id = refseq_id
    assertion_checks(resp.variation_descriptor, braf_v600e_nucleotide, "BRAF c.1799T>A")

    resp = await test_handler.normalize("NC_000007.13:g.140453136A>T")
    assertion_checks(
        resp.variation_descriptor, grch38_braf_genom_sub, "NC_000007.13:g.140453136A>T"
    )

    fixture_id = "normalize.variation:NC_000007.13%3Ag.140453136A%3ET"
    resp = await test_handler.normalize("7-140453136-A-T")  # 37
    assert resp.variation_descriptor.id == "normalize.variation:7-140453136-A-T"
    resp.variation_descriptor.id = fixture_id
    assertion_checks(
        resp.variation_descriptor, grch38_braf_genom_sub, "7-140453136-A-T"
    )

    resp = await test_handler.normalize("7-140753336-A-T")  # 38
    assert resp.variation_descriptor.id == "normalize.variation:7-140753336-A-T"
    resp.variation_descriptor.id = fixture_id
    assertion_checks(
        resp.variation_descriptor, grch38_braf_genom_sub, "7-140753336-A-T"
    )

    resp = await test_handler.normalize("BRAF V600E (g.140453136A>T)")
    assert (
        resp.variation_descriptor.id
        == "normalize.variation:BRAF%20V600E%20%28g.140453136A%3ET%29"
    )
    resp.variation_descriptor.id = refseq_id
    assertion_checks(
        resp.variation_descriptor, braf_v600e_nucleotide, "BRAF V600E (g.140453136A>T)"
    )

    resp = await test_handler.normalize("BRAF g.140453136A>T")
    assert resp.variation_descriptor.id == "normalize.variation:BRAF%20g.140453136A%3ET"
    resp.variation_descriptor.id = refseq_id
    assertion_checks(
        resp.variation_descriptor, braf_v600e_nucleotide, "BRAF g.140453136A>T"
    )

    # More than 1 gene (EGFR and EGFR-AS1)
    resp = await test_handler.normalize("NC_000007.13:g.55249071C>T")
    assertion_checks(
        resp.variation_descriptor, genomic_sub_grch38, "NC_000007.13:g.55249071C>T"
    )

    resp = await test_handler.normalize("EGFR g.55249071C>T")
    assert resp.variation_descriptor.id == "normalize.variation:EGFR%20g.55249071C%3ET"
    resp.variation_descriptor.id = "normalize.variation:NC_000007.13%3Ag.55249071C%3ET"
    assertion_checks(
        resp.variation_descriptor, genomic_substitution, "EGFR g.55249071C>T"
    )

    # MNV genomic substitution (CA009580)
    q = "5-112175770-GGAA-AGAA"
    resp = await test_handler.normalize(q)
    assertion_checks(resp.variation_descriptor, gnomad_vcf_genomic_sub_mnv, q)


@pytest.mark.asyncio
async def test_cdna_reference_agree(test_handler, cdna_reference_agree):
    """Test that cdna Reference Agree normalizes correctly."""
    resp = await test_handler.normalize("NM_004333.4:c.1799= ")
    assertion_checks(
        resp.variation_descriptor, cdna_reference_agree, "NM_004333.4:c.1799="
    )

    fixture_id = "normalize.variation:NM_004333.4%3Ac.1799%3D"

    resp = await test_handler.normalize("ENST00000288602.11:c.1799=")
    assert (
        resp.variation_descriptor.id
        == "normalize.variation:ENST00000288602.11%3Ac.1799%3D"
    )
    resp.variation_descriptor.id = fixture_id
    assertion_checks(
        resp.variation_descriptor, cdna_reference_agree, "ENST00000288602.11:c.1799="
    )

    resp = await test_handler.normalize("BRAF    c.1799=")
    assert resp.variation_descriptor.id == "normalize.variation:BRAF%20c.1799%3D"
    resp.variation_descriptor.id = fixture_id
    assertion_checks(resp.variation_descriptor, cdna_reference_agree, "BRAF    c.1799=")

    resp = await test_handler.normalize("  BRAF  V600E  c.1799=  ")
    assert (
        resp.variation_descriptor.id == "normalize.variation:BRAF%20V600E%20c.1799%3D"
    )
    resp.variation_descriptor.id = fixture_id
    assertion_checks(
        resp.variation_descriptor, cdna_reference_agree, "BRAF  V600E  c.1799="
    )


@pytest.mark.asyncio
async def test_genomic_reference_agree(
    test_handler, nc_000007_reference_agree, grch38_braf_genom_reference_agree
):
    """Test that genomic reference agree normalizes correctly."""
    resp = await test_handler.normalize("NC_000007.13:g.140453136=")
    assertion_checks(
        resp.variation_descriptor,
        grch38_braf_genom_reference_agree,
        "NC_000007.13:g.140453136=",
    )

    fixture_id = "normalize.variation:NC_000007.13%3Ag.140453136%3D"
    resp = await test_handler.normalize("7-140453136-A-A")
    assert resp.variation_descriptor.id == "normalize.variation:7-140453136-A-A"
    resp.variation_descriptor.id = fixture_id
    assertion_checks(
        resp.variation_descriptor, grch38_braf_genom_reference_agree, "7-140453136-A-A"
    )

    resp = await test_handler.normalize("7-140753336-A-A")
    assert resp.variation_descriptor.id == "normalize.variation:7-140753336-A-A"
    resp.variation_descriptor.id = fixture_id
    assertion_checks(
        resp.variation_descriptor, grch38_braf_genom_reference_agree, "7-140753336-A-A"
    )

    q = "7-140753336-ACT-ACT"
    resp = await test_handler.normalize(q)
    assertion_checks(
        resp.variation_descriptor, grch38_braf_genom_reference_agree, q, ignore_id=True
    )

    resp = await test_handler.normalize("BRAF g.140453136=")
    assert resp.variation_descriptor.id == "normalize.variation:BRAF%20g.140453136%3D"
    resp.variation_descriptor.id = fixture_id
    assertion_checks(
        resp.variation_descriptor, nc_000007_reference_agree, "BRAF g.140453136="
    )


@pytest.mark.asyncio
async def test_cdna_delins(test_handler, nm_004448_cdna_delins, nm_000551):
    """Test that cdna DelIns normalizes correctly."""
    resp = await test_handler.normalize("    NM_004448.4:c.2326_2327delinsCT    ")
    assertion_checks(
        resp.variation_descriptor,
        nm_004448_cdna_delins,
        "NM_004448.4:c.2326_2327delinsCT",
    )

    resp = await test_handler.normalize("NM_000551.3:c.615delinsAA")
    nm_000551.id = "normalize.variation:NM_000551.3%3Ac.615delinsAA"
    assertion_checks(resp.variation_descriptor, nm_000551, "NM_000551.3:c.615delinsAA")


@pytest.mark.asyncio
async def test_genomic_delins(
    test_handler,
    grch38_genomic_delins1,
    grch38_genomic_delins2,
    genomic_delins_gene,
    gnomad_vcf_genomic_delins1,
    gnomad_vcf_genomic_delins2,
    gnomad_vcf_genomic_delins3,
    gnomad_vcf_genomic_delins4,
    gnomad_vcf_genomic_delins5,
    genomic_del1_lse,
    genomic_del2_lse,
):
    """Test that Genomic DelIns normalizes correctly."""
    resp = await test_handler.normalize("NC_000007.13:g.140453135_140453136delinsAT")
    assertion_checks(
        resp.variation_descriptor,
        grch38_genomic_delins1,
        "NC_000007.13:g.140453135_140453136delinsAT",
    )

    resp = await test_handler.normalize("NC_000003.12:g.10149938delinsAA")
    assertion_checks(
        resp.variation_descriptor,
        grch38_genomic_delins2,
        "NC_000003.12:g.10149938delinsAA",
    )

    q = "3-10149938-C-AA"
    resp = await test_handler.normalize(q)
    assertion_checks(
        resp.variation_descriptor, grch38_genomic_delins2, q, ignore_id=True
    )

    q = "BRAF g.140453135_140453136delinsAT"
    resp = await test_handler.normalize(q)
    assertion_checks(resp.variation_descriptor, genomic_delins_gene, q)

    # https://reg.clinicalgenome.org/redmine/projects/registry/genboree_registry/
    # allele?hgvsOrDescriptor=NM_000249.3%3Ac.489_498delinsGAGGCTTT
    q = "3-37050340-AAAAGCTTTA-GAGGCTTT"
    resp = await test_handler.normalize(q)
    assertion_checks(resp.variation_descriptor, gnomad_vcf_genomic_delins1, q)

    q = "16-68846036-AG-TGAGTTT"
    resp = await test_handler.normalize(q)
    assertion_checks(resp.variation_descriptor, gnomad_vcf_genomic_delins2, q)

    # NC_000023.10:g.70350063_70350064delinsAGGCAGCGCATAAAGCGCATTCTCCG
    # NC_000023.10:g.70350063_70350064insGGCAGCGCATAAAGCGCATTCTCC
    # https://reg.clinicalgenome.org/redmine/projects/registry/genboree_registry/
    # allele?hgvsOrDescriptor=NC_000023.11%3Ag.71130213_71130214insGGCAGCGCATAAAGCGCATTCTCC noqa: E501
    q = "X-70350063-AG-AGGCAGCGCATAAAGCGCATTCTCCG"
    resp = await test_handler.normalize(q)
    assertion_checks(resp.variation_descriptor, gnomad_vcf_genomic_delins3, q)

    # CA523275412
    q = "1-55509715-AC-A"
    resp = await test_handler.normalize(q)
    assertion_checks(resp.variation_descriptor, gnomad_vcf_genomic_delins4, q)

    # CA497925643
    q = "17-7578455-CGCGG-CGCG"
    resp = await test_handler.normalize(q)
    assertion_checks(resp.variation_descriptor, gnomad_vcf_genomic_delins5, q)

    q = "3-10146594-AATGTTGACGGACAGCCTAT-A"
    resp = await test_handler.normalize(q)
    assertion_checks(resp.variation_descriptor, genomic_del2_lse, q, ignore_id=True)

    q = "3-10188278-AATGTTGACGGACAGCCTAT-A"
    resp = await test_handler.normalize(q)
    assertion_checks(resp.variation_descriptor, genomic_del2_lse, q, ignore_id=True)

    q = "3-10149810-CT-C"  # 38
    resp = await test_handler.normalize(q)
    assertion_checks(resp.variation_descriptor, genomic_del1_lse, q, ignore_id=True)

    # gnomad should always return lse even if provided other hgvs dup del mode option
    resp = await test_handler.normalize(q, HGVSDupDelModeOption.COPY_NUMBER_COUNT)
    assertion_checks(resp.variation_descriptor, genomic_del1_lse, q, ignore_id=True)


@pytest.mark.asyncio
async def test_protein_delins(test_handler, protein_delins):
    """Test that Amnio Acid DelIns normalizes correctly."""
    resp = await test_handler.normalize("NP_001333827.1:p.Leu747_Thr751delinsPro")
    assertion_checks(
        resp.variation_descriptor,
        protein_delins,
        "NP_001333827.1:p.Leu747_Thr751delinsPro",
    )

    resp = await test_handler.normalize("EGFR p.Leu747_Thr751delinsPro")
    assert (
        resp.variation_descriptor.id
        == "normalize.variation:EGFR%20p.Leu747_Thr751delinsPro"
    )
    resp.variation_descriptor.id = (
        "normalize.variation:NP_001333827.1%3Ap.Leu747_Thr751delinsPro"
    )
    assertion_checks(
        resp.variation_descriptor, protein_delins, "EGFR p.Leu747_Thr751delinsPro"
    )

    resp = await test_handler.normalize("EGFR Leu747_Thr751delinsPro")
    assert (
        resp.variation_descriptor.id
        == "normalize.variation:EGFR%20Leu747_Thr751delinsPro"
    )
    resp.variation_descriptor.id = (
        "normalize.variation:NP_001333827.1%3Ap.Leu747_Thr751delinsPro"
    )
    assertion_checks(
        resp.variation_descriptor, protein_delins, "EGFR Leu747_Thr751delinsPro"
    )

    resp = await test_handler.normalize("EGFR L747_T751delinsP")
    assert resp.variation_descriptor.id == "normalize.variation:EGFR%20L747_T751delinsP"
    resp.variation_descriptor.id = (
        "normalize.variation:NP_001333827.1%3Ap.Leu747_Thr751delinsPro"
    )
    assertion_checks(resp.variation_descriptor, protein_delins, "EGFR L747_T751delinsP")


@pytest.mark.asyncio
async def test_protein_deletion(test_handler, protein_deletion_np_range):
    """Test that Protein Deletion normalizes correctly."""
    resp = await test_handler.normalize("NP_004439.2:p.Leu755_Thr759del")
    assertion_checks(
        resp.variation_descriptor,
        protein_deletion_np_range,
        "NP_004439.2:p.Leu755_Thr759del",
    )

    resp = await test_handler.normalize("ERBB2 p.Leu755_Thr759del")
    assert (
        resp.variation_descriptor.id == "normalize.variation:ERBB2%20p.Leu755_Thr759del"
    )
    resp.variation_descriptor.id = (
        "normalize.variation:NP_004439.2%3Ap.Leu755_Thr759del"
    )
    assertion_checks(
        resp.variation_descriptor, protein_deletion_np_range, "ERBB2 p.Leu755_Thr759del"
    )

    resp = await test_handler.normalize("ERBB2 Leu755_Thr759del")
    assert (
        resp.variation_descriptor.id == "normalize.variation:ERBB2%20Leu755_Thr759del"
    )
    resp.variation_descriptor.id = (
        "normalize.variation:NP_004439.2%3Ap.Leu755_Thr759del"
    )
    assertion_checks(
        resp.variation_descriptor, protein_deletion_np_range, "ERBB2 Leu755_Thr759del"
    )

    resp1 = await test_handler.normalize("EGFR L747_T751del")
    resp2 = await test_handler.normalize("EGFR L747_T751delLREAT")
    assert (
        resp1.variation_descriptor.variation.id
        == resp2.variation_descriptor.variation.id
    )

    # incorrect deleted sequence
    resp = await test_handler.normalize("EGFR L747_T751delLREA")
    assert not resp.variation_descriptor


@pytest.mark.asyncio
async def test_cdna_deletion(test_handler, cdna_deletion):
    """Test that cdna deletion normalizes correctly."""
    # https://reg.clinicalgenome.org/redmine/projects/registry/genboree_registry/by_caid?caid=CA645372623  # noqa: E501
    q = "NM_004448.3:c.2264_2278delTGAGGGAAAACACAT"
    resp1 = await test_handler.normalize(q)
    assertion_checks(resp1.variation_descriptor, cdna_deletion, q)

    # incorrected deleted sequence
    resp = await test_handler.normalize("NM_004448.3:c.2264_2278delTGAGGGAAAACACTA")
    assert not resp.variation_descriptor

    resp2 = await test_handler.normalize("NM_004448.3:c.2264_2278del")
    assert (
        resp1.variation_descriptor.variation.id
        == resp2.variation_descriptor.variation.id
    )

    q = "ERBB2 c.2264_2278delTGAGGGAAAACACAT"
    resp = await test_handler.normalize(q)
    assert (
        resp.variation_descriptor.id
        == "normalize.variation:ERBB2%20c.2264_2278delTGAGGGAAAACACAT"
    )
    resp.variation_descriptor.id = (
        "normalize.variation:NM_004448.3%3Ac.2264_2278delTGAGGGAAAACACAT"
    )
    assertion_checks(resp.variation_descriptor, cdna_deletion, q)


@pytest.mark.asyncio
async def test_genomic_deletion(test_handler, genomic_deletion):
    """Test that genomic deletion normalizes correctly"""
    # CA915940709
    q = "NC_000003.12:g.10146527_10146528del"
    resp1 = await test_handler.normalize(q)
    assertion_checks(resp1.variation_descriptor, genomic_deletion, q)

    resp2 = await test_handler.normalize("NC_000003.12:g.10146527_10146528delCT")
    assert (
        resp2.variation_descriptor.variation.id
        == resp1.variation_descriptor.variation.id
    )

    resp3 = await test_handler.normalize("3-10146526-TCT-T")
    assert (
        resp3.variation_descriptor.variation.id
        == resp2.variation_descriptor.variation.id
    )

    # incorrect deleted sequence
    resp = await test_handler.normalize("NC_000003.12:g.10146527_10146528delCC")
    assert not resp.variation_descriptor


@pytest.mark.asyncio
async def test_protein_insertion(test_handler, protein_insertion):
    """Test that protein insertion normalizes correctly."""
    resp = await test_handler.normalize("NP_005219.2:p.Asp770_Asn771insGlyLeu")
    assertion_checks(
        resp.variation_descriptor,
        protein_insertion,
        "NP_005219.2:p.Asp770_Asn771insGlyLeu",
    )

    def change_resp(response):
        fixture_id = "normalize.variation:NP_005219.2%3Ap.Asp770_Asn771insGlyLeu"
        response.id = fixture_id

    resp = await test_handler.normalize("EGFR D770_N771insGL")
    assert resp.variation_descriptor.id == "normalize.variation:EGFR%20D770_N771insGL"
    change_resp(resp.variation_descriptor)
    assertion_checks(
        resp.variation_descriptor, protein_insertion, "EGFR D770_N771insGL"
    )

    resp = await test_handler.normalize("EGFR p.D770_N771insGL")
    assert resp.variation_descriptor.id == "normalize.variation:EGFR%20p.D770_N771insGL"
    change_resp(resp.variation_descriptor)
    assertion_checks(
        resp.variation_descriptor, protein_insertion, "EGFR p.D770_N771insGL"
    )

    resp = await test_handler.normalize("EGFR Asp770_Asn771insGlyLeu")
    assert (
        resp.variation_descriptor.id
        == "normalize.variation:EGFR%20Asp770_Asn771insGlyLeu"
    )
    change_resp(resp.variation_descriptor)
    assertion_checks(
        resp.variation_descriptor, protein_insertion, "EGFR Asp770_Asn771insGlyLeu"
    )

    resp = await test_handler.normalize("EGFR p.Asp770_Asn771insGlyLeu")
    assert (
        resp.variation_descriptor.id
        == "normalize.variation:EGFR%20p.Asp770_Asn771insGlyLeu"
    )
    change_resp(resp.variation_descriptor)
    assertion_checks(
        resp.variation_descriptor, protein_insertion, "EGFR p.Asp770_Asn771insGlyLeu"
    )


@pytest.mark.asyncio
async def test_cdna_insertion(test_handler, cdna_insertion):
    """Test that cdna insertion normalizes correctly."""
    resp = await test_handler.normalize("ENST00000331728.9:c.2049_2050insA")
    assertion_checks(
        resp.variation_descriptor, cdna_insertion, "ENST00000331728.9:c.2049_2050insA"
    )


@pytest.mark.asyncio
async def test_genomic_insertion(
    test_handler, genomic_insertion, grch38_genomic_insertion
):
    """Test that genomic insertion normalizes correctly."""
    resp = await test_handler.normalize(
        "NC_000017.10:g.37880993_37880994insGCTTACGTGATG"
    )
    assertion_checks(
        resp.variation_descriptor,
        grch38_genomic_insertion,
        "NC_000017.10:g.37880993_37880994insGCTTACGTGATG",
    )

    resp = await test_handler.normalize("ERBB2 g.37880993_37880994insGCTTACGTGATG")
    assert (
        resp.variation_descriptor.id
        == "normalize.variation:ERBB2%20g.37880993_37880994insGCTTACGTGATG"
    )
    resp.variation_descriptor.id = "normalize.variation:NC_000017.10%3Ag.37880993_37880994insGCTTACGTGATG"  # noqa: E501
    assertion_checks(
        resp.variation_descriptor,
        genomic_insertion,
        "ERBB2 g.37880993_37880994insGCTTACGTGATG",
    )

    q = "17-37880993-G-GGCTTACGTGATG"
    resp = await test_handler.normalize(q)
    assertion_checks(
        resp.variation_descriptor, grch38_genomic_insertion, q, ignore_id=True
    )


@pytest.mark.asyncio
async def test_amplification(test_handler, braf_amplification, prpf8_amplification):
    """Test that amplification normalizes correctly."""
    q = "BRAF Amplification"
    resp = await test_handler.normalize(q)
    assertion_checks(resp.variation_descriptor, braf_amplification, q)

    # Gene with > 1 sequence location
    q = "PRPF8 AMPLIFICATION"
    resp = await test_handler.normalize(q)
    assertion_checks(resp.variation_descriptor, prpf8_amplification, q)

    # Gene with no location. This should NOT return a variation
    resp = await test_handler.normalize("IFNR amplification")
    assert resp.variation_descriptor is None


@pytest.mark.asyncio
async def test_valid_queries(test_handler):
    """Test that valid queries don"t throw exceptions. Used for queries that
    revealed bugs in service.
    """
    assert await test_handler.normalize("CCND1 Y44D")

    resp = await test_handler.normalize("NC_000002.12:g.73448098_73448100delCTC")
    assert resp
    assert resp.variation_descriptor.variation.state.sequence == "CTC"
    assert (
        resp.variation_descriptor.variation.id
        == "ga4gh:VA.M0oqztcpIsGcwxeOfxR32Q9_2Xvpu4BE"
    )

    # Test ambiguous IUPAC code N
    for q in [
        "NC_000017.10:g.7572948_7572949insTTTTTTTTTNNNNN",
        "NC_000007.13:g.140453136A>N",
        "NC_000007.13:g.140453135_140453136delinsATN",
        "NM_007294.3:c.2902_2903insTCN",
        "NM_004333.4:c.1799T>N",
        "NM_001289937.1:c.2326_2327delinsCTN",
    ]:
        resp = await test_handler.normalize(q)
        assert resp.variation_descriptor, q


@pytest.mark.asyncio
async def test_no_matches(test_handler):
    """Test no matches work correctly."""
    queries = [
        "braf",  # no change
        "braf v600e",  # incorrect case
        "braf v600000932092039e",  # invalid pos
        "NP_000213.1:cp.Leu862=",  # cp is invalid
        "NP_000213.1:cp.Leu862",  # cp is invalid
        "BRAF V600E 33",  # not supported query type
        "NP_004324.2:p.Glu600Val",  # not valid ref
        "NP_004324.2:p.Glu600Gal",  # not valid ref
        "NP_004324.2839:p.Glu600Val",  # not valid accession
        "NP_004324.2:t.Glu600Val",  # t is invalid
        "this:c.54G>H",  # not a valid accession
        "NC_000007.13:g.4T<A",  # invalid hgvs format
        "test",
        "131",
        "braf Z600E",  # invalid ref
        "braf E600Z",  # invalid ref
        "Thr790Met",  # no gene/accession
        "p.Tyr365Ter",  # no accession
        "ERBB2 G776delinsVCZ",
        "NP005219.2:p.Glu746_Thr751delinsValAla",
        "NP_005219.2:p.Glu746Thr751delinsValAla",
        "EGFR L747_L474delinsP",
        "NP_005219.2:p.Glu746_Thr751delinssValAla",
        "EGFR delins",
        "NM_004333.4:c.1799_1800delTGinsAT",
        "NM_173851.3(SLC30A8):c.973C>T%20(p.Arg325Trp)",
        "NG_008212.3:g.5426_5445del",  # NG accessions not supported
        "NC_000010.11-87925523-C-G",  # invalid format
    ]
    for q in queries:
        resp = await test_handler.normalize(q, untranslatable_returns_text=True)
        assert resp.variation_descriptor.type == "VariationDescriptor", q
        assert resp.variation_descriptor.variation.type == "Text", q
        assert resp.variation_descriptor.label == q.strip(), q

    resp = await test_handler.normalize("clinvar:10")
    assert resp.variation_descriptor is None

    resp = await test_handler.normalize("   ")
    assert resp.variation_descriptor is None

    resp = await test_handler.normalize("")
    assert resp.variation_descriptor is None


@pytest.mark.asyncio
async def test_service_meta():
    """Test that service meta info populates correctly."""
    response = await normalize_get_response("BRAF v600e", "default")
    service_meta = response.service_meta_
    assert service_meta.name == "variation-normalizer"
    assert service_meta.version
    assert isinstance(service_meta.response_datetime, datetime)
    assert (
        service_meta.url == "https://github.com/cancervariants/variation-normalization"
    )

    response = await normalize_get_response("this-wont-normalize", "default")
    service_meta = response.service_meta_
    assert service_meta.name == "variation-normalizer"
    assert service_meta.version
    assert isinstance(service_meta.response_datetime, datetime)
    assert (
        service_meta.url == "https://github.com/cancervariants/variation-normalization"
    )

    response = await to_vrs_get_response("this-wont-normalize")
    service_meta = response.service_meta_
    assert service_meta.name == "variation-normalizer"
    assert service_meta.version
    assert isinstance(service_meta.response_datetime, datetime)
    assert (
        service_meta.url == "https://github.com/cancervariants/variation-normalization"
    )
