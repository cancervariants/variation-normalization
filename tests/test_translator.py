"""Module for testing translators"""
import pytest
import yaml

from tests import PROJECT_ROOT
from variation.hgvs_dup_del_mode import HGVSDupDelMode
from variation.translators import (
    Amplification,
    CdnaDeletion,
    CdnaDelIns,
    CdnaInsertion,
    CdnaReferenceAgree,
    CdnaSubstitution,
    GenomicDeletion,
    GenomicDeletionAmbiguous,
    GenomicDelIns,
    GenomicDuplication,
    GenomicDuplicationAmbiguous,
    GenomicInsertion,
    GenomicReferenceAgree,
    GenomicSubstitution,
    ProteinDeletion,
    ProteinDelIns,
    ProteinInsertion,
    ProteinReferenceAgree,
    ProteinStopGain,
    ProteinSubstitution,
)
from variation.validate import Validate
from variation.vrs_representation import VRSRepresentation


@pytest.fixture(scope="module")
def all_fixtures():
    """Create fixture for translators"""
    with open(f"{PROJECT_ROOT}/tests/fixtures/translators.yml") as stream:
        return yaml.safe_load(stream)


@pytest.fixture(scope="module")
def trans_params(test_cool_seq_tool):
    """Create fixture for translator params"""
    vrs_representation = VRSRepresentation(test_cool_seq_tool.seqrepo_access)
    hgvs_dup_del_mode = HGVSDupDelMode(test_cool_seq_tool.seqrepo_access)
    return [
        test_cool_seq_tool.seqrepo_access,
        test_cool_seq_tool.mane_transcript,
        test_cool_seq_tool.uta_db,
        vrs_representation,
        hgvs_dup_del_mode,
    ]


@pytest.fixture(scope="module")
def test_validator(val_params):
    """Create fixture for validate class"""
    return Validate(*val_params)


async def translator_checks(
    all_fixtures,
    test_tokenizer,
    test_classifier,
    test_validator,
    trans_params,
    fixture_name,
    translator_instance,
):
    """Ensure that fixtures exist for fixture name and that translator response matches
    expected
    """
    fixtures = all_fixtures.get(fixture_name, {"tests": []})
    assert fixtures["tests"], f"{fixture_name} has no tests"

    for x in fixtures["tests"]:
        query = x["query"]
        expected = x["variations"]

        tokens = test_tokenizer.perform(query, [])
        classification = test_classifier.perform(tokens)
        validation_summary = await test_validator.perform(classification)
        translations = []
        for vr in validation_summary.valid_results:
            translation_result = await translator_instance(*trans_params).translate(
                vr, []
            )
            vrs_variation = translation_result.vrs_variation
            if vrs_variation and vrs_variation not in translations:
                assert vrs_variation in expected, f"{query}: {vrs_variation['_id']}"
                translations.append(vrs_variation)

        assert len(translations) == len(expected), query


@pytest.mark.asyncio
async def test_protein_substitution(
    all_fixtures, trans_params, test_tokenizer, test_classifier, test_validator
):
    """Test that protein substitution validator works correctly"""
    translator_instance = ProteinSubstitution
    fixture_name = "protein_substitution"
    await translator_checks(
        all_fixtures,
        test_tokenizer,
        test_classifier,
        test_validator,
        trans_params,
        fixture_name,
        translator_instance,
    )


@pytest.mark.asyncio
async def test_cdna_substitution(
    all_fixtures, trans_params, test_tokenizer, test_classifier, test_validator
):
    """Test that cdna substitution validator works correctly"""
    translator_instance = CdnaSubstitution
    fixture_name = "cdna_substitution"
    await translator_checks(
        all_fixtures,
        test_tokenizer,
        test_classifier,
        test_validator,
        trans_params,
        fixture_name,
        translator_instance,
    )


@pytest.mark.asyncio
async def test_genomic_substitution(
    all_fixtures, trans_params, test_tokenizer, test_classifier, test_validator
):
    """Test that genomic substitution validator works correctly"""
    translator_instance = GenomicSubstitution
    fixture_name = "genomic_substitution"
    await translator_checks(
        all_fixtures,
        test_tokenizer,
        test_classifier,
        test_validator,
        trans_params,
        fixture_name,
        translator_instance,
    )


@pytest.mark.asyncio
async def test_protein_stop_gain(
    all_fixtures, trans_params, test_tokenizer, test_classifier, test_validator
):
    """Test that protein stop gain validator works correctly"""
    translator_instance = ProteinStopGain
    fixture_name = "protein_stop_gain"
    await translator_checks(
        all_fixtures,
        test_tokenizer,
        test_classifier,
        test_validator,
        trans_params,
        fixture_name,
        translator_instance,
    )


@pytest.mark.asyncio
async def test_protein_reference_agree(
    all_fixtures, trans_params, test_tokenizer, test_classifier, test_validator
):
    """Test that protein reference agree validator works correctly"""
    translator_instance = ProteinReferenceAgree
    fixture_name = "protein_reference_agree"
    await translator_checks(
        all_fixtures,
        test_tokenizer,
        test_classifier,
        test_validator,
        trans_params,
        fixture_name,
        translator_instance,
    )


@pytest.mark.asyncio
async def test_cdna_reference_agree(
    all_fixtures, trans_params, test_tokenizer, test_classifier, test_validator
):
    """Test that cdna reference agree validator works correctly"""
    translator_instance = CdnaReferenceAgree
    fixture_name = "cdna_reference_agree"
    await translator_checks(
        all_fixtures,
        test_tokenizer,
        test_classifier,
        test_validator,
        trans_params,
        fixture_name,
        translator_instance,
    )


@pytest.mark.asyncio
async def test_genomic_reference_agree(
    all_fixtures, trans_params, test_tokenizer, test_classifier, test_validator
):
    """Test that genomic reference agree validator works correctly"""
    translator_instance = GenomicReferenceAgree
    fixture_name = "genomic_reference_agree"
    await translator_checks(
        all_fixtures,
        test_tokenizer,
        test_classifier,
        test_validator,
        trans_params,
        fixture_name,
        translator_instance,
    )


@pytest.mark.asyncio
async def test_protein_delins(
    all_fixtures, trans_params, test_tokenizer, test_classifier, test_validator
):
    """Test that protein delins validator works correctly"""
    translator_instance = ProteinDelIns
    fixture_name = "protein_delins"
    await translator_checks(
        all_fixtures,
        test_tokenizer,
        test_classifier,
        test_validator,
        trans_params,
        fixture_name,
        translator_instance,
    )


@pytest.mark.asyncio
async def test_cdna_delins(
    all_fixtures, trans_params, test_tokenizer, test_classifier, test_validator
):
    """Test that cdna delins validator works correctly"""
    translator_instance = CdnaDelIns
    fixture_name = "cdna_delins"
    await translator_checks(
        all_fixtures,
        test_tokenizer,
        test_classifier,
        test_validator,
        trans_params,
        fixture_name,
        translator_instance,
    )


@pytest.mark.asyncio
async def test_genomic_delins(
    all_fixtures, trans_params, test_tokenizer, test_classifier, test_validator
):
    """Test that genomic delins validator works correctly"""
    translator_instance = GenomicDelIns
    fixture_name = "genomic_delins"
    await translator_checks(
        all_fixtures,
        test_tokenizer,
        test_classifier,
        test_validator,
        trans_params,
        fixture_name,
        translator_instance,
    )


@pytest.mark.asyncio
async def test_protein_deletion(
    all_fixtures, trans_params, test_tokenizer, test_classifier, test_validator
):
    """Test that protein deletion validator works correctly"""
    translator_instance = ProteinDeletion
    fixture_name = "protein_deletion"
    await translator_checks(
        all_fixtures,
        test_tokenizer,
        test_classifier,
        test_validator,
        trans_params,
        fixture_name,
        translator_instance,
    )


@pytest.mark.asyncio
async def test_cdna_deletion(
    all_fixtures, trans_params, test_tokenizer, test_classifier, test_validator
):
    """Test that protein deletion validator works correctly"""
    translator_instance = CdnaDeletion
    fixture_name = "cdna_deletion"
    await translator_checks(
        all_fixtures,
        test_tokenizer,
        test_classifier,
        test_validator,
        trans_params,
        fixture_name,
        translator_instance,
    )


@pytest.mark.asyncio
async def test_genomic_deletion(
    all_fixtures, trans_params, test_tokenizer, test_classifier, test_validator
):
    """Test that genomic deletion validator works correctly"""
    translator_instance = GenomicDeletion
    fixture_name = "genomic_deletion"
    await translator_checks(
        all_fixtures,
        test_tokenizer,
        test_classifier,
        test_validator,
        trans_params,
        fixture_name,
        translator_instance,
    )


@pytest.mark.asyncio
async def test_genomic_deletion_ambiguous(
    all_fixtures, trans_params, test_tokenizer, test_classifier, test_validator
):
    """Test that genomic deletion ambiguous validator works correctly"""
    translator_instance = GenomicDeletionAmbiguous
    fixture_name = "genomic_deletion_ambiguous"
    await translator_checks(
        all_fixtures,
        test_tokenizer,
        test_classifier,
        test_validator,
        trans_params,
        fixture_name,
        translator_instance,
    )


@pytest.mark.asyncio
async def test_protein_insertion(
    all_fixtures, trans_params, test_tokenizer, test_classifier, test_validator
):
    """Test that protein insertion validator works correctly"""
    translator_instance = ProteinInsertion
    fixture_name = "protein_insertion"
    await translator_checks(
        all_fixtures,
        test_tokenizer,
        test_classifier,
        test_validator,
        trans_params,
        fixture_name,
        translator_instance,
    )


@pytest.mark.asyncio
async def test_cdna_insertion(
    all_fixtures, trans_params, test_tokenizer, test_classifier, test_validator
):
    """Test that cdna insertion validator works correctly"""
    translator_instance = CdnaInsertion
    fixture_name = "cdna_insertion"
    await translator_checks(
        all_fixtures,
        test_tokenizer,
        test_classifier,
        test_validator,
        trans_params,
        fixture_name,
        translator_instance,
    )


@pytest.mark.asyncio
async def test_genomic_insertion(
    all_fixtures, trans_params, test_tokenizer, test_classifier, test_validator
):
    """Test that genomic insertion validator works correctly"""
    translator_instance = GenomicInsertion
    fixture_name = "genomic_insertion"
    await translator_checks(
        all_fixtures,
        test_tokenizer,
        test_classifier,
        test_validator,
        trans_params,
        fixture_name,
        translator_instance,
    )


@pytest.mark.asyncio
async def test_genomic_duplication(
    all_fixtures, trans_params, test_tokenizer, test_classifier, test_validator
):
    """Test that genomic duplication validator works correctly"""
    translator_instance = GenomicDuplication
    fixture_name = "genomic_duplication"
    await translator_checks(
        all_fixtures,
        test_tokenizer,
        test_classifier,
        test_validator,
        trans_params,
        fixture_name,
        translator_instance,
    )


@pytest.mark.asyncio
async def test_genomic_duplication_ambiguous(
    all_fixtures, trans_params, test_tokenizer, test_classifier, test_validator
):
    """Test that genomic duplication ambiguous validator works correctly"""
    translator_instance = GenomicDuplicationAmbiguous
    fixture_name = "genomic_duplication_ambiguous"
    await translator_checks(
        all_fixtures,
        test_tokenizer,
        test_classifier,
        test_validator,
        trans_params,
        fixture_name,
        translator_instance,
    )


@pytest.mark.asyncio
async def test_amplification(
    all_fixtures, trans_params, test_tokenizer, test_classifier, test_validator
):
    """Test that amplification validator works correctly"""
    translator_instance = Amplification
    fixture_name = "amplification"
    await translator_checks(
        all_fixtures,
        test_tokenizer,
        test_classifier,
        test_validator,
        trans_params,
        fixture_name,
        translator_instance,
    )
