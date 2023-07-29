"""Module for testing validators"""
import pytest
import yaml

from tests import PROJECT_ROOT
from variation.validators import (
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


@pytest.fixture(scope="module")
def all_fixtures():
    """Create fixture for validators"""
    with open(f"{PROJECT_ROOT}/tests/fixtures/validators.yml") as stream:
        return yaml.safe_load(stream)


async def validator_checks(
    all_fixtures,
    test_tokenizer,
    test_classifier,
    val_params,
    fixture_name,
    validator_instance,
):
    """Ensure that fixtures exist for fixture name and that validator response matches
    expected
    """
    labels = ["should_match", "should_not_match"]
    if fixture_name == "amplification":
        # Amplification are always valid
        labels = labels[:-1]
        fixtures = all_fixtures.get(fixture_name, {labels[0]: []})
    else:
        labels = ["should_match", "should_not_match"]
        fixtures = all_fixtures.get(fixture_name, {labels[0]: [], labels[1]: []})

    for label in labels:
        assert fixtures[label], f"{fixture_name} has no {label} queries"

        for x in fixtures[label]:
            query = x["query"]
            tokens = test_tokenizer.perform(query, [])
            classification = test_classifier.perform(tokens)

            try:
                validation_results = await validator_instance(*val_params).validate(
                    classification
                )
            except Exception as e:
                raise Exception(f"{e}: {query}")
            else:
                validator_instance
                is_valid = False
                for vr in validation_results:
                    if vr.is_valid:
                        is_valid = True
                        break

                assert is_valid if label == "should_match" else not is_valid, query


@pytest.mark.asyncio
async def test_protein_substitution(
    all_fixtures, test_tokenizer, test_classifier, val_params
):
    """Test that protein substitution validator works correctly"""
    fixture_name = "protein_substitution"
    validator_instance = ProteinSubstitution
    await validator_checks(
        all_fixtures,
        test_tokenizer,
        test_classifier,
        val_params,
        fixture_name,
        validator_instance,
    )


@pytest.mark.asyncio
async def test_cdna_substitution(
    all_fixtures, test_tokenizer, test_classifier, val_params
):
    """Test that cdna substitution validator works correctly"""
    fixture_name = "cdna_substitution"
    validator_instance = CdnaSubstitution
    await validator_checks(
        all_fixtures,
        test_tokenizer,
        test_classifier,
        val_params,
        fixture_name,
        validator_instance,
    )


@pytest.mark.asyncio
async def test_genomic_substitution(
    all_fixtures, test_tokenizer, test_classifier, val_params
):
    """Test that genomic substitution validator works correctly"""
    fixture_name = "genomic_substitution"
    validator_instance = GenomicSubstitution
    await validator_checks(
        all_fixtures,
        test_tokenizer,
        test_classifier,
        val_params,
        fixture_name,
        validator_instance,
    )


@pytest.mark.asyncio
async def test_protein_stop_gain(
    all_fixtures, test_tokenizer, test_classifier, val_params
):
    """Test that protein stop gain validator works correctly"""
    fixture_name = "protein_stop_gain"
    validator_instance = ProteinStopGain
    await validator_checks(
        all_fixtures,
        test_tokenizer,
        test_classifier,
        val_params,
        fixture_name,
        validator_instance,
    )


@pytest.mark.asyncio
async def test_protein_reference_agree(
    all_fixtures, test_tokenizer, test_classifier, val_params
):
    """Test that protein reference agree validator works correctly"""
    fixture_name = "protein_reference_agree"
    validator_instance = ProteinReferenceAgree
    await validator_checks(
        all_fixtures,
        test_tokenizer,
        test_classifier,
        val_params,
        fixture_name,
        validator_instance,
    )


@pytest.mark.asyncio
async def test_cdna_reference_agree(
    all_fixtures, test_tokenizer, test_classifier, val_params
):
    """Test that cdna reference agree validator works correctly"""
    fixture_name = "cdna_reference_agree"
    validator_instance = CdnaReferenceAgree
    await validator_checks(
        all_fixtures,
        test_tokenizer,
        test_classifier,
        val_params,
        fixture_name,
        validator_instance,
    )


@pytest.mark.asyncio
async def test_genomic_reference_agree(
    all_fixtures, test_tokenizer, test_classifier, val_params
):
    """Test that genomic reference agree validator works correctly"""
    fixture_name = "genomic_reference_agree"
    validator_instance = GenomicReferenceAgree
    await validator_checks(
        all_fixtures,
        test_tokenizer,
        test_classifier,
        val_params,
        fixture_name,
        validator_instance,
    )


@pytest.mark.asyncio
async def test_protein_delins(
    all_fixtures, test_tokenizer, test_classifier, val_params
):
    """Test that protein delins validator works correctly"""
    fixture_name = "protein_delins"
    validator_instance = ProteinDelIns
    await validator_checks(
        all_fixtures,
        test_tokenizer,
        test_classifier,
        val_params,
        fixture_name,
        validator_instance,
    )


@pytest.mark.asyncio
async def test_cdna_delins(all_fixtures, test_tokenizer, test_classifier, val_params):
    """Test that cdna delins validator works correctly"""
    fixture_name = "cdna_delins"
    validator_instance = CdnaDelIns
    await validator_checks(
        all_fixtures,
        test_tokenizer,
        test_classifier,
        val_params,
        fixture_name,
        validator_instance,
    )


@pytest.mark.asyncio
async def test_genomic_delins(
    all_fixtures, test_tokenizer, test_classifier, val_params
):
    """Test that genomic delins validator works correctly"""
    fixture_name = "genomic_delins"
    validator_instance = GenomicDelIns
    await validator_checks(
        all_fixtures,
        test_tokenizer,
        test_classifier,
        val_params,
        fixture_name,
        validator_instance,
    )


@pytest.mark.asyncio
async def test_protein_deletion(
    all_fixtures, test_tokenizer, test_classifier, val_params
):
    """Test that protein deletion validator works correctly"""
    fixture_name = "protein_deletion"
    validator_instance = ProteinDeletion
    await validator_checks(
        all_fixtures,
        test_tokenizer,
        test_classifier,
        val_params,
        fixture_name,
        validator_instance,
    )


@pytest.mark.asyncio
async def test_cdna_deletion(all_fixtures, test_tokenizer, test_classifier, val_params):
    """Test that cdna deletion validator works correctly"""
    fixture_name = "cdna_deletion"
    validator_instance = CdnaDeletion
    await validator_checks(
        all_fixtures,
        test_tokenizer,
        test_classifier,
        val_params,
        fixture_name,
        validator_instance,
    )


@pytest.mark.asyncio
async def test_genomic_deletion(
    all_fixtures, test_tokenizer, test_classifier, val_params
):
    """Test that genomic deletion validator works correctly"""
    fixture_name = "genomic_deletion"
    validator_instance = GenomicDeletion
    await validator_checks(
        all_fixtures,
        test_tokenizer,
        test_classifier,
        val_params,
        fixture_name,
        validator_instance,
    )


@pytest.mark.asyncio
async def test_genomic_deletion_ambiguous(
    all_fixtures, test_tokenizer, test_classifier, val_params
):
    """Test that genomic deletion ambiguous validator works correctly"""
    fixture_name = "genomic_deletion_ambiguous"
    validator_instance = GenomicDeletionAmbiguous
    await validator_checks(
        all_fixtures,
        test_tokenizer,
        test_classifier,
        val_params,
        fixture_name,
        validator_instance,
    )


@pytest.mark.asyncio
async def test_protein_insertion(
    all_fixtures, test_tokenizer, test_classifier, val_params
):
    """Test that protein insertion validator works correctly"""
    fixture_name = "protein_insertion"
    validator_instance = ProteinInsertion
    await validator_checks(
        all_fixtures,
        test_tokenizer,
        test_classifier,
        val_params,
        fixture_name,
        validator_instance,
    )


@pytest.mark.asyncio
async def test_cdna_insertion(
    all_fixtures, test_tokenizer, test_classifier, val_params
):
    """Test that cdna insertion validator works correctly"""
    fixture_name = "cdna_insertion"
    validator_instance = CdnaInsertion
    await validator_checks(
        all_fixtures,
        test_tokenizer,
        test_classifier,
        val_params,
        fixture_name,
        validator_instance,
    )


@pytest.mark.asyncio
async def test_genomic_insertion(
    all_fixtures, test_tokenizer, test_classifier, val_params
):
    """Test that genomic insertion validator works correctly"""
    fixture_name = "genomic_insertion"
    validator_instance = GenomicInsertion
    await validator_checks(
        all_fixtures,
        test_tokenizer,
        test_classifier,
        val_params,
        fixture_name,
        validator_instance,
    )


@pytest.mark.asyncio
async def test_genomic_duplication(
    all_fixtures, test_tokenizer, test_classifier, val_params
):
    """Test that genomic duplication validator works correctly"""
    fixture_name = "genomic_duplication"
    validator_instance = GenomicDuplication
    await validator_checks(
        all_fixtures,
        test_tokenizer,
        test_classifier,
        val_params,
        fixture_name,
        validator_instance,
    )


@pytest.mark.asyncio
async def test_genomic_duplication_ambiguous(
    all_fixtures, test_tokenizer, test_classifier, val_params
):
    """Test that genomic duplication ambiguous validator works correctly"""
    fixture_name = "genomic_duplication_ambiguous"
    validator_instance = GenomicDuplicationAmbiguous
    await validator_checks(
        all_fixtures,
        test_tokenizer,
        test_classifier,
        val_params,
        fixture_name,
        validator_instance,
    )


@pytest.mark.asyncio
async def test_amplification(all_fixtures, test_tokenizer, test_classifier, val_params):
    """Test that amplification validator works correctly"""
    fixture_name = "amplification"
    validator_instance = Amplification
    await validator_checks(
        all_fixtures,
        test_tokenizer,
        test_classifier,
        val_params,
        fixture_name,
        validator_instance,
    )
