"""Module for testing classifiers"""
import pytest
import yaml

from tests import PROJECT_ROOT
from variation.schemas.classification_response_schema import (
    AmplificationClassification,
    CdnaDeletionClassification,
    CdnaDelInsClassification,
    CdnaInsertionClassification,
    CdnaReferenceAgreeClassification,
    CdnaSubstitutionClassification,
    GenomicDeletionAmbiguousClassification,
    GenomicDeletionClassification,
    GenomicDelInsClassification,
    GenomicDuplicationAmbiguousClassification,
    GenomicDuplicationClassification,
    GenomicInsertionClassification,
    GenomicReferenceAgreeClassification,
    GenomicSubstitutionClassification,
    ProteinDeletionClassification,
    ProteinDelInsClassification,
    ProteinInsertionClassification,
    ProteinReferenceAgreeClassification,
    ProteinStopGainClassification,
    ProteinSubstitutionClassification,
)


@pytest.fixture(scope="module")
def all_fixtures():
    """Create fixture for classifiers"""
    with open(f"{PROJECT_ROOT}/tests/fixtures/classifiers.yml") as stream:
        return yaml.safe_load(stream)


def classifier_checks(
    all_fixtures, test_tokenizer, test_classifier, fixture_name, expected_classification
):
    """Ensure that fixtures exist for fixture name and that classifier response matches
    expected
    """
    fixtures = all_fixtures.get(
        fixture_name, {"should_match": [], "should_not_match": []}
    )

    for label in ["should_match", "should_not_match"]:
        assert fixtures[label], f"{fixture_name} has no {label} queries"

        for x in fixtures[label]:
            query = x["query"]
            tokens = test_tokenizer.perform(query, [])
            classification = test_classifier.perform(tokens)

            if label == "should_match":
                assert isinstance(classification, expected_classification), query
            else:
                assert classification is None, query


def test_amplification(all_fixtures, test_tokenizer, test_classifier):
    """Test that amplification classifier works"""
    fixture_name = "amplification"
    expected_classification = AmplificationClassification
    classifier_checks(
        all_fixtures,
        test_tokenizer,
        test_classifier,
        fixture_name,
        expected_classification,
    )


def test_protein_substitution(all_fixtures, test_tokenizer, test_classifier):
    """Test that protein substitution classifier works"""
    fixture_name = "protein_substitution"
    expected_classification = ProteinSubstitutionClassification
    classifier_checks(
        all_fixtures,
        test_tokenizer,
        test_classifier,
        fixture_name,
        expected_classification,
    )


def test_cdna_substitution(all_fixtures, test_tokenizer, test_classifier):
    """Test that cdna substitution classifier works"""
    fixture_name = "cdna_substitution"
    expected_classification = CdnaSubstitutionClassification
    classifier_checks(
        all_fixtures,
        test_tokenizer,
        test_classifier,
        fixture_name,
        expected_classification,
    )


def test_genomic_substitution(all_fixtures, test_tokenizer, test_classifier):
    """Test that genomic substitution classifier works"""
    fixture_name = "genomic_substitution"
    expected_classification = GenomicSubstitutionClassification
    classifier_checks(
        all_fixtures,
        test_tokenizer,
        test_classifier,
        fixture_name,
        expected_classification,
    )


def test_protein_stop_gain(all_fixtures, test_tokenizer, test_classifier):
    """Test that protein stop gain classifier works"""
    fixture_name = "protein_stop_gain"
    expected_classification = ProteinStopGainClassification
    classifier_checks(
        all_fixtures,
        test_tokenizer,
        test_classifier,
        fixture_name,
        expected_classification,
    )


def test_protein_reference_agree(all_fixtures, test_tokenizer, test_classifier):
    """Test that protein reference agree classifier works"""
    fixture_name = "protein_reference_agree"
    expected_classification = ProteinReferenceAgreeClassification
    classifier_checks(
        all_fixtures,
        test_tokenizer,
        test_classifier,
        fixture_name,
        expected_classification,
    )


def test_cdna_reference_agree(all_fixtures, test_tokenizer, test_classifier):
    """Test that cdna reference agree classifier works"""
    fixture_name = "cdna_reference_agree"
    expected_classification = CdnaReferenceAgreeClassification
    classifier_checks(
        all_fixtures,
        test_tokenizer,
        test_classifier,
        fixture_name,
        expected_classification,
    )


def test_genomic_reference_agree(all_fixtures, test_tokenizer, test_classifier):
    """Test that genomic reference agree classifier works"""
    fixture_name = "genomic_reference_agree"
    expected_classification = GenomicReferenceAgreeClassification
    classifier_checks(
        all_fixtures,
        test_tokenizer,
        test_classifier,
        fixture_name,
        expected_classification,
    )


def test_protein_delins(all_fixtures, test_tokenizer, test_classifier):
    """Test that protein delins classifier works"""
    fixture_name = "protein_delins"
    expected_classification = ProteinDelInsClassification
    classifier_checks(
        all_fixtures,
        test_tokenizer,
        test_classifier,
        fixture_name,
        expected_classification,
    )


def test_cdna_delins(all_fixtures, test_tokenizer, test_classifier):
    """Test that cdna delins classifier works"""
    fixture_name = "cdna_delins"
    expected_classification = CdnaDelInsClassification
    classifier_checks(
        all_fixtures,
        test_tokenizer,
        test_classifier,
        fixture_name,
        expected_classification,
    )


def test_genomic_delins(all_fixtures, test_tokenizer, test_classifier):
    """Test that genomic delins classifier works"""
    fixture_name = "genomic_delins"
    expected_classification = GenomicDelInsClassification
    classifier_checks(
        all_fixtures,
        test_tokenizer,
        test_classifier,
        fixture_name,
        expected_classification,
    )


def test_protein_deletion(all_fixtures, test_tokenizer, test_classifier):
    """Test that protein deletion classifier works"""
    fixture_name = "protein_deletion"
    expected_classification = ProteinDeletionClassification
    classifier_checks(
        all_fixtures,
        test_tokenizer,
        test_classifier,
        fixture_name,
        expected_classification,
    )


def test_cdna_deletion(all_fixtures, test_tokenizer, test_classifier):
    """Test that cdna deletion classifier works"""
    fixture_name = "cdna_deletion"
    expected_classification = CdnaDeletionClassification
    classifier_checks(
        all_fixtures,
        test_tokenizer,
        test_classifier,
        fixture_name,
        expected_classification,
    )


def test_genomic_deletion(all_fixtures, test_tokenizer, test_classifier):
    """Test that genomic deletion classifier works"""
    fixture_name = "genomic_deletion"
    expected_classification = GenomicDeletionClassification
    classifier_checks(
        all_fixtures,
        test_tokenizer,
        test_classifier,
        fixture_name,
        expected_classification,
    )


def test_genomic_deletion_ambiguous(all_fixtures, test_tokenizer, test_classifier):
    """Test that genomic deletion ambiguous classifier works"""
    fixture_name = "genomic_deletion_ambiguous"
    expected_classification = GenomicDeletionAmbiguousClassification
    classifier_checks(
        all_fixtures,
        test_tokenizer,
        test_classifier,
        fixture_name,
        expected_classification,
    )


def test_protein_insertion(all_fixtures, test_tokenizer, test_classifier):
    """Test that protein insertion classifier works"""
    fixture_name = "protein_insertion"
    expected_classification = ProteinInsertionClassification
    classifier_checks(
        all_fixtures,
        test_tokenizer,
        test_classifier,
        fixture_name,
        expected_classification,
    )


def test_cdna_insertion(all_fixtures, test_tokenizer, test_classifier):
    """Test that cdna insertion classifier works"""
    fixture_name = "cdna_insertion"
    expected_classification = CdnaInsertionClassification
    classifier_checks(
        all_fixtures,
        test_tokenizer,
        test_classifier,
        fixture_name,
        expected_classification,
    )


def test_genomic_insertion(all_fixtures, test_tokenizer, test_classifier):
    """Test that genomic insertion classifier works"""
    fixture_name = "genomic_insertion"
    expected_classification = GenomicInsertionClassification
    classifier_checks(
        all_fixtures,
        test_tokenizer,
        test_classifier,
        fixture_name,
        expected_classification,
    )


def test_genomic_duplication(all_fixtures, test_tokenizer, test_classifier):
    """Test that genomic duplication classifier works"""
    fixture_name = "genomic_duplication"
    expected_classification = GenomicDuplicationClassification
    classifier_checks(
        all_fixtures,
        test_tokenizer,
        test_classifier,
        fixture_name,
        expected_classification,
    )


def test_genomic_duplication_ambiguous(all_fixtures, test_tokenizer, test_classifier):
    """Test that genomic duplication ambiguous classifier works"""
    fixture_name = "genomic_duplication_ambiguous"
    expected_classification = GenomicDuplicationAmbiguousClassification
    classifier_checks(
        all_fixtures,
        test_tokenizer,
        test_classifier,
        fixture_name,
        expected_classification,
    )
