"""Module for testing tokenizers"""
import pytest
import yaml

from tests import PROJECT_ROOT
from variation.schemas.token_response_schema import (
    AmplificationToken,
    CdnaDeletionToken,
    CdnaDelInsToken,
    CdnaInsertionToken,
    CdnaReferenceAgreeToken,
    CdnaSubstitutionToken,
    GenomicDeletionAmbiguousToken,
    GenomicDeletionToken,
    GenomicDelInsToken,
    GenomicDuplicationAmbiguousToken,
    GenomicDuplicationToken,
    GenomicInsertionToken,
    GenomicReferenceAgreeToken,
    GenomicSubstitutionToken,
    ProteinDeletionToken,
    ProteinDelInsToken,
    ProteinInsertionToken,
    ProteinReferenceAgreeToken,
    ProteinStopGainToken,
    ProteinSubstitutionToken,
)
from variation.tokenizers import (
    CdnaDeletion,
    CdnaDelIns,
    CdnaGenomicReferenceAgree,
    CdnaInsertion,
    CdnaSubstitution,
    FreeTextCategorical,
    GenomicDeletion,
    GenomicDelIns,
    GenomicDuplication,
    GenomicInsertion,
    GenomicSubstitution,
    ProteinDeletion,
    ProteinDelIns,
    ProteinInsertion,
    ProteinReferenceAgree,
    ProteinSubstitution,
)


@pytest.fixture(scope="module")
def all_fixtures():
    """Create fixture for tokenizers"""
    with open(f"{PROJECT_ROOT}/tests/fixtures/tokenizers.yml") as stream:
        return yaml.safe_load(stream)


def tokenizer_checks(all_fixtures, fixture_name, tokenizer_instance, expected_token):
    """Ensure that fixtures exist for fixture name and that tokenizer response matches
    expected
    """
    labels = ["should_match", "should_not_match"]
    fixtures = all_fixtures.get(fixture_name, {labels[0]: [], labels[1]: []})

    for label in labels:
        assert fixtures[label], f"{fixture_name} has no {label} queries"

        for x in fixtures[label]:
            query = x["token"]
            token = tokenizer_instance().match(query)

            if label == "should_match":
                assert isinstance(token, expected_token), query
            else:
                assert not isinstance(token, expected_token), query


def test_amplification(all_fixtures):
    """Test that amplification tokenizer works"""
    fixture_name = "amplification"
    tokenizer_instance = FreeTextCategorical
    expected_token = AmplificationToken
    tokenizer_checks(all_fixtures, fixture_name, tokenizer_instance, expected_token)


def test_protein_substitution(all_fixtures):
    """Test that protein substitution tokenizer works"""
    fixture_name = "protein_substitution"
    tokenizer_instance = ProteinSubstitution
    expected_token = ProteinSubstitutionToken
    tokenizer_checks(all_fixtures, fixture_name, tokenizer_instance, expected_token)


def test_cdna_substitution(all_fixtures):
    """Test that cdna substitution tokenizer works"""
    fixture_name = "cdna_substitution"
    tokenizer_instance = CdnaSubstitution
    expected_token = CdnaSubstitutionToken
    tokenizer_checks(all_fixtures, fixture_name, tokenizer_instance, expected_token)


def test_genomic_substitution(all_fixtures):
    """Test that genomic substitution tokenizer works"""
    fixture_name = "genomic_substitution"
    tokenizer_instance = GenomicSubstitution
    expected_token = GenomicSubstitutionToken
    tokenizer_checks(all_fixtures, fixture_name, tokenizer_instance, expected_token)


def test_protein_stop_gain(all_fixtures):
    """Test that protein stop gain tokenizer works"""
    fixture_name = "protein_stop_gain"
    tokenizer_instance = ProteinSubstitution
    expected_token = ProteinStopGainToken
    tokenizer_checks(all_fixtures, fixture_name, tokenizer_instance, expected_token)


def test_protein_reference_agree(all_fixtures):
    """Test that protein reference agree tokenizer works"""
    fixture_name = "protein_reference_agree"
    tokenizer_instance = ProteinReferenceAgree
    expected_token = ProteinReferenceAgreeToken
    tokenizer_checks(all_fixtures, fixture_name, tokenizer_instance, expected_token)


def test_cdna_reference_agree(all_fixtures):
    """Test that cdna reference agree tokenizer works"""
    fixture_name = "cdna_reference_agree"
    tokenizer_instance = CdnaGenomicReferenceAgree
    expected_token = CdnaReferenceAgreeToken
    tokenizer_checks(all_fixtures, fixture_name, tokenizer_instance, expected_token)


def test_genomic_reference_agree(all_fixtures):
    """Test that genomic reference agree tokenizer works"""
    fixture_name = "genomic_reference_agree"
    tokenizer_instance = CdnaGenomicReferenceAgree
    expected_token = GenomicReferenceAgreeToken
    tokenizer_checks(all_fixtures, fixture_name, tokenizer_instance, expected_token)


def test_protein_delins(all_fixtures):
    """Test that protein delins tokenizer works"""
    fixture_name = "protein_delins"
    tokenizer_instance = ProteinDelIns
    expected_token = ProteinDelInsToken
    tokenizer_checks(all_fixtures, fixture_name, tokenizer_instance, expected_token)


def test_cdna_delins(all_fixtures):
    """Test that cdna delins tokenizer works"""
    fixture_name = "cdna_delins"
    tokenizer_instance = CdnaDelIns
    expected_token = CdnaDelInsToken
    tokenizer_checks(all_fixtures, fixture_name, tokenizer_instance, expected_token)


def test_genomic_delins(all_fixtures):
    """Test that genomic delins tokenizer works"""
    fixture_name = "genomic_delins"
    tokenizer_instance = GenomicDelIns
    expected_token = GenomicDelInsToken
    tokenizer_checks(all_fixtures, fixture_name, tokenizer_instance, expected_token)


def test_protein_deletion(all_fixtures):
    """Test that protein deletion tokenizer works"""
    fixture_name = "protein_deletion"
    tokenizer_instance = ProteinDeletion
    expected_token = ProteinDeletionToken
    tokenizer_checks(all_fixtures, fixture_name, tokenizer_instance, expected_token)


def test_cdna_deletion(all_fixtures):
    """Test that cdna deletion tokenizer works"""
    fixture_name = "cdna_deletion"
    tokenizer_instance = CdnaDeletion
    expected_token = CdnaDeletionToken
    tokenizer_checks(all_fixtures, fixture_name, tokenizer_instance, expected_token)


def test_genomic_deletion(all_fixtures):
    """Test that genomic deletion tokenizer works"""
    fixture_name = "genomic_deletion"
    tokenizer_instance = GenomicDeletion
    expected_token = GenomicDeletionToken
    tokenizer_checks(all_fixtures, fixture_name, tokenizer_instance, expected_token)


def test_genomic_deletion_ambiguous(all_fixtures):
    """Test that genomic deletion ambiguous tokenizer works"""
    fixture_name = "genomic_deletion_ambiguous"
    tokenizer_instance = GenomicDeletion
    expected_token = GenomicDeletionAmbiguousToken
    tokenizer_checks(all_fixtures, fixture_name, tokenizer_instance, expected_token)


def test_protein_insertion(all_fixtures):
    """Test that protein insertion tokenizer works"""
    fixture_name = "protein_insertion"
    tokenizer_instance = ProteinInsertion
    expected_token = ProteinInsertionToken
    tokenizer_checks(all_fixtures, fixture_name, tokenizer_instance, expected_token)


def test_cdna_insertion(all_fixtures):
    """Test that cdna insertion tokenizer works"""
    fixture_name = "cdna_insertion"
    tokenizer_instance = CdnaInsertion
    expected_token = CdnaInsertionToken
    tokenizer_checks(all_fixtures, fixture_name, tokenizer_instance, expected_token)


def test_genomic_insertion(all_fixtures):
    """Test that genomic insertion tokenizer works"""
    fixture_name = "genomic_insertion"
    tokenizer_instance = GenomicInsertion
    expected_token = GenomicInsertionToken
    tokenizer_checks(all_fixtures, fixture_name, tokenizer_instance, expected_token)


def test_genomic_duplication(all_fixtures):
    """Test that genomic duplication tokenizer works"""
    fixture_name = "genomic_duplication"
    tokenizer_instance = GenomicDuplication
    expected_token = GenomicDuplicationToken
    tokenizer_checks(all_fixtures, fixture_name, tokenizer_instance, expected_token)


def test_genomic_duplication_ambiguous(all_fixtures):
    """Test that genomic duplication ambiguous tokenizer works"""
    fixture_name = "genomic_duplication_ambiguous"
    tokenizer_instance = GenomicDuplication
    expected_token = GenomicDuplicationAmbiguousToken
    tokenizer_checks(all_fixtures, fixture_name, tokenizer_instance, expected_token)
