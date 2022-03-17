"""A module for testing gnomad VCF tokenization."""
import pytest

from variation.schemas.token_response_schema import ChromosomeToken, \
    GenomicDeletionToken, CoordinateType, GenomicSubstitutionToken, \
    GenomicInsertionToken, GenomicSilentMutationToken
from variation.tokenizers import GnomadVCF


@pytest.fixture(scope="module")
def tokenizer():
    """Create GnomadVCF Tokenizer test fixture."""
    return GnomadVCF()


def chr_assertion_checks(token, input_str, token_val, chromosome_val):
    """Chromosome token assertion checks."""
    assert isinstance(token, ChromosomeToken)
    assert token.token == token_val
    assert token.input_string == input_str
    assert token.chromosome == chromosome_val
    assert token.nomenclature == "gnomad_vcf"


def test_matches(tokenizer):
    """Test that valid input returns correct output"""
    # Genomic Substitution
    input_str = "7-140453136-A-T"
    tokens = tokenizer.match(input_str)
    assert tokens
    assert len(tokens) == 2
    chr_assertion_checks(tokens[0], input_str, "7", "chr7")
    token = tokens[1]
    assert isinstance(token, GenomicSubstitutionToken)
    assert token.position == 140453136
    assert token.ref_nucleotide == "A"
    assert token.new_nucleotide == "T"
    assert token.coordinate_type == CoordinateType.LINEAR_GENOMIC
    assert token.token_type == "GenomicSubstitution"
    assert token.alt_type == "substitution"
    assert token.so_id == "SO:0001483"
    assert token.molecule_context == "genomic"
    assert token.nomenclature == "gnomad_vcf"

    # Genomic Insertion
    input_str = "chromosome22-6548-G-GCTAG"
    tokens = tokenizer.match(input_str)
    assert tokens
    assert len(tokens) == 2
    chr_assertion_checks(tokens[0], input_str, "chromosome22", "chr22")
    token = tokens[1]
    assert isinstance(token, GenomicInsertionToken)
    assert token.start_pos_flank == 6548
    assert token.end_pos_flank == 6549
    assert token.inserted_sequence == "CTAG"
    assert token.inserted_sequence2 is None
    assert token.coordinate_type == CoordinateType.LINEAR_GENOMIC
    assert token.token_type == "GenomicInsertion"
    assert token.alt_type == "insertion"
    assert token.so_id == "SO:0000667"
    assert token.molecule_context == "genomic"
    assert token.nomenclature == "gnomad_vcf"

    # Genomic Deletion
    input_str = "x-1234-ac-a"
    tokens = tokenizer.match(input_str)
    assert tokens
    assert len(tokens) == 2
    chr_assertion_checks(tokens[0], input_str, "x", "chrX")
    token = tokens[1]
    assert isinstance(token, GenomicDeletionToken)
    assert token.start_pos_del == 1235
    assert token.end_pos_del == 1235
    assert token.coordinate_type == CoordinateType.LINEAR_GENOMIC
    assert token.deleted_sequence == "C"
    assert token.token_type == "GenomicDeletion"
    assert token.alt_type == "deletion"
    assert token.so_id == "SO:0000159"
    assert token.molecule_context == "genomic"
    assert token.nomenclature == "gnomad_vcf"

    # Genomic Silent Mutation
    input_str = "chrY-451-C-C"
    tokens = tokenizer.match(input_str)
    assert tokens
    assert len(tokens) == 2
    chr_assertion_checks(tokens[0], input_str, "chrY", "chrY")
    token = tokens[1]
    assert isinstance(token, GenomicSilentMutationToken)
    assert token.position == 451
    assert token.ref_nucleotide == "C"
    assert token.new_nucleotide == "="
    assert token.coordinate_type == CoordinateType.LINEAR_GENOMIC
    assert token.token_type == "GenomicSilentMutation"
    assert token.alt_type == "silent_mutation"
    assert token.so_id == "SO:0002073"
    assert token.molecule_context == "genomic"
    assert token.nomenclature == "gnomad_vcf"


def test_no_matches(tokenizer):
    """Test that invalid input returns no tokens"""
    assert tokenizer.match("") is None
    assert tokenizer.match("7-A-C-T") is None
    assert tokenizer.match("7-131-G-B") is None
    assert tokenizer.match("7-131-S-A") is None
    assert tokenizer.match("23-24-A-C") is None
