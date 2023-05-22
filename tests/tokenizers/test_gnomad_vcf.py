"""A module for testing gnomad VCF tokenization."""
import pytest

from variation.schemas.token_response_schema import (
    CoordinateType, GnomadVcfToken, TokenType
)
from variation.tokenizers import GnomadVCF


@pytest.fixture(scope="module")
def tokenizer():
    """Create GnomadVCF Tokenizer test fixture."""
    return GnomadVCF()


def test_matches(tokenizer):
    """Test that valid input returns correct output"""
    # Genomic Substitution
    input_str = "7-140453136-A-T"
    token = tokenizer.match(input_str)
    assert isinstance(token, GnomadVcfToken)
    assert token.token == input_str
    assert token.input_string == input_str
    assert token.chromosome == "7"
    assert token.pos == 140453136
    assert token.ref == "A"
    assert token.alt == "T"
    assert token.coordinate_type == CoordinateType.LINEAR_GENOMIC
    assert token.token_type == TokenType.GNOMAD_VCF

    # Genomic Insertion
    input_str = "chromosome22-6548-G-GCTAG"
    token = tokenizer.match(input_str)
    assert isinstance(token, GnomadVcfToken)
    assert token.token == "22-6548-G-GCTAG"
    assert token.input_string == input_str
    assert token.chromosome == "22"
    assert token.pos == 6548
    assert token.ref == "G"
    assert token.alt == "GCTAG"
    assert token.coordinate_type == CoordinateType.LINEAR_GENOMIC
    assert token.token_type == TokenType.GNOMAD_VCF

    # Genomic Deletion
    input_str = "x-1234-ac-a"
    token = tokenizer.match(input_str)
    assert isinstance(token, GnomadVcfToken)
    assert token.token == "X-1234-AC-A"
    assert token.input_string == input_str
    assert token.chromosome == "X"
    assert token.pos == 1234
    assert token.ref == "AC"
    assert token.alt == "A"
    assert token.coordinate_type == CoordinateType.LINEAR_GENOMIC
    assert token.token_type == TokenType.GNOMAD_VCF

    # Genomic Reference Agree
    input_str = "chrY-451-C-C"
    token = tokenizer.match(input_str)
    assert isinstance(token, GnomadVcfToken)
    assert token.token == "Y-451-C-C"
    assert token.input_string == input_str
    assert token.chromosome == "Y"
    assert token.pos == 451
    assert token.ref == "C"
    assert token.alt == "C"
    assert token.coordinate_type == CoordinateType.LINEAR_GENOMIC
    assert token.token_type == TokenType.GNOMAD_VCF


def test_no_matches(tokenizer):
    """Test that invalid input returns no tokens"""
    assert tokenizer.match("") is None
    assert tokenizer.match("7-A-C-T") is None
    assert tokenizer.match("7-131-G-B") is None
    assert tokenizer.match("7-131-S-A") is None
    assert tokenizer.match("23-24-A-C") is None
