"""Module for testing Codon Table class."""
import pytest

from variation.data_sources.codon_table import CodonTable
from variation.tokenizers.caches import AminoAcidCache


@pytest.fixture(scope="module")
def test_codon_table():
    """Build codon table test fixture."""
    return CodonTable(AminoAcidCache())


def test_get_codons(test_codon_table):
    """Test that _get_codons method works correctly."""
    ala_codons = ["GCA", "GCC", "GCG", "GCU"]
    assert test_codon_table.get_codons("A") == ala_codons
    assert test_codon_table.get_codons("a") == ala_codons
    assert test_codon_table.get_codons("Ala") == ala_codons
    assert test_codon_table.get_codons("ala") == ala_codons
    assert test_codon_table.get_codons("AlA") == ala_codons

    leu_codons = ["CUA", "CUC", "CUG", "CUU", "UUA", "UUG"]
    assert test_codon_table.get_codons("L") == leu_codons
    assert test_codon_table.get_codons("Leu") == leu_codons

    arg_codons = ["AGA", "AGG", "CGA", "CGC", "CGG", "CGU"]
    assert test_codon_table.get_codons("r") == arg_codons
    assert test_codon_table.get_codons("Arg") == arg_codons


def test_dna_to_rna(test_codon_table):
    """Test that dna_to_rna method works correctly."""
    resp = test_codon_table.dna_to_rna("GTA")
    assert resp == "CAU"

    resp = test_codon_table.dna_to_rna("AAGTGACA")
    assert resp == "UUCACUGU"
