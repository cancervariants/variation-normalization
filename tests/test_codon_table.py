"""Module for testing Codon Table class."""
import pytest

from variation.data_sources.codon_table import CodonTable


@pytest.fixture(scope="module")
def test_codon_table():
    """Build codon table test fixture."""
    return CodonTable()


def test_dna_to_rna(test_codon_table):
    """Test that dna_to_rna method works correctly."""
    resp = test_codon_table.dna_to_rna("GTA")
    assert resp == "CAU"

    resp = test_codon_table.dna_to_rna("AAGTGACA")
    assert resp == "UUCACUGU"
