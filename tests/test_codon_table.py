"""Module for testing Codon Table class."""
import pytest
from variation.data_sources.codon_table import CodonTable
from variation.tokenizers.caches import AminoAcidCache


@pytest.fixture(scope='module')
def test_codon_table():
    """Build codon table test fixture."""
    class TestCodonTable:

        def __init__(self):
            self.test_codon_table = CodonTable(AminoAcidCache())

        def _get_codons(self, amino_acid):
            return self.test_codon_table.get_codons(amino_acid)

    return TestCodonTable()


def test_get_codons(test_codon_table):
    """Test that _get_codons method works correctly."""
    ala_codons = ['GCA', 'GCC', 'GCG', 'GCT']
    assert test_codon_table._get_codons('A') == ala_codons
    assert test_codon_table._get_codons('a') == ala_codons
    assert test_codon_table._get_codons('Ala') == ala_codons
    assert test_codon_table._get_codons('ala') == ala_codons
    assert test_codon_table._get_codons('AlA') == ala_codons

    leu_codons = ['CTA', 'CTC', 'CTG', 'CTT', 'TTA', 'TTG']
    assert test_codon_table._get_codons('L') == leu_codons
    assert test_codon_table._get_codons('Leu') == leu_codons

    arg_codons = ['AGA', 'AGG', 'CGA', 'CGC', 'CGG', 'CGT']
    assert test_codon_table._get_codons('r') == arg_codons
    assert test_codon_table._get_codons('Arg') == arg_codons
