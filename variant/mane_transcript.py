"""Module for retrieving MANE transcript."""
from typing import Optional, Tuple, List
from variant.schemas.token_response_schema import ReferenceSequence
from variant.data_sources import CodonTable
from variant.data_sources.codon_table import MULTIPLE_CODONS
import hgvs.parser
import logging


logger = logging.getLogger('variant')
logger.setLevel(logging.DEBUG)


class MANETranscript:
    """Class for retrieving MANE transcripts."""

    def __init__(self, transcript_mappings, amino_acid_cache) -> None:
        """Initialize the MANETranscript class.

        :param TranscriptMappings transcript_mappings: Access to transcript
            accession mappings
        :param AminoAcidCache amino_acid_cache: Access to amino acid codes
            and conversions
        """
        self.hgvs_parser = hgvs.parser.Parser()
        self.transcript_mappings = transcript_mappings
        self.codon_table = CodonTable(amino_acid_cache)

    def p_to_c(self, transcript, token) -> Optional[str]:
        """Convert protein (p.) to annotation to cDNA (c.) annotation.

        :param str transcript: Transcript accession
        :param Token token: Classification token
        :return: HGVS expression at cDNA reference
        """
        if token.reference_sequence != ReferenceSequence.PROTEIN:
            logger.debug(f"{token} does not have a "
                         f"protein reference sequence.")

        if transcript.startswith('NP_'):
            ac = self.transcript_mappings.np_to_nm[transcript.split(':')[0]]
        elif transcript.startswith('ENSP'):
            ac = \
                self.transcript_mappings.ensp_to_enst[transcript.split(':')[0]]
        else:
            return None

        pos = token.position * 3 - 1
        hgvs_c = f"{ac}:c.{pos}"
        ref, alt = self._get_proteins(token)
        ref_aa_codons, alt_aa_codons = self._get_codons(ref, alt)
        change = self._c_change(ref_aa_codons, alt_aa_codons, ref, alt)

        if change:
            hgvs_c += change
        else:
            return None
        return hgvs_c

    def _get_proteins(self, token) -> Tuple[str, str]:
        """Return reference and new amino acid for a given token.

        :param Token token: Classification token
        :return: [Reference amino acid, New amino acid]
        """
        ref = None
        alt = None
        if token.token_type == 'AminoAcidSubstitution':
            ref = token.ref_protein
            alt = token.alt_protein
        return ref, alt

    def _get_codons(self, ref, alt) -> Tuple[List[str], List[str]]:
        """Get codons for given reference and new amino acids.

        :param str ref: Reference amino acid
        :param str alt: New amino acid
        :return: [Codons for reference amino acid, codons for new amino acid]
        """
        ref_aa_codons = self.codon_table.get_codons(ref)
        alt_aa_codons = self.codon_table.get_codons(alt)
        return ref_aa_codons, alt_aa_codons

    def _get_nuc_index(self, ref_aa_codons, alt_aa_codons) -> Tuple[str, int]:
        """Return nucleotide and index that match.

        :param list ref_aa_codons: Codons for reference amino acid
        :param list alt_aa_codons: Codons for new amino acid
        :return: [Matching nucleotide, Index where nucleotides match]
        """
        i = 0
        while i < 3:
            for ref_codon in ref_aa_codons:
                for alt_codon in alt_aa_codons:
                    if ref_codon[i] == alt_codon[i]:
                        return ref_codon[i], i
            i += 1

    def _c_change(self, ref_aa_codons, alt_aa_codons, ref, alt) -> str:
        """Return variant at the cDNA reference.

        :param list ref_aa_codons: Codons for reference amino acid
        :param list alt_aa_codons: Codons for new amino acid
        :param str ref: Reference amino acid
        :param str alt: New amino acid
        :return: Nucleotide change
        """
        if {ref, alt}.intersection(MULTIPLE_CODONS):
            match, i = self._get_nuc_index(ref_aa_codons, alt_aa_codons)
            ref_aa_codons = \
                [codon for codon in ref_aa_codons if codon[i] == match]
            alt_aa_codons = \
                [codon for codon in alt_aa_codons if codon[i] == match]

        i = 0
        while i < 3:
            for ref_codon in ref_aa_codons:
                for alt_codon in alt_aa_codons:
                    ref_nuc = ref_codon[i]
                    alt_nuc = alt_codon[i]
                    if ref_nuc != alt_nuc:
                        return f"{ref_nuc}>{alt_nuc}"
            i += 1
