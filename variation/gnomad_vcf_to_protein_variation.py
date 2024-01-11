"""Module for translating VCF-like to protein VRS Allele representation"""
from datetime import datetime
from typing import List, Tuple

from cool_seq_tool.handlers import SeqRepoAccess
from cool_seq_tool.mappers import ManeTranscript
from cool_seq_tool.schemas import Strand
from ga4gh.core import ga4gh_identify
from ga4gh.vrs import models, normalize

from variation.classify import Classify
from variation.schemas.classification_response_schema import (
    Nomenclature,
)
from variation.schemas.normalize_response_schema import (
    NormalizeService,
    ServiceMeta,
)
from variation.schemas.token_response_schema import AltType
from variation.schemas.validation_response_schema import ValidationResult
from variation.tokenize import Tokenize
from variation.translate import Translate
from variation.validate import Validate
from variation.version import __version__


class GnomadVcfToProteinError(Exception):
    """Custom exception for Gnomad VCF To Protein"""


def _get_char_match_count(
    min_length: int, ref: str, alt: str, trim_prefix: bool = True
) -> int:
    """Get the count of matched sequential prefixes

    :param min_length: Length of the shortest sequence (using ``ref`` or ``alt``)
    :param ref: Reference sequence
    :param alt: Alternate sequence
    :param trim_prefix: ``True`` if trimming prefixes. ``False`` if trimming suffixes
    :return: The number of sequential characters that were the same in ``ref`` and
        ``alt``
    """
    matched = 0
    num_seq = range(min_length) if trim_prefix else reversed(range(min_length))
    for i in num_seq:
        if alt[i] == ref[i]:
            matched += 1
        else:
            break
    return matched


def _trim_prefix_or_suffix(
    aa_ref: str, aa_alt: str, aa_start_pos: int = 0, trim_prefix: bool = True
) -> Tuple[str, str, int]:
    """Trim prefix or suffix matches

    :param aa_ref: Amino acid reference sequence
    :param aa_alt: Amino acid alternate sequence
    :param aa_start_pos: Amino acid start position. Only required when ``trim_prefix``
        is ``True``
    :param trim_prefix: ``True`` if trimming prefixes. ``False`` if trimming suffixes
    :return: Tuple containing trimmed ``aa_ref``, trimmed `aa_alt``, and
        ``aa_start_pos`` after trimming
    """
    if (aa_ref and aa_alt) and (aa_ref != aa_alt):
        aa_match = 0
        len_aa_ref = len(aa_ref)
        len_aa_alt = len(aa_alt)

        # Trim prefixes
        range_len = len_aa_ref if len_aa_ref < len_aa_alt else len_aa_alt
        aa_match = _get_char_match_count(
            range_len, aa_ref, aa_alt, trim_prefix=trim_prefix
        )
        if aa_match:
            aa_start_pos += aa_match
            aa_alt = aa_alt[aa_match:] if trim_prefix else aa_alt[:-aa_match]
            aa_ref = aa_ref[aa_match:] if trim_prefix else aa_ref[:-aa_ref]

    return aa_ref, aa_alt, aa_start_pos


RNA_CODON_TO_1AA = {
    "AUA": "I",
    "AUC": "I",
    "AUU": "I",
    "AUG": "M",
    "ACA": "T",
    "ACC": "T",
    "ACG": "T",
    "ACU": "T",
    "AAC": "N",
    "AAU": "N",
    "AAA": "K",
    "AAG": "K",
    "AGC": "S",
    "AGU": "S",
    "AGA": "R",
    "AGG": "R",
    "CUA": "L",
    "CUC": "L",
    "CUG": "L",
    "CUU": "L",
    "CCA": "P",
    "CCC": "P",
    "CCG": "P",
    "CCU": "P",
    "CAC": "H",
    "CAU": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGA": "R",
    "CGC": "R",
    "CGG": "R",
    "CGU": "R",
    "GUA": "V",
    "GUC": "V",
    "GUG": "V",
    "GUU": "V",
    "GCA": "A",
    "GCC": "A",
    "GCG": "A",
    "GCU": "A",
    "GAC": "D",
    "GAU": "D",
    "GAA": "E",
    "GAG": "E",
    "GGA": "G",
    "GGC": "G",
    "GGG": "G",
    "GGU": "G",
    "UCA": "S",
    "UCC": "S",
    "UCG": "S",
    "UCU": "S",
    "UUC": "F",
    "UUU": "F",
    "UUA": "L",
    "UUG": "L",
    "UAC": "Y",
    "UAU": "Y",
    "UAA": "*",
    "UAG": "*",
    "UGC": "C",
    "UGU": "C",
    "UGA": "*",
    "UGG": "W",
}


class GnomadVcfToProteinVariation:
    """Class for translating gnomAD VCF-like representation to VRS Allele protein
    representation
    """

    def __init__(
        self,
        seqrepo_access: SeqRepoAccess,
        tokenizer: Tokenize,
        classifier: Classify,
        validator: Validate,
        translator: Translate,
        mane_transcript: ManeTranscript,
    ) -> None:
        """Initialize the GnomadVcfToProteinVariation class

        :param seqrepo_access: Access to SeqRepo
        :param tokenizer: Tokenizer class for tokenizing
        :param classifier: Classifier class for classifying tokens
        :param validator: Validator class for validating valid inputs
        :param translator: Translating valid inputs
        :param mane_transcript: Access MANE Transcript information
        """
        self.seqrepo_access = seqrepo_access
        self.tokenizer = tokenizer
        self.classifier = classifier
        self.validator = validator
        self.translator = translator
        self.mane_transcript = mane_transcript

    async def _get_valid_result(
        self, vcf_query: str, warnings: List
    ) -> List[ValidationResult]:
        """Get gnomad vcf validation summary

        :param vcf_query: gnomad vcf input query
        :param warnings: List of warnings
        :raises GnomadVcfToProteinError: If no tokens, classifications, or valid results
            are found. Also if ``vcf_query`` is not a gnomAD VCF-like query.
        :return: List of valid results for a gnomad VCF query
        """
        tokens = self.tokenizer.perform(vcf_query, warnings)
        if not tokens:
            raise GnomadVcfToProteinError("No tokens found")

        classification = self.classifier.perform(tokens)
        if not classification:
            raise GnomadVcfToProteinError("No classification found")

        if classification.nomenclature != Nomenclature.GNOMAD_VCF:
            raise GnomadVcfToProteinError(
                f"{vcf_query} is not a gnomAD VCF-like query (`chr-pos-ref-alt`)"
            )

        validation_summary = await self.validator.perform(classification)
        valid_results = validation_summary.valid_results
        if valid_results:
            # Temporary work around until issue-490 complete
            valid_results.sort(
                key=lambda v: int(v.accession.split(".")[-1]),
                reverse=True,
            )
            return valid_results[0]
        else:
            raise GnomadVcfToProteinError(
                f"{vcf_query} is not a valid gnomad vcf query"
            )

    @staticmethod
    def _get_alt_type_and_prefix_match(
        len_g_ref: int, len_g_alt: int, g_ref: str, g_alt: str
    ) -> Tuple[AltType, int]:
        """Get genomic alteration type and number of prefixes match

        :param len_g_ref: Length of genomic reference sequence
        :param len_g_alt: Length of genomic alternate sequence
        :param g_ref: Genomic reference sequence
        :param g_alt: Genomic alternate sequence
        :return: Tuple containing genomic alteration type and the number of prefixes matched
        """
        if len_g_ref == len_g_alt:
            num_prefix_matched = _get_char_match_count(
                len_g_ref, g_ref, g_alt, trim_prefix=True
            )
            alt_type = AltType.SUBSTITUTION
        elif len_g_ref > len_g_alt:
            num_prefix_matched = _get_char_match_count(
                len_g_alt, g_ref, g_alt, trim_prefix=True
            )
            alt_type = (
                AltType.DELETION if num_prefix_matched == len_g_alt else AltType.DELINS
            )
        else:
            num_prefix_matched = _get_char_match_count(
                len_g_ref, g_ref, g_alt, trim_prefix=True
            )
            alt_type = (
                AltType.INSERTION if num_prefix_matched == len_g_ref else AltType.DELINS
            )

        return alt_type, num_prefix_matched

    def _get_genomic_pos_range(
        self,
        c_start_pos: int,
        c_end_pos: int,
        strand: Strand,
        g_start_pos: int,
        g_end_pos: int,
    ) -> Tuple[int, int, int]:
        """Get genomic positions to cover the range of codons

        :param c_start_pos: cDNA start position
        :param c_end_pos: cDNA end position
        :param strand: Strand
        :param g_start_pos: Original genomic start position for change
        :param g_end_pos: Original genomic end position for change
        :return: Tuple containing genomic start and end positions and the start index
            for the original position change
        """
        # Get cDNA reading frame
        start_reading_frame = self.mane_transcript._get_reading_frame(c_start_pos + 1)
        end_reading_frame = self.mane_transcript._get_reading_frame(c_end_pos)

        # Get genomic position range change
        # This ensures that there 3 nucleotides needed for codon
        strand = strand
        start_ix = start_reading_frame - 1
        if strand == Strand.NEGATIVE:
            new_g_end_pos = g_end_pos + (start_reading_frame - 1)
            new_g_start_pos = g_start_pos - (3 - end_reading_frame)
        else:
            new_g_start_pos = g_start_pos - (start_reading_frame - 1)
            new_g_end_pos = g_end_pos + (3 - end_reading_frame)
        return new_g_start_pos, new_g_end_pos, start_ix

    def _get_genomic_alt(
        self,
        g_ac: str,
        g_input_alt: str,
        g_end_pos: int,
        alt_type: AltType,
        genomic_start_ix: int,
        strand: Strand,
        ref: str,
    ) -> str:
        """Get entire genomic altered sequence

        :param g_ac: Genomic accession
        :param g_input_alt: Original alteration provided by VCF-like query
        :param g_end_pos: Genomic end position for codon
        :param alt_type: The type of alteration
        :param genomic_start_ix: The start index for the original genomic start position
        :param strand: Strand
        :param ref: The genomic reference sequence
        :return: The updated genomic alteration
        """
        if alt_type == AltType.DELETION:
            alt = ""
        else:
            alt = ref[:genomic_start_ix]

            if strand == Strand.POSITIVE:
                alt += g_input_alt
            else:
                alt += g_input_alt[::-1]

            if alt_type == AltType.SUBSTITUTION:
                alt += ref[len(alt) :]
            else:
                alt += ref[len(ref) - genomic_start_ix :]

            # We need to get the entire inserted sequence. It needs to be a factor of 3
            # since DNA (3 nuc) -> RNA (3 nuc) -> Protein (1 aa). The reason why we
            # DO NOT do this for insertions, is because we only want the provided
            # insertion sequence
            if alt_type != AltType.INSERTION:
                len_alt = len(alt)
                rem_alt = len_alt % 3
                if rem_alt != 0:
                    tmp_g_end_pos = g_end_pos + (3 - rem_alt)
                    tmp_ref, _ = self.seqrepo_access.get_reference_sequence(
                        g_ac, g_end_pos, tmp_g_end_pos
                    )
                    alt += tmp_ref
        return alt

    @staticmethod
    def _dna_to_aa(dna_seq: str, strand: Strand) -> str:
        """Get amino acid(s) from DNA sequence

        :param dna_seq: DNA sequence
        :raises ValueError: If DNA character is not supported
        :return: Amino acid(s)
        """
        # DNA -> RNA
        rna_seq = ""
        if strand == strand.NEGATIVE:
            # Since it's on the negative strand, we need to reverse
            for char in dna_seq:
                if char == "A":
                    rna_seq += "U"
                elif char == "T":
                    rna_seq += "A"
                elif char == "G":
                    rna_seq += "C"
                elif char == "C":
                    rna_seq += "G"
                else:
                    raise ValueError(f"{char} is not a supported nucleotide")
        else:
            # We only need to replace T/U for DNA->RNA
            rna_seq = dna_seq.replace("T", "U")

        # RNA -> 1 letter Amino Acid codes
        aa = ""
        for i in range(int(len(rna_seq) / 3)):
            aa += RNA_CODON_TO_1AA[rna_seq[3 * i : (3 * i) + 3]]
        return aa

    def _get_protein_representation(
        self, ga4gh_seq_id: str, aa_start_pos: int, aa_end_pos: int, aa_alt: str
    ) -> models.Allele:
        """Create VRS Allele for protein representation

        :param ga4gh_seq_id: GA4GH identifier for protein accession
        :param aa_start_pos: Protein start position (inter-residue coordinates)
        :param aa_end_pos: Protein end position (inter-residue coordinates)
        :param aa_alt: Protein alternate sequence
        :raises GnomadVcfToProteinError: If VRS-Python is unable to perform fully
            justified allele normalization
        :return: Normalized VRS Allele on the protein sequence
        """
        variation = models.Allele(
            location=models.SequenceLocation(
                sequenceReference=models.SequenceReference(
                    refgetAccession=ga4gh_seq_id[0].split("ga4gh:")[-1]
                ),
                start=aa_start_pos,
                end=aa_end_pos,
            ),
            state=models.LiteralSequenceExpression(sequence=aa_alt),
        )

        # Perform fully justified allele normalization
        try:
            variation = normalize(variation, self.seqrepo_access)
        except (KeyError, AttributeError) as e:
            raise GnomadVcfToProteinError(f"VRS-Python unable to normalize allele: {e}")

        # Add VRS digests for VRS Allele and VRS Sequence Location
        variation.id = ga4gh_identify(variation)
        variation.location.id = ga4gh_identify(variation.location)
        return variation

    async def gnomad_vcf_to_protein(self, vcf_query: str) -> NormalizeService:
        """Get protein consequence for gnomAD-VCF like expression
        Assumes input query uses GRCh38 representation

        :param vcf_query: gnomAD VCF-like expression (``chr-pos-ref-alt``) on the GRCh38
            assembly
        :return: Normalize Service containing protein VRS Allele, if translation was
            successful
        """
        variation = None
        warnings = []

        # First we need to validate the input query
        try:
            valid_result = await self._get_valid_result(vcf_query, warnings)
        except GnomadVcfToProteinError as e:
            warnings.append(str(e))
            return NormalizeService(
                variation_query=vcf_query,
                variation=variation,
                warnings=warnings,
                service_meta_=ServiceMeta(
                    version=__version__, response_datetime=datetime.now()
                ),
            )

        # Get relevant information from query
        token = valid_result.classification.matching_tokens[0]
        g_ac = valid_result.accession
        g_ref = token.ref
        g_alt = token.alt

        len_g_ref = len(g_ref)
        len_g_alt = len(g_alt)

        # Determine the type of alteration and the number of nucleotide prefixes matched
        alt_type, num_prefix_matched = self._get_alt_type_and_prefix_match(
            len_g_ref, len_g_alt, g_ref, g_alt
        )

        # Get genomic position change
        g_start_pos = token.pos
        if alt_type == AltType.SUBSTITUTION:
            g_end_pos = g_start_pos + (len_g_ref - (num_prefix_matched + 1))
        else:
            g_end_pos = g_start_pos + (len_g_ref - num_prefix_matched)

            if num_prefix_matched == 0:
                g_end_pos -= 1

        # Given genomic data, get associated cDNA and protein representation
        p_c_data = await self.mane_transcript.grch38_to_mane_c_p(
            g_ac, g_start_pos, g_end_pos, try_longest_compatible=True
        )
        if not p_c_data:
            warnings.append("Unable to get cDNA and protein representation")
            return NormalizeService(
                variation_query=vcf_query,
                variation=variation,
                warnings=warnings,
                service_meta_=ServiceMeta(
                    version=__version__, response_datetime=datetime.now()
                ),
            )

        p_data = p_c_data.protein
        c_data = p_c_data.cdna

        # Get GA4GH identifier for protein accession. This is used later, but we want
        # to fail fast
        p_ac = p_data.refseq or p_data.ensembl
        p_ga4gh_seq_id, w = self.seqrepo_access.translate_identifier(p_ac, "ga4gh")
        if w:
            warnings.append(w)
            return NormalizeService(
                variation_query=vcf_query,
                variation=variation,
                warnings=warnings,
                service_meta_=ServiceMeta(
                    version=__version__, response_datetime=datetime.now()
                ),
            )

        # Get genomic position range change
        # This ensures that there 3 nucleotides needed for codon
        strand = c_data.strand
        new_g_start_pos, new_g_end_pos, genomic_start_ix = self._get_genomic_pos_range(
            c_data.pos[0], c_data.pos[1], strand, g_start_pos, g_end_pos
        )

        # Get genomic reference sequence
        ref, w = self.seqrepo_access.get_reference_sequence(
            g_ac, new_g_start_pos, new_g_end_pos
        )
        if w:
            warnings.append(w)
            return NormalizeService(
                variation_query=vcf_query,
                variation=variation,
                warnings=warnings,
                service_meta_=ServiceMeta(
                    version=__version__, response_datetime=datetime.now()
                ),
            )

        if strand == Strand.NEGATIVE:
            ref = ref[::-1]

        # Get genomic altered sequence
        alt = self._get_genomic_alt(
            g_ac, g_alt, new_g_end_pos, alt_type, genomic_start_ix, strand, ref
        )

        # DNA -> RNA -> Protein (1 AA)
        aa_ref = self._dna_to_aa(ref, strand)
        aa_alt = self._dna_to_aa(alt, strand)

        # Trim AA prefixes / suffixes and update the protein start position accordingly
        aa_start_pos = p_data.pos[0]
        aa_ref, aa_alt, aa_start_pos = _trim_prefix_or_suffix(
            aa_ref, aa_alt, aa_start_pos=aa_start_pos, trim_prefix=True
        )
        aa_ref, aa_alt, _ = _trim_prefix_or_suffix(aa_ref, aa_alt, trim_prefix=False)

        # Get protein end position
        if alt_type == AltType.DELETION:
            aa_end_pos = aa_start_pos + (len(aa_ref) - 1)
        else:
            aa_end_pos = p_data.pos[1]

        # Create the protein VRS Allele
        try:
            variation = self._get_protein_representation(
                p_ga4gh_seq_id, aa_start_pos, aa_end_pos, aa_alt
            )
        except GnomadVcfToProteinError as e:
            warnings.append(str(e))

        return NormalizeService(
            variation_query=vcf_query,
            variation=variation,
            warnings=warnings,
            service_meta_=ServiceMeta(
                version=__version__, response_datetime=datetime.now()
            ),
        )
