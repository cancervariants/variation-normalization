"""Module for translating gnomAD-VCF to protein VRS Allele representation"""

import datetime
from typing import TYPE_CHECKING, Literal

from cool_seq_tool.handlers import SeqRepoAccess
from cool_seq_tool.mappers import ManeTranscript
from cool_seq_tool.schemas import CoordinateType, Strand
from ga4gh.core import ga4gh_identify
from ga4gh.core.models import MappableConcept
from ga4gh.vrs import models, normalize
from gene.query import QueryHandler as GeneQueryHandler
from gene.schemas import MatchType as GeneMatchType

from variation import __version__
from variation.classify import Classify
from variation.schemas.classification_response_schema import Nomenclature
from variation.schemas.gnomad_vcf_to_protein_schema import GnomadVcfToProteinService
from variation.schemas.normalize_response_schema import ServiceMeta
from variation.schemas.service_schema import ClinVarAssembly
from variation.schemas.token_response_schema import AltType, GnomadVcfToken
from variation.schemas.validation_response_schema import ValidationResult
from variation.tokenize import Tokenize
from variation.translate import Translate
from variation.utils import get_vrs_loc_seq
from variation.validate import Validate

if TYPE_CHECKING:
    from cool_seq_tool.mappers.mane_transcript import (
        CdnaRepresentation,
        DataRepresentation,
        GenomicRepresentation,
        ProteinAndCdnaRepresentation,
    )


class GnomadVcfToProteinError(Exception):
    """Custom exception for gnomAD VCF To Protein specific errors"""


def _get_char_match_count(
    min_length: int, ref: str, alt: str, trim_prefix: bool = True
) -> int:
    """Get the count of matched sequential prefixes

    :param min_length: Length of the shortest sequence (using `ref` or `alt`)
    :param ref: Reference sequence
    :param alt: Alternate sequence
    :param trim_prefix: `True` if trimming prefixes. `False` if trimming suffixes
    :return: The number of sequential characters that were the same in `ref` and `alt`
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
) -> tuple[str, str, int]:
    """Trim prefix or suffix matches

    :param aa_ref: Amino acid reference sequence
    :param aa_alt: Amino acid alternate sequence
    :param aa_start_pos: Amino acid start position. Only required when `trim_prefix` is
        `True`
    :param trim_prefix: `True` if trimming prefixes. `False` if trimming suffixes
    :return: Tuple containing trimmed `aa_ref`, trimmed aa_alt`, and `aa_start_pos`
        after trimming
    """
    if (aa_ref and aa_alt) and (aa_ref != aa_alt):
        aa_match = 0
        len_aa_ref = len(aa_ref)
        len_aa_alt = len(aa_alt)

        # Trim prefixes or suffixes
        range_len = min(len_aa_alt, len_aa_ref)
        aa_match = _get_char_match_count(
            range_len, aa_ref, aa_alt, trim_prefix=trim_prefix
        )
        if aa_match:
            aa_start_pos += aa_match
            aa_alt = aa_alt[aa_match:] if trim_prefix else aa_alt[:-aa_match]
            aa_ref = aa_ref[aa_match:] if trim_prefix else aa_ref[:-aa_match]

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
    """Class for translating gnomAD-VCF representation to VRS Allele protein
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
        gene_normalizer: GeneQueryHandler,
    ) -> None:
        """Initialize the GnomadVcfToProteinVariation class

        :param seqrepo_access: Access to SeqRepo
        :param tokenizer: Tokenizer class for tokenizing
        :param classifier: Classifier class for classifying tokens
        :param validator: Validator class for validating valid inputs
        :param translator: Translating valid inputs
        :param mane_transcript: Access MANE Transcript information
        :param gene_normalizer: Client for normalizing gene concepts
        """
        self.seqrepo_access = seqrepo_access
        self.tokenizer = tokenizer
        self.classifier = classifier
        self.validator = validator
        self.translator = translator
        self.mane_transcript = mane_transcript
        self.gene_normalizer = gene_normalizer

    async def _get_valid_result(
        self,
        vcf_query: str,
        warnings: list,
        input_assembly: Literal[ClinVarAssembly.GRCH37, ClinVarAssembly.GRCH38] | None,
    ) -> ValidationResult:
        """Perform validation steps

        :param vcf_query: gnomAD-VCF input query
        :param warnings: List of warnings. This
        :param input_assembly: Assembly used for `vcf_query`.
        :raises GnomadVcfToProteinError: If no tokens, classifications, or valid results
            are found. Also if `vcf_query` is not a gnomAD-VCF query.
        :return: List of valid results for a gnomAD-VCF query
        """
        tokens = self.tokenizer.perform(vcf_query, warnings)
        if not tokens:
            msg = "No tokens found"
            raise GnomadVcfToProteinError(msg)

        classification = self.classifier.perform(tokens)
        if not classification:
            msg = "No classification found"
            raise GnomadVcfToProteinError(msg)

        if classification.nomenclature != Nomenclature.GNOMAD_VCF:
            msg = f"{vcf_query} is not a gnomAD-VCF query (`chr-pos-ref-alt`)"
            raise GnomadVcfToProteinError(msg)

        validation_summary = await self.validator.perform(
            classification, input_assembly=input_assembly
        )
        valid_results = validation_summary.valid_results
        if valid_results:
            # Always use the latest assembly in valid results
            valid_results.sort(
                key=lambda v: int(v.accession.split(".")[-1]),
                reverse=True,
            )
            return valid_results[0]

        msg = f"{vcf_query} is not a valid gnomAD-VCF query"
        raise GnomadVcfToProteinError(msg)

    @staticmethod
    def _get_genomic_alt_type(
        len_g_ref: int, len_g_alt: int, g_ref: str, g_alt: str
    ) -> AltType:
        """Get genomic alteration type

        :param len_g_ref: Length of genomic reference sequence, `g_ref`
        :param len_g_alt: Length of genomic alternate sequence, `g_alt`
        :param g_ref: Genomic reference sequence
        :param g_alt: Genomic alternate sequence
        :return: The type of genomic variant (substitution, deletion, insertion)
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

        return alt_type

    def _get_codon_aligned_interval(
        self,
        c_start_pos: int,
        c_end_pos: int,
        strand: Strand,
        g_start_pos: int,
        g_end_pos: int,
    ) -> tuple[int, int, int]:
        """Get the codon-aligned interval for a given genomic change

        Given cDNA and the original genomic change interval, expand the genomic interval
        to include the complete codon(s) affected by the original genomic change (
            via gnomAD-VCF input)

        :param c_start_pos: cDNA start position (inter-residue coordinates)
        :param c_end_pos: cDNA end position (inter-residue coordinates)
        :param strand: Strand
        :param g_start_pos: Original genomic start position for change (residue
            coordinates)
        :param g_end_pos: Original genomic end position for change (residue coordinates)
        :return: Tuple containing codon aligned interval start position, codon aligned
            interval end position, and start index for original genomic change
        """
        # Get cDNA reading frames
        start_reading_frame = self.mane_transcript.get_reading_frame(c_start_pos + 1)
        end_reading_frame = self.mane_transcript.get_reading_frame(c_end_pos)

        # Get genomic position interval change
        # This ensures that there 3 nucleotides needed for codon
        genomic_start_ix = start_reading_frame - 1
        if strand == Strand.NEGATIVE:
            codon_aligned_interval_end = g_end_pos + (start_reading_frame - 1)
            codon_aligned_interval_start = g_start_pos - (3 - end_reading_frame)
        else:
            codon_aligned_interval_start = g_start_pos - (start_reading_frame - 1)
            codon_aligned_interval_end = g_end_pos + (3 - end_reading_frame)
        return (
            codon_aligned_interval_start,
            codon_aligned_interval_end,
            genomic_start_ix,
        )

    def _get_codon_aligned_alternate_sequence(
        self,
        g_ac: str,
        g_input_alt: str,
        len_g_ref: int,
        g_end_pos: int,
        alt_type: AltType,
        genomic_start_ix: int,
        strand: Strand,
        codon_aligned_ref_seq: str,
    ) -> str:
        """Build the genomic alteration sequence (or change) within a codon-aligned
        interval.

        The genomic reference sequence, `codon_aligned_ref_seq`, includes the complete
        codon(s) affected by the original genomic change (via gnomAD-VCF input). We
        apply the alteration to this codon-aligned sequence so translation reflects the
        correct in-frame protein sequence consequence.

        :param g_ac: Genomic RefSeq accession
        :param g_input_alt: Original alteration provided by the input gnomAD-VCF
        :param len_g_ref: Length of genomic reference sequence provided by the input
            gnomAD-VCF
        :param g_end_pos: Genomic end position of the codon-aligned interval
            (residue coordinates).
        :param alt_type: The type of alteration (substitution, deletion, insertion,
            delins)
        :param genomic_start_ix: The start index for the original genomic start position
            within `codon_aligned_ref_seq`
        :param strand: Strand
        :param codon_aligned_ref_seq: The codon-aligned genomic reference sequence
            fetched from SeqRepo
        :return: Genomic alteration sequence corresponding to the codon-aligned interval
        """
        # Codon-aligned interval structure:
        #
        #     codon_aligned_ref_seq =
        #         [ prefix ][ reference segment ][ suffix ]
        #                    ^
        #                    genomic_start_ix
        #
        # After applying the alteration:
        #
        #         [ prefix ][ input_alt ][ suffix ]
        #
        # - Substitution: reference segment is replaced
        # - Deletion: reference segment is removed
        # - Insertion: input_alt is inserted at the change position
        # - Delins: treated as replacement with length change

        input_alt = g_input_alt if strand == Strand.POSITIVE else g_input_alt[::-1]

        # Rebuild the codon-aligned sequence around the original genomic change.
        # Keep the sequence before the change, apply the alternate sequence (change),
        # and append the remaining sequence.
        prefix = codon_aligned_ref_seq[:genomic_start_ix]

        if alt_type == AltType.DELETION:
            suffix = codon_aligned_ref_seq[genomic_start_ix + len_g_ref :]
            alt = prefix + input_alt + suffix
        else:
            if alt_type == AltType.SUBSTITUTION:
                suffix = codon_aligned_ref_seq[genomic_start_ix + len_g_ref :]
            else:
                suffix = codon_aligned_ref_seq[
                    len(codon_aligned_ref_seq) - genomic_start_ix :
                ]

            alt = prefix + input_alt + suffix

            # Get the entire inserted sequence to ensure resulting DNA sequence length
            # is codon-aligned.
            # This needs to be a factor of 3 since
            # DNA (3 nucleotides) -> RNA (3 nucleotides) -> Protein (1 amino acid).
            # DO NOT do this for insertions because we only want the provided insertion
            # sequence
            if alt_type != AltType.INSERTION:
                remainder = len(alt) % 3
                if remainder:
                    tmp_g_end_pos = g_end_pos + (3 - remainder)
                    tmp_ref, _ = self.seqrepo_access.get_reference_sequence(
                        g_ac, g_end_pos, tmp_g_end_pos
                    )
                    alt += tmp_ref
        return alt

    @staticmethod
    def _dna_to_aa(dna_seq: str, strand: Strand) -> str:
        """Get amino acid(s) from DNA sequence

        :param dna_seq: DNA sequence
        :param strand: Strand
        :raises ValueError: If DNA character is not supported
        :return: Amino acid(s)
        """
        # DNA -> RNA
        rna_seq = ""
        if strand == strand.NEGATIVE:
            # Since it's on the negative strand, we need to get the nucleic acid complement
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
                    msg = f"{char} is not a supported nucleotide"
                    raise ValueError(msg)
        else:
            # We only need to replace T/U for DNA->RNA
            rna_seq = dna_seq.replace("T", "U")

        # RNA -> 1 letter Amino Acid codes
        aa = ""
        for i in range(int(len(rna_seq) / 3)):
            aa += RNA_CODON_TO_1AA[rna_seq[3 * i : (3 * i) + 3]]
        return aa

    def _get_protein_representation(
        self,
        ga4gh_seq_id: str,
        p_ac: str,
        aa_start_pos: int,
        aa_end_pos: int,
        aa_alt: str,
    ) -> models.Allele:
        """Create VRS Allele for protein representation

        :param ga4gh_seq_id: GA4GH identifier for protein accession (`ga4gh:SQ.`)
        :param p_ac: RefSeq or Ensembl protein accession
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
            msg = f"VRS-Python unable to normalize allele: {e}"
            raise GnomadVcfToProteinError(msg) from e

        loc_seq = get_vrs_loc_seq(
            self.seqrepo_access, p_ac, variation.location.start, variation.location.end
        )
        if loc_seq:
            variation.location.sequence = models.sequenceString(root=loc_seq)

        # Add VRS digests for VRS Allele and VRS Sequence Location
        variation.id = ga4gh_identify(variation)
        variation.location.id = ga4gh_identify(variation.location)
        return variation

    def _get_gene_context(self, gene: str) -> MappableConcept | None:
        """Get additional gene information from gene-normalizer

        :param gene: Gene symbol
        :return: Gene data from gene-normalizer if match found
        """
        gene_norm_resp = self.gene_normalizer.normalize(gene)
        return (
            gene_norm_resp.gene
            if gene_norm_resp.match_type != GeneMatchType.NO_MATCH
            else None
        )

    async def gnomad_vcf_to_protein(
        self,
        vcf_query: str,
        input_assembly: Literal[ClinVarAssembly.GRCH37, ClinVarAssembly.GRCH38]
        | None = None,
    ) -> GnomadVcfToProteinService:
        """Given genomic gnomAD-VCF expression, return associated protein consequence

        Genomic variant -> cDNA variant -> protein variant
        If `input_assembly` is GRCh37, will attempt to liftover to GRCh38

        :param vcf_query: gnomAD-VCF expression of the form `chr-pos-ref-alt`. For
            example, `7-140753336-A-T`. gnomAD-VCF uses 1-based inclusive coordinates.
        :param input_assembly: Assembly used for `vcf_query`. If not provided, will try
            to first validate against GRCh38 and then GRCh37
        :return: GnomadVcfToProteinService containing protein VRS Allele, if validation
            and translation was successful
        """
        variation = None
        warnings = []

        # Ensure `vcf_query` is valid (both syntax and reference sequence)
        try:
            valid_result: ValidationResult = await self._get_valid_result(
                vcf_query, warnings, input_assembly=input_assembly
            )
        except GnomadVcfToProteinError as e:
            warnings.append(str(e))
            return GnomadVcfToProteinService(
                variation_query=vcf_query,
                variation=variation,
                warnings=warnings,
                service_meta_=ServiceMeta(
                    version=__version__,
                    response_datetime=datetime.datetime.now(tz=datetime.UTC),
                ),
            )

        # Get relevant genomic information from input `vcf_query`
        token: GnomadVcfToken = valid_result.classification.matching_tokens[0]  # type: ignore
        g_ac: str = valid_result.accession  # type: ignore
        g_start_pos = token.pos  # residue (1-based), following gnomAD-VCF convention
        g_ref = token.ref
        g_alt = token.alt

        # For GRCh37 input assembly, we need to liftover to GRCh38 and update variables
        if input_assembly == ClinVarAssembly.GRCH37:
            grch38_rep: (
                GenomicRepresentation | None
            ) = await self.mane_transcript.g_to_grch38(
                ac=g_ac,
                start_pos=g_start_pos,
                end_pos=g_start_pos,
                coordinate_type=CoordinateType.RESIDUE,
            )
            if not grch38_rep:
                warnings.append(
                    f"Unable to liftover {vcf_query} to GRCh38 representation"
                )
                return GnomadVcfToProteinService(
                    variation_query=vcf_query,
                    variation=variation,
                    warnings=warnings,
                    service_meta_=ServiceMeta(
                        version=__version__,
                        response_datetime=datetime.datetime.now(tz=datetime.UTC),
                    ),
                )

            g_ac = grch38_rep.ac
            g_start_pos = grch38_rep.pos[0] + 1  # Change back to residue

        len_g_ref = len(g_ref)
        len_g_alt = len(g_alt)

        # Compute 1-based inclusive end position from reference sequence length
        # We subtract 1 because start already counts as the first base
        g_end_pos = g_start_pos + len_g_ref - 1

        # Determine the type of genomic alteration (or change) for `vcf_query`
        alt_type: AltType = self._get_genomic_alt_type(
            len_g_ref, len_g_alt, g_ref, g_alt
        )

        # Given genomic data, get associated cDNA and protein consequences
        p_c_data: (
            ProteinAndCdnaRepresentation | None
        ) = await self.mane_transcript.grch38_to_mane_c_p(
            g_ac,
            g_start_pos,
            g_end_pos,
            try_longest_compatible=True,
            coordinate_type=CoordinateType.RESIDUE,
        )
        if not p_c_data:
            warnings.append("Unable to get cDNA and protein representation")
            return GnomadVcfToProteinService(
                variation_query=vcf_query,
                variation=variation,
                warnings=warnings,
                service_meta_=ServiceMeta(
                    version=__version__,
                    response_datetime=datetime.datetime.now(tz=datetime.UTC),
                ),
            )

        # NOTE: These coordinates are inter-residue
        p_data: DataRepresentation = p_c_data.protein
        c_data: CdnaRepresentation = p_c_data.cdna

        # Get GA4GH identifier (`ga4gh:SQ.`) for protein accession.
        # This is used later, but we want to fail fast
        p_ac = p_data.refseq or p_data.ensembl
        p_ga4gh_seq_id, w = self.seqrepo_access.translate_identifier(p_ac, "ga4gh")
        if w:
            warnings.append(w)
            return GnomadVcfToProteinService(
                variation_query=vcf_query,
                variation=variation,
                warnings=warnings,
                service_meta_=ServiceMeta(
                    version=__version__,
                    response_datetime=datetime.datetime.now(tz=datetime.UTC),
                ),
            )

        # Get genomic position range change (NOTE: these are residue coordinates)
        # This ensures that there 3 nucleotides needed for codon
        strand = c_data.strand
        codon_aligned_interval_start, codon_aligned_interval_end, genomic_start_ix = (
            self._get_codon_aligned_interval(
                c_data.pos[0], c_data.pos[1], strand, g_start_pos, g_end_pos
            )
        )

        # Get genomic reference sequence for the codon-aligned interval
        codon_aligned_ref_seq, w = self.seqrepo_access.get_reference_sequence(
            g_ac,
            codon_aligned_interval_start,
            codon_aligned_interval_end,
            coordinate_type=CoordinateType.RESIDUE,
        )
        if w:
            warnings.append(w)
            return GnomadVcfToProteinService(
                variation_query=vcf_query,
                variation=variation,
                warnings=warnings,
                service_meta_=ServiceMeta(
                    version=__version__,
                    response_datetime=datetime.datetime.now(tz=datetime.UTC),
                ),
            )

        if strand == Strand.NEGATIVE:
            codon_aligned_ref_seq = codon_aligned_ref_seq[::-1]

        # Get genomic altered sequence within a codon-aligned interval
        codon_aligned_alt = self._get_codon_aligned_alternate_sequence(
            g_ac,
            g_alt,
            len_g_ref,
            codon_aligned_interval_end,
            alt_type,
            genomic_start_ix,
            strand,
            codon_aligned_ref_seq,
        )

        # DNA -> RNA -> Protein (1 AA)
        aa_ref = self._dna_to_aa(codon_aligned_ref_seq, strand)
        aa_alt = self._dna_to_aa(codon_aligned_alt, strand)

        # Trim AA prefixes / suffixes and update the protein start position accordingly
        # NOTE: These are inter-residue coordinates
        aa_start_pos, aa_end_pos = p_data.pos
        aa_ref, aa_alt, aa_start_pos = _trim_prefix_or_suffix(
            aa_ref, aa_alt, aa_start_pos=aa_start_pos, trim_prefix=True
        )
        aa_ref, aa_alt, _ = _trim_prefix_or_suffix(aa_ref, aa_alt, trim_prefix=False)

        # Construct the protein consequence (VRS Allele)
        try:
            variation = self._get_protein_representation(
                p_ga4gh_seq_id, p_ac, aa_start_pos, aa_end_pos, aa_alt
            )
        except GnomadVcfToProteinError as e:
            warnings.append(str(e))

        if p_data.gene and c_data.gene and p_data.gene != c_data.gene:
            warnings.append(
                f"Protein gene ({p_data.gene}) and cDNA gene ({c_data.gene}) mismatch"
            )
        gene = p_data.gene or c_data.gene
        gene_context = self._get_gene_context(gene) if gene else None

        return GnomadVcfToProteinService(
            variation_query=vcf_query,
            variation=variation,
            gene_context=gene_context,
            warnings=warnings,
            service_meta_=ServiceMeta(
                version=__version__,
                response_datetime=datetime.datetime.now(tz=datetime.UTC),
            ),
        )
