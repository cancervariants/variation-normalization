"""Module for going from gnomAD VCF to VRS variation on the protein coordinate"""
from copy import deepcopy
from datetime import datetime
from typing import List, Tuple
from urllib.parse import quote

from cool_seq_tool.handlers import SeqRepoAccess
from cool_seq_tool.mappers import ManeTranscript
from cool_seq_tool.schemas import Strand
from cool_seq_tool.sources import TranscriptMappings
from ga4gh.core import ga4gh_identify
from ga4gh.vrs import models, normalize
from ga4gh.vrsatile.pydantic.vrsatile_models import MoleculeContext
from gene.query import QueryHandler as GeneQueryHandler

from variation.classify import Classify
from variation.schemas.classification_response_schema import (
    Nomenclature,
)
from variation.schemas.normalize_response_schema import (
    NormalizeService,
    ServiceMeta,
)
from variation.schemas.service_schema import ClinVarAssembly
from variation.schemas.token_response_schema import AltType
from variation.schemas.translation_response_schema import TranslationResult
from variation.schemas.validation_response_schema import ValidationResult
from variation.to_vrsatile import ToVRSATILE
from variation.tokenize import Tokenize
from variation.translate import Translate
from variation.utils import no_variation_resp
from variation.validate import Validate
from variation.version import __version__


class GnomadVcfToProteinError(Exception):
    """Custom exception for Gnomad VCF To Protein"""


def _get_prefix_match_count(min_length: int, ref: str, alt: str) -> int:
    """Get the count of matched sequential prefixes

    :param min_length: Length of the shortest sequence (using ``ref`` or ``alt``)
    :param ref: Reference sequence
    :param alt: Alternate sequence
    :return: The number of sequential characters that were the same in ``ref`` and
        ``alt``
    """
    matched = 0
    for i in range(min_length):
        if alt[i] == ref[i]:
            matched += 1
        else:
            break
    return matched


def _trim_prefix_suffix(
    aa_ref: str, aa_alt: str, aa_start_pos: int = 0, trim_prefix: bool = True
) -> Tuple[str, str, int]:
    """Trim prefix or suffix matches

    :param aa_ref: Amino acid reference sequence
    :param aa_alt: Amino acid altered sequence
    :param aa_start_pos: Amino acid start position. Only required when ``trim_prefix``
        is ``True``
    :param trim_prefix: ``True`` if trimming prefixes. ``False`` if trimming suffixes
    :return: Tuple containing trimmed ``aa_ref``, trimmed `aa_alt``, and `
        `aa_start_pos`` after trimming
    """
    if (aa_ref and aa_alt) and (aa_ref != aa_alt):
        aa_match = 0
        len_aa_ref = len(aa_ref)
        len_aa_alt = len(aa_alt)

        # Trim prefixes
        if len_aa_ref < len_aa_alt:
            range_len = len_aa_ref
        else:
            range_len = len_aa_alt

        aa_match = _get_prefix_match_count(range_len, aa_ref, aa_alt)
        if aa_match:
            aa_start_pos += aa_match

            if trim_prefix:
                aa_alt = aa_alt[aa_match:]
                aa_ref = aa_ref[aa_match:]
            else:
                aa_alt = aa_alt[:-aa_match]
                aa_ref = aa_ref[:-aa_ref]

    return aa_ref, aa_alt, aa_start_pos


CODON_TABLE = {
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


class GnomadVcfToProteinVariation(ToVRSATILE):
    """Class for translating gnomAD VCF representation to protein representation"""

    def __init__(
        self,
        seqrepo_access: SeqRepoAccess,
        tokenizer: Tokenize,
        classifier: Classify,
        validator: Validate,
        translator: Translate,
        gene_normalizer: GeneQueryHandler,
        transcript_mappings: TranscriptMappings,
        mane_transcript: ManeTranscript,
    ) -> None:
        """Initialize the GnomadVcfToProteinVariation class

        :param seqrepo_access: Access to SeqRepo
        :param tokenizer: Tokenizer class for tokenizing
        :param classifier: Classifier class for classifying tokens
        :param validator: Validator class for validating valid inputs
        :param translator: Translating valid inputs
        :parm gene_normalizer: Client for normalizing gene concepts
        :param mane_transcript: Access MANE Transcript information
        """
        super().__init__(
            seqrepo_access,
            tokenizer,
            classifier,
            validator,
            translator,
            gene_normalizer,
            transcript_mappings,
        )
        self.mane_transcript = mane_transcript

    async def _get_valid_result(self, q: str, warnings: List) -> List[ValidationResult]:
        """Get gnomad vcf validation summary

        :param q: gnomad vcf input query
        :param warnings: List of warnings
        :return: List of valid results for a gnomad VCF query
        """
        tokens = self.tokenizer.perform(q.strip(), warnings)
        if not tokens:
            raise GnomadVcfToProteinError("Unable to get tokens")

        classification = self.classifier.perform(tokens)
        if not classification:
            raise GnomadVcfToProteinError("Unable to get classification")

        if classification.nomenclature != Nomenclature.GNOMAD_VCF:
            raise GnomadVcfToProteinError(f"{q} is not a supported gnomad vcf query")

        validation_summary = await self.validator.perform(
            classification, input_assembly=ClinVarAssembly.GRCH38
        )
        valid_results = validation_summary.valid_results
        if valid_results:
            len_valid_results = len(valid_results)
            if len_valid_results != 1:
                raise GnomadVcfToProteinError(
                    f"Expected 1 valid result, but found {len_valid_results}"
                )
            else:
                return valid_results[0]
        else:
            raise GnomadVcfToProteinError(f"{q} is not a valid gnomad vcf query")

    @staticmethod
    def _dna_to_aa(dna_seq: str, strand: Strand) -> str:
        """Get amino acid(s) from dna sequence

        :param dna_seq: DNA sequence
        :raises ValueError: If DNA character is not supported
        :return: Amino acid(s)
        """
        rna_seq = ""
        if strand == strand.NEGATIVE:
            # Since it's on the negative strand, we need to flip
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
                    raise ValueError
        else:
            # We only need to replace T/U for DNA->RNA
            rna_seq = dna_seq.replace("T", "U")

        # RNA -> Protein
        aa = ""
        for i in range(int(len(rna_seq) / 3)):
            aa += CODON_TABLE[rna_seq[3 * i : (3 * i) + 3]]
        return aa

    async def gnomad_vcf_to_protein(
        self, q: str, untranslatable_returns_text: bool = False
    ) -> NormalizeService:
        """Get MANE protein consequence for gnomad vcf (chr-pos-ref-alt).
        Assumes using GRCh38 coordinates

        :param q: gnomad vcf (chr-pos-ref-alt)
        :param untranslatable_returns_text: `True` return VRS Text Object when unable to
            translate or normalize query. `False` return `None` when unable to translate
            or normalize query.
        :return: Normalize Service containing variation descriptor and warnings
        """
        q = q.strip()
        vd = None
        warnings = []
        _id = f"normalize.variation:{quote(q)}"

        # First we need to validate the input query
        try:
            valid_result = await self._get_valid_result(q, warnings)
        except GnomadVcfToProteinError as e:
            warnings.append(str(e))
            vd, warnings = no_variation_resp(
                q, _id, warnings, untranslatable_returns_text
            )
            return NormalizeService(
                variation_query=q,
                variation_descriptor=vd,
                warnings=warnings,
                service_meta_=ServiceMeta(
                    version=__version__, response_datetime=datetime.now()
                ),
            )

        # Get relevant information from query
        classification = valid_result.classification
        token = classification.matching_tokens[0]
        g_ac = valid_result.accession
        g_ref = token.ref
        g_alt = token.alt

        len_g_ref = len(g_ref)
        len_g_alt = len(g_alt)

        # Determine the type of alteration and the number of nucleotide prefixes matched
        if len_g_ref == len_g_alt:
            num_prefix_matched = _get_prefix_match_count(len_g_ref, g_ref, g_alt)
            alt_type = AltType.SUBSTITUTION
        else:
            if len_g_ref > len_g_alt:
                num_prefix_matched = _get_prefix_match_count(len_g_alt, g_ref, g_alt)
                if num_prefix_matched == len_g_alt:
                    alt_type = AltType.DELETION
                else:
                    alt_type = AltType.DELINS
            else:
                num_prefix_matched = _get_prefix_match_count(len_g_ref, g_ref, g_alt)
                if num_prefix_matched == len_g_ref:
                    alt_type = AltType.INSERTION
                else:
                    alt_type = AltType.DELINS

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
            warnings.append("Unable to get MANE cDNA and protein representation")
            vd, warnings = no_variation_resp(
                q, _id, warnings, untranslatable_returns_text
            )
            return NormalizeService(
                variation_query=q,
                variation_descriptor=vd,
                warnings=warnings,
                service_meta_=ServiceMeta(
                    version=__version__, response_datetime=datetime.now()
                ),
            )

        p_data = p_c_data.protein
        c_data = p_c_data.cdna

        # Get reading frame
        start_reading_frame = self.mane_transcript._get_reading_frame(c_data.pos[0] + 1)
        end_reading_frame = self.mane_transcript._get_reading_frame(c_data.pos[1])

        # Get genomic position range change
        strand = c_data.strand
        start_ix = start_reading_frame - 1
        if strand == Strand.NEGATIVE:
            new_g_end_pos = g_end_pos + (start_reading_frame - 1)
            new_g_start_pos = g_start_pos - (3 - end_reading_frame)
        else:
            new_g_start_pos = g_start_pos - (start_reading_frame - 1)
            new_g_end_pos = g_end_pos + (3 - end_reading_frame)

        # Get reference sequence
        ref, _ = self.seqrepo_access.get_reference_sequence(
            g_ac, new_g_start_pos, new_g_end_pos
        )

        if strand == Strand.NEGATIVE:
            ref = ref[::-1]

        # Get altered sequence. Deletion will always be empty string
        if alt_type == AltType.DELETION:
            alt = ""
        else:
            if strand == Strand.POSITIVE:
                alt = ref[:start_ix] + g_alt
            else:
                alt = ref[:start_ix] + g_alt[::-1]

            if alt_type == AltType.SUBSTITUTION:
                alt += ref[len(alt) :]
            else:
                len_ref = len(ref)
                alt += ref[len_ref - start_ix :]

            # We need to get the entire inserted sequence. It needs to be a factor of 3
            # since DNA (3 nuc) -> RNA (3 nuc) -> Protein (1 aa). The reason why we
            # DO NOT do this for insertions, is because we only want the provided
            # insertion sequence and do not want to mess with it
            if alt_type != AltType.INSERTION:
                len_alt = len(alt)
                rem_alt = len_alt % 3
                if rem_alt != 0:
                    tmp_g_end_pos = new_g_end_pos + (3 - rem_alt)
                    tmp_ref, _ = self.seqrepo_access.get_reference_sequence(
                        g_ac, new_g_end_pos, tmp_g_end_pos
                    )
                    alt += tmp_ref

        # Get protein sequence
        aa_ref = self._dna_to_aa(ref, strand)
        aa_alt = self._dna_to_aa(alt, strand)

        # Get protein start position
        # We need to trim prefixes / suffixes and update the position accordingly
        aa_start_pos = p_data.pos[0]
        aa_ref, aa_alt, aa_start_pos = _trim_prefix_suffix(
            aa_ref, aa_alt, aa_start_pos=aa_start_pos, trim_prefix=True
        )
        aa_ref, aa_alt, _ = _trim_prefix_suffix(aa_ref, aa_alt, trim_prefix=False)

        seq_id = p_data.refseq or p_data.ensembl
        ga4gh_seq_id, w = self.seqrepo_access.translate_identifier(seq_id, "ga4gh")
        if w:
            warnings.append(w)
            vd, warnings = no_variation_resp(
                q, _id, warnings, untranslatable_returns_text
            )
            return NormalizeService(
                variation_query=q,
                variation_descriptor=vd,
                warnings=warnings,
                service_meta_=ServiceMeta(
                    version=__version__, response_datetime=datetime.now()
                ),
            )

        if alt_type == AltType.DELETION:
            aa_end_pos = aa_start_pos + (len(aa_ref) - 1)
        else:
            aa_end_pos = p_data.pos[1]

        a = models.Allele(
            location=models.SequenceLocation(
                sequence_id=ga4gh_seq_id[0],
                interval=models.SequenceInterval(
                    start=models.Number(value=aa_start_pos),
                    end=models.Number(value=aa_end_pos),
                ),
            ),
            state=models.LiteralSequenceExpression(sequence=aa_alt),
        )
        try:
            a = normalize(a, self.seqrepo_access)
        except (KeyError, AttributeError) as e:
            warnings.append(f"VRS-Python unable to normalize allele: {e}")
            vd, warnings = no_variation_resp(
                q, _id, warnings, untranslatable_returns_text
            )
            return NormalizeService(
                variation_query=q,
                variation_descriptor=vd,
                warnings=warnings,
                service_meta_=ServiceMeta(
                    version=__version__, response_datetime=datetime.now()
                ),
            )
        a._id = ga4gh_identify(a)
        a.location._id = ga4gh_identify(a.location)

        translation_result = TranslationResult(
            vrs_variation=a.as_dict(), validation_result=valid_result
        )
        tr_copy = deepcopy(translation_result)
        tr_copy.vrs_seq_loc_ac = p_data.refseq or p_data.ensembl
        tr_copy.vrs_seq_loc_ac_status = p_data.status

        vd, warnings = self.get_variation_descriptor(
            q,
            tr_copy,
            valid_result,
            _id,
            warnings,
            gene=p_data.gene,
        )
        vd.molecule_context = MoleculeContext.PROTEIN

        return NormalizeService(
            variation_query=q,
            variation_descriptor=vd,
            warnings=warnings,
            service_meta_=ServiceMeta(
                version=__version__, response_datetime=datetime.now()
            ),
        )
