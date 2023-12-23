"""Module for going from gnomAD VCF to VRS variation on the protein coordinate"""
from copy import deepcopy
from datetime import datetime
from typing import List
from urllib.parse import quote

from cool_seq_tool.handlers import SeqRepoAccess
from cool_seq_tool.mappers import MANETranscript
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
from variation.schemas.translation_response_schema import TranslationResult
from variation.schemas.validation_response_schema import ValidationResult
from variation.to_vrsatile import ToVRSATILE
from variation.tokenize import Tokenize
from variation.translate import Translate
from variation.validate import Validate
from variation.version import __version__


class GnomadVcfToProteinError(Exception):
    """Custom exception for Gnomad VCF To Protein"""


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
        mane_transcript: MANETranscript,
    ) -> None:
        """Initialize the GnomadVcfToProteinVariation class

        :param SeqRepoAccess seqrepo_access: Access to SeqRepo
        :param Tokenize tokenizer: Tokenizer class for tokenizing
        :param Classify classifier: Classifier class for classifying tokens
        :param Validate validator: Validator class for validating valid inputs
        :param Translate translator: Translating valid inputs
        :parm GeneQueryHandler gene_normalizer: Client for normalizing gene concepts
        :param MANETranscript mane_transcript: Access MANE Transcript information
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
            rna_seq = dna_seq.replace("T", "U")

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
        _id = f"normalize.variation:{quote(' '.join(q.split()))}"

        try:
            valid_result = await self._get_valid_result(q, warnings)
        except GnomadVcfToProteinError as e:
            warnings.append(str(e))
            return NormalizeService(
                variation_query=q,
                variation_descriptor=vd,
                warnings=warnings,
                service_meta_=ServiceMeta(
                    version=__version__, response_datetime=datetime.now()
                ),
            )

        classification = valid_result.classification
        token = classification.matching_tokens[0]
        g_ac = valid_result.accession
        g_ref = token.ref
        g_alt = token.alt
        g_start_pos = token.pos
        g_end_pos = g_start_pos + (len(g_ref) - 1)

        p_c_data = await self.mane_transcript.grch38_to_mane_c_p(
            g_ac, g_start_pos, g_end_pos, try_longest_compatible=True
        )
        if not p_c_data:
            warnings.append("Unable to get MANE cDNA and protein representation")
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

        start_reading_frame = self.mane_transcript._get_reading_frame(c_data.pos[0] + 1)
        end_reading_frame = self.mane_transcript._get_reading_frame(c_data.pos[1])

        strand = c_data.strand
        start_ix = start_reading_frame - 1
        if strand == Strand.NEGATIVE:
            new_g_end_pos = g_end_pos + (start_reading_frame - 1)
            new_g_start_pos = g_start_pos - (3 - end_reading_frame)
        else:
            new_g_start_pos = g_start_pos - (start_reading_frame - 1)
            new_g_end_pos = g_end_pos + (3 - end_reading_frame)

        ref, _ = self.seqrepo_access.get_reference_sequence(
            g_ac, new_g_start_pos, new_g_end_pos
        )

        if strand == Strand.NEGATIVE:
            ref = ref[::-1]

        if strand == Strand.POSITIVE:
            alt = ref[:start_ix] + g_alt
        else:
            alt = ref[:start_ix] + g_alt[::-1]
        alt += ref[len(alt) :]

        aa_alt = self._dna_to_aa(alt, strand)
        aa_ref = self._dna_to_aa(ref, strand)

        aa_start_pos = p_data.pos[0]
        if aa_ref != aa_alt:
            aa_match = 0
            for i in range(len(aa_ref)):
                if aa_ref[i] == aa_alt[i]:
                    aa_match += 1
                else:
                    break
            aa_start_pos += aa_match
            aa_alt = aa_alt[aa_match:]

        seq_id = p_data.refseq or p_data.ensembl
        ga4gh_seq_id, _ = self.seqrepo_access.translate_identifier(seq_id, "ga4gh")

        a = models.Allele(
            location=models.SequenceLocation(
                sequence_id=ga4gh_seq_id[0],
                interval=models.SequenceInterval(
                    start=models.Number(value=aa_start_pos),
                    end=models.Number(value=p_data.pos[1]),
                ),
            ),
            state=models.LiteralSequenceExpression(sequence=aa_alt),
        )
        a = normalize(a, self.seqrepo_access)
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
