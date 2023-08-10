"""Module for going from gnomAD VCF to VRS variation on the protein coordinate"""
from copy import deepcopy
from datetime import datetime
from typing import Dict, List, Optional, Tuple
from urllib.parse import quote

from cool_seq_tool.data_sources import (
    MANETranscript,
    MANETranscriptMappings,
    SeqRepoAccess,
    TranscriptMappings,
    UTADatabase,
)
from cool_seq_tool.schemas import ResidueMode
from ga4gh.vrsatile.pydantic.vrsatile_models import MoleculeContext
from gene.query import QueryHandler as GeneQueryHandler

from variation.classify import Classify
from variation.schemas.app_schemas import Endpoint
from variation.schemas.classification_response_schema import (
    ClassificationType,
    Nomenclature,
)
from variation.schemas.normalize_response_schema import (
    HGVSDupDelModeOption,
    NormalizeService,
    ServiceMeta,
)
from variation.schemas.token_response_schema import AltType, Token
from variation.schemas.translation_response_schema import TranslationResult
from variation.schemas.validation_response_schema import ValidationResult
from variation.to_vrsatile import ToVRSATILE
from variation.tokenize import Tokenize
from variation.translate import Translate
from variation.utils import no_variation_resp
from variation.validate import Validate
from variation.version import __version__

DNA_TO_RNA = {"T": "A", "A": "U", "G": "C", "C": "G"}

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


def dna_to_rna(dna_codon: str) -> str:
    """Convert DNA codon to RNA codon.

    :param str dna_codon: DNA codon
    :return: RNA codon
    """
    dna_codon_list = list(dna_codon)
    rna_codon = ""
    for char in dna_codon_list:
        rna_codon += DNA_TO_RNA[char]
    return rna_codon


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
        uta: UTADatabase,
        mane_transcript: MANETranscript,
        mane_transcript_mappings: MANETranscriptMappings,
    ) -> None:
        """Initialize the GnomadVcfToProteinVariation class

        :param SeqRepoAccess seqrepo_access: Access to SeqRepo
        :param Tokenize tokenizer: Tokenizer class for tokenizing
        :param Classify classifier: Classifier class for classifying tokens
        :param Validate validator: Validator class for validating valid inputs
        :param Translate translator: Translating valid inputs
        :parm GeneQueryHandler gene_normalizer: Client for normalizing gene concepts
        :param UTADatabase uta: Access to db containing alignment data
        :param MANETranscript mane_transcript: Access MANE Transcript
            information
        :param MANETranscriptMappings mane_transcript_mappings: Mappings for
            MANE Transcript data
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
        self.uta = uta
        self.mane_transcript = mane_transcript
        self.mane_transcript_mappings = mane_transcript_mappings

    async def _get_valid_results(
        self, q: str, warnings: List
    ) -> List[ValidationResult]:
        """Get gnomad vcf validation summary

        :param q: gnomad vcf input query
        :param warnings: List of warnings
        :return: List of valid results for a gnomad VCF query
        """
        tokens = self.tokenizer.perform(q.strip(), warnings)
        if not tokens:
            return None

        classification = self.classifier.perform(tokens)
        if not classification:
            return None

        if classification.nomenclature != Nomenclature.GNOMAD_VCF:
            warnings.append(f"{q} is not a supported gnomad vcf query")
            return None

        validation_summary = await self.validator.perform(classification)
        if validation_summary.valid_results:
            valid_results = validation_summary.valid_results
        else:
            warnings.append(f"{q} is not a valid gnomad vcf query")
            valid_results = []
        return valid_results

    def _get_refseq_alt_ac_from_variation(self, variation: Dict) -> str:
        """Get genomic ac from variation sequence_id

        :param Dict variation: VRS variation object
        :return: RefSeq genomic accession
        """
        # genomic ac should always be in 38
        alt_ac = variation["location"]["sequence_id"]
        aliases = self.seqrepo_access.sr.translate_identifier(
            alt_ac, target_namespaces="refseq"
        )
        return aliases[0].split("refseq:")[-1]

    def _update_gnomad_vcf_mane_c_pos(
        self,
        reading_frame: int,
        mane_c_ac: str,
        mane_c_pos_change: Tuple[int, int],
        coding_start_site: int,
        warnings: List,
    ) -> Optional[Tuple[int, int]]:
        """Return updated mane c position change for a gnomad vcf variation
        depending on reading frame base

        :param int reading_frame: reading frame base
        :param str mane_c_ac: Mane transcript accession
        :param Tuple[int, int] mane_c_pos_change: Mane transcript position
            change
        :param int coding_start_site: Coding start site
        :param List warnings: List of warnings
        :return: Mane c pos start and end
        """
        if reading_frame == 1:
            # first pos
            mane_c_pos_change = mane_c_pos_change[0], mane_c_pos_change[0] + 2
        elif reading_frame == 2:
            # middle pos
            mane_c_pos_change = mane_c_pos_change[0] - 1, mane_c_pos_change[0] + 1
        elif reading_frame == 3:
            # last pos
            mane_c_pos_change = mane_c_pos_change[0] - 2, mane_c_pos_change[0]

        if not self.mane_transcript._validate_index(
            mane_c_ac, mane_c_pos_change, coding_start_site
        ):
            warnings.append(
                f"{mane_c_pos_change} are not valid positions on "
                f"{mane_c_ac} with coding start site "
                f"{coding_start_site}"
            )
            return None
        return mane_c_pos_change

    def _get_gnomad_vcf_protein_alt(
        self,
        classification_token: Token,
        alt_type: AltType,
        reading_frame: int,
        strand: str,
        alt_ac: str,
        g_start_pos: int,
        g_end_pos: int,
    ) -> Optional[str]:
        """Return protein alteration that corresponds to gnomad VCF alteration

        :param classification_token: Classification token for query
        :param alt_type: Alteration type
        :param reading_frame: cDNA reading frame number (1, 2, 3)
        :param strand: Strand for query
        :param alt_ac: RefSeq genomic accession
        :param g_start_pos: Genomic start position
        :param g_end_pos: Genomic end position
        :return: Amino acid alteration (using 1-letter codes)
        """
        alt = None
        residue_mode = ResidueMode.INTER_RESIDUE
        if alt_type in {AltType.SUBSTITUTION, AltType.REFERENCE_AGREE}:
            alt_nuc = classification_token.matching_tokens[0].alt

            ref = None
            if reading_frame == 1:
                # first pos
                if strand == "-":
                    ref, _ = self.seqrepo_access.get_reference_sequence(
                        alt_ac,
                        g_start_pos - 2,
                        g_end_pos + 1,
                        residue_mode=residue_mode,
                    )
                    alt = alt_nuc + ref[1] + ref[0]
                else:
                    ref, _ = self.seqrepo_access.get_reference_sequence(
                        alt_ac, g_start_pos, g_end_pos + 3, residue_mode=residue_mode
                    )
                    alt = alt_nuc + ref[1] + ref[2]
            elif reading_frame == 2:
                # middle pos
                ref, _ = self.seqrepo_access.get_reference_sequence(
                    alt_ac, g_start_pos - 1, g_end_pos + 2, residue_mode=residue_mode
                )

                if strand == "-":
                    alt = ref[2] + alt_nuc + ref[0]
                else:
                    alt = ref[0] + alt_nuc + ref[2]
            elif reading_frame == 3:
                # last pos
                if strand == "-":
                    ref, _ = self.seqrepo_access.get_reference_sequence(
                        alt_ac, g_start_pos, g_end_pos + 3, residue_mode=residue_mode
                    )
                    alt = ref[2] + ref[1] + alt_nuc
                else:
                    ref, _ = self.seqrepo_access.get_reference_sequence(
                        alt_ac,
                        g_start_pos - 2,
                        g_end_pos + 1,
                        residue_mode=residue_mode,
                    )
                    alt = ref[0] + ref[1] + alt_nuc
            if alt and strand == "-":
                alt = dna_to_rna(alt)
            else:
                alt = alt.replace("T", "U")
        elif alt_type == AltType.INSERTION:
            alt = classification_token.inserted_sequence[1:].replace("T", "U")
            if strand == "-":
                alt = alt[::-1]
        else:
            return None

        if alt is None:
            return None
        else:
            if len(alt) % 3 != 0:
                return None

            aa_alt = ""
            for i in range(int(len(alt) / 3)):
                aa_alt += CODON_TABLE[alt[3 * i : (3 * i) + 3]]
            return aa_alt

    async def gnomad_vcf_to_protein(
        self, q: str, untranslatable_returns_text: bool = False
    ) -> NormalizeService:
        """Get MANE protein consequence for gnomad vcf (chr-pos-ref-alt).
        Assumes using GRCh38 coordinates

        :param str q: gnomad vcf (chr-pos-ref-alt)
        :param bool untranslatable_returns_text: `True` return VRS Text Object when
            unable to translate or normalize query. `False` return `None` when
            unable to translate or normalize query.
        :return: Normalize Service containing variation descriptor and warnings
        """
        q = q.strip()
        vd = None
        warnings = []
        _id = f"normalize.variation:{quote(' '.join(q.split()))}"

        valid_results = await self._get_valid_results(q, warnings)
        if valid_results:
            translations, warnings = await self.get_translations(
                valid_results,
                warnings,
                Endpoint.NORMALIZE,
                hgvs_dup_del_mode=HGVSDupDelModeOption.LITERAL_SEQ_EXPR,
            )

            if translations:
                translations.sort(
                    key=lambda t: (t.og_ac.split(".")[0], int(t.og_ac.split(".")[1])),
                    reverse=True,
                )

                all_warnings = set()
                checked_valid_results = []
                for translation in translations:
                    warnings = []
                    # all gnomad vcf will be alleles with a literal seq expression
                    variation = translation.vrs_variation
                    validation_result = translation.validation_result
                    classification_token = validation_result.classification

                    # We do not need to check the same variation that has the same
                    # classification
                    checked_tuple = (
                        variation["_id"],
                        translation.vrs_seq_loc_ac,
                        classification_token.classification_type.value,
                    )
                    if checked_tuple in checked_valid_results:
                        continue

                    checked_valid_results.append(checked_tuple)
                    alt_ac = self._get_refseq_alt_ac_from_variation(variation)

                    # 0-based
                    alt_type = None
                    g_start_pos = None
                    g_end_pos = None
                    if (
                        classification_token.classification_type
                        == ClassificationType.GENOMIC_DELINS
                    ):
                        g_start_pos = classification_token.pos0
                        g_end_pos = (
                            classification_token.pos1
                            if classification_token.pos1
                            else classification_token.pos0
                        )

                        # Right now, deletions and insertions are classified as delins
                        # Only support simple deletions and insertions
                        gnomad_vcf_token = classification_token.matching_tokens[0]
                        ref = gnomad_vcf_token.ref
                        alt = gnomad_vcf_token.alt

                        if ref[0] == alt[0]:
                            if len(alt) == 1:
                                alt_type = AltType.DELETION
                                g_start_pos += 1
                                g_end_pos += 1
                            elif len(ref) == 1:
                                alt_type = AltType.INSERTION
                            else:
                                alt_type = AltType.DELINS
                    elif classification_token.classification_type in {
                        ClassificationType.GENOMIC_SUBSTITUTION,
                        ClassificationType.GENOMIC_REFERENCE_AGREE,
                    }:
                        g_start_pos = classification_token.pos
                        g_end_pos = classification_token.pos
                        ref_seq, w = self.seqrepo_access.get_reference_sequence(
                            alt_ac, g_start_pos
                        )
                        if not ref_seq:
                            all_warnings.add(w)
                        else:
                            if ref_seq != classification_token.matching_tokens[0].ref:
                                all_warnings.add(
                                    f"Expected {classification_token.ref} but found "
                                    f"{ref_seq} on {alt_ac} at position {g_start_pos}"
                                )
                                continue

                        if (
                            classification_token.classification_type
                            == ClassificationType.GENOMIC_SUBSTITUTION
                        ):
                            alt_type = AltType.SUBSTITUTION
                        else:
                            alt_type = AltType.REFERENCE_AGREE
                    else:
                        all_warnings.add(
                            f"{classification_token.classification_type} classification_type not supported"  # noqa: E501
                        )
                        continue

                    mane_data = self.mane_transcript_mappings.get_mane_data_from_chr_pos(  # noqa: E501
                        alt_ac, g_start_pos, g_end_pos
                    )

                    mane_data_len = len(mane_data)
                    g_start_pos -= 1
                    g_end_pos -= 1

                    for i in range(mane_data_len):
                        current_mane_data = mane_data[i]
                        mane_c_ac = current_mane_data["RefSeq_nuc"]
                        mane_tx_genomic_data = await self.uta.get_mane_c_genomic_data(
                            mane_c_ac, alt_ac, g_start_pos, g_end_pos
                        )
                        if not mane_tx_genomic_data:
                            all_warnings.add(
                                f"Unable to get MANE data for {mane_c_ac} using "
                                f"{alt_ac} at positions {g_start_pos} to {g_end_pos}"
                            )
                            continue

                        coding_start_site = mane_tx_genomic_data["coding_start_site"]
                        mane_c_pos_change = self.mane_transcript.get_mane_c_pos_change(
                            mane_tx_genomic_data, coding_start_site
                        )

                        # We use 1-based
                        reading_frame = self.mane_transcript._get_reading_frame(
                            mane_c_pos_change[0] + 1
                        )
                        if classification_token.classification_type in {
                            ClassificationType.GENOMIC_SUBSTITUTION,
                            ClassificationType.GENOMIC_REFERENCE_AGREE,
                        }:
                            mane_c_pos_change = self._update_gnomad_vcf_mane_c_pos(
                                reading_frame,
                                mane_c_ac,
                                mane_c_pos_change,
                                coding_start_site,
                                warnings,
                            )
                            if mane_c_pos_change is None:
                                if len(warnings) > 0:
                                    all_warnings.add(warnings[0])
                                continue

                        mane_p = self.mane_transcript._get_mane_p(
                            current_mane_data,
                            (mane_c_pos_change[0] + 1, mane_c_pos_change[1] + 1),
                        )
                        if mane_p["pos"][0] > mane_p["pos"][1]:
                            mane_p["pos"] = (mane_p["pos"][1], mane_p["pos"][0])
                        p_ac = mane_p["refseq"]
                        aa_alt = self._get_gnomad_vcf_protein_alt(
                            classification_token,
                            alt_type,
                            reading_frame,
                            mane_tx_genomic_data["strand"],
                            alt_ac,
                            g_start_pos,
                            g_end_pos,
                        )
                        # Deletions don't have an aa_alt
                        if aa_alt or alt_type == AltType.DELETION:
                            # mane_p is 0-based, but to_vrs allele takes 1-based
                            variation = self.to_vrs_allele(
                                p_ac,
                                mane_p["pos"][0],
                                mane_p["pos"][1],
                                "p",
                                alt_type,
                                [],
                                alt=aa_alt,
                            )
                            if variation:
                                translation_result = TranslationResult(
                                    vrs_variation=variation,
                                    validation_result=validation_result,
                                )

                                tr_copy = deepcopy(translation_result)
                                tr_copy.vrs_seq_loc_ac = p_ac
                                tr_copy.vrs_seq_loc_ac_status = mane_p["status"]

                                vd, warnings = self.get_variation_descriptor(
                                    q,
                                    tr_copy,
                                    validation_result,
                                    _id,
                                    warnings,
                                    gene=current_mane_data["HGNC_ID"],
                                )
                                if not vd:
                                    continue

                                vd.molecule_context = MoleculeContext.PROTEIN

                                return NormalizeService(
                                    variation_query=q,
                                    variation_descriptor=vd,
                                    warnings=warnings,
                                    service_meta_=ServiceMeta(
                                        version=__version__,
                                        response_datetime=datetime.now(),
                                    ),
                                )
                        else:
                            all_warnings.add(
                                "Unable to get associated amino acid change"
                            )

                if all_warnings:
                    vd, warnings = no_variation_resp(
                        q, _id, list(all_warnings), untranslatable_returns_text
                    )
                else:
                    vd, warnings = no_variation_resp(
                        q,
                        _id,
                        [f"Unable to get protein variation for {q}"],
                        untranslatable_returns_text,
                    )
            else:
                vd, warnings = no_variation_resp(
                    q, _id, warnings, untranslatable_returns_text
                )
        else:
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
