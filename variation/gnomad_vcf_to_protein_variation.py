"""Module for going from gnomAD VCF to VRS variation on the protein coordinate"""
from typing import Optional, Tuple, List, Dict
from urllib.parse import quote
from datetime import datetime

from uta_tools.data_sources import SeqRepoAccess, UTADatabase, MANETranscript,\
    MANETranscriptMappings
from uta_tools.schemas import ResidueMode
from gene.query import QueryHandler as GeneQueryHandler
from ga4gh.vrs.dataproxy import SeqRepoDataProxy

from variation.classifiers.classify import Classify
from variation.data_sources.codon_table import CodonTable
from variation.hgvs_dup_del_mode import HGVSDupDelMode
from variation.to_vrsatile import ToVRSATILE
from variation.tokenizers.tokenize import Tokenize
from variation.translators.translate import Translate
from variation.utils import no_variation_entered, no_variation_resp
from variation.validators.validate import Validate
from variation.schemas.validation_response_schema import ValidationSummary
from variation.schemas.token_response_schema import Nomenclature, Token, \
    CoordinateType, SequenceOntology
from variation.schemas.app_schemas import Endpoint
from variation.schemas.normalize_response_schema\
    import HGVSDupDelMode as HGVSDupDelModeEnum, NormalizeService, ServiceMeta
from variation.version import __version__


class GnomadVcfToProteinVariation(ToVRSATILE):
    """Class for translating gnomAD VCF representation to protein representation"""

    def __init__(self, seqrepo_access: SeqRepoAccess, dp: SeqRepoDataProxy,
                 tokenizer: Tokenize, classifier: Classify, validator: Validate,
                 translator: Translate, hgvs_dup_del_mode: HGVSDupDelMode,
                 gene_normalizer: GeneQueryHandler,
                 uta: UTADatabase, mane_transcript: MANETranscript,
                 mane_transcript_mappings: MANETranscriptMappings,
                 codon_table: CodonTable) -> None:
        """Initialize the GnomadVcfToProteinVariation class

        :param SeqRepoAccess seqrepo_access: Access to SeqRepo via UTA Tools
        :param SeqRepoDataProxy dp: Access to SeqRepo via VRS Python
        :param Tokenize tokenizer: Tokenizer class for tokenizing
        :param Classify classifier: Classifier class for classifying tokens
        :param Validate validator: Validator class for validating valid inputs
        :param Translate translator: Translating valid inputs
        :param HGVSDupDelMode hgvs_dup_del_mode: Class for handling
            HGVS dup/del expressions
        :parm GeneQueryHandler gene_normalizer: Client for normalizing gene concepts
        :param UTADatabase uta: Access to db containing alignment data
        :param MANETranscript mane_transcript: Access MANE Transcript
            information
        :param MANETranscriptMappings mane_transcript_mappings: Mappings for
            MANE Transcript data
        :param CodonTable codon_table: Codon table data
        """
        super().__init__(seqrepo_access, dp, tokenizer, classifier, validator,
                         translator, hgvs_dup_del_mode, gene_normalizer)
        self.uta = uta
        self.mane_transcript = mane_transcript
        self.mane_transcript_mappings = mane_transcript_mappings
        self.codon_table = codon_table

    async def _get_gnomad_vcf_validations(
            self, q: str, warnings: List) -> Optional[ValidationSummary]:
        """Get gnomad vcf validation summary

        :param str q: Input query
        :param List warnings: List of warnings
        :return: ValidationSummary for a gnomad VCF query
        """
        tokens = self.tokenizer.perform(q.strip(), warnings)
        for t in tokens:
            if t.nomenclature != Nomenclature.GNOMAD_VCF:
                warnings.append(f"{q} is not a supported gnomad vcf query")
                return None
        classifications = self.classifier.perform(tokens)
        validation_summary = await self.validator.perform(
            classifications, Endpoint.NORMALIZE, warnings,
            hgvs_dup_del_mode=HGVSDupDelModeEnum.LITERAL_SEQ_EXPR
        )
        if not validation_summary:
            warnings.append(f"{q} is not a valid gnomad vcf query")
            return None
        else:
            return validation_summary

    def _get_refseq_alt_ac_from_variation(self, variation: Dict) -> str:
        """Get genomic ac from variation sequence_id

        :param Dict variation: VRS variation object
        :return: RefSeq genomic accession
        """
        # genomic ac should always be in 38
        alt_ac = variation["location"]["sequence_id"]
        aliases = self.seqrepo_access.seqrepo_client.translate_identifier(
            alt_ac, target_namespaces="refseq")
        return aliases[0].split("refseq:")[-1]

    def _update_gnomad_vcf_mane_c_pos(
            self, reading_frame: int, mane_c_ac: str,
            mane_c_pos_change: Tuple[int, int], coding_start_site: int,
            warnings: List) -> Optional[Tuple[int, int]]:
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
            mane_c_pos_change = \
                mane_c_pos_change[0], mane_c_pos_change[0] + 2
        elif reading_frame == 2:
            # middle pos
            mane_c_pos_change = \
                mane_c_pos_change[0] - 1, mane_c_pos_change[0] + 1
        elif reading_frame == 3:
            # last pos
            mane_c_pos_change = \
                mane_c_pos_change[0] - 2, mane_c_pos_change[0]

        if not self.mane_transcript._validate_index(
                mane_c_ac, mane_c_pos_change, coding_start_site):
            warnings.append(
                f"{mane_c_pos_change} are not valid positions on "
                f"{mane_c_ac} with coding start site "
                f"{coding_start_site}")
            return None
        return mane_c_pos_change

    def _get_gnomad_vcf_protein_alt(
            self, classification_token: Token, reading_frame: int, strand: str,
            alt_ac: str, g_start_pos: int, g_end_pos: int) -> Optional[str]:
        """Return protein alteration that corresponds to gnomad VCF alteration

        :param Token classification_token: Classification token for query
        :param int reading_frame: cDNA reading frame number (1, 2, 3)
        :param str strand: Strand for query
        :param str alt_ac: RefSeq genomic accession
        :param int g_start_pos: Genomic start position
        :param int g_end_pos: Genomic end position
        :return: Amino acid alteration (using 1-letter codes)
        """
        alt = None
        residue_mode = ResidueMode.INTER_RESIDUE
        classification_token.coordinate_type = CoordinateType.PROTEIN
        classification_token.molecule_context = "protein"
        if classification_token.alt_type in ["substitution",
                                             "silent_mutation"]:
            if classification_token.alt_type == "substitution":
                alt_nuc = classification_token.new_nucleotide
                classification_token.so_id = \
                    SequenceOntology.PROTEIN_SUBSTITUTION
            else:
                alt_nuc = classification_token.ref_nucleotide
                classification_token.so_id = SequenceOntology.SILENT_MUTATION

            ref = None
            if reading_frame == 1:
                # first pos
                if strand == "-":
                    ref, _ = self.seqrepo_access.get_reference_sequence(
                        alt_ac, g_start_pos - 2, g_end_pos + 1,
                        residue_mode=residue_mode)
                    alt = alt_nuc + ref[1] + ref[0]
                else:
                    ref, _ = self.seqrepo_access.get_reference_sequence(
                        alt_ac, g_start_pos, g_end_pos + 3,
                        residue_mode=residue_mode)
                    alt = alt_nuc + ref[1] + ref[2]
            elif reading_frame == 2:
                # middle pos
                ref, _ = self.seqrepo_access.get_reference_sequence(
                    alt_ac, g_start_pos - 1, g_end_pos + 2,
                    residue_mode=residue_mode)

                if strand == "-":
                    alt = ref[2] + alt_nuc + ref[0]
                else:
                    alt = ref[0] + alt_nuc + ref[2]
            elif reading_frame == 3:
                # last pos
                if strand == "-":
                    ref, _ = self.seqrepo_access.get_reference_sequence(
                        alt_ac, g_start_pos, g_end_pos + 3,
                        residue_mode=residue_mode)
                    alt = ref[2] + ref[1] + alt_nuc
                else:
                    ref, _ = self.seqrepo_access.get_reference_sequence(
                        alt_ac, g_start_pos - 2, g_end_pos + 1,
                        residue_mode=residue_mode)
                    alt = ref[0] + ref[1] + alt_nuc
            if alt and strand == "-":
                alt = self.codon_table.dna_to_rna(alt)
            else:
                alt = alt.replace("T", "U")
        elif classification_token.alt_type == "deletion":
            # There is no alt for a deletion
            classification_token.so_id = SequenceOntology.PROTEIN_DELETION
        elif classification_token.alt_type == "insertion":
            classification_token.so_id = SequenceOntology.PROTEIN_INSERTION
            alt = classification_token.inserted_sequence.replace("T", "U")
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
                aa_alt += self.codon_table.table[alt[3 * i:(3 * i) + 3]]
            return aa_alt

    async def gnomad_vcf_to_protein(
        self, q: str, untranslatable_returns_text: bool = False
    ) -> NormalizeService:
        """Get protein consequence for gnomad vcf

        :param str q: gnomad vcf (chr-pos-ref-alt)
        :param bool untranslatable_returns_text: `True` return VRS Text Object when
            unable to translate or normalize query. `False` return `None` when
            unable to translate or normalize query.
        :return: Normalize Service containing variation descriptor and warnings
        """
        q = q.strip()
        vd = None
        warnings = list()
        if q:
            _id = f"normalize.variation:{quote(' '.join(q.split()))}"
            warnings = list()
            valid_list = list()
            validations = await self._get_gnomad_vcf_validations(q, warnings)
            if validations:
                all_warnings = list()
                for valid_result in validations.valid_results:
                    warnings = list()
                    # all gnomad vcf will be alleles with a literal seq expression
                    variation = valid_result.variation
                    classification_token = valid_result.classification_token
                    alt_ac = self._get_refseq_alt_ac_from_variation(variation)

                    # 0-based
                    g_start_pos = None
                    g_end_pos = None
                    if classification_token.alt_type == "deletion":
                        g_start_pos = classification_token.start_pos_del
                        g_end_pos = classification_token.end_pos_del
                    elif classification_token.alt_type == "insertion":
                        g_start_pos = classification_token.start_pos_flank
                        g_end_pos = classification_token.end_pos_flank
                    elif classification_token.alt_type in ["silent_mutation",
                                                           "substitution"]:
                        g_start_pos = classification_token.position
                        g_end_pos = classification_token.position
                        ref_seq, w = self.seqrepo_access.get_reference_sequence(
                            alt_ac, g_start_pos)
                        if not ref_seq:
                            all_warnings.append(w)
                        else:
                            if ref_seq != classification_token.ref_nucleotide:
                                all_warnings.append(
                                    f"Expected {classification_token.ref_nucleotide}"
                                    f" but found {ref_seq} on {alt_ac} at position"
                                    f" {g_start_pos}"
                                )
                                continue
                    else:
                        all_warnings.append(
                            f"{classification_token.alt_type} alt_type not supported"
                        )
                        continue

                    g_start_pos -= 1
                    g_end_pos -= 1

                    transcripts = await self.uta.get_transcripts_from_genomic_pos(
                        alt_ac, g_start_pos)

                    if not transcripts:
                        all_warnings.append(f"Unable to get transcripts given {alt_ac}"
                                            f" and {g_start_pos + 1}")
                        continue
                    mane_data = self.mane_transcript_mappings.get_mane_from_transcripts(transcripts)  # noqa: E501
                    mane_data_len = len(mane_data)
                    for i in range(mane_data_len):
                        index = mane_data_len - i - 1
                        current_mane_data = mane_data[index]
                        mane_c_ac = current_mane_data["RefSeq_nuc"]
                        mane_tx_genomic_data = await self.uta.get_mane_c_genomic_data(
                            mane_c_ac, alt_ac, g_start_pos, g_end_pos)
                        if not mane_tx_genomic_data:
                            all_warnings.append("Unable to get mane transcript and "
                                                "genomic data")

                        coding_start_site = mane_tx_genomic_data["coding_start_site"]
                        mane_c_pos_change = \
                            self.mane_transcript.get_mane_c_pos_change(
                                mane_tx_genomic_data, coding_start_site)

                        # We use 1-based
                        reading_frame = self.mane_transcript._get_reading_frame(
                            mane_c_pos_change[0] + 1)
                        if classification_token.alt_type in ["substitution",
                                                             "silent_mutation"]:
                            mane_c_pos_change = self._update_gnomad_vcf_mane_c_pos(
                                reading_frame, mane_c_ac, mane_c_pos_change,
                                coding_start_site, warnings)
                            if mane_c_pos_change is None:
                                if len(warnings) > 0:
                                    all_warnings.append(warnings[0])
                                continue

                        mane_p = self.mane_transcript._get_mane_p(
                            current_mane_data,
                            (mane_c_pos_change[0] + 1, mane_c_pos_change[1] + 1)
                        )
                        if mane_p["pos"][0] > mane_p["pos"][1]:
                            mane_p["pos"] = (mane_p["pos"][1], mane_p["pos"][0])
                        p_ac = mane_p["refseq"]
                        valid_result.identifier = p_ac
                        aa_alt = self._get_gnomad_vcf_protein_alt(
                            classification_token, reading_frame,
                            mane_tx_genomic_data["strand"], alt_ac,
                            g_start_pos, g_end_pos)
                        if aa_alt or classification_token.alt_type == "deletion":
                            variation = self.to_vrs_allele(
                                p_ac, mane_p["pos"][0], mane_p["pos"][1], "p",
                                classification_token.alt_type, [], alt=aa_alt
                            )
                            if variation:
                                valid_list.append(
                                    self.get_variation_descriptor(
                                        q, variation, valid_result, _id, warnings,
                                        gene=current_mane_data["HGNC_ID"]))

                if valid_list:
                    vd, warnings = valid_list[0]
                else:
                    if all_warnings:
                        vd, warnings = no_variation_resp(q, _id, all_warnings,
                                                         untranslatable_returns_text)
                    else:
                        vd, warnings = no_variation_resp(
                            q, _id, [f"Unable to get protein variation for {q}"],
                            untranslatable_returns_text)
            else:
                vd, warnings = no_variation_resp(q, _id, warnings,
                                                 untranslatable_returns_text)
        else:
            vd, warnings = no_variation_entered()

        return NormalizeService(
            variation_query=q,
            variation_descriptor=vd,
            warnings=warnings,
            service_meta_=ServiceMeta(
                version=__version__,
                response_datetime=datetime.now()
            )
        )
