"""This module provides methods for handling queries."""
from typing import Tuple, Optional, List, Union, Dict
from urllib.parse import quote

from gene.query import QueryHandler as GeneQueryHandler
from ga4gh.vrs.dataproxy import SeqRepoDataProxy
from ga4gh.vrs.extras.translator import Translator
from ga4gh.core import ga4gh_identify
from ga4gh.vrs import models
from variation import SEQREPO_DATA_PATH, TRANSCRIPT_MAPPINGS_PATH, \
    REFSEQ_GENE_SYMBOL_PATH, AMINO_ACID_PATH, UTA_DB_URL, REFSEQ_MANE_PATH
from variation.schemas.token_response_schema import Nomenclature, Token, \
    ReferenceSequence, SequenceOntology
from variation.schemas.validation_response_schema import ValidationSummary
from variation.to_vrs import ToVRS
from variation.vrs import VRS
from variation.normalize import Normalize
from variation.classifiers import Classify
from variation.tokenizers import Tokenize
from variation.validators import Validate
from variation.translators import Translate
from variation.data_sources import SeqRepoAccess, TranscriptMappings, \
    UTA, MANETranscriptMappings, CodonTable
from variation.mane_transcript import MANETranscript
from variation.hgvs_dup_del_mode import HGVSDupDelMode
from variation.tokenizers import GeneSymbol
from variation.tokenizers.caches import AminoAcidCache
from ga4gh.vrsatile.pydantic.vrs_models import Text, Allele, CopyNumber, \
    Haplotype, VariationSet
from ga4gh.vrsatile.pydantic.vrsatile_models import VariationDescriptor
from variation.schemas.normalize_response_schema\
    import HGVSDupDelMode as HGVSDupDelModeEnum


class QueryHandler:
    """Class for handling queries."""

    def __init__(self,
                 dynamodb_url: str = "",
                 dynamodb_region: str = "us-east-2",
                 seqrepo_data_path: str = SEQREPO_DATA_PATH,
                 amino_acids_file_path: str = AMINO_ACID_PATH,
                 uta_db_url: str = UTA_DB_URL,
                 uta_db_pwd: Optional[str] = None) -> None:
        """Initialize QueryHandler instance.
        :param str dynamodb_url: URL to gene-normalizer database source.
        :param str dynamodb_region: AWS default region for gene-normalizer.
        :param str seqrepo_data_path: Path to seqrepo data directory
        :param str amino_acids_file_path: Path to amino acids file
        :param str uta_db_url: URL for UTA database
        :param Optional[str] uta_db_pwd: Password for UTA database user
        """
        self.amino_acid_cache = AminoAcidCache(
            amino_acids_file_path=amino_acids_file_path
        )
        self.gene_normalizer = GeneQueryHandler(db_url=dynamodb_url,
                                                db_region=dynamodb_region)
        self.seqrepo_access = SeqRepoAccess(
            seqrepo_data_path=seqrepo_data_path
        )
        self.codon_table = CodonTable(self.amino_acid_cache)
        self.uta = UTA(db_url=uta_db_url, db_pwd=uta_db_pwd)
        self.dp = SeqRepoDataProxy(self.seqrepo_access.seq_repo_client)
        self.hgvs_dup_del_mode = HGVSDupDelMode(self.seqrepo_access)
        self.vrs = VRS(self.dp, self.seqrepo_access)
        self.to_vrs_handler = self._init_to_vrs()
        self.normalize_handler = Normalize(
            self.seqrepo_access, self.uta, self.gene_normalizer
        )

    def _init_to_vrs(self,
                     transcript_file_path: str = TRANSCRIPT_MAPPINGS_PATH,
                     refseq_file_path: str = REFSEQ_GENE_SYMBOL_PATH,
                     mane_data_path: str = REFSEQ_MANE_PATH) -> ToVRS:
        """Return toVRS instance

        :param str transcript_file_path: Path to transcript mappings file
        :param str refseq_file_path: Path to refseq gene symbol file
        :param str mane_data_path: Path to refseq mane data file
        :return: toVRS instance
        """
        gene_symbol = GeneSymbol(self.gene_normalizer)
        tokenizer = Tokenize(self.amino_acid_cache, gene_symbol)
        classifier = Classify()
        transcript_mappings = TranscriptMappings(
            transcript_file_path=transcript_file_path,
            refseq_file_path=refseq_file_path
        )
        mane_transcript_mappings = MANETranscriptMappings(
            mane_data_path=mane_data_path
        )
        tlr = Translator(data_proxy=self.dp)
        mane_transcript = MANETranscript(
            self.seqrepo_access, transcript_mappings,
            mane_transcript_mappings, self.uta
        )
        validator = Validate(
            self.seqrepo_access, transcript_mappings, gene_symbol,
            mane_transcript, self.uta, self.dp, tlr,
            self.amino_acid_cache, self.gene_normalizer, self.vrs
        )
        translator = Translate()
        return ToVRS(
            tokenizer, classifier, self.seqrepo_access, transcript_mappings,
            gene_symbol, self.amino_acid_cache, self.uta,
            mane_transcript_mappings, mane_transcript, validator, translator,
            self.gene_normalizer, self.hgvs_dup_del_mode
        )

    def to_vrs(self, q: str)\
            -> Tuple[Optional[Union[List[Allele], List[CopyNumber],
                                    List[Text], List[Haplotype],
                                    List[VariationSet]]],
                     Optional[List[str]]]:
        """Return a VRS-like representation of all validated variations for a query.  # noqa: E501

        :param str q: The variation to translate
        :return: Validated VRS Variations and list of warnings
        """
        validations, warnings = \
            self.to_vrs_handler.get_validations(q)
        translations, warnings = \
            self.to_vrs_handler.get_translations(validations, warnings)

        if not translations:
            if q and q.strip():
                text = models.Text(definition=q)
                text._id = ga4gh_identify(text)
                translations = [Text(**text.as_dict())]
            else:
                translations = None
        return translations, warnings

    def normalize(
            self, q: str,
            hgvs_dup_del_mode: Optional[HGVSDupDelModeEnum] = HGVSDupDelModeEnum.DEFAULT  # noqa: E501
    ) -> Optional[VariationDescriptor]:
        """Return normalized Variation Descriptor for variation.

        :param q: Variation to normalize
        :param Optional[HGVSDupDelModeEnum] hgvs_dup_del_mode:
            Must be set when querying HGVS dup/del expressions.
            Must be: `default`, `cnv`, `repeated_seq_expr`, `literal_seq_expr`.
            This parameter determines how to interpret HGVS dup/del expressions
            in VRS.
        :return: Variation Descriptor for variation
        """
        validations, warnings = \
            self.to_vrs_handler.get_validations(
                q, normalize_endpoint=True,
                hgvs_dup_del_mode=hgvs_dup_del_mode
            )
        if not validations:
            self.normalize_handler.warnings = warnings
            return None
        return self.normalize_handler.normalize(q, validations,
                                                warnings)

    def _get_gnomad_vcf_validations(
            self, q: str, warnings: List) -> Optional[ValidationSummary]:
        """Get gnomad vcf validation summary

        :param str q: Input query
        :param List warnings: List of warnings
        :return: ValidationSummary for a gnomad VCF query
        """
        tokens = self.to_vrs_handler.tokenizer.perform(q.strip(), warnings)
        for t in tokens:
            if t.nomenclature != Nomenclature.GNOMAD_VCF:
                warnings.append(f"{q} is not a supported gnomad vcf query")
                return None
        classifications = self.to_vrs_handler.classifier.perform(tokens)
        validation_summary = self.to_vrs_handler.validator.perform(
            classifications, True, warnings,
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
        aliases = self.seqrepo_access.seq_repo_client.translate_identifier(
            alt_ac, target_namespaces="refseq")
        return aliases[0].split("refseq:")[-1]

    def _update_gnomad_vcf_mane_c_pos(
            self, reading_frame: int, mane_c_ac: str,
            mane_c_pos_change: Tuple[int, int], coding_start_site: int,
            warnings: List) -> Optional[Tuple[int, int]]:
        """Return updated mane c position change for a gnomad vcf variation"""
        if reading_frame == 1:
            # first pos
            mane_c_pos_change = \
                mane_c_pos_change[0], mane_c_pos_change[0] + 2
        elif reading_frame == 2:
            # middle pos
            mane_c_pos_change = \
                mane_c_pos_change[0] - 1, mane_c_pos_change[0] + 1
            pass
        elif reading_frame == 3:
            # last pos
            mane_c_pos_change = \
                mane_c_pos_change[0] - 2, mane_c_pos_change[0]

        if not self.to_vrs_handler.mane_transcript._validate_index(
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
        classification_token.reference_sequence = ReferenceSequence.PROTEIN
        classification_token.molecule_context = "protein"
        if classification_token.alt_type in ["substitution",
                                             "silent_mutation"]:
            if classification_token.alt_type == "substitution":
                alt_nuc = classification_token.new_nucleotide
                classification_token.so_id = \
                    SequenceOntology.AMINO_ACID_SUBSTITUTION
            else:
                alt_nuc = classification_token.ref_nucleotide
                classification_token.so_id = SequenceOntology.SILENT_MUTATION
            if reading_frame == 1:
                # first pos
                ref = self.seqrepo_access.get_sequence(
                    alt_ac, g_start_pos, g_end_pos + 2)
                alt = alt_nuc + ref[1] + ref[2]
            elif reading_frame == 2:
                # middle pos
                ref = self.seqrepo_access.get_sequence(
                    alt_ac, g_start_pos - 1, g_end_pos + 1)
                alt = ref[0] + alt_nuc + ref[2]
            elif reading_frame == 3:
                # last pos
                ref = self.seqrepo_access.get_sequence(
                    alt_ac, g_start_pos - 2, g_end_pos)
                alt = ref[0] + ref[1] + alt_nuc
            if alt:
                # TODO: Check if we only do this for substitutions
                alt = self.codon_table.dna_to_rna(alt)
        elif classification_token.alt_type == "deletion":
            classification_token.so_id = SequenceOntology.AMINO_ACID_DELETION
        elif classification_token.alt_type == "insertion":
            classification_token.so_id = SequenceOntology.AMINO_ACID_INSERTION
            alt = classification_token.inserted_sequence.replace("T", "U")

        if alt is None:
            return None
        else:
            if len(alt) % 3 != 0:
                return None

            if strand == "-":
                alt = alt[::-1]

            aa_alt = ""
            for i in range(int(len(alt) / 3)):
                aa_alt += self.codon_table.table[alt[3 * i:(3 * i) + 3]]
            return aa_alt

    def gnomad_vcf_to_protein(
            self, q: str) -> Tuple[VariationDescriptor, List]:
        """Get protein consequence for gnomad vcf

        :param str q: gnomad vcf (chr-pos-ref-alt)
        :return: Variation Descriptor, list of warnings
        """
        q = q.strip()
        _id = f"normalize.variation:{quote(' '.join(q.split()))}"
        warnings = list()
        valid_list = list()
        invalid_list = list()
        validations = self._get_gnomad_vcf_validations(q, warnings)
        if not validations:
            return self.normalize_handler.text_variation_resp(q, _id, warnings)

        for valid_result in validations.valid_results:
            warnings = list()
            # all gnomad vcf will be alleles with a literal seq expression
            variation = valid_result.variation
            classification_token = valid_result.classification_token
            alt_ac = self._get_refseq_alt_ac_from_variation(variation)

            # 1-based
            g_start_pos =\
                variation["location"]["interval"]["start"]["value"] + 1
            g_end_pos = variation["location"]["interval"]["end"]["value"]
            transcripts = self.uta.get_transcripts_from_genomic_pos(
                alt_ac, g_start_pos)
            if not transcripts:
                warnings.append(f"Unable to get transcripts given {alt_ac} "
                                f"and {g_start_pos}, {g_end_pos}")
                invalid_list.append(
                    self.normalize_handler.text_variation_resp(
                        q, _id, warnings)
                )
                continue
            mane_data = self.to_vrs_handler.mane_transcript_mappings.get_mane_from_transcripts(transcripts)  # noqa: E501

            mane_data_len = len(mane_data)
            for i in range(mane_data_len):
                index = mane_data_len - i - 1
                current_mane_data = mane_data[index]
                mane_c_ac = current_mane_data["RefSeq_nuc"]
                mane_tx_genomic_data = self.uta.get_mane_c_genomic_data(
                    mane_c_ac, alt_ac, g_start_pos, g_end_pos
                )
                if not mane_tx_genomic_data:
                    warnings.append("Unable to get mane transcript and "
                                    "genomic data")
                    invalid_list.append(
                        self.normalize_handler.text_variation_resp(
                            q, _id, warnings)
                    )
                    continue

                coding_start_site = mane_tx_genomic_data["coding_start_site"]
                mane_c_pos_change = \
                    self.to_vrs_handler.mane_transcript.get_mane_c_pos_change(
                        mane_tx_genomic_data, (g_start_pos, g_end_pos),
                        coding_start_site)
                if mane_c_pos_change[0] > mane_c_pos_change[1]:
                    mane_c_pos_change = (mane_c_pos_change[1],
                                         mane_c_pos_change[0])

                reading_frame = self.to_vrs_handler.mane_transcript.get_reading_frame(mane_c_pos_change[0])  # noqa: E501
                if classification_token.alt_type in ["substitution",
                                                     "silent_mutation"]:
                    mane_c_pos_change = self._update_gnomad_vcf_mane_c_pos(
                        reading_frame, mane_c_ac, mane_c_pos_change,
                        coding_start_site, warnings)
                    if mane_c_pos_change is None:
                        invalid_list.append(
                            self.normalize_handler.text_variation_resp(
                                q, _id, warnings)
                        )
                        continue

                mane_p = self.to_vrs_handler.mane_transcript._get_mane_p(
                    current_mane_data, mane_c_pos_change)
                p_ac = mane_p["refseq"]
                valid_result.identifier = p_ac
                aa_alt = self._get_gnomad_vcf_protein_alt(
                    classification_token, reading_frame,
                    mane_tx_genomic_data["strand"], alt_ac,
                    g_start_pos, g_end_pos)
                if aa_alt or classification_token.alt_type == "deletion":
                    variation = self.vrs.to_vrs_allele(
                        p_ac, mane_p["pos"][0], mane_p["pos"][1], "p",
                        classification_token.alt_type, [], alt=aa_alt
                    )
                    if variation:
                        valid_list.append(
                            self.normalize_handler.get_variation_descriptor(
                                variation, valid_result, _id, warnings,
                                gene=current_mane_data["HGNC_ID"]))

        if valid_list:
            return valid_list[0]
        else:
            if invalid_list:
                return invalid_list[0]
            else:
                return self.normalize_handler.text_variation_resp(
                    q, _id, [f"Unable to get protein variation for {q}"])
