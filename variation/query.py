"""This module provides methods for handling queries."""
from typing import Tuple, Optional, List, Union, Dict
from urllib.parse import quote
import copy
import json
from os import environ

import python_jsonschema_objects
from gene.query import QueryHandler as GeneQueryHandler
from ga4gh.vrs.dataproxy import SeqRepoDataProxy
from ga4gh.vrs.extras.translator import Translator
from ga4gh.core import sha512t24u, ga4gh_identify
from ga4gh.vrs import models
from ga4gh.vrsatile.pydantic.vrs_models import Text, Allele, AbsoluteCopyNumber, \
    Haplotype, VariationSet, RelativeCopyNumber
from ga4gh.vrsatile.pydantic.vrsatile_models import VariationDescriptor, \
    CanonicalVariation, ComplexVariation
from uta_tools import SEQREPO_DATA_PATH, TRANSCRIPT_MAPPINGS_PATH, \
    LRG_REFSEQGENE_PATH, MANE_SUMMARY_PATH, UTATools

from variation import AMINO_ACID_PATH, UTA_DB_URL
from variation.schemas.app_schemas import Endpoint
from variation.schemas.hgvs_to_copy_number_schema import VALID_RELATIVE_COPY_CLASS, \
    CopyNumberType, RelativeCopyClass
from variation.schemas.token_response_schema import Nomenclature, Token, \
    CoordinateType, SequenceOntology
from variation.schemas.validation_response_schema import ValidationSummary
from variation.to_vrs import ToVRS
from variation.vrs import VRS
from variation.normalize import Normalize
from variation.classifiers import Classify
from variation.tokenizers import Tokenize
from variation.validators import Validate
from variation.translators import Translate
from variation.data_sources import CodonTable
from variation.hgvs_dup_del_mode import HGVSDupDelMode
from variation.tokenizers import GeneSymbol
from variation.tokenizers.caches import AminoAcidCache
from variation.schemas.normalize_response_schema\
    import HGVSDupDelMode as HGVSDupDelModeEnum


class QueryHandler:
    """Class for handling queries."""

    def __init__(self,
                 dynamodb_url: str = "",
                 dynamodb_region: str = "us-east-2",
                 seqrepo_data_path: str = None,
                 amino_acids_file_path: str = AMINO_ACID_PATH,
                 transcript_file_path: str = None,
                 refseq_file_path: str = None,
                 mane_data_path: str = None,
                 uta_db_url: str = UTA_DB_URL,
                 uta_db_pwd: Optional[str] = None) -> None:
        """Initialize QueryHandler instance.
        :param str dynamodb_url: URL to gene-normalizer database source.
        :param str dynamodb_region: AWS default region for gene-normalizer.
        :param str seqrepo_data_path: Path to seqrepo data directory
        :param str amino_acids_file_path: Path to amino acids file
        :param str transcript_file_path: Path to transcript mappings file
        :param str refseq_file_path: Path to refseq gene symbol file
        :param str mane_data_path: Path to refseq mane data file
        :param str uta_db_url: URL for UTA database
        :param Optional[str] uta_db_pwd: Password for UTA database user
        """
        self.amino_acid_cache = AminoAcidCache(
            amino_acids_file_path=amino_acids_file_path
        )
        self.gene_normalizer = GeneQueryHandler(db_url=dynamodb_url,
                                                db_region=dynamodb_region)

        if not seqrepo_data_path:
            seqrepo_data_path = environ.get("SEQREPO_DATA_PATH", SEQREPO_DATA_PATH)
        if not transcript_file_path:
            transcript_file_path = environ.get("TRANSCRIPT_MAPPINGS_PATH",
                                               TRANSCRIPT_MAPPINGS_PATH)
        if not refseq_file_path:
            refseq_file_path = environ.get("LRG_REFSEQGENE_PATH", LRG_REFSEQGENE_PATH)
        if not mane_data_path:
            mane_data_path = environ.get("MANE_SUMMARY_PATH", MANE_SUMMARY_PATH)

        self.uta_tools = UTATools(seqrepo_data_path=seqrepo_data_path,
                                  transcript_file_path=transcript_file_path,
                                  lrg_refseqgene_path=refseq_file_path,
                                  mane_data_path=mane_data_path,
                                  db_url=uta_db_url, db_pwd=uta_db_pwd)
        self.seqrepo_access = self.uta_tools.seqrepo_access
        self.codon_table = CodonTable(self.amino_acid_cache)
        self.uta = self.uta_tools.uta_db
        self.dp = SeqRepoDataProxy(self.seqrepo_access.seqrepo_client)
        self.tlr = Translator(data_proxy=self.dp)
        self.hgvs_dup_del_mode = HGVSDupDelMode(self.seqrepo_access)
        self.vrs = VRS(self.dp, self.seqrepo_access)
        self.to_vrs_handler = self._init_to_vrs()
        self.normalize_handler = Normalize(
            self.seqrepo_access, self.uta, self.gene_normalizer
        )

    def _init_to_vrs(self) -> ToVRS:
        """Return toVRS instance"""
        gene_symbol = GeneSymbol(self.gene_normalizer)
        tokenizer = Tokenize(self.amino_acid_cache, gene_symbol)
        classifier = Classify()
        transcript_mappings = self.uta_tools.transcript_mappings
        mane_transcript_mappings = self.uta_tools.mane_transcript_mappings
        mane_transcript = self.uta_tools.mane_transcript
        validator = Validate(
            self.seqrepo_access, transcript_mappings, gene_symbol,
            mane_transcript, self.uta, self.dp, self.tlr,
            self.amino_acid_cache, self.gene_normalizer, self.vrs
        )
        translator = Translate()
        return ToVRS(
            tokenizer, classifier, self.seqrepo_access, transcript_mappings,
            gene_symbol, self.amino_acid_cache, self.uta,
            mane_transcript_mappings, mane_transcript, validator, translator,
            self.gene_normalizer, self.hgvs_dup_del_mode
        )

    async def to_vrs(self, q: str)\
            -> Tuple[Optional[Union[List[Allele], List[AbsoluteCopyNumber],
                                    List[Text], List[Haplotype],
                                    List[VariationSet]]],
                     Optional[List[str]]]:
        """Return a VRS-like representation of all validated variations for a query.  # noqa: E501

        :param str q: The variation to translate
        :return: Validated VRS Variations and list of warnings
        """
        validations, warnings = await self.to_vrs_handler.get_validations(q)
        translations, warnings = \
            self.to_vrs_handler.get_translations(validations, warnings)

        if not translations:
            if q and q.strip():
                text = models.Text(definition=q, type="Text")
                text._id = ga4gh_identify(text)
                translations = [Text(**text.as_dict())]
            else:
                translations = None
        return translations, warnings

    async def normalize(
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
        validations, warnings = await self.to_vrs_handler.get_validations(
            q, endpoint_name=Endpoint.NORMALIZE, hgvs_dup_del_mode=hgvs_dup_del_mode
        )
        if not validations:
            self.normalize_handler.warnings = warnings
            return None
        return self.normalize_handler.normalize(q, validations,
                                                warnings)

    async def _get_gnomad_vcf_validations(
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
        validation_summary = await self.to_vrs_handler.validator.perform(
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
                        alt_ac, g_start_pos - 1, g_end_pos + 2)
                    alt = alt_nuc + ref[1] + ref[0]
                else:
                    ref, _ = self.seqrepo_access.get_reference_sequence(
                        alt_ac, g_start_pos + 1, g_end_pos + 4)
                    alt = alt_nuc + ref[1] + ref[2]
            elif reading_frame == 2:
                # middle pos
                ref, _ = self.seqrepo_access.get_reference_sequence(
                    alt_ac, g_start_pos, g_end_pos + 3)
                alt = ref[0] + alt_nuc + ref[2]
            elif reading_frame == 3:
                # last pos
                if strand == "-":
                    ref, _ = self.seqrepo_access.get_reference_sequence(
                        alt_ac, g_start_pos + 1, g_end_pos + 4)
                    alt = ref[2] + ref[1] + alt_nuc
                else:
                    ref, _ = self.seqrepo_access.get_reference_sequence(
                        alt_ac, g_start_pos - 1, g_end_pos + 2)
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
            self, q: str) -> Tuple[Optional[VariationDescriptor], List]:
        """Get protein consequence for gnomad vcf

        :param str q: gnomad vcf (chr-pos-ref-alt)
        :return: Variation Descriptor, list of warnings
        """
        q = q.strip()
        if not q:
            return self.normalize_handler.no_variation_entered()
        _id = f"normalize.variation:{quote(' '.join(q.split()))}"
        warnings = list()
        valid_list = list()
        validations = await self._get_gnomad_vcf_validations(q, warnings)
        if not validations:
            return self.normalize_handler.text_variation_resp(q, _id, warnings)

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
            mane_data = self.to_vrs_handler.mane_transcript_mappings.get_mane_from_transcripts(transcripts)  # noqa: E501
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
                    self.to_vrs_handler.mane_transcript.get_mane_c_pos_change(
                        mane_tx_genomic_data, coding_start_site)

                # We use 1-based
                reading_frame = self.to_vrs_handler.mane_transcript.get_reading_frame(
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

                mane_p = self.to_vrs_handler.mane_transcript._get_mane_p(
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
                    variation = self.vrs.to_vrs_allele(
                        p_ac, mane_p["pos"][0], mane_p["pos"][1], "p",
                        classification_token.alt_type, [], alt=aa_alt
                    )
                    if variation:
                        valid_list.append(
                            self.normalize_handler.get_variation_descriptor(
                                q, variation, valid_result, _id, warnings,
                                gene=current_mane_data["HGNC_ID"]))

        if valid_list:
            return valid_list[0]
        else:
            if all_warnings:
                return self.normalize_handler.text_variation_resp(
                    q, _id, all_warnings)
            else:
                return self.normalize_handler.text_variation_resp(
                    q, _id, [f"Unable to get protein variation for {q}"])

    def canonical_spdi_to_categorical_variation(
            self, q: str, complement: bool = False
    ) -> Tuple[Optional[Union[CanonicalVariation, ComplexVariation]], List]:
        """Return categorical variation for canonical SPDI

        :param str q: Canonical SPDI
        :param bool complement: This field indicates that a categorical
            variation is defined to include (false) or exclude (true) variation
             concepts matching the categorical variation. This is equivalent to a
             logical NOT operation on the categorical variation properties.
        :return: Tuple containing CanonicalVariation and list of warnings
        """
        q = q.strip()
        if not q:
            return self.normalize_handler.no_variation_entered()
        warnings = list()

        variation = None
        try:
            variation = self.tlr.translate_from(q, fmt="spdi")
        except (ValueError, python_jsonschema_objects.validators.ValidationError) as e:
            warnings.append(f"vrs-python translator raised error: {e}")
        except KeyError as e:
            warnings.append(f"vrs-python translator raised error: "
                            f"seqrepo could not translate identifier {e}")
        if warnings or not variation:
            return None, warnings

        spdi_parts = q.split(":")
        ac = spdi_parts[0]
        start_pos = int(spdi_parts[1])  # inter-residue, 0 based
        deleted_seq = spdi_parts[2]
        end_pos = start_pos + len(deleted_seq)
        try:
            sequence = self.seqrepo_access.seqrepo_client.fetch(
                ac, start_pos, end=end_pos)
        except ValueError as e:
            warnings.append(str(e))
        else:
            if not sequence:
                warnings.append(f"Position, {start_pos}, does not exist on {ac}")
            else:
                if deleted_seq != sequence:
                    warnings.append(f"Expected to find reference sequence"
                                    f" {deleted_seq} but found {sequence} on {ac}")
        if warnings:
            return None, warnings

        variation.location._id = ga4gh_identify(variation.location)
        canonical_variation = {
            "type": "CanonicalVariation",
            "complement": complement,
            "variation": variation.as_dict()
        }

        cpy_canonical_variation = copy.deepcopy(canonical_variation)
        cpy_canonical_variation["variation"] = canonical_variation["variation"]["_id"].split(".")[-1]  # noqa: E501
        serialized = json.dumps(
            cpy_canonical_variation, sort_keys=True, separators=(",", ":"), indent=None
        ).encode("utf-8")
        digest = sha512t24u(serialized)
        # VCC = variation categorical canonical
        canonical_variation["_id"] = f"ga4gh:VCC.{digest}"
        canonical_variation = CanonicalVariation(**canonical_variation)
        return canonical_variation, warnings

    def _hgvs_to_cnv_resp(
        self, copy_number_type: CopyNumberType, hgvs_expr: str, do_liftover: bool,
        validations: Tuple[Optional[ValidationSummary], Optional[List[str]]],
        warnings: List[str]
    ) -> Tuple[Optional[Union[AbsoluteCopyNumber, RelativeCopyNumber, Text]], List[str]]:  # noqa: E501
        """Return copy number variation and warnings response

        :param CopyNumberType copy_number_type: The type of copy number variation
        :param str hgvs_expr: HGVS expression
        :param bool do_liftover: Whether or not to liftover to GRCh38 assembly
        :param Tuple[Optional[ValidationSummary], Optional[List[str]]]: Validation
            summary and warnings for hgvs_expr
        :param List[str] warnings: List of warnings
        :return: CopyNumberVariation and warnings
        """
        variation = None
        if do_liftover:
            valid_result = self.normalize_handler.get_valid_result(
                hgvs_expr, validations, [])
            if valid_result:
                variation = valid_result.variation
            else:
                warnings.append(f"Unable to translate {hgvs_expr} to "
                                f"copy number variation")
        else:
            translations, warnings = \
                self.to_vrs_handler.get_translations(validations, warnings)
            if translations:
                variation = translations[0]

        if not variation:
            if hgvs_expr and hgvs_expr.strip():
                text = models.Text(definition=hgvs_expr, type="Text")
                text._id = ga4gh_identify(text)
                variation = Text(**text.as_dict())
        else:
            if copy_number_type == CopyNumberType.ABSOLUTE:
                variation = AbsoluteCopyNumber(**variation)
            else:
                variation = RelativeCopyNumber(**variation)
        return variation, warnings

    async def hgvs_to_absolute_copy_number(
        self, hgvs_expr: str, baseline_copies: Optional[int] = None,
        do_liftover: bool = False
    ) -> Tuple[Optional[AbsoluteCopyNumber], List[str]]:
        """Given hgvs, return abolute copy number variation

        :param str hgvs_expr: HGVS expression
        :param Optional[int] baseline_copies: Baseline copies number
        :param bool do_liftover: Whether or not to liftover to GRCh38 assembly
        :return: Absolute Copy Number Variation and warnings
        """
        validations, warnings = await self.to_vrs_handler.get_validations(
            hgvs_expr, endpoint_name=Endpoint.HGVS_TO_ABSOLUTE_CN,
            hgvs_dup_del_mode=CopyNumberType.ABSOLUTE,
            baseline_copies=baseline_copies, do_liftover=do_liftover
        )
        return self._hgvs_to_cnv_resp(
            CopyNumberType.ABSOLUTE, hgvs_expr, do_liftover, validations, warnings)

    async def hgvs_to_relative_copy_number(
        self, hgvs_expr: str, relative_copy_class: RelativeCopyClass,
        do_liftover: bool = False
    ) -> Tuple[Optional[RelativeCopyNumber], List[str]]:
        """Given hgvs, return relative copy number variation

        :param str hgvs_expr: HGVS expression
        :param RelativeCopyClass relative_copy_class: The relative copy class
        :param bool do_liftover: Whether or not to liftover to GRCh38 assembly
        :return: Relative Copy Number Variation and warnings
        """
        if relative_copy_class and relative_copy_class.lower() not in VALID_RELATIVE_COPY_CLASS:  # noqa: E501
            return None, [f"{relative_copy_class} is not a valid relative copy class: "
                          f"{VALID_RELATIVE_COPY_CLASS}"]

        validations, warnings = await self.to_vrs_handler.get_validations(
            hgvs_expr, endpoint_name=Endpoint.HGVS_TO_RELATIVE_CN,
            hgvs_dup_del_mode=CopyNumberType.RELATIVE,
            relative_copy_class=relative_copy_class, do_liftover=do_liftover
        )

        return self._hgvs_to_cnv_resp(
            CopyNumberType.RELATIVE, hgvs_expr, do_liftover, validations, warnings)
