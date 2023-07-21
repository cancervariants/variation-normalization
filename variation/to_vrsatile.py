"""Module for representing VRSATILE objects"""
from typing import List, Optional, Tuple, Callable

from ga4gh.vrsatile.pydantic.vrs_models import VRSTypes
from ga4gh.vrsatile.pydantic.vrsatile_models import (
    VariationDescriptor, GeneDescriptor, MoleculeContext
)
from cool_seq_tool.data_sources import SeqRepoAccess, TranscriptMappings
from gene.query import QueryHandler as GeneQueryHandler

from variation.schemas.classification_response_schema import (
    ClassificationType, Nomenclature
)
from variation.to_vrs import ToVRS
from variation.tokenizers import Tokenize
from variation.classifiers import Classify
from variation.validators import Validate
from variation.translators import Translate
from variation.hgvs_dup_del_mode import HGVSDupDelMode
from variation.schemas.token_response_schema import GeneToken
from variation.schemas.validation_response_schema import ValidationResult
from variation.schemas.translation_response_schema import TranslationResult


class ToVRSATILE(ToVRS):
    """Class for representing VRSATILE objects"""

    def __init__(
        self, seqrepo_access: SeqRepoAccess, tokenizer: Tokenize, classifier: Classify,
        validator: Validate, translator: Translate, hgvs_dup_del_mode: HGVSDupDelMode,
        gene_normalizer: GeneQueryHandler, transcript_mappings: TranscriptMappings
    ) -> None:
        """Initialize the ToVRSATILE class.

        :param SeqRepoAccess seqrepo_access: Access to SeqRepo
        :param Tokenize tokenizer: Tokenizer class for tokenizing
        :param Classify classifier: Classifier class for classifying tokens
        :param Validate validator: Validator class for validating valid inputs
        :param Translate translator: Translating valid inputs
        :param HGVSDupDelMode hgvs_dup_del_mode: Class for handling
            HGVS dup/del expressions
        :param GeneQueryHandler gene_normalizer: Client for normalizing gene concepts
        """
        super().__init__(
            seqrepo_access, tokenizer, classifier, validator, translator,
            hgvs_dup_del_mode, gene_normalizer
        )
        self.transcript_mappings = transcript_mappings

        # AC to Gene mapping
        self._protein_ac_to_gene_mapping = [
            self.transcript_mappings.get_gene_symbol_from_ensembl_protein,
            self.transcript_mappings.get_gene_symbol_from_refeq_protein
        ]
        self._cdna_ac_to_gene_mapping = [
            self.transcript_mappings.get_gene_symbol_from_refseq_rna,
            self.transcript_mappings.get_gene_symbol_from_ensembl_transcript
        ]

    def get_variation_descriptor(
        self, label: str, translation_result: TranslationResult,
        valid_result: ValidationResult, _id: str, warnings: List,
        gene: Optional[str] = None
    ) -> Tuple[Optional[VariationDescriptor], List]:
        """Return variation descriptor and warnings

        :param str label: Initial input query
        :param Dict variation: VRS variation object
        :param ValidationResult valid_result: Valid result for query
        :param str _id: _id field for variation descriptor
        :param List warnings: List of warnings
        :param Optional[str] gene: Gene symbol
        :return: Variation descriptor, warnings
        """
        try:
            variation = translation_result.vrs_variation
        except AttributeError as e:
            warnings.append(str(e))
            return None, warnings

        variation_id = variation["_id"]

        classification_type = valid_result.classification.classification_type
        if classification_type not in {ClassificationType.GENOMIC_DELETION_AMBIGUOUS,
                                       ClassificationType.GENOMIC_DUPLICATION_AMBIGUOUS,
                                       ClassificationType.AMPLIFICATION}:
            variation_type = variation["type"]
            if variation_type in {VRSTypes.ALLELE.value,
                                  VRSTypes.COPY_NUMBER_COUNT.value,
                                  VRSTypes.COPY_NUMBER_CHANGE.value}:
                key = "location" if variation_type == VRSTypes.ALLELE else "subject"
                vrs_ref_allele_seq = self.get_ref_allele_seq(
                    variation[key], translation_result.vrs_seq_loc_ac
                )
        else:
            vrs_ref_allele_seq = None

        molecule_context = valid_result.classification.molecule_context

        if valid_result.classification.gene_token:
            gene_context = self.get_gene_descriptor(
                gene_token=valid_result.classification.gene_token
            )
        else:
            if gene:
                gene_context = self.get_gene_descriptor(gene=gene)
            else:
                gene_context = None
                if valid_result.classification.nomenclature == Nomenclature.HGVS:
                    ac = valid_result.classification.ac
                    gene_token = self._get_hgvs_gene_context(ac, molecule_context)
                    if gene_token:
                        gene_context = self.get_gene_descriptor(gene_token=gene_token)

        return VariationDescriptor(
            id=_id,
            label=label,
            variation_id=variation_id,
            variation=variation,
            molecule_context=molecule_context,
            structural_type=valid_result.classification.so_id,
            vrs_ref_allele_seq=vrs_ref_allele_seq if vrs_ref_allele_seq else None,
            gene_context=gene_context
        ), warnings

    def _get_hgvs_gene_context(
        self, accession: str, molecule_context: MoleculeContext
    ) -> Optional[GeneToken]:
        def _get_gene_token(ac: str, mappings: List[Callable]) -> Optional[GeneToken]:
            gene_token = None
            for mapping in mappings:
                gene_symbol = mapping(ac)
                if gene_symbol:
                    gene_token = self.tokenizer.gene_symbol.match(gene_symbol)
                    if gene_token:
                        break
            return gene_token

        if molecule_context == MoleculeContext.PROTEIN:
            gene_token = _get_gene_token(accession, self._protein_ac_to_gene_mapping)
        elif molecule_context == MoleculeContext.TRANSCRIPT:
            gene_token = _get_gene_token(accession, self._cdna_ac_to_gene_mapping)
        else:
            gene_token = None
        return gene_token

    def get_gene_descriptor(
        self, gene_token: Optional[GeneToken] = None, gene: Optional[str] = None
    ) -> Optional[GeneDescriptor]:
        """Return a GA4GH Gene Descriptor using Gene Normalization.

        :param Optional[GeneToken] gene_token: A gene token
        :param Optional[str] gene: Gene query
        :return: A gene descriptor for a given gene if a record exists in
            gene-normalizer.
        """
        if gene_token:
            if gene_token.gene_descriptor:
                return gene_token.gene_descriptor

            gene_query = gene_token.matched_value
        elif gene:
            gene_query = gene
        else:
            return None

        response = self.gene_normalizer.normalize(gene_query)
        if response.gene_descriptor:
            return response.gene_descriptor
        else:
            return None
