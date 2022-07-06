"""Module for representing VRSATILE objects"""
from typing import Dict, List, Optional, Tuple

from ga4gh.vrsatile.pydantic.vrs_models import VRSTypes
from ga4gh.vrsatile.pydantic.vrsatile_models import VariationDescriptor, GeneDescriptor
from uta_tools.data_sources import SeqRepoAccess
from ga4gh.vrs.dataproxy import SeqRepoDataProxy
from gene.query import QueryHandler as GeneQueryHandler

from variation.classifiers.classify import Classify
from variation.hgvs_dup_del_mode import HGVSDupDelMode
from variation.to_vrs import ToVRS
from variation.schemas.token_response_schema import GeneMatchToken
from variation.schemas.validation_response_schema import ValidationResult
from variation.tokenizers.tokenize import Tokenize
from variation.translators.translate import Translate
from variation.validators.validate import Validate


class ToVRSATILE(ToVRS):
    """Class for represnting VRSATILE objects"""

    def __init__(self, seqrepo_access: SeqRepoAccess, dp: SeqRepoDataProxy,
                 tokenizer: Tokenize, classifier: Classify, validator: Validate,
                 translator: Translate, hgvs_dup_del_mode: HGVSDupDelMode,
                 gene_normalizer: GeneQueryHandler) -> None:
        """Initialize the ToVRSATILE class

        :param SeqRepoAccess seqrepo_access: Access to SeqRepo via UTA Tools
        :param SeqRepoDataProxy dp: Access to SeqRepo via VRS Python
        :param Tokenize tokenizer: Tokenizer class for tokenizing
        :param Classify classifier: Classifier class for classifying tokens
        :param Validate validator: Validator class for validating valid inputs
        :param Translate translator: Translating valid inputs
        :param HGVSDupDelMode hgvs_dup_del_mode: Class for handling
            HGVS dup/del expressions
        :parm GeneQueryHandler gene_normalizer: Client for normalizing gene concepts
        """
        super().__init__(seqrepo_access, dp, tokenizer, classifier, validator,
                         translator, hgvs_dup_del_mode)
        self.gene_normalizer = gene_normalizer

    def get_variation_descriptor(
            self, label: str, variation: Dict, valid_result: ValidationResult,
            _id: str, warnings: List, gene: Optional[str] = None
    ) -> Tuple[VariationDescriptor, List]:
        """Return variation descriptor and warnings

        :param str label: Initial input query
        :param Dict variation: VRS variation object
        :param ValidationResult valid_result: Valid result for query
        :param str _id: _id field for variation descriptor
        :param List warnings: List of warnings
        :param Optional[str] gene: Gene symbol
        :return: Variation descriptor, warnings
        """
        variation_id = variation["_id"]
        identifier = valid_result.identifier
        token_type = valid_result.classification_token.token_type.lower()

        vrs_ref_allele_seq = None
        if "uncertain" in token_type:
            warnings = ["Ambiguous regions cannot be normalized"]
        elif "range" not in token_type:
            if variation["type"] == VRSTypes.ALLELE.value:
                vrs_ref_allele_seq = self.get_ref_allele_seq(
                    variation, identifier
                )
            elif variation["type"] in [VRSTypes.ABSOLUTE_COPY_NUMBER.value,
                                       VRSTypes.RELATIVE_COPY_NUMBER.value]:
                loc = {
                    "location": variation["subject"]
                }
                vrs_ref_allele_seq = self.get_ref_allele_seq(loc, identifier)

        if valid_result.gene_tokens:
            gene_token = valid_result.gene_tokens[0]
            gene_context = self.get_gene_descriptor(gene_token=gene_token)
        else:
            if gene:
                gene_context = self.get_gene_descriptor(gene=gene)
            else:
                gene_context = None

        return VariationDescriptor(
            id=_id,
            label=label,
            variation_id=variation_id,
            variation=variation,
            molecule_context=valid_result.classification_token.molecule_context,
            structural_type=valid_result.classification_token.so_id,
            vrs_ref_allele_seq=vrs_ref_allele_seq if vrs_ref_allele_seq else None,
            gene_context=gene_context,
            type="VariationDescriptor"
        ), warnings

    def get_gene_descriptor(
            self, gene_token: Optional[GeneMatchToken] = None,
            gene: Optional[str] = None) -> Optional[GeneDescriptor]:
        """Return a GA4GH Gene Descriptor using Gene Normalization.

        :param Optional[GeneMatchToken] gene_token: A gene token
        :param Optional[str] gene: Gene query
        :return: A gene descriptor for a given gene if a record exists in
            gene-normalizer.
        """
        if gene_token:
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
