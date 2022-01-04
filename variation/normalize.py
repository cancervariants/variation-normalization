"""Module for Variation Normalization."""
from typing import Optional, List, Tuple, Dict
from ga4gh.vrsatile.pydantic.vrs_models import Text
from ga4gh.vrsatile.pydantic.vrsatile_models import VariationDescriptor, \
    GeneDescriptor
from ga4gh.vrs import models
from ga4gh.core import ga4gh_identify
from variation.data_sources import SeqRepoAccess, UTA
from urllib.parse import quote
from variation import logger
from gene.query import QueryHandler as GeneQueryHandler
from variation.schemas.token_response_schema import GeneMatchToken, Token
from variation.schemas.validation_response_schema import ValidationSummary, \
    ValidationResult


class Normalize:
    """The Normalize class used to normalize a given variation."""

    def __init__(self, seqrepo_access: SeqRepoAccess, uta: UTA,
                 gene_normalizer: GeneQueryHandler) -> None:
        """Initialize Normalize class.

        :param SeqRepoAccess seqrepo_access: Access to SeqRepo data queries
        :param UTA uta: Access to UTA database and queries
        :parm QueryHandler gene_normalizer: Access to gene-normalizer queries
        """
        self.seqrepo_access = seqrepo_access
        self.uta = uta
        self.warnings = list()
        self._gene_norm_cache = dict()
        self.gene_normalizer = gene_normalizer

    @staticmethod
    def get_valid_result(q: str, toVRSATILE_endpoint: bool,
                         validations: ValidationSummary,
                         warnings: List) -> ValidationResult:
        """Get valid result from ValidationSummary

        :param str q: Query string
        :param ValidationSummary validations: Validation summary for query
        :param List warnings: List of warnings
        :return: Valid Validation Result
        """
        # For now, only use first valid result
        valid_result = None
        listOf_valid_res = list()
        listOf_possible_ac = list()
        for r in validations.valid_results:
            if not toVRSATILE_endpoint:
                if r.is_mane_transcript and r.variation:
                    valid_result = r
                    break
            else:
                if r.is_valid and r.variation:
                    listOf_valid_res.append(r)
                    listOf_possible_ac.append(r.possible_ac[:])
        if not valid_result: # TODO: if not listOf_valid_res?
            warning = f"Unable to find MANE Transcript for {q}."
            logger.warning(warning)
            warnings.append(warning)
            valid_result = validations.valid_results[0]
        return valid_result, listOf_valid_res, listOf_possible_ac

    def normalize(self, q: str, toVRSATILE_endpoint: bool,
                  validations: ValidationSummary,
                  warnings: List) -> Optional[VariationDescriptor]:
        """Normalize a given variation.

        :param str q: The variation to normalize
        :param ValidationSummary validations: Invalid and valid results
        :param List warnings: List of warnings
        :return: An variation descriptor for a valid result if one exists.
            Else, None.
        """
        if not q:
            resp, warnings = self._no_variation_entered()
        else:
            _id = f"normalize.variation:{quote(' '.join(q.strip().split()))}"
            if len(validations.valid_results) > 0:
                valid_result, listOf_valid_res, listOf_possible_ac = \
                    self.get_valid_result(
                        q, toVRSATILE_endpoint, validations, warnings)
                resp, warnings = self.get_variation_descriptor(
                    valid_result.variation, valid_result, _id, warnings)

                # get list of toVRSATILE objects
                if toVRSATILE_endpoint:
                    listOf_resp, listOf_warnings = list(), list()
                    for val_res in listOf_valid_res:
                        toVRSATILE_resp, toVRSATILE_warnings = \
                            self.get_variation_descriptor(
                                val_res.variation, val_res,\
                                _id, toVRSATILE_warnings)
                        listOf_resp.append(toVRSATILE_resp)
                        listOf_warnings.append(toVRSATILE_warnings)
            else:
                if not q.strip():
                    resp, warnings = self._no_variation_entered()
                else:
                    resp, warnings = self.text_variation_resp(q, _id, warnings)

        self.warnings = warnings
        if toVRSATILE_endpoint:
            return listOf_resp, listOf_possible_ac, listOf_warnings
        return resp

    @staticmethod
    def text_variation_resp(
            q: str, _id: str,
            warnings: List) -> Tuple[VariationDescriptor, List]:
        """Return text variation for queries that could not be normalized

        :param str q: query
        :param str _id: _id field for variation descriptor
        :param List warnings: List of warnings
        :return: Variation descriptor, warnings
        """
        warning = f"Unable to normalize {q}"
        text = models.Text(definition=q)
        text._id = ga4gh_identify(text)
        resp = VariationDescriptor(
            id=_id,
            variation=Text(**text.as_dict())
        )
        if not warnings:
            warnings.append(warning)
            logger.warning(warning)
        return resp, warnings

    def get_variation_descriptor(
            self, variation: Dict, valid_result: ValidationResult,
            _id: str, warnings: List, gene: Optional[str] = None
    ) -> Tuple[VariationDescriptor, List]:
        """Return variation descriptor and warnings

        :param Dict variation: VRS variation object
        :param ValidationResult valid_result: Valid result for query
        :param str _id: _id field for variation descriptor
        :param List warnings: List of warnings
        :param Optional[str] gene: Gene symbol
        :return: Variation descriptor, warnings
        """
        variation_id = variation['_id']
        identifier = valid_result.identifier
        token_type = \
            valid_result.classification_token.token_type.lower()

        vrs_ref_allele_seq = None
        if 'uncertain' in token_type:
            warnings = ['Ambiguous regions cannot be normalized']
        elif 'range' not in token_type:
            if variation['type'] == 'Allele':
                vrs_ref_allele_seq = self.get_ref_allele_seq(
                    variation, identifier
                )
            elif variation['type'] == 'CopyNumber':
                vrs_ref_allele_seq = self.get_ref_allele_seq(
                    variation['subject'], identifier
                )

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
            variation_id=variation_id,
            variation=variation,
            molecule_context=valid_result.classification_token.molecule_context,  # noqa: E501
            structural_type=valid_result.classification_token.so_id,
            vrs_ref_allele_seq=vrs_ref_allele_seq if vrs_ref_allele_seq else None,  # noqa: E501
            gene_context=gene_context
        ), warnings

    def _no_variation_entered(self) -> Tuple[None, List[str]]:
        """Return response when no variation queried.

        :return: None, list of warnings
        """
        warnings = ["No variation was entered to normalize"]
        logger.warning(warnings)
        return None, warnings

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

        if gene_query in self._gene_norm_cache:
            return self._gene_norm_cache[gene_query]
        else:
            response = self.gene_normalizer.normalize(gene_query)
            if response.gene_descriptor:
                gene_descriptor = response.gene_descriptor
                self._gene_norm_cache[gene_query] = gene_descriptor
                return gene_descriptor
            return None

    def get_ref_allele_seq(self, allele: Dict,
                           identifier: str) -> Optional[str]:
        """Return ref allele seq for transcript.

        :param Dict allele: VRS Allele object
        :param str identifier: Identifier for allele
        :return: Ref seq allele
        """
        start = None
        end = None
        interval = allele['location']['interval']
        ival_type = interval['type']
        if ival_type == 'SequenceInterval':
            if interval['start']['type'] == 'Number':
                start = interval['start']['value'] + 1
                end = interval['end']['value']

        if start is None and end is None:
            return None

        return self.seqrepo_access.get_sequence(identifier, start, end)

    def _is_token_type(self, valid_result_tokens: List,
                       token_type: str) -> bool:
        """Return whether or not token_type is in valid_result_tokens.

        :param List valid_result_tokens: Valid token matches
        :param str token_type: The token's type
        :return: Whether or not token_type is in valid_result_tokens
        """
        for t in valid_result_tokens:
            if t.token_type == token_type:
                return True
        return False

    def _get_instance_type_token(self, valid_result_tokens: List,
                                 instance_type: Token) -> Optional[Token]:
        """Return the tokens for a given instance type.

        :param List valid_result_tokens: A list of valid tokens for the input
            string
        :param Token instance_type: The instance type to check
        :return: Token for a given instance type
        """
        for t in valid_result_tokens:
            if isinstance(t, instance_type):
                return t
        return None
