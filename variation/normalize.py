"""Module for Variation Normalization."""
from typing import Optional, List, Tuple
from ga4gh.vrsatile.pydantic.vrsatile_model import VariationDescriptor
from ga4gh.vrsatile.pydantic.vrs_model import Text
from variation.data_sources import SeqRepoAccess, UTA
from urllib.parse import quote
from variation import logger
from gene.query import QueryHandler as GeneQueryHandler


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

    def normalize(self, q, validations, warnings):
        """Normalize a given variation.

        :param str q: The variation to normalize
        :param ValidationSummary validations: Invalid and valid results
        :param list warnings: List of warnings
        :return: An variation descriptor for a valid result if one exists.
            Else, None.
        """
        if not q:
            resp, warnings = self._no_variation_entered()
        else:
            _id = f"normalize.variation:{quote(' '.join(q.strip().split()))}"
            if len(validations.valid_results) > 0:
                # For now, only use first valid result
                valid_result = None
                for r in validations.valid_results:
                    if r.is_mane_transcript and r.variation:
                        valid_result = r
                        break
                if not valid_result:
                    warning = f"Unable to find MANE Select Transcript for {q}."
                    logger.warning(warning)
                    warnings.append(warning)
                    valid_result = validations.valid_results[0]

                variation = valid_result.variation

                variation_id = variation.pop('_id')
                identifier = valid_result.identifier

                if variation['type'] == 'Allele':
                    vrs_ref_allele_seq = self.get_ref_allele_seq(
                        variation, identifier
                    )
                elif variation['type'] == 'CopyNumber':
                    vrs_ref_allele_seq = self.get_ref_allele_seq(
                        variation['subject'], identifier
                    )
                else:
                    vrs_ref_allele_seq = None

                if valid_result.gene_tokens:
                    gene_token = valid_result.gene_tokens[0]
                    gene_context = self.get_gene_descriptor(gene_token)
                else:
                    gene_context = None

                if 'Uncertain' in valid_result.classification_token.token_type:
                    warnings = ['Ambiguous regions cannot be normalized']

                resp = VariationDescriptor(
                    id=_id,
                    variation_id=variation_id,
                    variation=variation,
                    molecule_context=valid_result.classification_token.molecule_context,  # noqa: E501
                    structural_type=valid_result.classification_token.so_id,
                    vrs_ref_allele_seq=vrs_ref_allele_seq,
                    gene_context=gene_context
                )
            else:
                if not q.strip():
                    resp, warnings = self._no_variation_entered()
                else:
                    warning = f"Unable to normalize {q}"
                    resp = VariationDescriptor(
                        id=_id,
                        variation=Text(definition=q)
                    )
                    if not warnings:
                        warnings.append(warning)
                    logger.warning(warning)
        self.warnings = warnings
        return resp

    def _no_variation_entered(self) -> Tuple[None, List[str]]:
        """Return response when no variation queried.

        :return: None, list of warnings
        """
        warnings = ["No variation was entered to normalize"]
        logger.warning(warnings)
        return None, warnings

    def get_gene_descriptor(self, gene_token):
        """Return a GA4GH Gene Descriptor using Gene Normalization.

        :param GeneMatchToken gene_token: A gene token
        :return: A gene descriptor for a given gene if a record exists in
            gene-normalizer.
        """
        gene_symbol = gene_token.matched_value
        if gene_symbol in self._gene_norm_cache:
            return self._gene_norm_cache[gene_symbol]
        else:
            response = self.gene_normalizer.normalize(gene_symbol)
            if 'gene_descriptor' in response and response['gene_descriptor']:
                gene_descriptor = response['gene_descriptor']
                self._gene_norm_cache[gene_symbol] = gene_descriptor
                return gene_descriptor
            return None

    def get_ref_allele_seq(self, allele, identifier) -> Optional[str]:
        """Return ref allele seq for transcript.

        :param dict allele: VRS Allele object
        :param str identifier: Identifier for allele
        :return: Ref seq allele
        """
        start = None
        end = None
        interval = allele['location']['interval']
        ival_type = interval['type']
        if ival_type == 'SimpleInterval':
            if interval['start'] != interval['end']:
                start = interval['start'] + 1
                end = interval['end']
        elif ival_type == 'SequenceInterval':
            if interval['start']['type'] == 'Number':
                start = interval['start']['value'] + 1
                end = interval['end']['value']

        if start is None and end is None:
            return None

        return self.seqrepo_access.get_sequence(identifier, start, end)

    def _is_token_type(self, valid_result_tokens, token_type) -> bool:
        """Return whether or not token_type is in valid_result_tokens.

        :param list valid_result_tokens: Valid token matches
        :param str token_type: The token's type
        :return: Whether or not token_type is in valid_result_tokens
        """
        for t in valid_result_tokens:
            if t.token_type == token_type:
                return True
        return False

    def _get_instance_type_token(self, valid_result_tokens, instance_type):
        """Return the tokens for a given instance type.

        :param list valid_result_tokens: A list of valid tokens for the input
            string
        :param Token instance_type: The instance type to check
        :return: Token for a given instance type
        """
        for t in valid_result_tokens:
            if isinstance(t, instance_type):
                return t
        return None
