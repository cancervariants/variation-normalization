"""Module for Variation Normalization."""
from typing import Optional, List, Tuple
from variation import GENE_NORMALIZER
from variation.schemas.ga4gh_vrsatile import VariationDescriptor
from variation.schemas.ga4gh_vrs import Text
from variation.data_sources import SeqRepoAccess, UTA
from urllib.parse import quote
from variation import logger


class Normalize:
    """The Normalize class used to normalize a given variation."""

    def __init__(self, seqrepo_access: SeqRepoAccess, uta: UTA):
        """Initialize Normalize class."""
        self.seqrepo_access = seqrepo_access
        self.uta = uta
        self.warnings = list()

    def normalize(self, q, validations, warnings):
        """Normalize a given variation.

        :param str q: The variation to normalize
        :param ValidationSummary validations: Invalid and valid results
        :param list warnings: List of warnings
        :return: An allele descriptor for a valid result if one exists. Else,
            None.
        """
        if not q:
            resp, warnings = self._no_variation_entered()
        else:
            _id = f"normalize.variation:{quote(' '.join(q.strip().split()))}"
            if len(validations.valid_results) > 0:
                # For now, only use first valid result
                valid_result = None
                for r in validations.valid_results:
                    if r.is_mane_transcript and r.allele:
                        valid_result = r
                        break
                if not valid_result:
                    warning = f"Unable to find MANE Select Transcript for {q}."
                    logger.warning(warning)
                    warnings.append(warning)
                    valid_result = validations.valid_results[0]

                allele = valid_result.allele
                allele_id = allele.pop('_id')
                identifier = valid_result.identifier
                vrs_ref_allele_seq = self.get_ref_allele_seq(
                    allele, identifier
                )

                if valid_result.gene_tokens:
                    gene_token = valid_result.gene_tokens[0]
                    gene_context = self.get_gene_descriptor(gene_token)
                else:
                    gene_context = None

                resp = VariationDescriptor(
                    id=_id,
                    value_id=allele_id,
                    value=allele,
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
                        value=Text(definition=q)
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
        response = GENE_NORMALIZER.normalize(gene_symbol)
        if 'gene_descriptor' in response and response['gene_descriptor']:
            return response['gene_descriptor']
        return None

    def get_ref_allele_seq(self, allele, identifier) -> Optional[str]:
        """Return ref allele seq for transcript.

        :param dict allele: VRS Allele object
        :param str identifier: Identifier for allele
        :return: Ref seq allele
        """
        interval = allele['location']['interval']
        if interval['start'] != interval['end']:
            start = interval['start'] + 1
            end = interval['end']
        else:
            return None

        if start and end:
            refseq_list = list()
            while start <= end:
                refseq_list.append(self.seqrepo_access.sequence_at_position(
                    identifier, start
                ))
                start += 1
            try:
                return ''.join(refseq_list)
            except TypeError:
                pass
        return None

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
