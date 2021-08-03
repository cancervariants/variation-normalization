"""Module for Variation Normalization."""
from typing import Optional
from variation import GENE_NORMALIZER
from variation.schemas.ga4gh_vod import Gene, VariationDescriptor, \
    GeneDescriptor
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
            ref_allele_seq = self.get_ref_allele_seq(
                allele, identifier
            )

            if valid_result.gene_tokens:
                gene_token = valid_result.gene_tokens[0]
                gene_context = self.get_gene_descriptor(gene_token)
            else:
                gene_context = None

            variation_descriptor = VariationDescriptor(
                id=f"normalize.variation:{quote(' '.join(q.strip().split()))}",
                value_id=allele_id,
                value=allele,
                molecule_context=valid_result.classification_token.molecule_context,  # noqa: E501
                structural_type=valid_result.classification_token.so_id,
                ref_allele_seq=ref_allele_seq,
                gene_context=gene_context
            )
        else:
            variation_descriptor = None
            warning = f"Unable to normalize {q}."
            if not warnings:
                warnings.append(warning)
            logger.warning(warning)

        self.warnings = warnings
        return variation_descriptor

    def get_gene_descriptor(self, gene_token):
        """Return a GA4GH Gene Descriptor using Gene Normalization.

        :param GeneMatchToken gene_token: A gene token
        :return: A gene descriptor for a given gene if a record exists in
            gene-normalizer.
        """
        gene_symbol = gene_token.matched_value
        response = GENE_NORMALIZER.search_sources(gene_symbol, incl='hgnc')
        if response['source_matches'][0]['records']:
            record = response['source_matches'][0]['records'][0]
            record_location = record.locations[0] if record.locations else None

            return GeneDescriptor(
                id=f"normalize.gene:{quote(' '.join(gene_symbol.strip().split()))}",  # noqa: E501
                label=gene_symbol,
                value=Gene(id=record.concept_id),
                xrefs=record.xrefs if record.xrefs else [],
                alternate_labels=[record.label] + record.aliases + record.previous_symbols,  # noqa: E501
                extensions=self.get_extensions(record, record_location)
            )
        return None

    def get_extensions(self, record, record_location):
        """Return a list of ga4gh extensions.

        :param gene.schemas.Gene record: The record from the normalization
            service
        :param gene.schemas.ChromosomeLocation record_location: The record's
            location
        :return: List of extensions providing additional information
        """
        extensions = list()
        self.add_extension(extensions, 'strand', record.strand)
        self.add_extension(extensions, 'symbol_status', record.symbol_status)
        self.add_extension(extensions, 'associated_with',
                           record.associated_with)
        self.add_extension(extensions, 'chromosome_location',
                           record_location.dict(by_alias=True))
        return extensions

    def add_extension(self, extensions, name, value):
        """Add extension to list of extensions.

        :param list extensions: List of ga4gh extensions
        :param str name: name of extension
        :param any value: value of extension
        """
        if value:
            extensions.append({
                'type': 'Extension',
                'name': name,
                'value': value
            })

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
