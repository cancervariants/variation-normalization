"""Module for Variant Normalization."""
from typing import Optional
from variant.schemas.token_response_schema import PolypeptideSequenceVariant
from variant import GENE_NORMALIZER
from variant.schemas.ga4gh_vod import Gene, VariationDescriptor, GeneDescriptor
from variant.data_sources import SeqRepoAccess
from urllib.parse import quote
from variant import logger


class Normalize:
    """The Normalize class used to normalize a given variant."""

    def __init__(self):
        """Initialize Normalize class."""
        self.seqrepo_access = SeqRepoAccess()
        self.warnings = list()

    def normalize(self, q, validations, amino_acid_cache):
        """Normalize a given variant.

        :param str q: The variant to normalize
        :param ValidationSummary validations: Invalid and valid results
        :param AminoAcidCache amino_acid_cache: Amino Acid Code and Conversion
        :return: An allele descriptor for a valid result if one exists. Else,
            None.
        """
        warnings = list()
        if len(validations.valid_results) > 0:
            # For now, only use first valid result
            valid_result = None
            label = None
            for r in validations.valid_results:
                if r.mane_transcript:
                    valid_result = r
                    label = valid_result.mane_transcript.strip()
                    break
            if not valid_result:
                warning = f"Unable to find MANE Select Transcript for {q}."
                logger.warning(warning)
                warnings.append(warning)
                valid_result = validations.valid_results[0]
                label = ' '.join(q.strip().split())

            valid_result_tokens = valid_result.classification.all_tokens
            allele = valid_result.allele
            allele_id = allele.pop('_id')
            variant_token = None
            molecule_context = None
            structural_type = None

            for token in valid_result_tokens:
                try:
                    molecule_context = token.molecule_context
                    structural_type = token.so_id
                    variant_token = token
                except AttributeError:
                    continue

            ref_allele_seq = self._get_ref_allele_seq(
                variant_token, amino_acid_cache, label, allele
            )

            if valid_result.gene_tokens:
                gene_token = valid_result.gene_tokens[0]
                gene_context = self.get_gene_descriptor(gene_token)
            else:
                gene_context = None

            variation_descriptor = VariationDescriptor(
                id=f"normalize.variant:{quote(' '.join(q.strip().split()))}",
                value_id=allele_id,
                label=label,
                value=allele,
                molecule_context=molecule_context,
                structural_type=structural_type,
                ref_allele_seq=ref_allele_seq,
                gene_context=gene_context
            )
        else:
            variation_descriptor = None
            warning = f"Unable to normalize {q}."
            logger.warning(warning)
            warnings.append(warning)

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
                xrefs=record.other_identifiers,
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
        if record.strand:
            self.add_extension(extensions, 'strand', record.strand)
        if record.symbol_status:
            self.add_extension(extensions, 'symbol_status',
                               record.symbol_status)
        self.add_extension(extensions, 'associated_with', record.xrefs)
        self.add_extension(extensions, 'chromosome_location',
                           record_location.dict(by_alias=True))
        return extensions

    def add_extension(self, extensions, name, value):
        """Add extension to list of extensions.

        :param list extensions: List of ga4gh extensions
        :param str name: name of extension
        :param str value: value of extension
        """
        extensions.append({
            'type': 'Extension',
            'name': name,
            'value': value
        })

    def _get_ref_allele_seq(self, variant_token, amino_acid_cache,
                            label, allele) -> str:
        """Return ref_allele_seq for a variant.

        :return: ref_allele_seq
        """
        if variant_token.token_type == PolypeptideSequenceVariant:
            # convert 3 letter to 1 letter amino acid code
            if len(variant_token.ref_protein) == 3:
                for one, three in \
                        amino_acid_cache.amino_acid_code_conversion.items():
                    if three == variant_token.ref_protein:
                        variant_token.ref_protein = one
            ref_allele_seq = variant_token.ref_protein
        else:
            if variant_token.token_type in ['CodingDNASubstitution',
                                            'GenomicSubstitution']:
                ref_allele_seq = variant_token.ref_nucleotide
            else:
                ref_allele_seq = self.get_delins_ref_allele_seq(allele, label)
        return ref_allele_seq

    def get_delins_ref_allele_seq(self, allele, label) -> Optional[str]:
        """Return ref allele seq for transcript.

        :param dict allele: VRS Allele object
        :param str label: Transcript label
        :return: Ref seq allele
        """
        label = label.split(':')[0]
        interval = allele['location']['interval']
        start = interval['start'] + 1
        end = interval['end']

        if start and end:
            refseq_list = list()
            while start <= end:
                refseq_list.append(self.seqrepo_access.sequence_at_position(
                    label, start
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
