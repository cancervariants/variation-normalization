"""Module for Variant Normalization."""
from variant.schemas.token_response_schema import PolypeptideSequenceVariant,\
    SingleNucleotideVariant, DelIns
from variant.schemas.ga4gh_vod import Gene, VariationDescriptor, GeneDescriptor
from variant.data_sources import SeqRepoAccess
from gene.query import QueryHandler as GeneQueryHandler
from urllib.parse import quote
from variant import logger


class Normalize:
    """The Normalize class used to normalize a given variant."""

    def __init__(self):
        """Initialize Normalize class."""
        self.gene_query_handler = GeneQueryHandler()
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
            allele_id = allele['_id']
            del allele['_id']
            molecule_context, structural_type, ref_allele_seq = \
                self._get_molecule_context_structural_type_ref_allele_seq(
                    valid_result_tokens, amino_acid_cache, label, allele)

            if valid_result.gene_tokens:
                gene_token = valid_result.gene_tokens[0]
                gene_context = self.get_gene_descriptor(gene_token)
            else:
                # TODO: Find gene context for genomic substitution
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
        response = self.gene_query_handler.search_sources(gene_symbol,
                                                          incl='hgnc')
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

    def _get_molecule_context_structural_type_ref_allele_seq(self,
                                                             valid_result_tokens,  # noqa: E501
                                                             amino_acid_cache,
                                                             label, allele):
        """Return context for a token.

        :return: (molecule_context, structural_type, ref_allele_seq)
        """
        polypeptide_sequence_variant_token = \
            self._get_instance_type_token(valid_result_tokens,
                                          PolypeptideSequenceVariant)
        dna_sequence_variant_token = \
            self._get_instance_type_token(valid_result_tokens,
                                          SingleNucleotideVariant)

        delins_token = self._get_instance_type_token(valid_result_tokens,
                                                     DelIns)

        if polypeptide_sequence_variant_token and not \
                dna_sequence_variant_token:
            molecule_context = 'protein'

            if self._is_token_type(valid_result_tokens,
                                   'AminoAcidSubstitution'):
                structural_type = 'SO:0001606'
            elif self._is_token_type(valid_result_tokens,
                                     'PolypeptideTruncation'):
                structural_type = 'SO:0001617'
            elif self._is_token_type(valid_result_tokens,
                                     'SilentMutation'):
                structural_type = 'SO:0001017'
            else:
                structural_type = None

            # convert 3 letter to 1 letter amino acid code
            if len(polypeptide_sequence_variant_token.ref_protein) == 3:
                for one, three in \
                        amino_acid_cache.amino_acid_code_conversion.items():
                    if three == polypeptide_sequence_variant_token.ref_protein:
                        polypeptide_sequence_variant_token.ref_protein = one
            ref_allele_seq = polypeptide_sequence_variant_token.ref_protein
        elif dna_sequence_variant_token:
            molecule_context = 'genomic'
            if self._is_token_type(valid_result_tokens,
                                   'CodingDNASubstitution') or \
                    self._is_token_type(valid_result_tokens,
                                        'GenomicSubstitution'):
                structural_type = 'SO:0001483'
            else:
                structural_type = None
            ref_allele_seq = dna_sequence_variant_token.ref_nucleotide
        elif delins_token:
            if delins_token.reference_sequence in ['c', 'g']:
                molecule_context = 'genomic'
            else:
                # TODO
                molecule_context = None
            structural_type = 'SO:1000032'
            ref_allele_seq = self.get_delins_ref_allele_seq(allele,
                                                            label)
        else:
            molecule_context = None
            structural_type = None
            ref_allele_seq = None
        return molecule_context, structural_type, ref_allele_seq

    def get_delins_ref_allele_seq(self, allele, label):
        """Return ref allele seq for transcript.

        :param dict allele: VRS Allele object
        :param str label: Transcript label
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

    def _is_token_type(self, valid_result_tokens, token_type):
        """Return whether or not token_type is in valid_result_tokens."""
        for t in valid_result_tokens:
            if t.token_type == token_type:
                return True
        return False

    def _get_instance_type_token(self, valid_result_tokens, instance_type):
        """Return the tokens for a given instance type.

        :return: A list of tokens
        """
        for t in valid_result_tokens:
            if isinstance(t, instance_type):
                return t
        return None
