"""Base class for variation normalization"""
from typing import Tuple, Optional, List, Dict
from ga4gh.vrsatile.pydantic.vrs_models import Text
from ga4gh.vrsatile.pydantic.vrsatile_models import VariationDescriptor, \
    GeneDescriptor
from variation.hgvs_dup_del_mode import HGVSDupDelMode
from variation.schemas.token_response_schema import GeneMatchToken
from variation.schemas.validation_response_schema import ValidationResult
from variation.classifiers import Classify
from variation.tokenizers import Tokenize
from variation.validators import Validate
from variation.translators import Translate
from variation.data_sources import SeqRepoAccess, TranscriptMappings, \
    UTA, MANETranscriptMappings
from variation.mane_transcript import MANETranscript
from variation.tokenizers import GeneSymbol
from variation.tokenizers.caches import AminoAcidCache
from gene.query import QueryHandler as GeneQueryHandler
from ga4gh.vrs import models
from ga4gh.core import ga4gh_identify
from variation import logger


class Base:
    """Base class"""

    def __init__(self, tokenizer: Tokenize, classifier: Classify,
                 seqrepo_access: SeqRepoAccess,
                 transcript_mappings: TranscriptMappings,
                 gene_symbol: GeneSymbol, amino_acid_cache: AminoAcidCache,
                 uta: UTA, mane_transcript_mappings: MANETranscriptMappings,
                 mane_transcript: MANETranscript, validator: Validate,
                 translator: Translate,
                 gene_normalizer: GeneQueryHandler,
                 hgvs_dup_del_mode: HGVSDupDelMode) -> None:
        """Initialize the ToVRS class.

        :param Tokenize tokenizer: Tokenizer class for tokenizing
        :param Classify classifier: Classifier class for classifying tokens
        :param SeqRepoAccess seqrepo_access: Access to SeqRepo
        :param TranscriptMappings transcript_mappings: Transcript mappings
            data class
        :param GeneSymbol gene_symbol: Class for identifying gene symbols
        :param AminoAcidCache amino_acid_cache: Amino Acid data class
        :param UTA uta: UTA DB and queries
        :param MANETranscriptMappings mane_transcript_mappings: Class for
            getting mane transcript data from gene
        :param MANETranscript mane_transcript: Mane transcript data class
        :param Validate validator: Validator class for validating valid inputs
        :param Translate translator: Translating valid inputs
        :param GeneQueryHandler gene_normalizer: Gene normalizer access
        :param HGVSDupDelMode hgvs_dup_del_mode: Class for handling
            HGVS dup/del expressions
        """
        self.tokenizer = tokenizer
        self.classifier = classifier
        self.seqrepo_access = seqrepo_access
        self.transcript_mappings = transcript_mappings
        self.gene_symbol = gene_symbol
        self.amino_acid_cache = amino_acid_cache
        self.uta = uta
        self.mane_transcript_mappings = mane_transcript_mappings
        self.mane_transcript = mane_transcript
        self.validator = validator
        self.translator = translator
        self.gene_normalizer = gene_normalizer
        self.hgvs_dup_del_mode = hgvs_dup_del_mode
        self._gene_norm_cache = dict()
        self.warnings = list()

    @staticmethod
    def no_variation_entered() -> Tuple[None, List[str]]:
        """Return response when no variation queried.

        :return: None, list of warnings
        """
        warnings = ["No variation was entered to normalize"]
        logger.warning(warnings)
        return None, warnings

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
