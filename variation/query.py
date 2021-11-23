"""This module provides methods for handling queries."""
from typing import Tuple, Optional, List, Union
from gene.query import QueryHandler as GeneQueryHandler
from ga4gh.vrs.dataproxy import SeqRepoDataProxy
from ga4gh.vrs.extras.translator import Translator
from ga4gh.core import ga4gh_identify
from ga4gh.vrs import models
from variation import SEQREPO_DATA_PATH, TRANSCRIPT_MAPPINGS_PATH, \
    REFSEQ_GENE_SYMBOL_PATH, AMINO_ACID_PATH, UTA_DB_URL, REFSEQ_MANE_PATH
from variation.to_vrs import ToVRS
from variation.normalize import Normalize
from variation.classifiers import Classify
from variation.tokenizers import Tokenize
from variation.validators import Validate
from variation.translators import Translate
from variation.data_sources import SeqRepoAccess, TranscriptMappings, \
    UTA, MANETranscriptMappings
from variation.mane_transcript import MANETranscript
from variation.hgvs_dup_del_mode import HGVSDupDelMode
from variation.tokenizers import GeneSymbol
from variation.tokenizers.caches import AminoAcidCache
from ga4gh.vrsatile.pydantic.vrs_models import Text, Allele, CopyNumber, \
    Haplotype, VariationSet
from ga4gh.vrsatile.pydantic.vrsatile_models import VariationDescriptor
from variation.schemas.normalize_response_schema\
    import HGVSDupDelMode as HGVSDupDelModeEnum


class QueryHandler:
    """Class for handling queries."""

    def __init__(self,
                 dynamodb_url: str = '',
                 dynamodb_region: str = 'us-east-2',
                 seqrepo_data_path: str = SEQREPO_DATA_PATH,
                 uta_db_url: str = UTA_DB_URL,
                 uta_db_pwd: Optional[str] = None) -> None:
        """Initialize QueryHandler instance.
        :param str dynamodb_url: URL to gene-normalizer database source.
        :param str dynamodb_region: AWS default region for gene-normalizer.
        :param str seqrepo_data_path: Path to seqrepo data directory
        :param str uta_db_url: URL for UTA database
        :param Optional[str] uta_db_pwd: Password for UTA database user
        """
        self.gene_normalizer = GeneQueryHandler(db_url=dynamodb_url,
                                                db_region=dynamodb_region)
        self.seqrepo_access = SeqRepoAccess(
            seqrepo_data_path=seqrepo_data_path
        )
        self.uta = UTA(db_url=uta_db_url, db_pwd=uta_db_pwd)
        self.dp = SeqRepoDataProxy(self.seqrepo_access.seq_repo_client)
        self.hgvs_dup_del_mode = HGVSDupDelMode(self.seqrepo_access)
        self.to_vrs_handler = self._init_to_vrs()
        self.normalize_handler = Normalize(
            self.seqrepo_access, self.uta, self.gene_normalizer
        )

    def _init_to_vrs(self,
                     transcript_file_path: str = TRANSCRIPT_MAPPINGS_PATH,
                     refseq_file_path: str = REFSEQ_GENE_SYMBOL_PATH,
                     amino_acids_file_path: str = AMINO_ACID_PATH,
                     mane_data_path: str = REFSEQ_MANE_PATH) -> ToVRS:
        """Return toVRS instance

        :param str transcript_file_path: Path to transcript mappings file
        :param str refseq_file_path: Path to refseq gene symbol file
        :param str amino_acids_file_path: Path to amino acids file
        :param str mane_data_path: Path to refseq mane data file
        :return: toVRS instance
        """
        gene_symbol = GeneSymbol(self.gene_normalizer)
        amino_acid_cache = AminoAcidCache(
            amino_acids_file_path=amino_acids_file_path
        )
        tokenizer = Tokenize(amino_acid_cache, gene_symbol)
        classifier = Classify()
        transcript_mappings = TranscriptMappings(
            transcript_file_path=transcript_file_path,
            refseq_file_path=refseq_file_path
        )
        mane_transcript_mappings = MANETranscriptMappings(
            mane_data_path=mane_data_path
        )

        tlr = Translator(data_proxy=self.dp)
        mane_transcript = MANETranscript(
            self.seqrepo_access, transcript_mappings,
            mane_transcript_mappings, self.uta
        )
        validator = Validate(
            self.seqrepo_access, transcript_mappings, gene_symbol,
            mane_transcript, self.uta, self.dp, tlr,
            amino_acid_cache, self.gene_normalizer
        )
        translator = Translate()
        return ToVRS(
            tokenizer, classifier, self.seqrepo_access, transcript_mappings,
            gene_symbol, amino_acid_cache, self.uta, mane_transcript_mappings,
            mane_transcript, validator, translator, self.gene_normalizer
        )

    def to_vrs(self, q: str)\
            -> Tuple[Optional[Union[List[Allele], List[CopyNumber],
                                    List[Text], List[Haplotype],
                                    List[VariationSet]]],
                     Optional[List[str]]]:
        """Return a VRS-like representation of all validated variations for a query.  # noqa: E501

        :param str q: The variation to translate
        :return: Validated VRS Variations and list of warnings
        """
        validations, warnings = \
            self.to_vrs_handler.get_validations(q)
        translations, warnings = \
            self.to_vrs_handler.get_translations(validations, warnings)

        if not translations:
            if q and q.strip():
                text = models.Text(definition=q)
                text._id = ga4gh_identify(text)
                translations = [Text(**text.as_dict())]
            else:
                translations = None
        return translations, warnings

    def normalize(
            self, q: str,
            hgvs_dup_del_mode: HGVSDupDelModeEnum = HGVSDupDelModeEnum.DEFAULT
    ) -> Optional[VariationDescriptor]:
        """Return normalized Variation Descriptor for variation.

        :param q: Variation to normalize
        :param HGVSDupDelModeEnum hgvs_dup_del_mode: Must be: `default`, `cnv`,
            `repeated_seq_expr`, `literal_seq_expr`.
            This parameter determines how to interpret HGVS dup/del expressions
            in VRS.
        :return: Variation Descriptor for variation
        """
        if hgvs_dup_del_mode:
            hgvs_dup_del_mode = hgvs_dup_del_mode.strip().lower()
            if not self.hgvs_dup_del_mode.is_valid_mode(hgvs_dup_del_mode):
                self.normalize_handler.warnings = \
                    [f"hgvs_dup_del_mode must be one of: "
                     f"{self.hgvs_dup_del_mode.valid_modes}"]
                return None
            else:
                validations, warnings = \
                    self.to_vrs_handler.get_validations(
                        q, normalize_endpoint=True,
                        hgvs_dup_del_mode=hgvs_dup_del_mode
                    )
                return self.normalize_handler.normalize(q, validations,
                                                        warnings)
        else:
            self.normalize_handler.warnings = \
                ["hgvs_dup_del_mode cannot be None"]
            return None
