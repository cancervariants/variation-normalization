"""A module for testing validator classes."""
import yaml
from ga4gh.vrs.dataproxy import SeqRepoDataProxy
from ga4gh.vrs.extras.translator import Translator
from gene.query import QueryHandler as GeneQueryHandler
from uta_tools.data_sources import TranscriptMappings, SeqRepoAccess, \
    MANETranscriptMappings, UTADatabase, MANETranscript

from tests import PROJECT_ROOT
from variation.schemas.app_schemas import Endpoint
from variation.vrs_representation import VRSRepresentation
from variation.tokenizers import Tokenize, GeneSymbol
from variation.tokenizers.caches import AminoAcidCache


class ValidatorBase:
    """The validator base class."""

    @classmethod
    def setUpClass(cls):
        """Set up the test cases."""
        with open(f"{PROJECT_ROOT}/tests/fixtures/validators.yml") as stream:
            cls.all_fixtures = yaml.safe_load(stream)
        amino_acid_cache = AminoAcidCache()
        gene_normalizer = GeneQueryHandler()
        gene_symbol = GeneSymbol(gene_normalizer)
        cls.tokenizer = Tokenize(amino_acid_cache, gene_symbol)

        seqrepo_access = SeqRepoAccess()
        transcript_mappings = TranscriptMappings()
        uta = UTADatabase()
        dp = SeqRepoDataProxy(seqrepo_access.seqrepo_client)
        tlr = Translator(data_proxy=dp)
        mane_transcript = MANETranscript(
            seqrepo_access, transcript_mappings, MANETranscriptMappings(), uta,
            gene_normalizer)
        vrs = VRSRepresentation(dp, seqrepo_access)

        cls.aa_params = [
            seqrepo_access, transcript_mappings, gene_symbol,
            mane_transcript, uta, tlr, gene_normalizer,
            vrs, amino_acid_cache
        ]
        cls.params = cls.aa_params[:-1]

    def classifier_instance(self):
        """Check that the classifier_instance method is implemented."""
        raise NotImplementedError()

    def validator_instance(self, **kwargs):
        """Check that the validator_instance method is implemented."""
        raise NotImplementedError()

    def fixture_name(self):
        """Check that the fixture_name method is implemented."""
        raise NotImplementedError()

    def set_up(self):
        """Initialize fixtures, classifier + validator"""
        self.fixtures = self.all_fixtures.get(
            self.fixture_name(),
            {"should_match": [], "should_not_match": []}
        )
        self.classifier = self.classifier_instance()
        self.validator = self.validator_instance()

    async def test_matches(self):
        """Test that validator matches correctly."""
        self.set_up()
        for x in self.fixtures["should_match"]:
            tokens = self.tokenizer.perform(x["query"], [])
            classification = self.classifier.match(tokens)
            validation_results = await self.validator.validate(
                classification, hgvs_dup_del_mode="default",
                endpoint_name=Endpoint.TO_VRS
            )
            is_valid = False
            for vr in validation_results:
                if vr.is_valid:
                    is_valid = True
                    break
            self.assertTrue(is_valid, msg=x)
            self.assertIsNotNone(validation_results, msg=x)

    async def test_not_matches(self):
        """Test that validator matches correctly."""
        self.set_up()
        for x in self.fixtures["should_not_match"]:
            tokens = self.tokenizer.perform(x["query"], [])
            classification = self.classifier.match(tokens)
            validation_results = await self.validator.validate(
                classification, endpoint_name=Endpoint.TO_VRS
            )
            is_valid = False
            for vr in validation_results:
                if vr.is_valid:
                    is_valid = True
                    break
            self.assertFalse(is_valid, msg=x)
