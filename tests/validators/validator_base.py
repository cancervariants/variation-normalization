"""A module for testing validator classes."""
import yaml
from biocommons.seqrepo import SeqRepo
from ga4gh.vrs.extras.translator import Translator
from gene.database.dynamodb import DynamoDbDatabase
from gene.query import QueryHandler as GeneQueryHandler
from cool_seq_tool import SEQREPO_ROOT_DIR
from cool_seq_tool.data_sources import TranscriptMappings, SeqRepoAccess, \
    MANETranscriptMappings, UTADatabase, MANETranscript

from tests import PROJECT_ROOT
from variation.schemas.app_schemas import Endpoint
from variation.vrs_representation import VRSRepresentation
from variation.tokenizers import Tokenize, GeneSymbol


class ValidatorBase:
    """The validator base class."""

    @classmethod
    def setUpClass(cls):
        """Set up the test cases."""
        with open(f"{PROJECT_ROOT}/tests/fixtures/validators.yml") as stream:
            cls.all_fixtures = yaml.safe_load(stream)
        gene_normalizer = GeneQueryHandler(DynamoDbDatabase())
        gene_symbol = GeneSymbol(gene_normalizer)
        cls.tokenizer = Tokenize(gene_symbol)

        sr = SeqRepo(root_dir=SEQREPO_ROOT_DIR)
        seqrepo_access = SeqRepoAccess(sr)
        transcript_mappings = TranscriptMappings()
        uta = UTADatabase()
        tlr = Translator(data_proxy=seqrepo_access)
        mane_transcript = MANETranscript(
            seqrepo_access, transcript_mappings, MANETranscriptMappings(), uta,
            gene_normalizer)
        vrs = VRSRepresentation(seqrepo_access)

        cls.aa_params = [
            seqrepo_access, transcript_mappings, gene_symbol,
            mane_transcript, uta, tlr, gene_normalizer, vrs
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
