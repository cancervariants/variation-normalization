"""A module for testing translator classes."""
import yaml
from biocommons.seqrepo import SeqRepo
from gene.database.dynamodb import DynamoDbDatabase
from gene.query import QueryHandler as GeneQueryHandler
from cool_seq_tool.paths import SEQREPO_ROOT_DIR
from cool_seq_tool.data_sources import TranscriptMappings, SeqRepoAccess, \
    MANETranscriptMappings, UTADatabase, MANETranscript

from tests import PROJECT_ROOT
from variation.schemas.app_schemas import Endpoint
from variation.vrs_representation import VRSRepresentation
from variation.tokenizers import Tokenize, GeneSymbol
from variation.hgvs_dup_del_mode import HGVSDupDelMode


class TranslatorBase:
    """The translator base class."""

    @classmethod
    def setUpClass(cls):
        """Set up the test cases."""
        with open(f"{PROJECT_ROOT}/tests/fixtures/translators.yml") as stream:
            cls.all_fixtures = yaml.safe_load(stream)
        gene_normalizer = GeneQueryHandler(DynamoDbDatabase())
        gene_symbol = GeneSymbol(gene_normalizer)
        cls.tokenizer = Tokenize(gene_symbol)

        sr = SeqRepo(root_dir=SEQREPO_ROOT_DIR)
        seqrepo_access = SeqRepoAccess(sr)
        transcript_mappings = TranscriptMappings()
        uta = UTADatabase()
        mane_transcript = MANETranscript(
            seqrepo_access, transcript_mappings, MANETranscriptMappings(), uta,
            gene_normalizer)
        vrs = VRSRepresentation(seqrepo_access)
        hgvs_dup_del_mode = HGVSDupDelMode(seqrepo_access)

        cls.aa_params = [
            seqrepo_access, mane_transcript, uta, vrs, hgvs_dup_del_mode
        ]
        cls.params = cls.aa_params[:-1]

    def classifier_instance(self):
        """Check that the classifier_instance method is implemented."""
        raise NotImplementedError()

    def validator_instance(self):
        """Check that the validator_instance method is implemented."""
        raise NotImplementedError()

    def translator_instance(self):
        """Check that the translator_instance method is implemented."""
        raise NotImplementedError()

    def fixture_name(self):
        """Check that the fixture_name method is implemented."""
        raise NotImplementedError()

    def set_up(self):
        """Initialize fixtures, classifier, validator + translator"""
        self.fixtures = self.all_fixtures.get(
            self.fixture_name(),
            {"tests": []}
        )
        self.classifier = self.classifier_instance()
        self.validator = self.validator_instance()
        self.translator = self.translator_instance()

    async def test_translator(self):
        """Test that translator matches correctly."""
        self.set_up()
        for x in self.fixtures["tests"]:
            tokens = self.tokenizer.perform(x["query"], [])
            classification = self.classifier.match(tokens)
            validation_results = await self.validator.validate(
                classification, endpoint_name=Endpoint.TO_VRS
            )
            num_valid = 0
            found = list()
            for vr in validation_results:
                if vr.is_valid:
                    variation = await self.translator.translate(vr)
                    if variation["type"] == "Allele":
                        variation["location"] = variation["location"]
                        if "id" in variation["location"].keys():
                            del variation["location"]["id"]
                    elif variation["type"] == "CopyNumberCount":
                        variation["subject"] = variation["subject"]
                        if "id" in variation["subject"]["location"].keys():
                            del variation["subject"]["location"]["id"]

                    if variation not in found:
                        found.append(variation)
                        num_valid += 1
                        self.assertIn(variation, x["variations"],
                                      msg=x["query"])
            self.assertEqual(len(x["variations"]), num_valid, msg=x["query"])
