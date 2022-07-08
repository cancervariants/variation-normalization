"""A module for testing translator classes."""
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


class TranslatorBase:
    """The translator base class."""

    @classmethod
    def setUpClass(cls):
        """Set up the test cases."""
        with open(f"{PROJECT_ROOT}/tests/fixtures/translators.yml") as stream:
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
                    variation = self.translator.translate(vr)
                    if variation["type"] == "Allele":
                        variation["location"] = variation["location"]
                        if "id" in variation["location"].keys():
                            del variation["location"]["id"]
                    elif variation["type"] == "AbsoluteCopyNumber":
                        variation["subject"] = variation["subject"]
                        if "id" in variation["subject"]["location"].keys():
                            del variation["subject"]["location"]["id"]

                    if variation not in found:
                        found.append(variation)
                        num_valid += 1
                        self.assertIn(variation, x["variations"],
                                      msg=x["query"])
            self.assertEqual(len(x["variations"]), num_valid, msg=x["query"])
