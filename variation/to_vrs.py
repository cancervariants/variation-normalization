"""Module for to_vrs endpoint."""
from typing import Tuple, Optional, List, Union, Dict
from urllib.parse import unquote
from datetime import datetime

from ga4gh.vrsatile.pydantic.vrs_models import Allele, Haplotype, CopyNumberCount,\
    VariationSet, Text
from ga4gh.vrs import models
from ga4gh.core import ga4gh_identify
from cool_seq_tool.schemas import ResidueMode
from cool_seq_tool.data_sources import SeqRepoAccess
from gene.query import QueryHandler as GeneQueryHandler

from variation.schemas.normalize_response_schema\
    import HGVSDupDelMode as HGVSDupDelModeEnum, ServiceMeta
from variation.hgvs_dup_del_mode import HGVSDupDelMode
from variation.schemas.app_schemas import Endpoint
from variation.schemas.hgvs_to_copy_number_schema import CopyChange
from variation.schemas.to_vrs_response_schema import ToVRSService
from variation.schemas.validation_response_schema import ValidationSummary
from variation.classifiers import Classify
from variation.tokenizers import Tokenize
from variation.validators import Validate
from variation.translators import Translate
from variation.vrs_representation import VRSRepresentation
from variation.version import __version__


class ToVRS(VRSRepresentation):
    """The class for translating variation strings to VRS representations."""

    def __init__(self, seqrepo_access: SeqRepoAccess,
                 tokenizer: Tokenize, classifier: Classify, validator: Validate,
                 translator: Translate, hgvs_dup_del_mode: HGVSDupDelMode,
                 gene_normalizer: GeneQueryHandler) -> None:
        """Initialize the ToVrsAndVrsatile class.

        :param SeqRepoAccess seqrepo_access: Access to SeqRepo
        :param Tokenize tokenizer: Tokenizer class for tokenizing
        :param Classify classifier: Classifier class for classifying tokens
        :param Validate validator: Validator class for validating valid inputs
        :param Translate translator: Translating valid inputs
        :param HGVSDupDelMode hgvs_dup_del_mode: Class for handling
            HGVS dup/del expressions
        :param GeneQueryHandler gene_normalizer: Client for normalizing gene concepts
        """
        super().__init__(seqrepo_access)
        self.tokenizer = tokenizer
        self.classifier = classifier
        self.validator = validator
        self.translator = translator
        self.hgvs_dup_del_mode = hgvs_dup_del_mode
        self.gene_normalizer = gene_normalizer

    async def get_translations(
        self,
        validation_summary: ValidationSummary,
        warnings: List,
        endpoint_name: Optional[Endpoint] = None,
        hgvs_dup_del_mode: HGVSDupDelModeEnum = HGVSDupDelModeEnum.DEFAULT,
        baseline_copies: Optional[int] = None,
        copy_change: Optional[CopyChange] = None,
        do_liftover: bool = False
    ) -> Tuple[
        Optional[
            Union[
                List[Allele],
                List[CopyNumberCount],
                List[Text], List[Haplotype],
                List[VariationSet]
            ]
        ],
        Optional[List[str]]
    ]:
        translations = []
        if not validation_summary.valid_results:
            warnings.append("No valid results found")
        else:
            for valid_result in validation_summary.valid_results:
                result = await self.translator.perform(
                    valid_result, warnings, endpoint_name=endpoint_name,
                    hgvs_dup_del_mode=hgvs_dup_del_mode,
                    baseline_copies=baseline_copies, copy_change=copy_change,
                    do_liftover=do_liftover
                )
                if result and result not in translations:
                    translations.append(result)

        if not translations and not warnings:
            warnings.append("Unable to translate variation")

        return translations, warnings

    def get_ref_allele_seq(self, location: Dict, ac: str) -> Optional[str]:
        """Return ref allele seq for transcript.

        :param Dict location: VRS Location object
        :param str identifier: Identifier for allele
        :return: Ref seq allele
        """
        start = None
        end = None
        interval = location["interval"]
        ival_type = interval["type"]

        if ival_type == "SequenceInterval":
            if interval["start"]["type"] == "Number":
                start = interval["start"]["value"]
                end = interval["end"]["value"]

                if start == end:
                    return None

        if start is None and end is None:
            return None

        ref, _ = self.seqrepo_access.get_reference_sequence(
            ac, start, end, residue_mode=ResidueMode.INTER_RESIDUE
        )

        return ref

    async def to_vrs(
        self, q: str, untranslatable_returns_text: bool = False
    ) -> ToVRSService:
        """Return a VRS-like representation of all validated variations for a query.

        :param str q: The variation to translate (HGVS, gnomAD VCF, or free text) on
            GRCh37 or GRCh38 assembly
        :param bool untranslatable_returns_text: `True` return VRS Text Object when
            unable to translate or normalize query. `False` returns empty list when
            unable to translate or normalize query.
        :return: ToVRSService containing VRS variations and warnings
        """
        warnings = []
        variations = []
        params = {
            "search_term": q,
            "variations": variations,
            "service_meta_": ServiceMeta(
                version=__version__,
                response_datetime=datetime.now()
            ),
            "warnings": warnings
        }

        tokens = self.tokenizer.perform(unquote(q.strip()), warnings)
        if warnings:
            params["warnings"] = warnings
            return ToVRSService(**params)

        classification = self.classifier.perform(tokens)
        if not classification:
            params["warnings"] = [f"Unable to find classification for: {q}"]
            return ToVRSService(**params)

        validations = await self.validator.perform(classification, warnings)

        if validations:
            translations, warnings = await self.get_translations(
                validations, warnings, endpoint_name=Endpoint.TO_VRS,
                hgvs_dup_del_mode=HGVSDupDelModeEnum.DEFAULT, do_liftover=False
            )
        else:
            translations = []

        if not translations:
            if untranslatable_returns_text and q and q.strip():
                text = models.Text(definition=q, type="Text")
                text._id = ga4gh_identify(text)
                variations = [Text(**text.as_dict())]
            else:
                variations = []
        else:
            variations = [tr.vrs_variation for tr in translations]

        params["warnings"] = warnings
        params["variations"] = variations
        return ToVRSService(**params)
