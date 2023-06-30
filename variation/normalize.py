"""Module for Variation Normalization."""
from typing import Optional
from urllib.parse import quote, unquote
from datetime import datetime

from gene.query import QueryHandler as GeneQueryHandler
from cool_seq_tool.data_sources import SeqRepoAccess, UTADatabase

from variation.classifiers.classify import Classify
from variation.to_vrsatile import ToVRSATILE
from variation.hgvs_dup_del_mode import HGVSDupDelMode
from variation.tokenizers.tokenize import Tokenize
from variation.translators.translate import Translate
from variation.utils import no_variation_resp, get_hgvs_dup_del_mode
from variation.validators.validate import Validate
from variation.schemas.app_schemas import Endpoint
from variation.schemas.normalize_response_schema\
    import HGVSDupDelMode as HGVSDupDelModeEnum, NormalizeService, ServiceMeta
from variation.schemas.hgvs_to_copy_number_schema import CopyChange
from variation.version import __version__


class Normalize(ToVRSATILE):
    """The Normalize class used to normalize a given variation."""

    def __init__(self, seqrepo_access: SeqRepoAccess,
                 tokenizer: Tokenize, classifier: Classify, validator: Validate,
                 translator: Translate, hgvs_dup_del_mode: HGVSDupDelMode,
                 gene_normalizer: GeneQueryHandler, uta: UTADatabase) -> None:
        """Initialize Normalize class.

        :param SeqRepoAccess seqrepo_access: Access to SeqRepo
        :param Tokenize tokenizer: Tokenizer class for tokenizing
        :param Classify classifier: Classifier class for classifying tokens
        :param Validate validator: Validator class for validating valid inputs
        :param Translate translator: Translating valid inputs
        :param HGVSDupDelMode hgvs_dup_del_mode: Class for handling
            HGVS dup/del expressions
        :parm GeneQueryHandler gene_normalizer: Client for normalizing gene concepts
        :param UTADatabase uta: Access to db containing alignment data
        """
        super().__init__(seqrepo_access, tokenizer, classifier, validator,
                         translator, hgvs_dup_del_mode, gene_normalizer)
        self.uta = uta

    async def normalize(
        self, q: str,
        hgvs_dup_del_mode: Optional[HGVSDupDelModeEnum] = HGVSDupDelModeEnum.DEFAULT,
        baseline_copies: Optional[int] = None,
        copy_change: Optional[CopyChange] = None,
        untranslatable_returns_text: bool = False
    ) -> NormalizeService:
        """Normalize a given variation.

        :param str q: HGVS, gnomAD VCF or Free Text description on GRCh37 or GRCh38
            assembly
        :param Optional[HGVSDupDelModeEnum] hgvs_dup_del_mode:
            Must be set when querying HGVS dup/del expressions.
            Must be: `default`, `copy_number_count`, `copy_number_change`,
            `repeated_seq_expr`, `literal_seq_expr`. This parameter determines how to
            interpret HGVS dup/del expressions in VRS.
        :param Optional[int] baseline_copies: Baseline copies for HGVS duplications and
            deletions
        :param Optional[CopyChange] copy_change: The copy change
            for HGVS duplications and deletions represented as Copy Number Change
            Variation.
        :param bool untranslatable_returns_text: `True` return VRS Text Object when
            unable to translate or normalize query. `False` return `None` when
            unable to translate or normalize query.
        :return: NormalizeService with variation descriptor and warnings
        """
        label = q.strip()
        _id = f"normalize.variation:{quote(' '.join(label.split()))}"
        vd = None
        warnings = list()
        params = {
            "variation_query": q,
            "variation_descriptor": vd,
            "warnings": warnings,
            "service_meta_": ServiceMeta(
                version=__version__,
                response_datetime=datetime.now()
            )
        }

        tokens = self.tokenizer.perform(unquote(q.strip()), warnings)
        if warnings:
            vd, warnings = no_variation_resp(
                label, _id, warnings, untranslatable_returns_text
            )
            params["variation_descriptor"] = vd
            params["warnings"] = warnings
            return NormalizeService(**params)

        hgvs_dup_del_mode, warning = get_hgvs_dup_del_mode(
            tokens, Endpoint.NORMALIZE, hgvs_dup_del_mode=hgvs_dup_del_mode,
            baseline_copies=baseline_copies
        )
        if warning:
            warnings.append(warning)
            vd, warnings = no_variation_resp(
                label, _id, warnings, untranslatable_returns_text
            )
            params["variation_descriptor"] = vd
            params["warnings"] = warnings
            return NormalizeService(**params)

        classifications = self.classifier.perform(tokens)
        if not classifications:
            warnings.append(f"Unable to find classification for: {q}")
            vd, warnings = no_variation_resp(
                label, _id, warnings, untranslatable_returns_text
            )
            params["variation_descriptor"] = vd
            params["warnings"] = warnings
            return NormalizeService(**params)

        validation_summary = await self.validator.perform(classifications)
        if not validation_summary:
            vd, warnings = no_variation_resp(
                label, _id, validation_summary.warnings, untranslatable_returns_text
            )
            params["variation_descriptor"] = vd
            params["warnings"] = warnings
            return NormalizeService(**params)

        if len(validation_summary.valid_results) > 0:
            translations, warnings = await self.get_translations(
                validation_summary, warnings, endpoint_name=Endpoint.NORMALIZE,
                hgvs_dup_del_mode=hgvs_dup_del_mode,
                baseline_copies=baseline_copies, copy_change=copy_change,
                do_liftover=True
            )
            if translations:
                translation_result = translations[0]
                valid_result = validation_summary.valid_results[0]
                vd, warnings = self.get_variation_descriptor(
                    label, translation_result, valid_result, _id, warnings
                )

                if not vd:
                    vd, warnings = no_variation_resp(
                        label, _id, warnings, untranslatable_returns_text
                    )
            else:
                vd, warnings = no_variation_resp(
                    label, _id, warnings, untranslatable_returns_text
                )
        else:
            vd, warnings = no_variation_resp(
                label, _id, warnings, untranslatable_returns_text
            )

        return NormalizeService(
            variation_query=q,
            variation_descriptor=vd,
            warnings=warnings,
            service_meta_=ServiceMeta(
                version=__version__,
                response_datetime=datetime.now()
            )
        )
