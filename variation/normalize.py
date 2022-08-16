"""Module for Variation Normalization."""
from typing import Optional
from urllib.parse import quote
from datetime import datetime

from gene.query import QueryHandler as GeneQueryHandler
from uta_tools.data_sources import SeqRepoAccess, UTADatabase
from ga4gh.vrs.dataproxy import SeqRepoDataProxy

from variation.classifiers.classify import Classify
from variation.to_vrsatile import ToVRSATILE
from variation.hgvs_dup_del_mode import HGVSDupDelMode
from variation.tokenizers.tokenize import Tokenize
from variation.translators.translate import Translate
from variation.utils import get_mane_valid_result, no_variation_entered, \
    no_variation_resp
from variation.validators.validate import Validate
from variation.schemas.app_schemas import Endpoint
from variation.schemas.normalize_response_schema\
    import HGVSDupDelMode as HGVSDupDelModeEnum, NormalizeService, ServiceMeta
from variation.schemas.hgvs_to_copy_number_schema import RelativeCopyClass
from variation.version import __version__


class Normalize(ToVRSATILE):
    """The Normalize class used to normalize a given variation."""

    def __init__(self, seqrepo_access: SeqRepoAccess, dp: SeqRepoDataProxy,
                 tokenizer: Tokenize, classifier: Classify, validator: Validate,
                 translator: Translate, hgvs_dup_del_mode: HGVSDupDelMode,
                 gene_normalizer: GeneQueryHandler, uta: UTADatabase) -> None:
        """Initialize Normalize class.

        :param SeqRepoAccess seqrepo_access: Access to SeqRepo via UTA Tools
        :param SeqRepoDataProxy dp: Access to SeqRepo via VRS Python
        :param Tokenize tokenizer: Tokenizer class for tokenizing
        :param Classify classifier: Classifier class for classifying tokens
        :param Validate validator: Validator class for validating valid inputs
        :param Translate translator: Translating valid inputs
        :param HGVSDupDelMode hgvs_dup_del_mode: Class for handling
            HGVS dup/del expressions
        :parm GeneQueryHandler gene_normalizer: Client for normalizing gene concepts
        :param UTADatabase uta: Access to db containing alignment data
        """
        super().__init__(seqrepo_access, dp, tokenizer, classifier, validator,
                         translator, hgvs_dup_del_mode, gene_normalizer)
        self.uta = uta

    async def normalize(
        self, q: str,
        hgvs_dup_del_mode: Optional[HGVSDupDelModeEnum] = HGVSDupDelModeEnum.DEFAULT,
        baseline_copies: Optional[int] = None,
        relative_copy_class: Optional[RelativeCopyClass] = None,
        untranslatable_returns_text: bool = False
    ) -> NormalizeService:
        """Normalize a given variation.

        :param str q: The variation to normalize
        :param Optional[HGVSDupDelModeEnum] hgvs_dup_del_mode:
            Must be set when querying HGVS dup/del expressions.
            Must be: `default`, `absolute_cnv`, `relative_cnv`, `repeated_seq_expr`,
            `literal_seq_expr`. This parameter determines how to interpret HGVS dup/del
            expressions in VRS.
        :param Optional[int] baseline_copies: Baseline copies for HGVS duplications and
            deletions
        :param Optional[RelativeCopyClass] relative_copy_class: The relative copy class
            for HGVS duplications and deletions represented as Relative Copy Number
            Variation.
        :param bool untranslatable_returns_text: `True` return VRS Text Object when
            unable to translate or normalize query. `False` return `None` when
            unable to translate or normalize query.
        :return: NormalizeService with variation descriptor and warnings
        """
        vd = None
        warnings = list()
        if not q:
            vd, warnings = no_variation_entered()
        else:
            validations, warnings = await self.get_validations(
                q, endpoint_name=Endpoint.NORMALIZE,
                hgvs_dup_del_mode=hgvs_dup_del_mode,
                baseline_copies=baseline_copies,
                relative_copy_class=relative_copy_class)

            if validations:
                label = q.strip()
                _id = f"normalize.variation:{quote(' '.join(label.split()))}"
                if len(validations.valid_results) > 0:
                    valid_result = get_mane_valid_result(q, validations, warnings)
                    vd, warnings = self.get_variation_descriptor(
                        label, valid_result.variation, valid_result, _id, warnings)
                else:
                    if not label:
                        vd, warnings = no_variation_entered()
                    else:
                        vd, warnings = no_variation_resp(label, _id, warnings,
                                                         untranslatable_returns_text)

        return NormalizeService(
            variation_query=q,
            variation_descriptor=vd,
            warnings=warnings,
            service_meta_=ServiceMeta(
                version=__version__,
                response_datetime=datetime.now()
            )
        )
