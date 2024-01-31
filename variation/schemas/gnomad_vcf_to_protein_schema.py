"""Module for gnomad vcf to protein response schema"""
from typing import Optional

from ga4gh.core import core_models

from variation.schemas.normalize_response_schema import NormalizeService


class GnomadVcfToProteinService(NormalizeService):
    """Define response for gnomad vcf to protein service"""

    gene_context: Optional[core_models.Gene] = None
