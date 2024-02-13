"""Module for gnomad vcf to protein response schema"""
from typing import Optional

from ga4gh.core import core_models
from pydantic import StrictStr

from variation.schemas.normalize_response_schema import NormalizeService


class GnomadVcfToProteinService(NormalizeService):
    """Define response for gnomad vcf to protein service"""

    gene_context: Optional[core_models.Gene] = None
    vrs_ref_allele_seq: Optional[StrictStr] = None
