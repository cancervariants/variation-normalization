"""Module for gnomad vcf to protein response schema"""

from ga4gh.core.models import MappableConcept

from variation.schemas.normalize_response_schema import NormalizeService


class GnomadVcfToProteinService(NormalizeService):
    """Define response for gnomad vcf to protein service"""

    gene_context: MappableConcept | None = None
