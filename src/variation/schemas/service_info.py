"""GA4GH service-info response models."""

from typing import Literal

from ga4gh.vrs import VRS_VERSION
from pydantic import BaseModel, ConfigDict, Field

from variation import __version__


class ServiceType(BaseModel):
    """Identify the API exposed by this service."""

    group: Literal["org.cancervariants"] = "org.cancervariants"
    artifact: Literal["Variation Normalizer API"] = "Variation Normalizer API"
    version: str = __version__


class ServiceOrganization(BaseModel):
    """Identify the organization providing the service."""

    name: Literal["Variant Interpretation for Cancer Consortium"] = (
        "Variant Interpretation for Cancer Consortium"
    )
    url: Literal["https://cancervariants.org"] = "https://cancervariants.org"


class SpecMetadata(BaseModel):
    """Describe the standards used by the service."""

    vrs_version: str = Field(default=VRS_VERSION, alias="vrsVersion")


class DataVersions(BaseModel):
    """Identify the reference data releases used by the service."""

    uta_schema: Literal["uta_20241220"] = Field(
        default="uta_20241220", alias="utaSchema"
    )
    seqrepo_version: Literal["2024-12-20"] = Field(
        default="2024-12-20", alias="seqrepoVersion"
    )


class ServiceInfo(BaseModel):
    """Describe this Variation Normalizer instance."""

    model_config = ConfigDict(populate_by_name=True)

    id: Literal["org.cancervariants.variation_normalizer"] = (
        "org.cancervariants.variation_normalizer"
    )
    name: Literal["Variation Normalizer"] = "Variation Normalizer"
    type: ServiceType
    organization: ServiceOrganization
    version: str = __version__
    description: Literal["Normalize variation descriptions to GA4GH VRS objects."] = (
        "Normalize variation descriptions to GA4GH VRS objects."
    )
    documentation_url: Literal[
        "https://github.com/cancervariants/variation-normalization"
    ] = Field(
        default="https://github.com/cancervariants/variation-normalization",
        alias="documentationUrl",
    )
    spec_metadata: SpecMetadata = Field(
        default_factory=SpecMetadata, alias="specMetadata"
    )
    data_versions: DataVersions = Field(
        default_factory=DataVersions, alias="dataVersions"
    )
