"""GA4GH service-info response models."""

from pydantic import BaseModel, ConfigDict, Field


class ServiceType(BaseModel):
    """Identify the API exposed by this service."""

    group: str
    artifact: str
    version: str


class Organization(BaseModel):
    """Identify the organization providing the service."""

    name: str
    url: str


class SpecMetadata(BaseModel):
    """Describe the standards used by the service."""

    vrs_version: str = Field(alias="vrsVersion")


class DataVersions(BaseModel):
    """Identify the reference data releases used by the service."""

    uta_schema: str = Field(alias="utaSchema")
    seqrepo_version: str = Field(alias="seqrepoVersion")


class ServiceInfo(BaseModel):
    """Describe this Variation Normalizer instance."""

    model_config = ConfigDict(populate_by_name=True)

    id: str
    name: str
    type: ServiceType
    organization: Organization
    version: str
    description: str
    documentation_url: str = Field(alias="documentationUrl")
    spec_metadata: SpecMetadata = Field(alias="specMetadata")
    data_versions: DataVersions = Field(alias="dataVersions")
