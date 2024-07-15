"""The Variation Normalization package."""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("variation-normalizer")
except PackageNotFoundError:
    __version__ = "unknown"
finally:
    del version, PackageNotFoundError
