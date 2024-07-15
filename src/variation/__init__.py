"""The Variation Normalization package."""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("{{ cookiecutter.project_slug }}")
except PackageNotFoundError:
    __version__ = "unknown"
finally:
    del version, PackageNotFoundError
