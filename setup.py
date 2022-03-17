"""Module for package and distribution."""
import setuptools

exec(open("variation/version.py").read())
setuptools.setup(version=__version__)  # noqa: F821
