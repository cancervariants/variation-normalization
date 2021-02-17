"""Module for package setup."""
from setuptools import setup
from variant.__version__ import __version__, __authors__, __author_email__, \
    __description__, __url__

setup(
    name='variant-normalization',
    version=__version__,
    packages=['variant'],
    url=__url__,
    license='MIT',
    author=', '.join(__authors__),
    author_email=__author_email__,
    description=__description__,
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3'
    ],
    install_requires=[],  # TODO: ADD
    extras_require={
        'test': [
            'pytest'
        ]
    },
    python_requires='>=3.7'
)
