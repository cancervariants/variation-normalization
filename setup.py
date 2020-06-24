from setuptools import setup
from varlexapp.__version__ import __version__, __authors__, __author_email__, \
    __description__, __url__

setup(
    name='varlex',
    version=__version__,
    packages=['varlexapp'],
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
    install_requires=[
        'apispec~=2.0.2',
        'apispec-oneofschema~=2.1.0',
        'apispec-webframeworks~=0.4.0',
        'Flask~=1.1.1',
        'Flask-Cors~=3.0.8'
    ],
    extras_require={
        'test': [
            'pytest'
        ]
    },
    python_requires='>=3.7'
)
