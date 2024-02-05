Variation Normalizer |version|
==============================

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.5894937.svg
   :alt: DOI
   :target: https://doi.org/10.5281/zenodo.5894937

.. image:: https://img.shields.io/pypi/v/variation-normalizer.svg
   :alt: PyPI version
   :target: https://pypi.python.org/pypi/variation-normalizer

.. image:: https://img.shields.io/pypi/l/variation-normalizer.svg
   :alt: License
   :target: https://github.com/cancervariants/variation-normalization/blob/main/LICENSE

.. image:: https://img.shields.io/pypi/pyversions/variation-normalizer?color=gr
   :alt: PyPI - supported Python versions

.. image:: https://github.com/cancervariants/variation-normalization/actions/workflows/checks.yaml/badge.svg
   :alt: tests status
   :target: https://github.com/cancervariants/variation-normalization/actions/workflows/checks.yaml

Variation Normalization works by using four main steps: tokenization, classification, validation, and translation. During tokenization, we split strings on whitespace and parse to determine the type of token. During classification, we specify the order of tokens a classification can have. We then do validation checks such as ensuring references for a nucleotide or amino acid matches the expected value and validating a position exists on the given transcript. During translation, we return a VRS Allele object.

Variation Normalization is limited to the following types of variants:

* HGVS expressions and text representations (ex: BRAF V600E):

  * protein (p.): substitution, deletion, insertion, deletion-insertion

  * coding DNA (c.): substitution, deletion, insertion, deletion-insertion

  * genomic (g.): substitution, deletion, ambiguous deletion, insertion, deletion-insertion, duplication

* gnomAD-style VCF (chr-pos-ref-alt, ex: 7-140753336-A-T)

  * genomic (g.): substitution, deletion, insertion

Variation Normalizer accepts input from GRCh37 or GRCh8 assemblies.

We are working towards adding more types of variations, coordinates, and representations.

.. toctree::
   :hidden:
   :maxdepth: 2

    Installation<install>
    API Reference<reference/index>
    Changelog<changelog>
    Contributing<contributing>
    License<license>
