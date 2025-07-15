# Variation Normalization

[![image](https://img.shields.io/pypi/v/variation-normalizer.svg)](https://pypi.python.org/pypi/variation-normalizer) [![image](https://img.shields.io/pypi/l/variation-normalizer.svg)](https://pypi.python.org/pypi/variation-normalizer) [![image](https://img.shields.io/pypi/pyversions/variation-normalizer.svg)](https://pypi.python.org/pypi/variation-normalizer) [![Actions status](https://github.com/cancervariants/variation-normalization/actions/workflows/checks.yaml/badge.svg)](https://github.com/cancervariants/variation-normalization/actions/checks.yaml)[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5894937.svg)](https://doi.org/10.5281/zenodo.5894937)

<!-- description -->
The Variation Normalizer parses and translates free-text descriptions of genomic variations into computable objects conforming to the [Variation Representation Specification (VRS)](https://vrs.ga4gh.org/en/latest), enabling consistent and accurate variant harmonization across a diversity of genomic knowledge resources.
<!-- /description -->

---

[Live OpenAPI endpoint](https://normalize.cancervariants.org/variation)

---

## Installation

Install from [PyPI](https://pypi.org/project/variation-normalizer):

```shell
python3 -m pip install variation-normalizer
```

---

| variation-normalization branch | variation-normalizer version | gene-normalizer version | VRS version |
| ---- | --- | ---- | --- |
| [main](https://github.com/cancervariants/variation-normalization/tree/main) | >=0.14.Z | >=0.9.Z | [2.0](https://github.com/ga4gh/vrs/tree/2.0) |

## About

Variation Normalization works by using four main steps: tokenization, classification, validation, and translation. During tokenization, we split strings on whitespace and parse to determine the type of token. During classification, we specify the order of tokens a classification can have. We then do validation checks such as ensuring references for a nucleotide or amino acid matches the expected value and validating a position exists on the given transcript. During translation, we return a VRS Allele object.

Variation Normalization is limited to the following types of variants:

* HGVS expressions and text representations (ex: `BRAF V600E`):
  * **protein (p.)**: substitution, deletion, insertion, deletion-insertion
  * **coding DNA (c.)**: substitution, deletion, insertion, deletion-insertion
  * **genomic (g.)**: substitution, deletion, ambiguous deletion, insertion, deletion-insertion, duplication
* gnomAD-style VCF (chr-pos-ref-alt, ex: `7-140753336-A-T`)
  * **genomic (g.)**: substitution, deletion, insertion

Variation Normalizer accepts input from GRCh37 or GRCh8 assemblies.

We are working towards adding more types of variations, coordinates, and representations.

### VRS Versioning

The variation-normalization repo depends on VRS models, and therefore each variation-normalizer package on PyPI uses a particular version of VRS. The correspondences between packages may be summarized as:

| variation-normalization branch | variation-normalizer version | gene-normalizer version | VRS version |
| ---- | --- | ---- | --- |
| [main](https://github.com/cancervariants/variation-normalization/tree/main) | >=0.14.Z | >=0.9.Z | [2.0](https://github.com/ga4gh/vrs/tree/2.0) |

### Previous VRS Versioning

The correspondences between the packages that are **no longer maintained** may be summarized as:

| variation-normalization branch | variation-normalizer version | gene-normalizer version | VRS version |
| ---- | --- | ---- | --- |
| [vrs-1.3](https://github.com/cancervariants/variation-normalization/tree/vrs-1.3) | 0.6.Z | 0.1.Z | [1.3](https://github.com/ga4gh/vrs/tree/1.3) |

### Available Endpoints

#### `/to_vrs`

Returns a list of validated VRS [Variations](https://vrs.ga4gh.org/en/stable/terms_and_model.html#variation).

#### `/normalize`

Returns a VRS Variation aligned to the prioritized transcript. The Variation Normalizer relies on [**C**ommon **O**perations **O**n **L**ots-of **Seq**uences Tool (cool-seq-tool)](https://github.com/GenomicMedLab/cool-seq-tool) for retrieving the prioritized transcript data. More information on the transcript selection algorithm can be found [here](https://github.com/GenomicMedLab/cool-seq-tool/blob/main/docs/TranscriptSelectionPriority.md).

If a genomic variation query _is_ given a gene (E.g. `BRAF g.140753336A>T`), the associated cDNA representation will be returned. This is because the gene provides additional strand context. If a genomic variation query is _not_ given a gene, the GRCh38 representation will be returned.

## Development

Clone the repo:

```shell
git clone https://github.com/cancervariants/variation-normalization.git
cd variation-normalization
```

For a development install, we recommend using Pipenv. See the
[pipenv docs](https://pipenv-fork.readthedocs.io/en/latest/#install-pipenv-today)
for direction on installing pipenv in your compute environment.

Once installed, from the project root dir, just run:

```shell
pipenv shell
pipenv update && pipenv install --dev
```

### Required resources

Variation Normalization relies on some local data caches which you will need to set up.
We provide instructions on how to setup your development environment using Docker.

* [SeqRepo](https://github.com/biocommons/biocommons.seqrepo): You must setup SeqRepo
locally following [these steps](#seqrepo).
* [Gene Normalizer](https://github.com/cancervariants/gene-normalization): The Variation
Normalizer uses Gene Normalizer to get normalized gene concept information.
* [Universal Transcript Archive (UTA)](https://github.com/biocommons/uta): The Variation
Normalizer uses [**C**ommon **O**perations **O**n **L**ots-of **Seq**uences Tool (cool-seq-tool)](https://github.com/GenomicMedLab/cool-seq-tool) which uses UTA as the underlying PostgreSQL database.

#### SeqRepo

Variation Normalization relies on [seqrepo](https://github.com/biocommons/biocommons.seqrepo), which you must download yourself.

Variation Normalizer uses seqrepo to retrieve sequences at given positions on a transcript.

From the _root_ directory:

```shell
pip install seqrepo
sudo mkdir /usr/local/share/seqrepo
sudo chown $USER /usr/local/share/seqrepo
seqrepo pull -i 2024-12-20/  # Replace with latest version using `seqrepo list-remote-instances` if outdated
```

If you get an error similar to the one below:

```shell
PermissionError: [Error 13] Permission denied: '/usr/local/share/seqrepo/2024-12-20/._fkuefgd' -> '/usr/local/share/seqrepo/2024-12-20/'
```

You will want to do the following:\
(_Might not be .\_fkuefgd, so replace with your error message path_)

```shell
sudo mv /usr/local/share/seqrepo/2024-12-20._fkuefgd /usr/local/share/seqrepo/2024-12-20
exit
```

Use the `SEQREPO_ROOT_DIR` environment variable to set the path of an already existing SeqRepo directory. The default is `/usr/local/share/seqrepo/latest`.

## Docker Installation (Preferred)

We recommend installing the Variation Normalizer using Docker.

### Requirements

* [Docker](https://docs.docker.com/get-started/get-docker/)

### Build, (re)create, and start containers

```shell
docker volume create --name=uta_vol
docker compose up
```

> [!IMPORTANT]
> This assumes you have a local [SeqRepo](https://github.com/biocommons/biocommons.seqrepo)
installed at `/usr/local/share/seqrepo/2024-12-20`. If you have it installed elsewhere,
please update the `SEQREPO_ROOT_DIR` environment variable in
[compose.yaml](./compose.yaml).\
> If you're using Docker Desktop, you'll want to go to Settings -> Resources -> File sharing
and add `/usr/local/share/seqrepo` under the `Virtual file shares` section. Otherwise,
you will get the following error:
`OSError: Unable to open SeqRepo directory /usr/local/share/seqrepo/2024-12-20`.

> [!TIP]
> If you want a clean slate, run `docker compose down -v` to remove containers and
> volumes, then `docker compose up --build` to rebuild and start fresh containers.

Point your browser to <http://localhost:8001/variation/>.

### Code QC

Code style is managed by [Ruff](https://docs.astral.sh/ruff/) and checked prior to commit.

To perform formatting and check style:

```shell
python3 -m ruff format . && python3 -m ruff check --fix .
```

We use [pre-commit](https://pre-commit.com/#usage) to run conformance tests.

This ensures:

* Style correctness
* No large files
* AWS credentials are present
* Private key is present

Pre-commit *must* be installed before your first commit. Use the following command:

```commandline
pre-commit install
```

### Testing

From the _root_ directory of the repository:

```shell
pytest tests/
```
