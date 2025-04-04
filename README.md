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

Variation Normalization relies on some local data caches which you will need to set up. It uses pipenv to manage its environment, which you will also need to install.

#### Gene Normalizer

Variation Normalization relies on data from [Gene Normalization](https://github.com/cancervariants/gene-normalization). You must load all sources _and_ merged concepts.

You must also have Gene Normalization's DynamoDB running in a separate terminal for the application to work.

For more information about the gene-normalizer and how to load the database, visit the [README](https://github.com/cancervariants/gene-normalization/blob/main/README.md).

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

#### UTA

Variation Normalizer also uses [**C**ommon **O**perations **O**n **L**ots-of **Seq**uences Tool (cool-seq-tool)](https://github.com/GenomicMedLab/cool-seq-tool) which uses [UTA](https://github.com/biocommons/uta) as the underlying PostgreSQL database.

We provide two options for installing UTA:

1. [Using Docker](#installing-uta-via-docker): This is the preferred way
2. [Locally](#installing-uta-locally)

##### Installing UTA via Docker

For this, you will need to install Docker. We recommend using
[Docker Desktop](https://docs.docker.com/desktop/).

Once Docker is running, from the root of the directory, run the following:

```shell
docker volume create --name=uta_vol
docker compose up
```

This should start the following container:

* [uta](https://github.com/biocommons/uta): a database of transcripts and alignments (localhost:5432)

Check that the container is running:

```shell
$ docker ps
CONTAINER ID        IMAGE                                    //  NAMES
a40576b8cf1f        biocommons/uta:uta_20241220              //  variation-normalization-uta-1
```

Depending on your network and host, the _first_ run is likely to take 5-15
minutes in order to download and install data. Subsequent startups should be
nearly instantaneous.

You can test UTA and seqrepo installations like so:

```shell
$ psql -XAt postgres://anonymous@localhost/uta -c 'select count(*) from uta_20241220.transcript'
329090
```

##### Installing UTA Locally

_The following commands will likely need modification appropriate for the installation environment._

1. Install [PostgreSQL](https://www.postgresql.org/)
2. Create user and database.

    ```shell
    createuser -U postgres uta_admin
    createuser -U postgres anonymous
    createdb -U postgres -O uta_admin uta
    ```

3. To install locally:

```shell
export UTA_VERSION=uta_20241220.pgd.gz
curl -O http://dl.biocommons.org/uta/$UTA_VERSION
gzip -cdq ${UTA_VERSION} | grep -v "^REFRESH MATERIALIZED VIEW" | psql -h localhost -U uta_admin --echo-errors --single-transaction -v ON_ERROR_STOP=1 -d uta -p 5432
```

If you have trouble installing UTA, you can visit [these two READMEs](https://github.com/ga4gh/vrs-python/tree/main/docs/setup_help).

##### Connecting to the UTA database

To connect to the UTA database, you can use the default url (`postgresql://uta_admin@localhost:5432/uta/uta_20241220`). If you do not wish to use the default, you must set the environment variable `UTA_DB_URL` which has the format of `driver://user:pass@host:port/database/schema`.

## Starting the Variation Normalization Service Locally

`gene-normalizer`s dynamodb and the `uta` database must be running.

To start the service, run the following:

```shell
uvicorn variation.main:app --reload
```

Next, view the OpenAPI docs on your local machine:
<http://127.0.0.1:8000/variation>

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
