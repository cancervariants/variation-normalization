[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5894937.svg)](https://doi.org/10.5281/zenodo.5894937)

# Variation Normalization

Services and guidelines for normalizing variation terms into [VRS](https://vrs.ga4gh.org/en/latest) and [VRSATILE](https://vrsatile.readthedocs.io/en/latest/) compatible representations.

Public OpenAPI endpoint: https://normalize.cancervariants.org/variation

Installing with pip:

```commandline
pip install variation-normalizer
```

## About
Variation Normalization works by using four main steps: tokenization, classification, validation, and translation. During tokenization, we split strings on whitespace and parse to determine the type of token. During classification, we specify the order of tokens a classification can have. We then do validation checks such as ensuring references for a nucleotide or amino acid matches the expected value and validating a position exists on the given transcript. During translation, we return a VRS Allele object.

Variation Normalization is limited to the following types of variants:
* HGVS expressions and text representations (ex: `BRAF V600E`):
  * **protein (p.)**: substitution, deletion, insertion, deletion-insertion
  * **coding DNA (c.)**: substitution, deletion, insertion, deletion-insertion
  * **genomic (g.)**: substitution, deletion, ambiguous deletion, insertion, deletion-insertion, duplication
* gnomAD-style VCF (chr-pos-ref-alt, ex: `7-140753336-A-T`)
  * **genomic (g.)**: substitution, deletion, insertion

We are working towards adding more types of variations, coordinates, and representations.

### Endpoints
#### /to_vrs
The `/to_vrs` endpoint returns a list of validated VRS [Variations](https://vrs.ga4gh.org/en/1.2.0/terms_and_model.html#variation).

The `/normalize` endpoint returns a [Variation Descriptor](https://vrsatile.readthedocs.io/en/latest/value_object_descriptor/vod_index.html#variation-descriptor) containing the MANE Transcript, if one is found. If a genomic query is not given a gene, `normalize` will return its GRCh38 representation.

The steps for retrieving MANE Transcript data is as follows:
1. Map starting annotation layer to genomic
2. Liftover to preferred GRCh38\
    *We only support lifting over from GRCh37.*
3. Select preferred compatible annotation
    1. MANE Select
    2. MANE Plus Clinical
    3. Longest Compatible Remaining Transcript
4. Map back to starting annotation layer

## Backend Services

Variation Normalization relies on some local data caches which you will need to set up. It uses pipenv to manage its environment, which you will also need to install.

Once pipenv is installed:
```commandline
pipenv shell
pipenv lock
pipenv sync
```

### Gene Normalizer

Variation Normalization relies on data from [Gene Normalization](https://github.com/cancervariants/gene-normalization). You must load all sources _and_ merged concepts.

You must also have Gene Normalization's DynamoDB running for the application to work.

For more information about the gene-normalizer, visit the [README](https://github.com/cancervariants/gene-normalization/blob/main/README.md).

### SeqRepo
Variation Normalization relies on [seqrepo](https://github.com/biocommons/biocommons.seqrepo), which you must download yourself.

Variation Normalizer uses seqrepo to retrieve sequences at given positions on a transcript.

From the _root_ directory:
```
pip install seqrepo
sudo mkdir /usr/local/share/seqrepo
sudo chown $USER /usr/local/share/seqrepo
seqrepo pull -i 2021-01-29
```

### UTA
Variation Normalizer also uses [uta](https://github.com/biocommons/uta).

Variation Normalizer uses UTA to retrieve MANE Transcript data.

_The following commands will likely need modification appropriate for the installation environment._
1. Install [PostgreSQL](https://www.postgresql.org/)
2. Create user and database.

    ```
    $ createuser -U postgres uta_admin
    $ createuser -U postgres anonymous
    $ createdb -U postgres -O uta_admin uta
    ```

3. To install locally, from the _variation/data_ directory:
```
export UTA_VERSION=uta_20210129.pgd.gz
curl -O http://dl.biocommons.org/uta/$UTA_VERSION
gzip -cdq ${UTA_VERSION} | grep -v "^REFRESH MATERIALIZED VIEW" | psql -h localhost -U uta_admin --echo-errors --single-transaction -v ON_ERROR_STOP=1 -d uta -p 5433
```

To connect to the UTA database, you can use the default url (`postgresql://uta_admin@localhost:5433/uta/uta_20210129`). If you use the default url, you must either set the password using environment variable `UTA_PASSWORD` or setting the parameter `db_pwd` in the UTA class.

If you do not wish to use the default, you must set the environment variable `UTA_DB_URL` which has the format of `driver://user:pass@host/database/schema`.

### PyLiftover

Variation Normalizer uses [PyLiftover](https://github.com/konstantint/pyliftover) to convert GRCh37 coordinates to GRCh38 coordinates.

## Starting the Variation Normalization Service Locally
`gene-normalizer`s dynamodb and the `uta` database must be running.

To start the service, run the following:

```commandline
uvicorn variation.main:app --reload
```

Next, view the OpenAPI docs on your local machine:
http://127.0.0.1:8000/variation

### Init coding style tests
Code style is managed by [flake8](https://github.com/PyCQA/flake8) and checked prior to commit.

We use [pre-commit](https://pre-commit.com/#usage) to run conformance tests.

This ensures:

* Check code style
* Check for added large files
* Detect AWS Credentials
* Detect Private Key

Before first commit run:

```
pre-commit install
```

### Testing
From the _root_ directory of the repository:
```
pytest tests/
```
