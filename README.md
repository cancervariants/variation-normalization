# Variant Normalization
Services and guidelines for normalizing variant terms

## Backend Services
Variant Normalization relies on some local data caches which you will need to set up. It uses pipenv to manage its environment, which you will also need to install.

### Installation
From the _root_ directory of the repository:
```
pipenv sync
mkdir -p variant/data/seqrepo
curl ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/non_alt_loci_set.txt > variant/data/gene_symbols.txt
seqrepo --root-directory variant/data/seqrepo pull
cd variant/data/seqrepo
chmod -R u+w variant/data/seqrepo/<DATE>
ln -s variant/data/seqrepo/<DATE> latest
```

#### Installing Dependencies
From the _root_ directory of the repository:
```
$ docker volume create --name=uta_vol
$ docker volume create --name=seqrepo_vol
$ docker-compose -f docker-compose.yml up
```

This should start three containers:

  * [seqrepo](https://github.com/biocommons/seqrepo): a non-redundant archive of sequences
  * [seqrepo-rest-service](https://github.com/biocommons/seqrepo-rest-service): a REST service on seqrepo (localhost:5000)
  * [uta](https://github.com/biocommons/uta): a database of transcripts and alignments (localhost:5432)

The seqrepo container will exit as soon as the data are downloaded.

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

### Starting the Variant Normalization Service
From the _root_ directory of the repository:
```
uvicorn main:app --reload
```
Next, view the OpenAPI docs on your local machine:
http://127.0.0.1:8000/variant