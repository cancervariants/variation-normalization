# varlex
[![Build Status](https://travis-ci.org/cancervariants/varlex.svg?branch=master)](https://travis-ci.org/cancervariants/varlex)

Repository for the Variant Lexicon normalization service


## Backend Services
VarLex relies on some local data caches which you will need to set up. It uses pipenv to manage its environment, which you will also need to install.

### Installation
From the _root_ directory of the repository:
```
pipenv sync
mkdir -p varlexapp/data/seqrepo
curl ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/non_alt_loci_set.txt > varlexapp/data/gene_symbols.txt
seqrepo --root-directory varlexapp/data/seqrepo pull
cd varlexapp/data/seqrepo
chmod -R u+w varlexapp/data/seqrepo/<DATE>
ln -s varlexapp/data/seqrepo/<DATE> latest
```

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

### Starting the VarLex Service
From the _root_ directory of the repository:
```
uvicorn main:app --reload
```
Next, view the OpenAPI docs on your local machine:
http://127.0.0.1:8000/variant