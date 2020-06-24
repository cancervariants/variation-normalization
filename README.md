# varlex
[![Build Status](https://travis-ci.org/cancervariants/varlex.svg?branch=master)](https://travis-ci.org/cancervariants/varlex)

Repository for the Variant Lexicon normalization service


## Backend Services

The varlex backend is a simple flask app, but it does rely on some local data caches which you will need to set up. It uses pipenv to manage its environment, which you will also need to install.

### Installation
From the _root_ directory of the repository:
```
pipenv sync
mkdir -p varlexapp/data/seqrepo
curl ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/non_alt_loci_set.txt > varlexapp/data/gene_symbols.txt
seqrepo --root-directory varlexapp/data/seqrepo pull
```

### Running
```
export FLASK_ENV=development
export FLASK_APP=varlexapp
flask run
```

### Testing
From the _root_ directory of the repository:
```
pytest tests/
```

## Frontend

Varlex includes a simple fronted that can be used for basic exploration of its features. It requires `yarn` wich can be installed with `npm install yarn` or via `brew`.

### Installation/Running
From the _frontend_ directory of the repository:
```
yarn start
```
