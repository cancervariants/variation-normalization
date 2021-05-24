# Variant Normalization
Services and guidelines for normalizing variant terms

## Backend Services
Variant Normalization relies on some local data caches which you will need to set up. It uses pipenv to manage its environment, which you will also need to install.

### Installation
Variant Normalization relies on [seqrepo](https://github.com/biocommons/biocommons.seqrepo), which you must download yourself.

From the _variant_ directory of the repository:
```
pipenv sync
pip install seqrepo
mkdir -p data/seqrepo
seqrepo -r data/seqrepo pull -i 2021-01-29
sudo chmod -R u+w data/seqrepo
cd data/seqrepo
seqrepo_date_dir=$(ls -d */)
sudo mv $seqrepo_date_dir latest
```

Variant Normalizer also uses [uta](https://github.com/biocommons/uta).

To install:
```
uta_v=uta_20180821
docker pull biocommons/uta:$uta_v
export UTA_DB_URL=postgresql://anonymous@localhost:5432/uta/uta_20180821
docker-compose -f docker-compose.yml up
```

### Data
Variant Normalization uses [Ensembl BioMart](http://www.ensembl.org/biomart/martview) to retrieve `variant/data/transcript_mappings.tsv`. We currently use `Human Genes (GRCh38.p13)` for the dataset and the following attributes we use are: Gene stable ID, Gene stable ID version, Transcript stable ID, Transcript stable ID version, Protein stable ID, Protein stable ID version, RefSeq match transcript (MANE Select), Gene name. 

![image](biomart.png)

### Setting up Gene Normalizer
Variant Normalization `normalize` endpoint relies on data from Gene Normalization. To install:
```shell script
pip install gene-normalizer
```

To setup, follow the instructions from the [Gene Normalization README](https://github.com/cancervariants/gene-normalization). 

You must have the Gene Normalizer DynamoDB running for the variant `normalize` endpoint to work.

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

`gene-normalizer`s dynamodb must be running and run the following:
```
docker-compose -f docker-compose.yml up
```

From the _root_ directory of the repository:
```
uvicorn variant.main:app --reload
```
Next, view the OpenAPI docs on your local machine:
http://127.0.0.1:8000/variant