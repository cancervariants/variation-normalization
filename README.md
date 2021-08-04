# Variation Normalization
Services and guidelines for normalizing variation terms

## Backend Services
Variation Normalization relies on some local data caches which you will need to set up. It uses pipenv to manage its environment, which you will also need to install.

### Installation
Variation Normalization relies on [seqrepo](https://github.com/biocommons/biocommons.seqrepo), which you must download yourself.

From the _root_ directory:
```
pipenv shell
pipenv lock
pipenv sync
cd variation
pip install seqrepo
mkdir -p data/seqrepo
seqrepo -r data/seqrepo pull -i 2021-01-29
sudo chmod -R u+w data/seqrepo
cd data/seqrepo
seqrepo_date_dir=$(ls -d */)
sudo mv $seqrepo_date_dir latest
```

Variation Normalizer also uses [uta](https://github.com/biocommons/uta).

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

### Data
Variation Normalization uses [Ensembl BioMart](http://www.ensembl.org/biomart/martview) to retrieve `variation/data/transcript_mappings.tsv`. We currently use `Human Genes (GRCh38.p13)` for the dataset and the following attributes we use are: Gene stable ID, Gene stable ID version, Transcript stable ID, Transcript stable ID version, Protein stable ID, Protein stable ID version, RefSeq match transcript (MANE Select), Gene name. 

![image](biomart.png)

### Setting up Gene Normalizer
Variation Normalization `normalize` endpoint relies on data from Gene Normalization. To install:
```shell script
pip install gene-normalizer
```

To setup, follow the instructions from the [Gene Normalization README](https://github.com/cancervariants/gene-normalization). 

You must have the Gene Normalizer DynamoDB running for the variation `normalize` endpoint to work.

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

### Starting the Variation Normalization Service

`gene-normalizer`s dynamodb must be running and run the following:
```
docker-compose -f docker-compose.yml up
```

From the _root_ directory of the repository:
```
uvicorn variation.main:app --reload
```
Next, view the OpenAPI docs on your local machine:
http://127.0.0.1:8000/variation
