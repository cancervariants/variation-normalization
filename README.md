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

## Normalization

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
| [main](https://github.com/cancervariants/variation-normalization/tree/main) | 0.6.X | 0.1.X | [1.X.X](https://github.com/ga4gh/vrs) |
| [staging](https://github.com/cancervariants/variation-normalization/tree/staging) | 0.8.X | 0.3.X | [2.0-alpha](https://github.com/ga4gh/vrs/tree/2.0-alpha) |


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
seqrepo pull -i 2021-01-29  # Replace with latest version using `seqrepo list-remote-instances` if outdated
```

If you get an error similar to the one below:

```shell
PermissionError: [Error 13] Permission denied: '/usr/local/share/seqrepo/2021-01-29._fkuefgd' -> '/usr/local/share/seqrepo/2021-01-29'
```

You will want to do the following:\
(*Might not be ._fkuefgd, so replace with your error message path*)

```shell
sudo mv /usr/local/share/seqrepo/2021-01-29._fkuefgd /usr/local/share/seqrepo/2021-01-29
exit
```

Use the `SEQREPO_ROOT_DIR` environment variable to set the path of an already existing SeqRepo directory. The default is `/usr/local/share/seqrepo/latest`.

#### UTA

Variation Normalizer also uses [**C**ommon **O**perations **O**n **L**ots-of **Seq**uences Tool (cool-seq-tool)](https://github.com/GenomicMedLab/cool-seq-tool) which uses [UTA](https://github.com/biocommons/uta) as the underlying PostgreSQL database.

_The following commands will likely need modification appropriate for the installation environment._

1. Install [PostgreSQL](https://www.postgresql.org/)
2. Create user and database.

    ```shell
    createuser -U postgres uta_admin
    createuser -U postgres anonymous
    createdb -U postgres -O uta_admin uta
    ```

3. To install locally, from the _variation/data_ directory:

```shell
export UTA_VERSION=uta_20210129.pgd.gz
curl -O http://dl.biocommons.org/uta/$UTA_VERSION
gzip -cdq ${UTA_VERSION} | grep -v "^REFRESH MATERIALIZED VIEW" | psql -h localhost -U uta_admin --echo-errors --single-transaction -v ON_ERROR_STOP=1 -d uta -p 5433
```

##### UTA Installation Issues

If you have trouble installing UTA, you can visit [these two READMEs](https://github.com/ga4gh/vrs-python/tree/main/docs/setup_help).

##### Connecting to the UTA database

To connect to the UTA database, you can use the default url (`postgresql://uta_admin@localhost:5433/uta/uta_20210129`). If you do not wish to use the default, you must set the environment variable `UTA_DB_URL` which has the format of `driver://user:pass@host:port/database/schema`.

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

Docker Setup:
This Section deals with setting up of Variation Normalizer via docker. For detailed instructions on variation Normalizer and its developer setup , Please refer to the variation normalizer Home page : https://github.com/cancervariants/variation-normalization/tree/main

The Variation Normalizer depends upon several Modules , therefore its recommended to setup docker containers of these Modules before starting the Variation Normalizer container.Please ensure the target machine( where the Variation Normalizer is to be deployed as docker isntalled in it.Otherwise docker commands wont work.
To Create Docker network , Please type following command.
command : docker network create <"name of the network> for e.g we have used "tulip-net"
Please follow below steps for Docker Setup of Variation Normalizer and its dependant containers.

SeqRepo Variation Normalizer depends on SeqRepo database. We need to create docker image for Seqrepo. It is recomended to start first with this image as volume attached to Seqrepo takes time to download and its size is depending upon the version is 10 GB +.
a.) Pull the image from Docker Hub Repository by typing following command in terminal.
Command : docker pull biocommons/seqrepo
b.) This will initiate display the output something like this:
Using default tag: latest latest: Pulling from biocommons/seqrepo 125a6e411906: Pull complete
4da135235d92: Pull complete
abfb8a2bf499: Pull complete
c987b6c75b9d: Pull complete
6cafe4b33812: Pull complete
03f7d4217df5: Pull complete
Digest: sha256:0390108e54c500f72afe5b187ecfb1eb9ef14f21fdc0af18e819660e7c9430c4 Status: Downloaded newer image for biocommons/seqrepo:latest docker.io/biocommons/seqrepo:latest c.) Once the image is downloaded , Start the container with the command : docker run -net <"name of the network> --name seqrepo biocommons/seqrepo The Name of the network is the network name which was created above. Running the above command will start downloading the sequences file required by Variation Normalizer. By default the volume of this container is sharable. Other containers which are on same network can access it by appending this the docker command: --volumes-from seqrepo where seqrepo is the name of the container.For efficiency , the container can be run in daemon mode or seperate terminal so that other tasks can be performed in parallel.

UTA
The Postgres UTA instance is another dependancy required for Variation Normalizer. To setup Container for UTA postgres Db instance. Follow the following steps:
a.) Pull the image from Docker Hub Repository by typing following command in terminal.
Set the uta_v env variable by typing command uta_v=<"name of the version>. For eg uta_v=uta_20210129b.
Command : docker pull biocommons/uta:$uta_v
b.) Once the image is downnloaded, Start the container with the command :
docker run 
-d 
-e POSTGRES_PASSWORD=some-password-that-you-make-up 
-v /tmp:/tmp 
-v uta_vol:/var/lib/postgresql/data 
--name $uta_v 
--net=<"name of the network> \
biocommons/uta:$uta_v

Dynamo db
The AWS provides docker image for the local instance. The Dynamo DB even though as a local instance requires AWS username and AWS password. We can provide dummy values for these environment variables. These variables have been initialized in the docker file.
a.) Pull the image from Docker Hub Repository and Start the container with the command in terminal.
Command : docker run --net tulip-net -d --name dynamodb -p 8001:8001 amazon/dynamodb-local:1.18.0 -jar DynamoDBLocal.jar -port 8001

Variation Normalizer
There is no image hosted on Docker hub for the Variation Normalizer. Hence we need to build image for Variation Normalizer from the docker File. The Docker File is already there in the repo.
a.) To build the image from the docker file. Run the command from the directory where Docker File is there.
command : docker build -t variation-normalization .
b.) Once the image is created, Start the container with the command :
command : docker run --net <network name> --name variationnormalizer -p 8000:80 --volumes-from seqrepo <"image name">:<"tag name">