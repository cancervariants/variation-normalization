# Installation

## Docker (recommended)

Use [Docker](https://docs.docker.com/get-started/get-docker/) with
[Docker Compose](https://docs.docker.com/compose/) to run the published Variation
Normalizer API.

### Get the Compose files

Clone the repository to get `compose.yaml` and `uta-setup.sql`:

```shell
git clone https://github.com/cancervariants/variation-normalization.git
cd variation-normalization
```

Run the commands below from the repository directory. The Compose file uses the
published API image from [GitHub Container Registry](https://github.com/cancervariants/variation-normalization/pkgs/container/variation-normalizer-api).

!!! tip
    If you only need these two files, download them directly:

    ```shell
    curl -LO https://raw.githubusercontent.com/cancervariants/variation-normalization/main/compose.yaml
    curl -LO https://raw.githubusercontent.com/cancervariants/variation-normalization/main/uta-setup.sql
    ```

### Prepare the data

1. Download the [SeqRepo 2024-12-20 snapshot](#seqrepo) to
   `/usr/local/share/seqrepo/2024-12-20`.
2. Download `uta_20241220.pgd.gz` from <https://dl.biocommons.org/uta/>
   using a web browser and place it beside `compose.yaml`.

If SeqRepo is stored elsewhere, update both the host volume mount and the
container's `SEQREPO_ROOT_DIR` in `compose.yaml`. The environment variable must
point to the snapshot directory inside the container.

For Docker Desktop, open **Settings > Resources > File sharing** and add
`/usr/local/share/seqrepo` to the shared directories. Without access to the
snapshot, the API reports `OSError: Unable to open SeqRepo directory`.

### Start the service

The first startup loads UTA into the persistent `uta_vol` volume (this may take a few minutes).

```shell
docker volume create uta_vol
docker compose -f compose.yaml up -d
```

Compose reports the status of all three containers.

```shell
[+] up 3/3
 ✔ Container variation-normalization-uta-1                 Running
 ✔ Container variation-normalization-gene-dynamodb-local-1 Running
 ✔ Container variation-normalization-api-1                 Started
```

Open the API interface at <http://localhost:8001/variation>.

## Python library

With Python 3.11 or newer, install the library from
[PyPI](https://pypi.org/project/variation-normalizer):

```shell
python3 -m pip install variation-normalizer
```

The library requires SeqRepo, UTA, and Gene Normalizer data.

## Required resources

### SeqRepo

[SeqRepo](https://github.com/biocommons/biocommons.seqrepo) supplies reference
sequences. To download the `2024-12-20` snapshot:

```shell
python3 -m pip install biocommons.seqrepo
sudo mkdir -p /usr/local/share/seqrepo
sudo chown "$USER" /usr/local/share/seqrepo
seqrepo pull -i 2024-12-20/
```

When using the Python library, set `SEQREPO_ROOT_DIR` to the snapshot directory.
The default is `/usr/local/share/seqrepo/latest`.

### UTA

The [Universal Transcript Archive (UTA)](https://github.com/biocommons/uta)
supplies transcript alignments through Cool-Seq-Tool.

To connect a Python process on the host to the Compose database:

```shell
export UTA_DB_URL=postgresql://anonymous@localhost:5432/uta/uta_20241220
```

For a separate database deployment, follow the
[UTA installation guide](https://github.com/biocommons/uta#installing-uta-locally).

### Gene Normalizer

[Gene Normalizer](https://github.com/cancervariants/gene-normalization) supplies
normalized gene concepts. The Compose setup includes its database image.
To connect a Python process on the host to that local database:

```shell
export GENE_NORM_DB_URL=http://localhost:8000
export AWS_ACCESS_KEY_ID=dummy
export AWS_SECRET_ACCESS_KEY=dummy
export AWS_DEFAULT_REGION=us-east-2
```

For a separate database deployment, follow the
[Gene Normalizer installation guide](https://gene-normalizer.readthedocs.io/latest/install.html).
