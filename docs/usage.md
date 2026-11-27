# Usage

## HTTP API

### Set the service URL

[Start the local API with Docker](installation.md#docker-recommended), then set
`NORMALIZER_URL` in the terminal where you'll run the examples. They use `curl`
to send requests and Python 3 to format the responses.

```shell
export NORMALIZER_URL="http://localhost:8001/variation"
```

### Translate a variant with `to_vrs`

Use `to_vrs` to list valid VRS representations of a variant. For example,
`BRAF V600E` can match several protein sequences. The endpoint uses fully
justified allele normalization and keeps the original assembly and transcript
for each match.

```shell
curl --fail-with-body --get "${NORMALIZER_URL%/}/to_vrs" \
  --data-urlencode 'q=BRAF V600E' \
  | python3 -m json.tool --indent 2
```

Expected response excerpt:

```jsonc
{
  "search_term": "BRAF V600E",
  "warnings": [],
  "variations": [
    // Protein sequence: NP_001341538.1:p.Val600Glu
    {
      "id": "ga4gh:VA.GGJOybg6mckctDlXxLY1kHZQ6dbB0U75",
      "location": {
        "id": "ga4gh:SL.0Y2ZW1zB9mf0qGesJx1kDYtwfB-67Gnf",
        "digest": "0Y2ZW1zB9mf0qGesJx1kDYtwfB-67Gnf",
        "type": "SequenceLocation",
        "start": 599,
        "end": 600,
        "sequenceReference": {
          "type": "SequenceReference",
          "refgetAccession": "SQ.ZJwurRo2HLY018wghYjDKSfIlEH0Y8At"
        }
      },
      "state": {
        "type": "LiteralSequenceExpression",
        "sequence": "E"
      },
      "type": "Allele",
    },
    // NP_001365403.1:p.Val600Glu
    {
      "id": "ga4gh:VA.jSy0uhhLefGH3396djPYcJeSVyDvRYGc",
      "type": "Allele",
      "digest": "jSy0uhhLefGH3396djPYcJeSVyDvRYGc",
      "location": {
        "id": "ga4gh:SL.CowRxWqyJfqVwlfs5YdswnrnZDQL5QCi",
        "type": "SequenceLocation",
        "digest": "CowRxWqyJfqVwlfs5YdswnrnZDQL5QCi",
        "sequenceReference": {
          "type": "SequenceReference",
          "refgetAccession": "SQ.lKdPZpuT-VNvRuKDjsUItNgutfWYgWQd"
        },
        "start": 599,
        "end": 600,
        "sequence": "V"
      },
      "state": {
        "type": "LiteralSequenceExpression",
        "sequence": "E"
      }
    }
  ]
}
```

### Normalize a variant with `normalize`

Use `normalize` when you need one VRS representation. For protein and coding
DNA queries, it prefers the MANE reference transcript when available. It
normalizes alleles and uses [cool-seq-tool](https://github.com/GenomicMedLab/cool-seq-tool)
to map genomic coordinates to GRCh38.

```shell
curl --fail-with-body --get "${NORMALIZER_URL%/}/normalize" \
  --data-urlencode 'q=BRAF V600E' \
  | python3 -m json.tool --indent 2
```

Expected response excerpt:

```json
{
  "variation_query": "BRAF V600E",
  "warnings": [],
  "variation": {
    "id": "ga4gh:VA.j4XnsLZcdzDIYa5pvvXM7t1wn9OITr0L",
    "type": "Allele",
    "digest": "j4XnsLZcdzDIYa5pvvXM7t1wn9OITr0L",
    "location": {
      "id": "ga4gh:SL.t-3DrWALhgLdXHsupI-e-M00aL3HgK3y",
      "type": "SequenceLocation",
      "digest": "t-3DrWALhgLdXHsupI-e-M00aL3HgK3y",
      "sequenceReference": {
        "type": "SequenceReference",
        "refgetAccession": "SQ.cQvw4UsHHRRlogxbWCB8W-mKD4AraM9y"
      },
      "start": 599,
      "end": 600,
      "sequence": "V"
    },
    "state": {
      "type": "LiteralSequenceExpression",
      "sequence": "E"
    }
  }
}
```

### Explore other endpoints

Open the [local interactive API documentation](http://localhost:8001/variation)
to try an endpoint and inspect its parameters and response.

## Python API

Configure the [required resources](installation.md#required-resources) before creating a query handler.
This example uses [Docker services](installation.md#docker-recommended).

```python
from os import environ

environ["SEQREPO_ROOT_DIR"] = "/usr/local/share/seqrepo/2024-12-20"
environ["UTA_DB_URL"] = "postgresql://anonymous@localhost:5432/uta/uta_20241220"
environ["AWS_ACCESS_KEY_ID"] = "dummy"
environ["AWS_SECRET_ACCESS_KEY"] = "dummy"
environ["AWS_DEFAULT_REGION"] = "us-east-2"
environ["GENE_NORM_DB_URL"] = "http://127.0.0.1:8000"

import asyncio

from variation.query import QueryHandler


async def main():
    handler = QueryHandler()
    result = await handler.normalize_handler.normalize("BRAF V600E")
    print(result.model_dump_json(exclude_none=True, indent=2))


asyncio.run(main())
```
