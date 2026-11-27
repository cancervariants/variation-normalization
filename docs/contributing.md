# Contributing to Variation Normalization

## Before you start

1. Check the open issues and pull requests before starting work.
2. For a substantial change, share your approach in a new or existing issue or start a GitHub discussion.
3. Keep each pull request focused on one change and include tests or
   documentation when they are relevant.

## Development branches

If you are not a repository maintainer, first [fork the repository on
GitHub](https://github.com/cancervariants/variation-normalization/fork).

Create your branch from `staging`.

```shell
git clone https://github.com/<your-user>/variation-normalization.git
cd variation-normalization
git remote add upstream https://github.com/cancervariants/variation-normalization.git
git fetch upstream staging
git switch -c <short-description> upstream/staging
```

Open feature and bug-fix pull requests against `staging`. The `main` branch holds released code.

## Docker development setup

Use [compose-dev.yaml](https://github.com/cancervariants/variation-normalization/blob/main/compose-dev.yaml)
to build the API from your checkout. It also starts UTA and the Gene Normalizer
database. The [installation guide](installation.md#docker-recommended) describes
the published image setup in `compose.yaml`.

After creating your contribution branch, prepare the
[SeqRepo snapshot and UTA dump](installation.md#prepare-the-data). Run these
commands from the repository root:

```shell
export VERSION=$(git describe --tags --abbrev=0)
docker volume create uta_vol
docker compose -f compose-dev.yaml up --build -d
```

Follow the logs during the initial UTA load:

```shell
docker compose -f compose-dev.yaml logs -f
```

Once the databases are initialized, open <http://localhost:8001/variation>.
If the API started before UTA finished loading, restart it:

```shell
docker compose -f compose-dev.yaml restart api
```

The image contains a copy of your source. Rebuild after changing Python code:

```shell
docker compose -f compose-dev.yaml up --build -d api
```

Stop the development services with:

```shell
docker compose -f compose-dev.yaml down
```

Use one Compose setup at a time: both files use the same service ports and
external `uta_vol` volume. The volume retains UTA data across shutdowns.

## Python development environment

To run tests and code checks on the host, use Python 3.11 or newer and install
the development and test dependencies:

```shell
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[dev,tests]'
```

The full test suite needs the [required resources](installation.md#required-resources).
With the development containers running, configure the host's
[UTA connection](installation.md#uta),
[Gene Normalizer connection](installation.md#gene-normalizer), and
`SEQREPO_ROOT_DIR` pointing to the local snapshot directory.

## Run checks locally

Run the test suite from the repository root:

```shell
pytest tests/
```

Format and lint the code with the project's configured Ruff version:

```shell
python3 -m ruff format .
python3 -m ruff check --fix .
```

Install the pre-commit hooks before making the first commit, then run them on
all files when needed:

```shell
prek install
prek run --all-files
```

CI installs dependencies and runs Ruff and pre-commit. If a missing resource
prevents a local check, mention it in the pull request and list the checks you ran.

## Pull request checklist

Before opening a pull request:

- [ ] The branch is based on the latest `staging` branch.
- [ ] The change is limited to the issue or goal being addressed.
- [ ] Tests cover changed behavior, or the pull request explains why tests
      are not applicable.
- [ ] `pytest tests/` passes when the required resources are available.
- [ ] Ruff and pre-commit checks pass.
- [ ] Documentation and examples are updated when behavior or setup changes.
- [ ] The pull request description summarizes the change and lists the
      validation commands and results.
- [ ] The relevant issue is linked in the pull request description.

## Dependency changes

Keep dependency declarations consistent:

- Add or update published package dependencies in `pyproject.toml`.
- Regenerate `requirements.txt` for the Elastic Beanstalk deployment with
  `uv pip compile pyproject.toml -o requirements.txt --no-annotate` when
  production dependencies change.
- Keep the fallback `Pipfile` in sync with production dependency changes.
- Development and test dependencies belong in the optional dependency groups
  in `pyproject.toml`.

## Reporting security issues

Report vulnerabilities privately through the repository's security reporting process.

## Documentation

Install the documentation tools and preview the site:

```shell
python3 -m pip install '.[docs]'
zensical serve
```

Check the site before submitting a change:

```shell
zensical build --clean --strict
```

Edit pages in `docs/`. API reference pages in `docs/api/` use `mkdocstrings`
directives to render Python signatures and Sphinx-style docstrings from
`src/variation`. Add a Markdown page and a navigation entry in `zensical.toml`
when documenting another module.
