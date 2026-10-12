# Contributing to Variation Normalization

Thank you for helping improve Variation Normalization. Contributions that
improve normalization behavior, tests, documentation, and developer tooling
are welcome.

## Before you start

1. Check the open issues and pull requests before starting work.
2. For a non-trivial change, open or comment on an issue so the proposed
   approach can be discussed.
3. Keep each pull request focused on one change and include tests or
   documentation when they are relevant.

## Development branches

Create contribution branches from `staging`, not `main`:

```shell
git clone https://github.com/<your-user>/variation-normalization.git
cd variation-normalization
git remote add upstream https://github.com/cancervariants/variation-normalization.git
git fetch upstream staging
git switch -c <short-description> upstream/staging
```

Pull requests for normal development changes should target the `staging`
branch. The `main` branch is used for released code and should not be used as
the target for ordinary feature or bug-fix pull requests.

## Set up a development environment

Variation Normalization supports Python 3.11 or newer. The project uses
Pipenv for the development environment:

```shell
python3 -m pip install --user pipenv
pipenv shell
pipenv update && pipenv install --dev
```

The CI lint job uses the equivalent `uv` workflow:

```shell
uv sync --extra dev
```

### Required local resources

The full test suite requires local sequence and transcript resources:

- Set up [SeqRepo](https://github.com/biocommons/biocommons.seqrepo) as
  described in the [README](README.md#seqrepo).
- Make the UTA database available as described in the
  [README](README.md#uta).
- Alternatively, use the Docker setup described in the
  [README](README.md#docker-installation-preferred).

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
pre-commit install
pre-commit run --all-files
```

The CI checks also run dependency installation, Ruff formatting/lint checks,
and the pre-commit hooks. If a check cannot be run locally because a required
resource is unavailable, describe that limitation and the checks you did run
in the pull request.

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

Please do not disclose security vulnerabilities in a public issue. Use the
repository's private security reporting process instead.
