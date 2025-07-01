# A simple container for variation-normalizer service.
# Runs service on port 80.
# Health checks service up every 5m.
FROM python:3.12-slim

WORKDIR /app

ARG VERSION

ENV SETUPTOOLS_SCM_PRETEND_VERSION_FOR_VARIATION_NORMALIZER=$VERSION

COPY src/ ./src/
COPY pyproject.toml .

RUN apt-get update && \
    apt-get install -y libpq-dev gcc && \
    rm -rf /var/lib/apt/lists/*

RUN pip install --upgrade pip setuptools setuptools_scm
RUN pip install '.'

EXPOSE 80
HEALTHCHECK --interval=5m --timeout=3s \
    CMD curl -f http://localhost/variation || exit 1

CMD ["uvicorn", "variation.main:app", "--port", "80", "--host", "0.0.0.0"]
