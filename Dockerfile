# A simple container for variant-service.
# Runs service on port 80.
# Healthchecks service up every 5m.

FROM python:3.12
RUN apt update ; apt install -y rsync
COPY . /app
WORKDIR /app
RUN pip install --upgrade pip
RUN pip install -e '.[docker]'
EXPOSE 80
HEALTHCHECK --interval=5m --timeout=3s \
    CMD curl -f http://localhost/variation || exit 1

CMD uvicorn variation.main:app  --port 80 --host 0.0.0.0
