# A simple container for variant-service.
# Runs service on port 80.
# Healthchecks service up every 5m.  

FROM python:3.9
RUN apt update ; apt install -y rsync
RUN pip install pipenv uvicorn[standard]
ENV SEQREPO_ROOT_DIR=/usr/local/share/seqrepo/2021-01-29
#ENV GENE_NORM_DB_URL=http://localhost:8001
ENV GENE_NORM_DB_URL=http://dynamodb:8001
ENV AWS_ACCESS_KEY_ID = 'DUMMYIDEXAMPLE'
ENV AWS_SECRET_ACCESS_KEY = 'DUMMYEXAMPLEKEY'
ENV AWS_DEFAULT_REGION = 'us-west-2'
COPY . /app
WORKDIR /app
RUN if [ ! -f "Pipfile.lock" ] ; then pipenv lock ; else echo Pipfile.lock exists ; fi
RUN pipenv sync
EXPOSE 80
HEALTHCHECK --interval=5m --timeout=3s \
    CMD curl -f http://localhost/variation || exit 1

CMD cd src pipenv run uvicorn variation.main:app --log-level debug --port 80 --host 0.0.0.0
#CMD pipenv run uvicorn variation.main:app --reload
