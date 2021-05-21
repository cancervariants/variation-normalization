"""The Variant Normalization package."""
from .version import __version__  # noqa: F401
from pathlib import Path
import os
import logging
from ftplib import FTP
from os import environ
from gene.query import QueryHandler as GeneQueryHandler


APP_ROOT = Path(__file__).resolve().parents[0]

if 'VARIANT_NORM_EB_PROD' in os.environ:
    # Elastic beanstalk
    LOG_FN = '/tmp/sample-app.log'
    environ['GENE_NORM_EB_PROD'] = "true"
else:
    LOG_FN = f'{APP_ROOT}/variant/variant.log'

logging.basicConfig(
    filename='variant.log',
    format='[%(asctime)s] - %(name)s - %(levelname)s : %(message)s')
logger = logging.getLogger('variant')
logger.setLevel(logging.DEBUG)


# Default DynamoDB url is http://localhost:8002
# To use a different connection, set `GENE_NORM_DB_URL`
GENE_NORMALIZER = GeneQueryHandler()


def data_download(path, domain, dir, fn):
    """Download files using FTP.

    :param Path path: The path to the file
    :param str domain: The domain of the FTP site
    :param str dir: The directory that the file is located in
    :param str fn: The file name to download
    """
    if not Path(path).exists():
        with FTP(domain) as ftp:
            ftp.login()
            ftp.cwd(dir)
            with open(path, 'wb') as fp:
                ftp.retrbinary(f'RETR {fn}', fp.write)


SEQREPO_DATA_PATH = f"{APP_ROOT}/data/seqrepo/latest"
TRANSCRIPT_MAPPINGS_PATH = f"{APP_ROOT}/data/transcript_mapping.tsv"
AMINO_ACID_PATH = f"{APP_ROOT}/data/amino_acids.csv"
HGNC_GENE_SYMBOL_PATH = f"{APP_ROOT}/data/hgnc_gene_symbols.txt"
data_download(HGNC_GENE_SYMBOL_PATH, 'ftp.ebi.ac.uk',
              'pub/databases/genenames/new/tsv/', 'hgnc_complete_set.txt')
REFSEQ_GENE_SYMBOL_PATH = f"{APP_ROOT}/data/refseq_gene_symbols.txt"
data_download(REFSEQ_GENE_SYMBOL_PATH, 'ftp.ncbi.nih.gov',
              'refseq/H_sapiens/RefSeqGene/', 'LRG_RefSeqGene')
SEQREPO_REST_SERVICE_URL = "http://localhost:5000/seqrepo"
