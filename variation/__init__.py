"""The Variation Normalization package."""
from .version import __version__  # noqa: F401
from pathlib import Path
import logging
from ftplib import FTP
from os import environ, remove
import gzip
import shutil
from gene.query import QueryHandler as GeneQueryHandler


APP_ROOT = Path(__file__).resolve().parents[0]

if 'VARIATION_NORM_EB_PROD' in os.environ:
    # Elastic beanstalk
    LOG_FN = '/tmp/sample-app.log'
    environ['GENE_NORM_EB_PROD'] = "true"
else:
    LOG_FN = f'{APP_ROOT}/variation/variation.log'

logging.basicConfig(
    filename='variation.log',
    format='[%(asctime)s] - %(name)s - %(levelname)s : %(message)s')
logger = logging.getLogger('variation')
logger.setLevel(logging.DEBUG)
logger.handlers = []


# Default DynamoDB url is http://localhost:8000
# To use a different connection, set `GENE_NORM_DB_URL`
from gene.query import QueryHandler as GeneQueryHandler  # noqa: E402
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
            with open(f"{path}.gz", 'wb') as fp:
                ftp.retrbinary(f'RETR {fn}', fp.write)
            if fn.endswith('.gz'):
                with gzip.open(f"{path}.gz", 'rb') as f_in:
                    with open(path, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                remove(f"{path}.gz")


SEQREPO_DATA_PATH = f"{APP_ROOT}/data/seqrepo/latest"
TRANSCRIPT_MAPPINGS_PATH = f"{APP_ROOT}/data/transcript_mapping.tsv"
AMINO_ACID_PATH = f"{APP_ROOT}/data/amino_acids.csv"
HGNC_GENE_SYMBOL_PATH = f"{APP_ROOT}/data/hgnc_gene_symbols.txt"
data_download(HGNC_GENE_SYMBOL_PATH, 'ftp.ebi.ac.uk',
              'pub/databases/genenames/new/tsv/', 'hgnc_complete_set.txt')
REFSEQ_GENE_SYMBOL_PATH = f"{APP_ROOT}/data/refseq_gene_symbols.txt"
data_download(REFSEQ_GENE_SYMBOL_PATH, 'ftp.ncbi.nih.gov',
              'refseq/H_sapiens/RefSeqGene/', 'LRG_RefSeqGene')
REFSEQ_MANE_PATH = f"{APP_ROOT}/data/MANE.GRCh38.v0.93.summary.txt"
data_download(REFSEQ_MANE_PATH, 'ftp.ncbi.nlm.nih.gov',
              'refseq/MANE/MANE_human/release_0.93/',
              'MANE.GRCh38.v0.93.summary.txt.gz')
UTA_DB_URL = 'postgresql://anonymous@localhost:5432/uta/uta_20180821'