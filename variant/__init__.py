"""The Variant Normalization package."""
from pathlib import Path
import os
import logging
from ftplib import FTP

__version__ = "0.1.12"

APP_ROOT = Path(__file__).resolve().parents[0]

if 'VARIANT_NORM_EB_PROD' in os.environ:
    # Elastic beanstalk
    LOG_FN = '/tmp/sample-app.log'
else:
    LOG_FN = f'{APP_ROOT}/variant/variant.log'

logger = logging.getLogger('variant')
if Path(LOG_FN).exists():
    fhandler = logging.FileHandler(filename=LOG_FN)
    formatter = logging.Formatter('[%(asctime)s] %(name)s - %(levelname)s : %(message)s')  # noqa: E501
    fhandler.setFormatter(formatter)
    logger.addHandler(fhandler)
logger.setLevel(logging.DEBUG)


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
