"""The Variant Normalization package."""
from pathlib import Path
import os
import logging

__version__ = "0.1.5"

APP_ROOT = Path(__file__).resolve().parents[0]

if 'DEV' in os.environ:
    LOG_FN = './variant/variant.log'
else:
    # Elastic beanstalk
    LOG_FN = '/tmp/sample-app.log'

logger = logging.getLogger('variant')
fhandler = logging.FileHandler(filename=LOG_FN)
formatter = logging.Formatter('[%(asctime)s] %(name)s - %(levelname)s : %(message)s')  # noqa: E501
fhandler.setFormatter(formatter)
logger.addHandler(fhandler)
logger.setLevel(logging.DEBUG)

SEQREPO_DATA_PATH = f"{APP_ROOT}/data/seqrepo/latest"
TRANSCRIPT_MAPPINGS_PATH = f"{APP_ROOT}/data/transcript_mapping.tsv"
AMINO_ACID_PATH = f"{APP_ROOT}/data/amino_acids.csv"
GENE_SYMBOL_PATH = f"{APP_ROOT}/data/gene_symbols.txt"
