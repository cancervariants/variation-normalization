"""The Variant Normalization package."""
from pathlib import Path
import os
import logging

__version__ = "0.1.11"

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

SEQREPO_DATA_PATH = f"{APP_ROOT}/data/seqrepo/2021-01-29"
TRANSCRIPT_MAPPINGS_PATH = f"{APP_ROOT}/data/transcript_mapping.tsv"
AMINO_ACID_PATH = f"{APP_ROOT}/data/amino_acids.csv"
GENE_SYMBOL_PATH = f"{APP_ROOT}/data/gene_symbols.txt"
