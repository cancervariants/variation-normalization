"""The Variation Normalization package."""
from pathlib import Path
import logging
from ftplib import FTP
from os import environ, remove
import gzip
import shutil

from .version import __version__  # noqa: F401


APP_ROOT = Path(__file__).resolve().parents[0]

if "VARIATION_NORM_EB_PROD" in environ:
    LOG_FN = "/tmp/variation.log"
else:
    LOG_FN = "variation.log"

logging.basicConfig(
    filename=LOG_FN,
    format="[%(asctime)s] - %(name)s - %(levelname)s : %(message)s")
logger = logging.getLogger("variation")
logger.setLevel(logging.DEBUG)
logger.handlers = []

logging.getLogger("cool_seq_tool").setLevel(logging.INFO)
logging.getLogger("boto3").setLevel(logging.INFO)
logging.getLogger("botocore").setLevel(logging.INFO)
logging.getLogger("urllib3").setLevel(logging.INFO)
logging.getLogger("python_jsonschema_objects").setLevel(logging.INFO)
logging.getLogger("hgvs.parser").setLevel(logging.INFO)
logging.getLogger("biocommons.seqrepo.seqaliasdb.seqaliasdb").setLevel(logging.INFO)
logging.getLogger("biocommons.seqrepo.fastadir.fastadir").setLevel(logging.INFO)
logging.getLogger("asyncio").setLevel(logging.INFO)

if "VARIATION_NORM_EB_PROD" in environ:
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    logger.addHandler(ch)


UTA_DB_URL = environ.get("UTA_DB_URL",
                         "postgresql://uta_admin@localhost:5433/uta/uta_20210129")
