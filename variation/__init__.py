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

logging.getLogger("uta_tools").setLevel(logging.INFO)
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


def data_download(path: Path, domain: str, dir: str, fn: str) -> None:
    """Download files using FTP.

    :param Path path: The path to the file
    :param str domain: The domain of the FTP site
    :param str dir: The directory that the file is located in
    :param str fn: The file name to download
    """
    if "VARIATION_NORM_EB_PROD" not in environ and not Path(path).exists():
        with FTP(domain) as ftp:
            ftp.login()
            ftp.cwd(dir)
            if fn.endswith(".gz") and not path.endswith(".gz"):
                path += ".gz"
            with open(f"{path}", "wb") as fp:
                ftp.retrbinary(f"RETR {fn}", fp.write)
            if fn.endswith(".gz"):
                with gzip.open(f"{path}", "rb") as f_in:
                    dest = path[:-3]
                    with open(dest, "wb") as f_out:
                        shutil.copyfileobj(f_in, f_out)
                remove(f"{path}")


AMINO_ACID_PATH = environ.get("AMINO_ACID_PATH", f"{APP_ROOT}/data/amino_acids.csv")
UTA_DB_URL = environ.get("UTA_DB_URL",
                         "postgresql://uta_admin@localhost:5433/uta/uta_20210129")
