"""Provide functions and variables for logging setup.

We may need to set up common logging protocols for a few different entry points when
running the Variation Normalizer as an app, but we shouldn't dictate anything to
downstream users. Functions in this module should only be called from locations in code
that are being executed directly by users, not from anywhere that might be imported as
a library.
"""

import logging
import os


def _quiet_upstream_libs() -> None:
    """Turn off debug logging for chatty upstream library loggers."""
    for lib in (
        "boto3",
        "botocore",
        "urllib3",
        "hgvs.parser",
        "biocommons.seqrepo.seqaliasdb.seqaliasdb",
        "biocommons.seqrepo.fastadir.fastadir",
        "asyncio",
    ):
        logging.getLogger(lib).setLevel(logging.INFO)


def configure_logging(
    log_level: int = logging.DEBUG, quiet_upstream: bool = True
) -> None:
    """Configure logging.

    :param log_level: global log level to set
    :param quiet_upstream: if True, turn off debug logging for a selection of libraries
    """
    log_filename = (
        "/tmp/variation.log"  # noqa: S108
        if "VARIATION_NORM_EB_PROD" in os.environ
        else "variation.log"
    )
    if quiet_upstream:
        _quiet_upstream_libs()
    logging.basicConfig(
        filename=log_filename,
        format="[%(asctime)s] - %(name)s - %(levelname)s : %(message)s",
    )
    logger = logging.getLogger("variation")
    logger.setLevel(log_level)

    if "METAKB_NORM_EB_PROD" in os.environ:
        # force debug logging in production server
        logger.handlers = []
        handler = logging.StreamHandler()
        handler.setLevel(logging.DEBUG)
        logger.addHandler(handler)
