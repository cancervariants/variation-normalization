"""Module for accessing UTA database."""
from typing import Dict
import psycopg2
from six.moves.urllib import parse as urlparse
import logging
from variant import UTA_DB_URL

logger = logging.getLogger('variant')
logger.setLevel(logging.DEBUG)


class UTA:
    """Class for accessing UTA database."""

    def __init__(self, db_url=UTA_DB_URL) -> None:
        """Initialize UTA class.

        :param str db_url: UTA DB url
        """
        self.db_url = db_url
        self.conn = psycopg2.connect(**self._get_conn_args())
        self.conn.autocommit = True

    def _get_conn_args(self) -> Dict:
        """Return connection arguments.

        :return: A dictionary containing db credentials
        """
        url = _parse_url(self.db_url)
        return dict(
            host=url.hostname,
            port=url.port,
            database=url.database,
            user=url.username,
            password=url.password,
            application_name='variant',
        )


class ParseResult(urlparse.ParseResult):
    """Subclass of url.ParseResult that adds database and schema methods,
    and provides stringification.
    Source: https://github.com/biocommons/hgvs
    """

    def __new__(cls, pr):
        """Create new instance."""
        return super(ParseResult, cls).__new__(cls, *pr)

    @property
    def database(self):
        """Create database property."""
        path_elems = self.path.split("/")
        return path_elems[1] if len(path_elems) > 1 else None

    @property
    def schema(self):
        """Create schema property."""
        path_elems = self.path.split("/")
        return path_elems[2] if len(path_elems) > 2 else None


def _parse_url(db_url) -> ParseResult:
    """Parse database connection urls into components.
    Source: https://github.com/biocommons/hgvs

    :param str db_url: UTA DB url
    :return: Parsed Result containing URL components
    """
    return ParseResult(urlparse.urlparse(db_url))
