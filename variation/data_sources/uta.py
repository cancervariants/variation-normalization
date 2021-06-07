"""Module for accessing UTA database."""
from typing import Dict, Optional, Tuple
import psycopg2
import psycopg2.extras
from six.moves.urllib import parse as urlparse
import logging
from variation import UTA_DB_URL
from pyliftover import LiftOver

logger = logging.getLogger('variation')
logger.setLevel(logging.DEBUG)


GRCH_TO_HG = {
    'GRCh36': 'hg18',
    'GRCh37': 'hg19',
    'GRCh38': 'hg38'
}


class UTA:
    """Class for accessing UTA database."""

    def __init__(self, db_url=UTA_DB_URL) -> None:
        """Initialize UTA class.

        :param str db_url: UTA DB url
        """
        self.db_url = db_url
        self.conn = psycopg2.connect(**self._get_conn_args())
        self.conn.autocommit = True
        self.cursor = \
            self.conn.cursor(cursor_factory=psycopg2.extras.DictCursor)

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
            application_name='variation',
        )

    def get_alt_tx_data(self, ac, pos) -> \
            Optional[Tuple[str, str, Tuple[int, int]]]:
        """Get altered transcript data given transcript data.

        :param str ac: cDNA transcript
        :param tuple pos: [cDNA pos start, cDNA pos end]
        :return: [Gene, NC accession,
            [Altered transcript start, Altered transcript end]]
        """
        query = (
            f"""
            SELECT *
            FROM tx_exon_aln_v
            WHERE tx_ac='{ac}'
            AND alt_ac LIKE 'NC_00%'
            AND {pos[0]} BETWEEN tx_start_i AND tx_end_i
            AND {pos[1]} BETWEEN tx_start_i AND tx_end_i
            ORDER BY ord, alt_ac;
            """
        )
        self.cursor.execute(query)
        results = self.cursor.fetchall()
        if not results:
            logger.warning(f"Unable to find transcript alignment for {ac}")
            return None
        result = results[-1]
        gene = result[0]
        nc_accession = result[2]
        tx_pos_range = result[6], result[7]
        alt_pos_range = result[8], result[9]

        if (tx_pos_range[1] - tx_pos_range[0]) != \
                (alt_pos_range[1] - alt_pos_range[0]):
            logger.warning(f"{nc_accession} tx_pos_range {tx_pos_range} "
                           f"is not the same length as alt_pos_range "
                           f"{alt_pos_range}.")
            return None

        return gene, nc_accession, alt_pos_range

    def liftover_to_38(self, nc_accession, alt_pos_range) -> Tuple[int, int]:
        """Liftover NC accession to GRCh38 version.

        :param str nc_accession: NC Accession
        :param tuple alt_pos_range:
            [Altered transcript pos start, Altered transcript pos end]
        :return: [Altered transcript pos start, Altered transcript pos end]
            for GRCh38
        """
        query = (
            f"""
            SELECT descr
            FROM _seq_anno_most_recent
            WHERE ac = '{nc_accession}'
            """
        )
        self.cursor.execute(query)
        result = self.cursor.fetchone()
        if result and result[0]:
            descr = result[0].split(',')
            chromosome = f"chr{descr[0].split()[-1]}"
            assembly = f"GRCh{descr[1].split('.')[0].split('GRCh')[-1]}"

            # Get most recent assembly version position
            lo = LiftOver(GRCH_TO_HG[assembly], 'hg38')
            alt_pos_range = \
                (lo.convert_coordinate(chromosome, alt_pos_range[0]),
                 lo.convert_coordinate(chromosome, alt_pos_range[1]))

        return alt_pos_range


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
