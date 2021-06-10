"""Module for accessing UTA database."""
from typing import Dict, Optional, Tuple, List
import psycopg2
import psycopg2.extras
from six.moves.urllib import parse as urlparse
import logging
from variation import UTA_DB_URL
from pyliftover import LiftOver
from os import environ

logger = logging.getLogger('variation')
logger.setLevel(logging.DEBUG)


GRCH_TO_HG = {
    'GRCh36': 'hg18',
    'GRCh37': 'hg19',
    'GRCh38': 'hg38'
}


class UTA:
    """Class for accessing UTA database."""

    def __init__(self, db_url=UTA_DB_URL, db_pwd=None) -> None:
        """Initialize UTA class.

        :param str db_url: UTA DB url
        """
        self.db_url = self._update_db_url(db_pwd, db_url)
        self.url = _parse_url(self.db_url)
        self.schema = self.url.schema
        self.args = self._get_conn_args()
        self.conn = psycopg2.connect(**self.args)
        self.conn.autocommit = True
        self.cursor = \
            self.conn.cursor(cursor_factory=psycopg2.extras.DictCursor)

    @staticmethod
    def _update_db_url(db_pwd, db_url) -> Optional[str]:
        """Return new db_url containing password.

        :param str db_pwd: uta_admin's user password
        :param str db_url: PostgreSQL db url
        :return: PostgreSQL db url with password included
        """
        if not db_pwd and 'UTA_PASSWORD' not in environ:
            raise Exception('Environment variable UTA_PASSWORD '
                            'or `db_pwd` param must be set')
        else:
            uta_password_in_environ = 'UTA_PASSWORD' in environ
            db_url = db_url.split('@')
            if uta_password_in_environ and db_pwd:
                if db_pwd != environ['UTA_PASSWORD']:
                    raise Exception('If both environment variable UTA_'
                                    'PASSWORD and param db_pwd is set, '
                                    'they must both be the same')
            else:
                if uta_password_in_environ and not db_pwd:
                    db_pwd = environ['UTA_PASSWORD']
            return f"{db_url[0]}:{db_pwd}@{db_url[1]}"

    def _get_conn_args(self) -> Dict:
        """Return connection arguments.

        :return: A dictionary containing db credentials
        """
        return dict(
            host=self.url.hostname,
            port=self.url.port,
            database=self.url.database,
            user=self.url.username,
            password=self.url.password,
            application_name='variation',
        )

    def get_alt_tx_data(self, ac, pos) -> \
            Optional[Tuple[str, str, Tuple[int, int], str]]:
        """Get altered transcript data given transcript data.

        :param str ac: cDNA transcript
        :param tuple pos: [cDNA pos start, cDNA pos end]
        :return: [Gene, NC accession,
            [Genomic change pos start, Genomic change pos end]]
        """
        query = (
            f"""
            SELECT *
            FROM {self.schema}.tx_exon_aln_v
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
        if result[4] == -1:
            strand = '-'
        else:
            strand = '+'
        tx_pos_range = result[6], result[7]
        alt_pos_range = result[8], result[9]

        if (tx_pos_range[1] - tx_pos_range[0]) != \
                (alt_pos_range[1] - alt_pos_range[0]):
            logger.warning(f"{nc_accession} tx_pos_range {tx_pos_range} "
                           f"is not the same length as alt_pos_range "
                           f"{alt_pos_range}.")
            return None

        tx_pos_change = pos[0] - tx_pos_range[0], tx_pos_range[1] - pos[1]
        alt_pos = (alt_pos_range[0] + tx_pos_change[0],
                   alt_pos_range[1] - tx_pos_change[1])

        return gene, nc_accession, alt_pos, strand

    def liftover_to_38(self, nc_accession, alt_pos_range, strand=None) \
            -> Optional[Tuple[Tuple[int, int], str]]:
        """Liftover NC accession to GRCh38 version.

        :param str nc_accession: NC Accession
        :param tuple alt_pos_range:
            [Altered transcript pos start, Altered transcript pos end]
        :param str strand: Strand
        :return: [Altered transcript pos start, Altered transcript pos end]
            for GRCh38
        """
        query = (
            f"""
            SELECT descr
            FROM {self.schema}._seq_anno_most_recent
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
            liftover_start_i = \
                self._get_liftover(lo, chromosome, alt_pos_range[0])
            liftover_end_i = \
                self._get_liftover(lo, chromosome, alt_pos_range[1])
            if liftover_start_i is None or liftover_end_i is None:
                return None

            alt_pos_range = liftover_start_i[1], liftover_end_i[1]
            if liftover_start_i[2] != liftover_end_i[2]:
                logger.warning("Strand must be the same.")
            strand = liftover_start_i[2]

        return alt_pos_range, strand

    @staticmethod
    def _get_liftover(lo, chromosome, pos):
        liftover = lo.convert_coordinate(chromosome, pos)
        if liftover is None or len(liftover) == 0:
            logger.warning(f"{pos} does not exist on {chromosome}")
            return None
        else:
            return liftover[0]

    def p_to_c_ac(self, p_ac) -> List[str]:
        """Return c. accession from p. accession.

        :param str p_ac: Protein accession
        :return: List of rows containing c. accessions that are associated
            with the given p. accession.
        """
        query = (
            f"""
            SELECT *
            FROM {self.schema}.associated_accessions
            WHERE pro_ac = '{p_ac}'
            """
        )
        self.cursor.execute(query)
        return self.cursor.fetchall()


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
