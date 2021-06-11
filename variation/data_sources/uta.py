"""Module for accessing UTA database."""
from typing import Dict, Optional, List, Tuple
import psycopg2
import psycopg2.extras
from six.moves.urllib import parse as urlparse
import logging
from variation import UTA_DB_URL
from pyliftover import LiftOver
from os import environ

logger = logging.getLogger('variation')
logger.setLevel(logging.DEBUG)


# Assembly mappings
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
        :param str db_pwd: UTA user uta_admin's password
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

    def get_mane_tx_c_data(self, mane_c_ac, nc_ac, genomic_change_range)\
            -> List:
        """Get MANE Transcript c. data.

        :param str mane_c_ac: MANE Transcript c. accession
        :param str nc_ac: NC accession
        :param tuple[int, int] genomic_change_range: [genomic change start
            position, genomic change end position]
        :return: Data about MANE Transcript on c. coordinate
        """
        query = (
            f"""
            SELECT *
            FROM {self.schema}.tx_exon_aln_v
            WHERE tx_ac = '{mane_c_ac}'
            AND alt_ac = '{nc_ac}'
            AND {genomic_change_range[0]} BETWEEN alt_start_i and alt_end_i
            AND {genomic_change_range[1]} BETWEEN alt_start_i and alt_end_i
            ORDER BY ord, alt_ac;
            """
        )
        self.cursor.execute(query)
        result = self.cursor.fetchall()
        if len(result) > 1:
            logger.debug(f"Found more than one match for tx_ac {mane_c_ac} "
                         f"and alt_ac {nc_ac}")
        return result

    def get_alt_tx_data(self, ac, pos) -> Dict:
        """Get altered transcript data given transcript data.

        :param str ac: cDNA transcript
        :param tuple pos: [cDNA pos start, cDNA pos end]
        :return: Gene, Transcript accession and position change,
            Altered transcript accession and position change, Strand
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

        return dict(
            gene=gene,
            tx_ac=ac,
            tx_pos_range=tx_pos_range,
            alt_ac=nc_accession,
            alt_pos_range=alt_pos,
            pos_change=tx_pos_change,
            strand=strand
        )

    def liftover_to_38(self, alt_tx_data) -> None:
        """Liftover alt_tx_data to hg38 assembly.

        :param dict alt_tx_data: Dictionary containing gene, nc_accession,
            alt_pos, and strand
        """
        query = (
            f"""
            SELECT descr
            FROM {self.schema}._seq_anno_most_recent
            WHERE ac = '{alt_tx_data['alt_ac']}'
            """
        )
        self.cursor.execute(query)
        result = self.cursor.fetchone()
        if result and result[0]:
            query = (
                f"""
                SELECT DISTINCT alt_ac
                FROM {self.schema}.tx_exon_aln_v
                WHERE tx_ac = '{alt_tx_data['tx_ac']}'
                """
            )
            self.cursor.execute(query)
            nc_acs = self.cursor.fetchall()
            if len(nc_acs) == 1:
                logger.warning(f"UTA does not have GRCh38 assembly for "
                               f"{alt_tx_data['alt_ac'].split('.')[0]}")
                return None

            descr = result[0].split(',')
            chromosome = f"chr{descr[0].split()[-1]}"
            assembly = f"GRCh{descr[1].split('.')[0].split('GRCh')[-1]}"

            # Get most recent assembly version position
            lo = LiftOver(GRCH_TO_HG[assembly], 'hg38')
            liftover_start_i = \
                self._get_liftover(lo, chromosome,
                                   alt_tx_data['alt_pos_range'][0])
            liftover_end_i = \
                self._get_liftover(lo, chromosome,
                                   alt_tx_data['alt_pos_range'][1])
            if liftover_start_i is None or liftover_end_i is None:
                return None
            alt_tx_data['alt_pos_range'] = \
                liftover_start_i[1], liftover_end_i[1]

            # Change alt_ac to most recent
            query = (
                f"""
                SELECT *
                FROM {self.schema}.seq_anno
                WHERE ac LIKE '{alt_tx_data['alt_ac'].split('.')[0]}%'
                ORDER BY ac
                """
            )
            self.cursor.execute(query)
            nc_acs = self.cursor.fetchall()
            alt_tx_data['alt_ac'] = nc_acs[-1][3]

    @staticmethod
    def _get_liftover(lo, chromosome, pos) -> Optional[Tuple]:
        """Get new genome assembly data for a position on a chromosome.

        :param str LiftOver lo: LiftOver object to convert coordinates
        :param str chromosome: The chromosome number
        :param int pos: Position on the chromosome
        :return: [Target chromosome, target position, target strand,
            conversion_chain_score] for hg38 assembly
        """
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
