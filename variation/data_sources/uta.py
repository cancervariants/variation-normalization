"""Module for transcript alignments and genome  assembly liftover.
Variation Normalization only supports liftover from GRCh37 assembly.
"""
from typing import Dict, Optional, List, Tuple
import psycopg2
import psycopg2.extras
from pydantic.types import StrictBool
from six.moves.urllib import parse as urlparse
import logging
from variation import UTA_DB_URL
from pyliftover import LiftOver
from os import environ
import pandas as pd
from urllib.parse import quote, unquote
import boto3


logger = logging.getLogger('variation')
logger.setLevel(logging.DEBUG)


class UTA:
    """Class for accessing UTA database."""

    def __init__(self, db_url=UTA_DB_URL, db_pwd=None) -> None:
        """Initialize UTA class.

        :param str db_url: UTA DB url
        :param str db_pwd: UTA user uta_admin's password
        """
        self.args = self._get_conn_args(('VARIATION_NORM_EB_PROD' in environ),
                                        db_pwd, db_url)
        self.conn = psycopg2.connect(**self.args)
        self.conn.autocommit = True
        self.cursor = \
            self.conn.cursor(cursor_factory=psycopg2.extras.DictCursor)
        self._create_genomic_table()
        self.liftover = LiftOver('hg19', 'hg38')

    def _update_db_url(self, db_pwd, db_url) -> Optional[str]:
        """Return new db_url containing password.

        :param str db_pwd: uta_admin's user password
        :param str db_url: PostgreSQL db url
        :return: PostgreSQL db url with password included
        """
        if 'UTA_DB_URL' in environ:
            return environ['UTA_DB_URL']

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

    def _get_conn_args(self, is_prod, db_pwd, db_url) -> Dict:
        """Return connection arguments.

        :return: A dictionary containing db credentials
        """
        if not is_prod:
            db_url = self._update_db_url(db_pwd, db_url)
            db_url = self._url_encode_password(db_url)
            url = ParseResult(urlparse.urlparse(db_url))
            self.schema = url.schema
            return self._get_args(url.hostname, url.port, url.database,
                                  url.username, unquote(url.password))
        else:
            self.schema = environ['UTA_SCHEMA']
            region = 'us-east-2'
            client = boto3.client('rds', region_name=region)
            token = client.generate_db_auth_token(
                DBHostname=environ['UTA_HOST'], Port=environ['UTA_PORT'],
                DBUsername=environ['UTA_USER'], Region=region
            )
            return self._get_args(environ['UTA_HOST'],
                                  int(environ['UTA_PORT']),
                                  environ['UTA_DATABASE'], environ['UTA_USER'],
                                  token)

    def _url_encode_password(self, db_url) -> str:
        """Update DB URL to url encode password.

        :param str db_url: URL for PostgreSQL Database
        :return: DB URL encoded
        """
        original_pwd = db_url.split('//')[-1].split('@')[0].split(':')[-1]
        return db_url.replace(original_pwd, quote(original_pwd))

    def _get_args(self, host, port, database, user, password) -> Dict:
        """Return db credentials.

        :param str host: DB Host name
        :param int port: DB port number
        :param str database: DB name
        :param str user: DB user name
        :param str password: DB password for user
        """
        return dict(
            host=host,
            port=port,
            database=database,
            user=user,
            password=password,
            application_name='variation'
        )

    def _create_genomic_table(self):
        """Create table containing genomic accession information."""
        check_table_exists = (
            f"""
            SELECT EXISTS (
               SELECT FROM information_schema.tables
               WHERE table_schema = '{self.schema}'
               AND table_name = 'genomic'
            );
            """
        )
        self.cursor.execute(check_table_exists)
        genomic_table_exists = self.cursor.fetchall()[0][0]
        if not genomic_table_exists:
            create_genomic_table = (
                f"""
                CREATE TABLE {self.schema}.genomic AS
                    SELECT t.hgnc, aes.alt_ac, aes.alt_aln_method,
                        aes.alt_strand, ae.start_i AS alt_start_i,
                        ae.end_i AS alt_end_i
                    FROM ((((({self.schema}.transcript t
                        JOIN {self.schema}.exon_set tes ON (((t.ac = tes.tx_ac)
                            AND (tes.alt_aln_method = 'transcript'::text))))
                        JOIN {self.schema}.exon_set aes ON (((t.ac = aes.tx_ac)
                            AND (aes.alt_aln_method <> 'transcript'::text))))
                        JOIN {self.schema}.exon te ON
                            ((tes.exon_set_id = te.exon_set_id)))
                        JOIN {self.schema}.exon ae ON
                            (((aes.exon_set_id = ae.exon_set_id)
                            AND (te.ord = ae.ord))))
                        LEFT JOIN {self.schema}.exon_aln ea ON
                            (((te.exon_id = ea.tx_exon_id) AND
                            (ae.exon_id = ea.alt_exon_id))));
                """
            )
            self.cursor.execute(create_genomic_table)

            indexes = [
                f"""CREATE INDEX alt_pos_index ON {self.schema}.genomic (alt_ac, alt_start_i, alt_end_i);""",  # noqa: E501
                f"""CREATE INDEX gene_alt_index ON {self.schema}.genomic (hgnc, alt_ac);"""  # noqa: E501
                f"""CREATE INDEX alt_ac_index ON {self.schema}.genomic (alt_ac);"""  # noqa: E501
            ]
            for create_index in indexes:
                self.cursor.execute(create_index)

    def get_cds_start_end(self, ac) \
            -> Optional[Dict[int, int]]:
        """Get coding start and end site

        :param str ac: Accession
        :return: [Coding start site, Coding end site]
        """
        if ac.startswith('ENS'):
            ac = ac.split('.')[0]
        query = (
            f"""
            SELECT cds_start_i, cds_end_i
            FROM {self.schema}.transcript
            WHERE ac='{ac}'
            """
        )
        self.cursor.execute(query)
        cds_start_end = self.cursor.fetchone()
        if cds_start_end and cds_start_end[0] is not None \
                and cds_start_end[1] is not None:
            return cds_start_end
        else:
            logger.warning(f"Unable to get coding start/end site for "
                           f"accession: {ac}")
            return None

    def get_newest_assembly_ac(self, ac) -> Optional[List]:
        """Find most recent genomic accession version.

        :param str ac: Genomic accession
        :return: List of most recent genomic accession version
        """
        query = (
            f"""
            SELECT ac
            FROM {self.schema}._seq_anno_most_recent
            WHERE ac LIKE '{ac.split('.')[0]}%'
            AND ((descr IS NULL) OR (descr = ''))
            ORDER BY ac
            """
        )
        self.cursor.execute(query)
        return self.cursor.fetchall()

    def validate_genomic_ac(self, ac) -> StrictBool:
        """Return whether or not genomic accession exists.

        :param str ac: Genomic accession
        :return: `True` if genomic accession exists. `False` otherwise.
        """
        query = (
            f"""
            SELECT EXISTS(
                SELECT ac
                FROM {self.schema}._seq_anno_most_recent
                WHERE ac = '{ac}'
            )
            """
        )
        self.cursor.execute(query)
        return self.cursor.fetchone()[0]

    def get_ac_descr(self, ac) -> Optional[str]:
        """Return accession description.
        Typically description exists if not GRCh38 assembly.

        :param str ac: Accession
        :return: Description containing assembly and chromosome
        """
        query = (
            f"""
            SELECT descr
            FROM {self.schema}._seq_anno_most_recent
            WHERE ac = '{ac}'
            """
        )
        self.cursor.execute(query)
        result = self.cursor.fetchone()
        if not result:
            logger.warning(f"Accession {ac} does not have a description")
            return None
        else:
            return result[0]

    def get_tx_exon_aln_v_data(self, ac, start_pos, end_pos, alt_ac=None,
                               use_tx_pos=True) -> Optional[List]:
        """Return queried data from tx_exon_aln_v table.

        :param str ac: Accession
        :param int start_pos: Start position change
        :param int end_pos: End position change
        :param str alt_ac: NC accession
        :param bool use_tx_pos: `True` if querying on transcript position.
            `False` if querying on genomic position.
        :return: tx_exon_aln_v data
        """
        if end_pos is None:
            end_pos = start_pos

        if ac.startswith('ENST'):
            temp_ac = ac.split('.')[0]
            aln_method = f"AND alt_aln_method='genebuild'"  # noqa: F541
        else:
            temp_ac = ac
            aln_method = f"AND alt_aln_method='splign'"  # noqa: F541

        if alt_ac:
            alt_ac_q = f"AND alt_ac = '{alt_ac}'"
        else:
            alt_ac_q = f"AND alt_ac LIKE 'NC_00%'"  # noqa: F541

        if use_tx_pos:
            pos_q = f"""tx_start_i AND tx_end_i"""  # noqa: F541
        else:
            pos_q = f"""alt_start_i AND alt_end_i"""  # noqa: F541

        query = (
            f"""
            SELECT hgnc, tx_ac, tx_start_i, tx_end_i, alt_ac, alt_start_i,
                alt_end_i, alt_strand, alt_aln_method, tx_exon_id, alt_exon_id
            FROM {self.schema}.tx_exon_aln_v
            WHERE tx_ac='{temp_ac}'
            {alt_ac_q}
            {aln_method}
            AND {start_pos} BETWEEN {pos_q}
            AND {end_pos} BETWEEN {pos_q}
            ORDER BY alt_ac;
            """
        )
        self.cursor.execute(query)
        results = self.cursor.fetchall()
        if not results:
            logger.warning(f"Unable to find transcript alignment for query: "
                           f"{query}")
            return None
        if alt_ac and not use_tx_pos:
            if len(results) > 1:
                logger.debug(f"Found more than one match for tx_ac {temp_ac} "
                             f"and alt_ac = {alt_ac}")

        return results

    def data_from_result(self, result) -> Optional[Dict]:
        """Return data found from result.

        :param list result: Data from tx_exon_aln_v table
        :return: Gene, strand, and position ranges for tx and alt_ac
        """
        gene = result[0]
        if result[7] == -1:
            strand = '-'
        else:
            strand = '+'
        tx_pos_range = result[2], result[3]
        alt_pos_range = result[5], result[6]
        alt_aln_method = result[8]
        tx_exon_id = result[9]
        alt_exon_id = result[10]

        if (tx_pos_range[1] - tx_pos_range[0]) != \
                (alt_pos_range[1] - alt_pos_range[0]):
            logger.warning(f"tx_pos_range {tx_pos_range} "
                           f"is not the same length as alt_pos_range "
                           f"{alt_pos_range}.")
            return None

        return dict(
            gene=gene,
            strand=strand,
            tx_pos_range=tx_pos_range,
            alt_pos_range=alt_pos_range,
            alt_aln_method=alt_aln_method,
            tx_exon_id=tx_exon_id,
            alt_exon_id=alt_exon_id,
        )

    def get_mane_c_genomic_data(self, ac, alt_ac, start_pos, end_pos):
        """Get MANE Transcript and genomic data.
        Used when going from g -> MANE c

        :param str ac: MANE Transcript accession
        :param str alt_ac: NC Accession
        :param int start_pos: Genomic start position change
        :param int end_pos: Genomic end position change
        """
        results = self.get_tx_exon_aln_v_data(
            ac, start_pos, end_pos, alt_ac=alt_ac, use_tx_pos=False
        )
        if not results:
            return None
        result = results[0]

        data = self.data_from_result(result)
        if not data:
            return None

        coding_start_site = self.get_cds_start_end(ac)
        if coding_start_site is None:
            logger.warning(f"Accession {ac} not found in UTA")
            return None

        data['tx_ac'] = result[1]
        data['alt_ac'] = result[4]
        data['coding_start_site'] = coding_start_site[0]
        data['coding_end_site'] = coding_start_site[1]
        data['alt_pos_change'] = (
            start_pos - data['alt_pos_range'][0],
            data['alt_pos_range'][1] - end_pos
        )
        data['alt_pos_change_range'] = (
            data['alt_pos_range'][0] + data['alt_pos_change'][0],
            data['alt_pos_range'][1] - data['alt_pos_change'][1]
        )
        return data

    def get_genomic_tx_data(self, ac, pos) -> Optional[Dict]:
        """Get transcript mapping to genomic data.
        Used when going from c -> g

        :param str ac: cDNA transcript
        :param tuple pos: [cDNA pos start, cDNA pos end]
        :return: Gene, Transcript accession and position change,
            Altered transcript accession and position change, Strand
        """
        results = self.get_tx_exon_aln_v_data(ac, pos[0], pos[1])
        if not results:
            return None

        result = results[-1]
        data = self.data_from_result(result)
        if not data:
            return None
        data['tx_ac'] = ac
        data['alt_ac'] = result[4]
        data['pos_change'] = (
            pos[0] - data['tx_pos_range'][0],
            data['tx_pos_range'][1] - pos[1]
        )
        data['alt_pos_change_range'] = (
            data['alt_pos_range'][0] + data['pos_change'][0],
            data['alt_pos_range'][1] - data['pos_change'][1]
        )
        return data

    def get_ac_from_gene(self, gene) -> Optional[List[str]]:
        """Return genomic accession for a gene.

        :param str gene: Gene symbol
        :return: List of genomic accessions
        """
        query = (
            f"""
            SELECT DISTINCT alt_ac
            FROM {self.schema}.genomic
            WHERE hgnc = '{gene}'
            AND alt_ac LIKE 'NC_00%'
            ORDER BY alt_ac DESC
            """
        )
        self.cursor.execute(query)
        results = self.cursor.fetchall()
        if not results:
            return []
        return [item for sublist in results for item in sublist]

    def get_gene_from_ac(self, ac, start_pos, end_pos) -> Optional[str]:
        """Get transcripts from NC accession and positions.

        :param str ac: NC Accession
        :param int start_pos: Start position change
        :param int end_pos: End position change
        :return: HGNC gene symbol
        """
        if end_pos is None:
            end_pos = start_pos
        query = (
            f"""
            SELECT DISTINCT hgnc
            FROM {self.schema}.genomic
            WHERE alt_ac = '{ac}'
            AND {start_pos} BETWEEN alt_start_i AND alt_end_i
            AND {end_pos} BETWEEN alt_start_i AND alt_end_i
            """
        )
        self.cursor.execute(query)
        results = self.cursor.fetchall()
        if not results:
            logger.warning(f"Unable to find gene between {start_pos} and"
                           f" {end_pos} on {ac}")
            return None
        else:
            if len(results) > 1:
                logger.info(f"Found more than one gene between "
                            f"{start_pos} and {end_pos} on {ac}")

        return results

    def get_transcripts_from_gene(self, gene, start_pos, end_pos)\
            -> pd.core.frame.DataFrame:
        """Get transcripts on c coordinate associated to a gene.

        :param str gene: Gene symbol
        :param int start_pos: Start position change on c. coordinate
        :param int end_pos: End position change on c. coordinate
        :return: Data Frame containing transcripts associated with a gene.
            Transcripts are ordered by most recent NC accession, then by
            descending transcript length.
        """
        query = (
            f"""
            SELECT AA.pro_ac, AA.tx_ac, ALIGN.alt_ac, T.cds_start_i
            FROM {self.schema}.associated_accessions as AA
            JOIN {self.schema}.transcript as T ON T.ac = AA.tx_ac
            JOIN {self.schema}.tx_exon_aln_v as ALIGN ON T.ac = ALIGN.tx_ac
            WHERE T.hgnc = '{gene}'
            AND ALIGN.alt_ac LIKE 'NC_00%'
            AND ALIGN.alt_aln_method = 'splign'
            AND {start_pos} + T.cds_start_i
                BETWEEN ALIGN.tx_start_i AND ALIGN.tx_end_i
            AND {end_pos} + T.cds_start_i
                BETWEEN ALIGN.tx_start_i AND ALIGN.tx_end_i
            ORDER BY ALIGN.alt_ac, ALIGN.tx_end_i - ALIGN.tx_start_i DESC;
            """
        )
        return pd.read_sql(query, self.conn)

    def get_chr_assembly(self, ac) -> Optional[Tuple[str, str]]:
        """Get chromosome and assembly for NC accession.

        :param str ac: NC accession
        :return: Chromosome and Assembly accession is on
        """
        descr = self.get_ac_descr(ac)
        if not descr:
            # Already GRCh38 Assembly
            return None
        descr = descr.split(',')
        chromosome = f"chr{descr[0].split()[-1]}"
        assembly = f"GRCh{descr[1].split('.')[0].split('GRCh')[-1]}"

        if assembly not in ['GRCh37', 'GRCh38']:
            logger.warning(f"Assembly not supported: {assembly}. "
                           f"Only GRCh37 and GRCh38 are supported.")
            return None

        return chromosome, assembly

    def liftover_to_38(self, genomic_tx_data) -> None:
        """Liftover genomic_tx_data to hg38 assembly.

        :param dict genomic_tx_data: Dictionary containing gene, nc_accession,
            alt_pos, and strand
        """
        descr = self.get_chr_assembly(genomic_tx_data['alt_ac'])
        if descr is None:
            return None
        chromosome, assembly = descr

        query = (
            f"""
            SELECT DISTINCT alt_ac
            FROM {self.schema}.tx_exon_aln_v
            WHERE tx_ac = '{genomic_tx_data['tx_ac']}'
            """
        )
        self.cursor.execute(query)
        nc_acs = self.cursor.fetchall()
        if len(nc_acs) == 1:
            logger.warning(f"UTA does not have GRCh38 assembly for "
                           f"{genomic_tx_data['alt_ac'].split('.')[0]}")
            return None

        # Get most recent assembly version position
        # Liftover range
        self._set_liftover(
            genomic_tx_data, 'alt_pos_range', chromosome
        )

        # Liftover changes range
        self._set_liftover(
            genomic_tx_data, 'alt_pos_change_range', chromosome
        )

        # Change alt_ac to most recent
        query = (
            f"""
            SELECT alt_ac
            FROM {self.schema}.genomic
            WHERE alt_ac LIKE '{genomic_tx_data['alt_ac'].split('.')[0]}%'
            ORDER BY alt_ac;
            """
        )
        self.cursor.execute(query)
        nc_acs = self.cursor.fetchall()
        genomic_tx_data['alt_ac'] = nc_acs[-1][0]

    def get_liftover(self, chromosome, pos) -> Optional[Tuple]:
        """Get new genome assembly data for a position on a chromosome.

        :param LiftOver lo: LiftOver object to convert coordinates
        :param str chromosome: The chromosome number
        :param int pos: Position on the chromosome
        :return: [Target chromosome, target position, target strand,
            conversion_chain_score] for hg38 assembly
        """
        liftover = self.liftover.convert_coordinate(chromosome, pos)
        if liftover is None or len(liftover) == 0:
            logger.warning(f"{pos} does not exist on {chromosome}")
            return None
        else:
            return liftover[0]

    def _set_liftover(self, genomic_tx_data, key, chromosome) -> None:
        """Update genomic_tx_data to have hg38 coordinates.

        :param dict genomic_tx_data:
        :param str key: Key to access coordinate positions
        :param LiftOver lo: LiftOver object to convert coordinates
        :param str chromosome: Chromosome
        """
        liftover_start_i = self.get_liftover(chromosome,
                                             genomic_tx_data[key][0])
        if liftover_start_i is None:
            logger.warning(f"Unable to liftover position "
                           f"{genomic_tx_data[key][0]} on {chromosome}")
            return None

        liftover_end_i = self.get_liftover(chromosome,
                                           genomic_tx_data[key][1])
        if liftover_end_i is None:
            logger.warning(f"Unable to liftover position "
                           f"{genomic_tx_data[key][1]} on {chromosome}")
            return None

        genomic_tx_data[key] = liftover_start_i[1], liftover_end_i[1]

    def p_to_c_ac(self, p_ac) -> List[str]:
        """Return c. accession from p. accession.

        :param str p_ac: Protein accession
        :return: List of rows containing c. accessions that are associated
            with the given p. accession.
        """
        query = (
            f"""
            SELECT tx_ac
            FROM {self.schema}.associated_accessions
            WHERE pro_ac = '{p_ac}'
            ORDER BY tx_ac;
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
