"""Module for retrieving MANE transcript."""
from typing import Optional
import hgvs.parser
import hgvs.dataproviders.uta
import hgvs.assemblymapper
import logging

logger = logging.getLogger('variant')
logger.setLevel(logging.DEBUG)


class MANETranscript:
    """Class for retrieving MANE transcripts."""

    def __init__(self):
        """Initialize the MANETranscript class."""
        self.hgvs_parser = hgvs.parser.Parser()
        self.hgvs_data_providers = hgvs.dataproviders.uta.connect()
        self.assembly_mapper = hgvs.assemblymapper.AssemblyMapper(
            self.hgvs_data_providers,
            assembly_name='GRCh38',
            alt_aln_method='splign',
            replace_reference=True
        )

    def protein_to_transcript(self, token) -> Optional[str]:
        """Convert protein to annotation to transcript annotation.
        p. -> c.

        :param Token token: Classification token
        :return:
        """
        pos = token.position * 3 - 1  # noqa: F841
