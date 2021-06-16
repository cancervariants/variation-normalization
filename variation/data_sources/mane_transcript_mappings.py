"""The module for MANE Transcript mappings."""
from typing import Dict, Optional, List
from variation import REFSEQ_MANE_PATH
import pandas as pd
import logging

logger = logging.getLogger('variation')
logger.setLevel(logging.DEBUG)


class MANETranscriptMappings:
    """The MANE Transcript mappings class."""

    def __init__(self, mane_data_path=REFSEQ_MANE_PATH) -> None:
        """Initialize the MANE Transcript mappings class."""
        self.mane_data_path = mane_data_path
        self.df = self._load_mane_transcript_data()

    def _load_mane_transcript_data(self) -> pd.core.frame.DataFrame:
        """Read RefSeq MANE data file into DataFrame.

        :return: DataFrame containing RefSeq MANE Transcript data
        """
        return pd.read_csv(self.mane_data_path, delimiter='\t')

    def get_gene_mane_data(self, gene_symbol) -> Optional[List[Dict]]:
        """Return MANE Transcript data for a gene.

        :param str gene_symbol: HGNC Gene Symbol
        :return: MANE Transcript data (Transcript accessions,
            gene, and location information)
        """
        data = self.df.loc[self.df['symbol'] == gene_symbol.upper()]

        if len(data) == 0:
            logger.warning(f"Unable to get MANE Transcript data for gene: "
                           f"{gene_symbol}")
            return None

        return data.to_dict('records')
