"""The module for loading MANE Transcript mappings to genes."""
from typing import Dict, Optional, List
from variation import REFSEQ_MANE_PATH
import pandas as pd
import logging

logger = logging.getLogger('variation')
logger.setLevel(logging.DEBUG)


class MANETranscriptMappings:
    """The MANE Transcript mappings class."""

    def __init__(self, mane_data_path=REFSEQ_MANE_PATH) -> None:
        """Initialize the MANE Transcript mappings class.

        :param str mane_data_path: Path to RefSeq MANE summary data
        """
        self.mane_data_path = mane_data_path
        self.df = self._load_mane_transcript_data()

    def _load_mane_transcript_data(self) -> pd.core.frame.DataFrame:
        """Load RefSeq MANE data file into DataFrame.

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

        # Ordering: MANE Plus Clinical (If it exists), MANE Select
        data = data.sort_values('MANE_status')
        return data.to_dict('records')
