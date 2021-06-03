"""The module for MANE Transcript mappings."""
from typing import Dict, Optional
from variant import REFSEQ_MANE_PATH
import pandas as pd
import logging

logger = logging.getLogger('variant')
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

    def get_gene_mane_data(self, gene_symbol) -> Optional[Dict]:
        """Return MANE Transcript data for a gene.

        :param str gene_symbol: HGNC Gene Symbol
        :return: MANE Transcript data for gene
        """
        data = self.df.loc[self.df['symbol'] == gene_symbol.upper()]

        len_data = len(data)
        if len_data == 0:
            logger.warning(f"Unable to get MANE Transcript data for gene: "
                           f"{gene_symbol}")
            return None
        elif len_data > 1:
            # Contains both MANE Select and MANE Plus
            data = data.loc[data['MANE_status'] == 'MANE Select']

        return data.to_dict('records')[0]
