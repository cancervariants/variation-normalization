"""A module for caching gene symbols."""
from typing import Set, Dict, Iterable
from csv import DictReader
from variation import HGNC_GENE_SYMBOL_PATH
import re


class GeneSymbolCache:
    """The gene symbol cache class."""

    def __init__(self,
                 gene_file_path: str = HGNC_GENE_SYMBOL_PATH) \
            -> None:
        """Initialize the gene symbol cache class."""
        self.__load_caches(gene_file_path)

    def __load_caches(self, gene_file_path: str) -> None:
        """Load gene symbol cache."""
        self.gene_symbols: Set[str] = set()
        self.gene_ids: Dict[str, str] = {}
        self.gene_aliases: Dict[str, str] = {}
        self.previous_identifiers: Dict[str, str] = {}

        with open(gene_file_path, 'r') as f:
            reader = DictReader(f, delimiter='\t')
            for row in reader:
                symbol = row['symbol'].upper()
                self.gene_symbols.add(symbol)

                if row['entrez_id']:
                    curie_string = f"ncbigene:{row['entrez_id']}".upper()
                    self.gene_ids[curie_string] = symbol
                if row['ensembl_gene_id']:
                    curie_string = f"ensembl:{row['ensembl_gene_id']}".upper()
                    self.gene_ids[curie_string] = symbol
                if row['hgnc_id']:
                    # already in curie format
                    self.gene_ids[row['hgnc_id'].upper()] = symbol

                for x in self.__process_field(row['name']):
                    self.gene_aliases[x] = symbol
                for x in self.__process_field(row['alias_name']):
                    self.gene_aliases[x] = symbol
                for x in self.__process_field(row['alias_symbol']):
                    self.gene_aliases[x] = symbol

                for x in self.__process_field(row['prev_name']):
                    self.previous_identifiers[x] = symbol
                for x in self.__process_field(row['prev_symbol']):
                    self.previous_identifiers[x] = symbol

    def __process_field(self, field: str) -> Iterable[str]:
        """Process a field."""
        if field:
            return map(lambda x: x.upper(), re.sub('"', '', field).split('|'))
        else:
            return []
