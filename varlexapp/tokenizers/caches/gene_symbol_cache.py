from csv import DictReader

import re

class GeneSymbolCache:
    def __init__(self, gene_file_path):
        self.__load_caches(gene_file_path)

    def __load_caches(self, gene_file_path):
        self.gene_symbols = set()
        self.gene_ids = {}
        self.gene_aliases = {}

        self.previous_identifiers = {}

        with open(gene_file_path, 'r') as f:
            reader = DictReader(f, delimiter='\t')
            for row in reader:
                symbol = row['symbol'].upper()
                self.gene_symbols.add(symbol)

                if row['entrez_id']:
                    self.gene_ids[row['entrez_id'].upper()] = symbol
                if row['ensembl_gene_id']:
                    self.gene_ids[row['ensembl_gene_id'].upper()] = symbol
                if row['hgnc_id']:
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

    def __process_field(self, field):
        if field:
            return map(lambda x: x.upper(), re.sub('"', '', field).split('|'))
        else:
           return []

