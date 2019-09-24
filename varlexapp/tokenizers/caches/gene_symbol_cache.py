from csv import DictReader

class GeneSymbolCache:
    def __init__(self, gene_file_path):
        self.__gene_symbol_cache = self.__load_gene_symbols(gene_file_path)

    def __contains__(self, item):
        return item in self.__gene_symbol_cache

    def __load_gene_symbols(self, gene_file_path):
        gene_symbols = set()
        with open(gene_file_path, 'r') as f:
            reader = DictReader(f, delimiter='\t')
            for row in reader:
                gene_symbols.add(row['symbol'])
        return gene_symbols
