import csv

from typing import Dict, List, Optional

from ..models import LookupType

class TranscriptMappings:
    def __init__(self, transcript_file_path: str) -> None:
        self.file_path = transcript_file_path
        self.protein_transcripts_for_gene_symbol: Dict[str, List[str]] = {}
        self.genomic_transcripts_for_gene_symbol: Dict[str, List[str]] = {}
        with open(transcript_file_path) as file:
            reader = csv.DictReader(file, delimiter="\t")
            for row in reader:
                gene = row['Gene name']
                if gene:
                    protein_transcript = row['Protein stable ID']
                    if protein_transcript:
                        self.protein_transcripts_for_gene_symbol \
                            .setdefault(gene, []) \
                            .append(protein_transcript)
                    genomic_transcript = row['Transcript stable ID']
                    if genomic_transcript:
                        self.genomic_transcripts_for_gene_symbol \
                            .setdefault(gene, []) \
                            .append(genomic_transcript)


    def protein_transcripts(self, identifier: str, lookup_type: LookupType) -> Optional[List[str]]:
        if lookup_type == LookupType.GENE_SYMBOL:
            return self.protein_transcripts_for_gene_symbol.get(identifier)
        else:
            return None

    def genomic_transcripts(self, identifier: str, lookup_type: LookupType) -> Optional[List[str]]:
        if lookup_type == LookupType.GENE_SYMBOL:
            return self.genomic_transcripts_for_gene_symbol.get(identifier)
        else:
            return None


