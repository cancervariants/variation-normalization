"""The module for Transcript Mappings."""
import csv
from typing import Dict, List, Optional
from variant.schemas.validation_response_schema import LookupType
from variant import TRANSCRIPT_MAPPINGS_PATH


class TranscriptMappings:
    """The transcript mappings class."""

    def __init__(self, transcript_file_path=TRANSCRIPT_MAPPINGS_PATH) -> None:
        """Initialize the transcript mappings class."""
        self.file_path = transcript_file_path
        self.protein_transcripts_for_gene_symbol: Dict[str, List[str]] = {}
        self.genomic_transcripts_for_gene_symbol: Dict[str, List[str]] = {}
        self.protein_stable_id_to_gene_symbol: Dict[str, str] = {}
        with open(transcript_file_path) as file:
            reader = csv.DictReader(file, delimiter="\t")
            for row in reader:
                gene = row['Gene name']
                if gene:
                    protein_transcript = row['Protein stable ID version']
                    if protein_transcript:
                        self.protein_transcripts_for_gene_symbol \
                            .setdefault(gene, []) \
                            .append(protein_transcript)
                        self.protein_stable_id_to_gene_symbol[
                            protein_transcript] = gene
                    genomic_transcript = row['Transcript stable ID version']
                    if genomic_transcript:
                        self.genomic_transcripts_for_gene_symbol \
                            .setdefault(gene, []) \
                            .append(genomic_transcript)

    def protein_transcripts(self, identifier: str,
                            lookup_type: LookupType) -> Optional[List[str]]:
        """Return the protein transcripts for a gene symbol."""
        if lookup_type == LookupType.GENE_SYMBOL:
            return self.protein_transcripts_for_gene_symbol.get(identifier)
        else:
            return None

    def genomic_transcripts(self, identifier: str,
                            lookup_type: LookupType) -> Optional[List[str]]:
        """Return the genomic transcripts for a gene symbol."""
        if lookup_type == LookupType.GENE_SYMBOL:
            return self.genomic_transcripts_for_gene_symbol.get(identifier)
        else:
            return None

    def refseq_gene_symbol(self, refseq: str):
        """Return the gene symbol for a given reference sequence."""
        return self.protein_stable_id_to_gene_symbol.get(refseq)
