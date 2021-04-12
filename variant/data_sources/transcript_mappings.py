"""The module for Transcript Mappings."""
import csv
from typing import Dict, List, Optional
from variant.schemas.validation_response_schema import LookupType
from variant import TRANSCRIPT_MAPPINGS_PATH, REFSEQ_GENE_SYMBOL_PATH


class TranscriptMappings:
    """The transcript mappings class."""

    def __init__(self, transcript_file_path=TRANSCRIPT_MAPPINGS_PATH,
                 refseq_file_path=REFSEQ_GENE_SYMBOL_PATH) -> None:
        """Initialize the transcript mappings class."""
        self.file_path = transcript_file_path
        self.protein_transcripts_versioned_for_gene_symbol: \
            Dict[str, List[str]] = {}
        self.protein_stable_id_version_to_gene_symbol: Dict[str, str] = {}
        self.protein_transcripts_for_gene_symbol: Dict[str, List[str]] = {}
        self.protein_stable_id_to_gene_symbol: Dict[str, str] = {}
        self.genomic_transcripts_for_gene_symbol: Dict[str, List[str]] = {}
        self.refseq_protein_transcripts_for_gene_symbol:\
            Dict[str, List[str]] = {}
        self.refseq_protein_transcript_to_gene_symbol: Dict[str, str] = {}
        with open(transcript_file_path) as file:
            reader = csv.DictReader(file, delimiter="\t")
            for row in reader:
                gene = row['Gene name']
                if gene:
                    versioned_protein_transcript = \
                        row['Protein stable ID version']
                    if versioned_protein_transcript:
                        self.protein_transcripts_versioned_for_gene_symbol \
                            .setdefault(gene, []) \
                            .append(versioned_protein_transcript)
                        self.protein_stable_id_version_to_gene_symbol[
                            versioned_protein_transcript] = gene
                    protein_transcript = row['Protein stable ID']
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
        with open(refseq_file_path) as file:
            reader = csv.DictReader(file, delimiter="\t")
            for row in reader:
                gene = row['Symbol']
                if gene:
                    refseq_transcript = row['Protein']
                    if refseq_transcript:
                        self.refseq_protein_transcripts_for_gene_symbol.\
                            setdefault(gene, []).\
                            append(refseq_transcript)
                        self.refseq_protein_transcript_to_gene_symbol[
                            refseq_transcript] = gene

    def protein_transcripts(self, identifier: str,
                            lookup_type: LookupType) -> Optional[List[str]]:
        """Return the versioned protein transcripts for a gene symbol."""
        protein_transcripts = list()
        if lookup_type == LookupType.GENE_SYMBOL:
            protein_transcripts += \
                self.protein_transcripts_versioned_for_gene_symbol.get(
                    identifier)
            protein_transcripts += \
                self.protein_transcripts_for_gene_symbol.get(identifier)
            protein_transcripts += \
                self.refseq_protein_transcripts_for_gene_symbol.get(identifier)
            return list(set(protein_transcripts))
        else:
            return None

    def genomic_transcripts(self, identifier: str,
                            lookup_type: LookupType) -> Optional[List[str]]:
        """Return the versioned genomic transcripts for a gene symbol."""
        if lookup_type == LookupType.GENE_SYMBOL:
            return self.genomic_transcripts_for_gene_symbol.get(identifier)
        else:
            return None

    def ensembl_gene_symbol(self, q: str, versioned=True):
        """Return the gene symbol for a Ensembl Protein."""
        if versioned:
            return self.protein_stable_id_version_to_gene_symbol.get(q)
        else:
            return self.protein_stable_id_to_gene_symbol.get(q)

    def refseq_gene_symbol(self, q: str):
        """Return the gene symbol for a Refseq Protein."""
        return self.refseq_protein_transcript_to_gene_symbol.get(q)
