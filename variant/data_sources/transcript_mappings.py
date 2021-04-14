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
        # ENSP <-> Gene Symbol
        self.ensembl_protein_version_for_gene_symbol: Dict[str, List[str]] = {}
        self.ensembl_protein_version_to_gene_symbol: Dict[str, str] = {}
        self.ensembl_protein_for_gene_symbol: Dict[str, List[str]] = {}
        self.ensembl_protein_to_gene_symbol: Dict[str, str] = {}

        # Gene Symbol -> ENST
        self.ensembl_transcript_version_for_gene_symbol: \
            Dict[str, List[str]] = {}

        # NP_ <-> Gene Symbol
        self.refseq_protein_for_gene_symbol: Dict[str, List[str]] = {}
        self.refseq_protein_to_gene_symbol: Dict[str, str] = {}

        # NM_ <-> Gene Symbol
        self.refseq_rna_for_gene_symbol: Dict[str, List[str]] = {}
        self.refseq_rna_to_gene_symbol: Dict[str, str] = {}

        self._load_transcript_mappings_data(transcript_file_path)
        self._load_refseq_gene_symbol_data(refseq_file_path)

    def _load_transcript_mappings_data(self, transcript_file_path):
        """Load transcript mappings file to dictionaries."""
        with open(transcript_file_path) as file:
            reader = csv.DictReader(file, delimiter="\t")
            for row in reader:
                gene = row['Gene name']
                if gene:
                    versioned_protein_transcript = \
                        row['Protein stable ID version']
                    if versioned_protein_transcript:
                        self.ensembl_protein_version_for_gene_symbol \
                            .setdefault(gene, []) \
                            .append(versioned_protein_transcript)
                        self.ensembl_protein_version_to_gene_symbol[
                            versioned_protein_transcript] = gene
                    protein_transcript = row['Protein stable ID']
                    if protein_transcript:
                        self.ensembl_protein_for_gene_symbol \
                            .setdefault(gene, []) \
                            .append(protein_transcript)
                        self.ensembl_protein_to_gene_symbol[
                            protein_transcript] = gene
                    versioned_genomic_transcript = \
                        row['Transcript stable ID version']
                    if versioned_genomic_transcript:
                        self.ensembl_transcript_version_for_gene_symbol \
                            .setdefault(gene, []) \
                            .append(versioned_genomic_transcript)

    def _load_refseq_gene_symbol_data(self, refseq_file_path):
        """Load data from RefSeq Gene Symbol file to dictionaries."""
        with open(refseq_file_path) as file:
            reader = csv.DictReader(file, delimiter="\t")
            for row in reader:
                gene = row['Symbol']
                if gene:
                    refseq_transcript = row['Protein']
                    if refseq_transcript:
                        self.refseq_protein_for_gene_symbol.\
                            setdefault(gene, []).\
                            append(refseq_transcript)
                        self.refseq_protein_to_gene_symbol[
                            refseq_transcript] = gene
                    rna_trasncript = row['RNA']
                    if rna_trasncript:
                        self.refseq_rna_for_gene_symbol.\
                            setdefault(gene, []).\
                            append(rna_trasncript)
                        self.refseq_rna_to_gene_symbol[
                            rna_trasncript] = gene

    def protein_transcripts(self, identifier: str,
                            lookup_type: LookupType) -> Optional[List[str]]:
        """Return a list of protein transcripts for a gene symbol."""
        protein_transcripts = list()
        if lookup_type == LookupType.GENE_SYMBOL:
            protein_transcripts += \
                self.ensembl_protein_version_for_gene_symbol.get(
                    identifier)
            protein_transcripts += \
                self.ensembl_protein_for_gene_symbol.get(identifier)
            protein_transcripts += \
                self.refseq_protein_for_gene_symbol.get(identifier)
            return list(set(protein_transcripts))
        else:
            return None

    def genomic_transcripts(self, identifier: str,
                            lookup_type: LookupType) -> Optional[List[str]]:
        """Return genomic transcripts for a gene symbol."""
        if lookup_type == LookupType.GENE_SYMBOL:
            return self.ensembl_transcript_version_for_gene_symbol.get(
                identifier)
        else:
            return None

    def get_gene_symbol_from_ensembl_protein(self, q: str, versioned=True):
        """Return the gene symbol for a Ensembl Protein."""
        if versioned:
            return self.ensembl_protein_version_to_gene_symbol.get(q)
        else:
            return self.ensembl_protein_to_gene_symbol.get(q)

    def get_gene_symbol_from_refeq_protein(self, q: str):
        """Return the gene symbol for a Refseq Protein."""
        return self.refseq_protein_to_gene_symbol.get(q)

    def get_gene_symbol_from_refseq_rna(self, q: str) -> Optional[List[str]]:
        """Return gene symbol for a Refseq RNA Transcript."""
        return self.refseq_rna_to_gene_symbol.get(q)
