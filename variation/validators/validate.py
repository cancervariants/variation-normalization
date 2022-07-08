"""Module for Validation."""
from typing import List, Optional

from ga4gh.vrsatile.pydantic.vrs_models import RelativeCopyClass
from ga4gh.vrs.extras.translator import Translator
from gene.query import QueryHandler as GeneQueryHandler
from uta_tools.data_sources import TranscriptMappings, SeqRepoAccess, UTADatabase, \
    MANETranscript

from variation.schemas.normalize_response_schema\
    import HGVSDupDelMode as HGVSDupDelModeEnum
from variation.vrs_representation import VRSRepresentation
from variation.schemas.app_schemas import Endpoint
from variation.schemas.validation_response_schema import ValidationSummary
from variation.schemas.classification_response_schema import Classification
from variation.tokenizers import GeneSymbol
from variation.tokenizers.caches import AminoAcidCache
from .protein_substitution import ProteinSubstitution
from .polypeptide_truncation import PolypeptideTruncation
from .silent_mutation import SilentMutation
from .coding_dna_substitution import CodingDNASubstitution
from .coding_dna_silent_mutation import CodingDNASilentMutation
from .genomic_silent_mutation import GenomicSilentMutation
from .genomic_substitution import GenomicSubstitution
from .protein_delins import ProteinDelIns
from .coding_dna_delins import CodingDNADelIns
from .genomic_delins import GenomicDelIns
from .protein_deletion import ProteinDeletion
from .coding_dna_deletion import CodingDNADeletion
from .genomic_deletion import GenomicDeletion
from .protein_insertion import ProteinInsertion
from .coding_dna_insertion import CodingDNAInsertion
from .genomic_insertion import GenomicInsertion
from .genomic_uncertain_deletion import GenomicUncertainDeletion
from .genomic_duplication import GenomicDuplication
from .genomic_deletion_range import GenomicDeletionRange


class Validate:
    """The validation class."""

    def __init__(self, seqrepo_access: SeqRepoAccess,
                 transcript_mappings: TranscriptMappings,
                 gene_symbol: GeneSymbol,
                 mane_transcript: MANETranscript,
                 uta: UTADatabase, tlr: Translator,
                 amino_acid_cache: AminoAcidCache,
                 gene_normalizer: GeneQueryHandler, vrs: VRSRepresentation) -> None:
        """Initialize the validate class.

        :param SeqRepoAccess seqrepo_access: Access to SeqRepo data
        :param TranscriptMappings transcript_mappings: Access to transcript
            mappings
        :param GeneSymbol gene_symbol: Gene symbol tokenizer
        :param MANETranscript mane_transcript: Access MANE Transcript
            information
        :param UTADatabase uta: Access to UTA queries
        :param Translator tlr: Class for translating nomenclatures to and from VRS
        :param AminoAcidCache amino_acid_cache: Amino Acid codes and conversions
        :param GeneQueryHandler gene_normalizer: Access to gene-normalizer
        :param VRSRepresentation vrs: Class for representing VRS objects
        :param amino_acid_cache: Amino Acid codes and conversions
        """
        params = [
            seqrepo_access, transcript_mappings, gene_symbol,
            mane_transcript, uta, tlr, gene_normalizer, vrs
        ]
        protein_params = params[:]
        protein_params.append(amino_acid_cache)
        self.validators = [
            ProteinSubstitution(*protein_params),
            PolypeptideTruncation(*protein_params),
            SilentMutation(*protein_params),
            CodingDNASubstitution(*params),
            GenomicSubstitution(*params),
            CodingDNASilentMutation(*params),
            GenomicSilentMutation(*params),
            ProteinDelIns(*protein_params),
            CodingDNADelIns(*params),
            GenomicDelIns(*params),
            ProteinDeletion(*protein_params),
            CodingDNADeletion(*params),
            GenomicDeletion(*params),
            ProteinInsertion(*protein_params),
            CodingDNAInsertion(*params),
            GenomicInsertion(*params),
            GenomicDeletionRange(*params),
            GenomicUncertainDeletion(*params),
            GenomicDuplication(*params)
        ]

    async def perform(
            self, classifications: List[Classification],
            endpoint_name: Optional[Endpoint] = None, warnings: List = None,
            hgvs_dup_del_mode: HGVSDupDelModeEnum = HGVSDupDelModeEnum.DEFAULT,
            baseline_copies: Optional[int] = None,
            relative_copy_class: Optional[RelativeCopyClass] = None,
            do_liftover: bool = False
    ) -> ValidationSummary:
        """Validate a list of classifications.

        :param List classifications: List of classifications
        :param Optional[Endpoint] endpoint_name: Then name of the endpoint being used
        :param List warnings: List of warnings
        :param HGVSDupDelModeEnum hgvs_dup_del_mode: Must be: `default`, `absolute_cnv`,
            `relative_cnv`, `repeated_seq_expr`, `literal_seq_expr`. This parameter
            determines how to represent HGVS dup/del expressions as VRS objects.
        :param Optional[int] baseline_copies: Baseline copies number
        :param Optional[RelativeCopyClass] relative_copy_class: The relative copy class
        :param bool do_liftover: Whether or not to liftover to GRCh38 assembly
        :return: ValidationSummary containing valid and invalid results
        """
        valid_possibilities = list()
        invalid_possibilities = list()
        if not warnings:
            warnings = list()

        found_valid_result = False
        invalid_classifications = set()
        for classification in classifications:
            for validator in self.validators:
                if validator.validates_classification_type(
                        classification.classification_type):
                    results = await validator.validate(
                        classification, hgvs_dup_del_mode=hgvs_dup_del_mode,
                        endpoint_name=endpoint_name, baseline_copies=baseline_copies,
                        relative_copy_class=relative_copy_class,
                        do_liftover=do_liftover)
                    for res in results:
                        if res.is_valid:
                            found_valid_result = True
                            valid_possibilities.append(res)
                        else:
                            invalid_possibilities.append(res)
                            invalid_classifications.add(
                                classification.classification_type.value)

                if found_valid_result:
                    break

        if not found_valid_result and not warnings:
            warnings.append(f"Unable to find valid result for classifications: "
                            f"{invalid_classifications}")

        return ValidationSummary(
            valid_results=valid_possibilities,
            invalid_results=invalid_possibilities,
            warnings=warnings
        )
