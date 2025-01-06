"""Module for Genomic Substitution Translation."""

from cool_seq_tool.schemas import AnnotationLayer, CoordinateType, Strand
from ga4gh.vrs import models

from variation.schemas.app_schemas import Endpoint
from variation.schemas.classification_response_schema import (
    CdnaSubstitutionClassification,
    ClassificationType,
    GenomicSubstitutionClassification,
)
from variation.schemas.normalize_response_schema import HGVSDupDelModeOption
from variation.schemas.token_response_schema import AltType
from variation.schemas.translation_response_schema import (
    TranslationResult,
    VrsSeqLocAcStatus,
)
from variation.schemas.validation_response_schema import ValidationResult
from variation.translators.translator import Translator


class GenomicSubstitution(Translator):
    """The Genomic Substitution Translator class."""

    def can_translate(self, classification_type: ClassificationType) -> bool:
        """Determine if it's possible to translate a classification.

        :param classification_type: Classification type found
        :return: `True` if `classification_type` matches translator's classification
            type. Otherwise, `False`
        """
        return classification_type == ClassificationType.GENOMIC_SUBSTITUTION

    async def translate(
        self,
        validation_result: ValidationResult,
        warnings: list[str],
        endpoint_name: Endpoint | None = None,
        hgvs_dup_del_mode: HGVSDupDelModeOption = HGVSDupDelModeOption.DEFAULT,
        baseline_copies: int | None = None,
        copy_change: models.CopyChange | None = None,
        do_liftover: bool = False,
    ) -> TranslationResult | None:
        """Translate validation result to VRS representation

        :param validation_result: Validation result for a classification
        :param endpoint_name: Name of endpoint that is being used
        :param hgvs_dup_del_mode: Mode to use for interpreting HGVS duplications and
            deletions
        :param baseline_copies: The baseline copies for a copy number count variation
        :param copy_change: The change for a copy number change variation
        :param do_liftover: Whether or not to liftover to GRCh38 assembly
        :return: Translation result if translation was successful. If translation was
            not successful, `None`
        """
        errors = []

        # First will translate valid result to VRS Allele
        classification: GenomicSubstitutionClassification = (
            validation_result.classification
        )
        vrs_allele = None
        vrs_seq_loc_ac = None
        vrs_seq_loc_ac_status = VrsSeqLocAcStatus.NA

        if endpoint_name == Endpoint.NORMALIZE:
            gene = (
                classification.gene_token.token if classification.gene_token else None
            )
            mane = await self.mane_transcript.get_mane_transcript(
                validation_result.accession,
                classification.pos,
                classification.pos,
                AnnotationLayer.GENOMIC,
                try_longest_compatible=True,
                coordinate_type=CoordinateType.RESIDUE,
                gene=gene,
            )

            if mane:
                vrs_seq_loc_ac_status = mane.status

                if gene:
                    if mane.strand == Strand.NEGATIVE:
                        ref_rev = classification.ref[::-1]
                        alt_rev = classification.alt[::-1]

                        complements = {"A": "T", "T": "A", "C": "G", "G": "C"}

                        ref = ""
                        alt = ""

                        for nt in ref_rev:
                            ref += complements[nt]
                        for nt in alt_rev:
                            alt += complements[nt]
                    else:
                        ref = classification.ref
                        alt = classification.alt

                    classification = CdnaSubstitutionClassification(
                        matching_tokens=classification.matching_tokens,
                        nomenclature=classification.nomenclature,
                        gene_token=classification.gene_token,
                        pos=mane.pos[0] + 1,  # 1-based for classification
                        ref=ref,
                        alt=alt,
                    )
                    vrs_seq_loc_ac = mane.refseq
                    coord_type = AnnotationLayer.CDNA
                    validation_result.classification = classification
                    extensions = None
                else:
                    vrs_seq_loc_ac = mane.ac
                    coord_type = AnnotationLayer.GENOMIC
                    extensions = self._mane_gene_extensions(mane.mane_genes)

                vrs_allele = self.vrs.to_vrs_allele(
                    vrs_seq_loc_ac,
                    mane.pos[0],
                    mane.pos[1],
                    coord_type,
                    AltType.SUBSTITUTION,
                    errors,
                    alt=classification.alt,
                    cds_start=mane.coding_start_site if gene else None,
                    coordinate_type=CoordinateType.INTER_RESIDUE,
                    extensions=extensions,
                )
        else:
            vrs_seq_loc_ac = validation_result.accession
            vrs_allele = self.vrs.to_vrs_allele(
                vrs_seq_loc_ac,
                classification.pos,
                classification.pos,
                AnnotationLayer.GENOMIC,
                AltType.SUBSTITUTION,
                errors,
                alt=classification.alt,
                coordinate_type=CoordinateType.RESIDUE,
            )

        if vrs_allele and vrs_seq_loc_ac:
            return TranslationResult(
                vrs_variation=vrs_allele,
                vrs_seq_loc_ac=vrs_seq_loc_ac,
                vrs_seq_loc_ac_status=vrs_seq_loc_ac_status,
                og_ac=validation_result.accession,
                validation_result=validation_result,
            )

        return None
