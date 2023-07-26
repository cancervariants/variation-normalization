"""Module for testing validators"""
import pytest

from variation.validators import (
    ProteinSubstitution, CdnaSubstitution, GenomicSubstitution, ProteinStopGain,
    ProteinReferenceAgree, CdnaReferenceAgree, GenomicReferenceAgree, ProteinDelIns,
    CdnaDelIns, GenomicDelIns, ProteinDeletion, CdnaDeletion, GenomicDeletion,
    GenomicDeletionAmbiguous, ProteinInsertion, CdnaInsertion, GenomicInsertion,
    GenomicDuplication, GenomicDuplicationAmbiguous, Amplification
)


async def validator_checks(
    test_tokenizer, test_classifier, query, params, validator_instance,
    expect_valid=True
):
    assert isinstance(query, str)
    tokens = test_tokenizer.perform(query, [])
    classification = test_classifier.perform(tokens)

    try:
        validation_results = await validator_instance(*params).validate(
            classification
        )
    except Exception as e:
        raise Exception(f"{e}: {query}")
    else:
        validator_instance
        is_valid = False
        for vr in validation_results:
            if vr.is_valid:
                is_valid = True
                break

        assert is_valid if expect_valid else not is_valid, query


@pytest.mark.asyncio
async def test_protein_substitution(test_tokenizer, test_classifier, val_params):
    """Test that protein substitution validator works correctly"""
    validator_instance = ProteinSubstitution
    for q in [
        "BRAF V600E",
        "NP_004324.2:p.Val600Glu",
        "NP_005219.2:p.Thr790Met",
        "EGFR Leu858Arg"
    ]:
        await validator_checks(
            test_tokenizer, test_classifier, q, val_params, validator_instance
        )

    for q in [
        "NP_004324.2:p.Val600000000000Glu",
        "NP_004324.2:p.Glu600Val",
        "NP_005148.2:p.Leu2733Gln",
        "NP_000000000542.1:p.Val66Gly",
        "BRAF V9999999999999999999999999999999E"
    ]:
        await validator_checks(
            test_tokenizer, test_classifier, q, val_params, validator_instance,
            expect_valid=False
        )


@pytest.mark.asyncio
async def test_cdna_substitution(test_tokenizer, test_classifier, val_params):
    """Test that cdna substitution validator works correctly"""
    validator_instance = CdnaSubstitution
    for q in [
        "NM_004333.4:c.1799T>A",
        "ENST00000288602.10:c.1799T>A",
        "BRAF (c.1799T>A)",
        "BRAF c.1799T>A",
        "BRAF V600E c.1799T>A"
    ]:
        await validator_checks(
            test_tokenizer, test_classifier, q, val_params, validator_instance
        )

    for q in [
        "BRAF c.18000000000000T>A",
        "NM_004333.4:c.17699T>A"
    ]:
        await validator_checks(
            test_tokenizer, test_classifier, q, val_params, validator_instance,
            expect_valid=False
        )


@pytest.mark.asyncio
async def test_genomic_substitution(test_tokenizer, test_classifier, val_params):
    """Test that genomic substitution validator works correctly"""
    validator_instance = GenomicSubstitution
    for q in [
        "NC_000007.13:g.140453136A>T",
        "NC_000007.13:g.55259515T>G",
        "7-140453136-A-T",
        "7-55259515-T-G"
    ]:
        await validator_checks(
            test_tokenizer, test_classifier, q, val_params, validator_instance
        )

    for q in [
        "NC_000007.13:g.1436A>T",
        "NC_000007.13:g.4T>A",
        "7-140453136-G-T"
    ]:
        await validator_checks(
            test_tokenizer, test_classifier, q, val_params, validator_instance,
            expect_valid=False
        )


@pytest.mark.asyncio
async def test_protein_stop_gain(test_tokenizer, test_classifier, val_params):
    """Test that protein stop gain validator works correctly"""
    validator_instance = ProteinStopGain
    for q in [
        "NP_060842.3:p.Tyr365Ter",
        "NP_000542.1:p.Tyr185Ter",
        "NP_000542.1:p.Tyr185*"
    ]:
        await validator_checks(
            test_tokenizer, test_classifier, q, val_params, validator_instance
        )

    await validator_checks(
        test_tokenizer, test_classifier, "NP_060842.3:p.Tyr3650000000000Ter",
        val_params, validator_instance, expect_valid=False
    )


@pytest.mark.asyncio
async def test_protein_reference_agree(test_tokenizer, test_classifier, val_params):
    """Test that protein reference agree validator works correctly"""
    validator_instance = ProteinReferenceAgree
    for q in [
        "NP_000542.1:p.Pro61=",
        "NP_000918.2:p.Ile1145="
    ]:
        await validator_checks(
            test_tokenizer, test_classifier, q, val_params, validator_instance
        )

    await validator_checks(
        test_tokenizer, test_classifier, "NP_000542.1:p.Pro62=", val_params,
        validator_instance, expect_valid=False
    )


@pytest.mark.asyncio
async def test_cdna_reference_agree(test_tokenizer, test_classifier, val_params):
    """Test that cdna reference agree validator works correctly"""
    validator_instance = CdnaReferenceAgree
    for q in [
        "NM_004006.2:c.123=",
        "NM_004333.4:c.1799=",
        "ENST00000288602.11:c.1799=",
        "BRAF c.1799=",
        "BRAF V600E c.1799="
    ]:
        await validator_checks(
            test_tokenizer, test_classifier, q, val_params, validator_instance
        )

    for q in [
        "NM_004006.2:c.13994=",
        "BRAF c.18000000000000=",
        "NM_000412.5:c.1930="  # pos out of index
    ]:
        await validator_checks(
            test_tokenizer, test_classifier, q, val_params, validator_instance,
            expect_valid=False
        )


@pytest.mark.asyncio
async def test_genomic_reference_agree(test_tokenizer, test_classifier, val_params):
    """Test that genomic reference agree validator works correctly"""
    validator_instance = GenomicReferenceAgree
    for q in [
        "NC_000007.13:g.140453136=",
        "NC_000007.13:g.55259515=",
        "7-140453136-A-A",
        "7-55259515-T-T"
    ]:
        await validator_checks(
            test_tokenizer, test_classifier, q, val_params, validator_instance
        )

    for q in [
        "NC_000007.13:g.159138664=",
        "7-140453136-C-C"
    ]:
        await validator_checks(
            test_tokenizer, test_classifier, q, val_params, validator_instance,
            expect_valid=False
        )


@pytest.mark.asyncio
async def test_protein_delins(test_tokenizer, test_classifier, val_params):
    """Test that protein delins validator works correctly"""
    validator_instance = ProteinDelIns
    for q in [
        "NP_001333827.1:p.Leu747_Thr751delinsPro",
        "NP_000542.1:p.Gln96_Pro97delinsHis",
        "NP_005219.2:p.Glu746_Thr751delinsValAla",
        "ERBB2 G776delinsVC",
        "KIT P577_W582delinsPYD"
    ]:
        await validator_checks(
            test_tokenizer, test_classifier, q, val_params, validator_instance
        )

    for q in [
        "NP_001333827.1:p.Cys747_Thr751delinsPro",
        "NP_001333827.1:p.Leu747_Pro751delinsPro"
    ]:
        await validator_checks(
            test_tokenizer, test_classifier, q, val_params, validator_instance,
            expect_valid=False
        )


@pytest.mark.asyncio
async def test_cdna_delins(test_tokenizer, test_classifier, val_params):
    """Test that cdna delins validator works correctly"""
    validator_instance = CdnaDelIns
    for q in [
        "NM_001289937.1:c.2326_2327delinsCT",
        "NM_000551.3:c.615delinsAA",
        "ENST00000440973.5:c.1607_1608delinsAG",
        "ENST00000318560.5:c.1423_1424delinsGT",
        "ENST00000256474.2:c.364_365delinsAT",
        "NM_000551.3:c.615delinsAA"
    ]:
        await validator_checks(
            test_tokenizer, test_classifier, q, val_params, validator_instance
        )

    for q in [
        "NM_005228:c.2237_2253delinsTTGCT",
        "ENST00000277541.6:c.7330479587395delinsACA",
        "NM_000551.3:c.4561delinsAA",
        "NM_000551.3:c.4561_4562delinsAA",
        "NM_000551.3:c.4560_4561delinsAA",
        "NM_001289937.1:c.2327_2326delinsCT",
        "NM_000551.3:c.4559delinsAA"  # pos out of index
    ]:
        await validator_checks(
            test_tokenizer, test_classifier, q, val_params, validator_instance,
            expect_valid=False
        )


@pytest.mark.asyncio
async def test_genomic_delins(test_tokenizer, test_classifier, val_params):
    """Test that genomic delins validator works correctly"""
    validator_instance = GenomicDelIns
    for q in [
        "NC_000007.13:g.140453135_140453136delinsAT",
        "NC_000007.13:g.159138662delinsAT",
        "NC_000023.11:g.32386323delinsGA",
        "NC_000003.12:g.10149938delinsAA"
    ]:
        await validator_checks(
            test_tokenizer, test_classifier, q, val_params, validator_instance
        )

    for q in [
        "NC_000023.21:g.32386323delinsGA",
        "NC_000007.13:g.159138664delinsAT",
        "NC_000007.13:g.159138663_159138664delinsAT",
        "NC_000023.11:g.3238646549879323delinsGA",
        "NC_000007.13:g.140453136_140453134delinsAT"
    ]:
        await validator_checks(
            test_tokenizer, test_classifier, q, val_params, validator_instance,
            expect_valid=False
        )


@pytest.mark.asyncio
async def test_protein_deletion(test_tokenizer, test_classifier, val_params):
    """Test that protein deletion validator works correctly"""
    validator_instance = ProteinDeletion
    for q in [
        "NP_003997.1:p.Lys23_Val25del",
        "NP_000542.1:p.Glu186del",
        "NP_000542.1:p.Arg82_Val84del",
        "ENSP00000256474.2:p.Phe76del",
        "KIT D419del",
        "KIT E554_V559del",
        "EGFR L747_T751del",
        "EGFR L747_T751delLREAT"
    ]:
        await validator_checks(
            test_tokenizer, test_classifier, q, val_params, validator_instance
        )

    for q in [
        "EGFR L747_T751delLREATS",
        "KIT V419del"
    ]:
        await validator_checks(
            test_tokenizer, test_classifier, q, val_params, validator_instance,
            expect_valid=False
        )


@pytest.mark.asyncio
async def test_cdna_deletion(test_tokenizer, test_classifier, val_params):
    """Test that cdna deletion validator works correctly"""
    validator_instance = CdnaDeletion
    for q in [
        "ENST00000269571.9:c.2263_2277del",
        "NM_004448.3:c.2263_2277delTTGAGGGAAAACACA",
        "ERBB2 c.2263_2277delTTGAGGGAAAACACA",
        "NM_004448.3:c.2263_2277del",
        "NM_000535.6:c.2117delA",
        "NM_000535.6:c.2117del"
    ]:
        await validator_checks(
            test_tokenizer, test_classifier, q, val_params, validator_instance
        )

    for q in [
        "NM_000535.6:c.21174568delT",
        "NM_000535.6:c.21145457delA",
        "ENST00000269571.9:c.2277_2263del"
    ]:
        await validator_checks(
            test_tokenizer, test_classifier, q, val_params, validator_instance,
            expect_valid=False
        )


@pytest.mark.asyncio
async def test_genomic_deletion(test_tokenizer, test_classifier, val_params):
    """Test that genomic deletion validator works correctly"""
    validator_instance = GenomicDeletion
    for q in [
        "NC_000003.11:g.10188279_10188297del",
        "NC_000003.11:g.10191486_10191487delAG",
        "NC_000003.12:g.10146527_10146528del",
        "NC_000003.11:g.10191495delT",
        "VHL g.10188279_10188297del",
        "16-2138199-GTGAG-G"
    ]:
        await validator_checks(
            test_tokenizer, test_classifier, q, val_params, validator_instance
        )

    for q in [
        "NC_000003.11:g.10191454654654654495delT",
        "NC_000003.11:g.10188297_10188279del"

    ]:
        await validator_checks(
            test_tokenizer, test_classifier, q, val_params, validator_instance,
            expect_valid=False
        )


@pytest.mark.asyncio
async def test_genomic_deletion_ambiguous(test_tokenizer, test_classifier, val_params):
    """Test that genomic deletion validator works correctly"""
    validator_instance = GenomicDeletionAmbiguous
    for q in [
        "NC_000023.11:g.(?_155980375)_(156013167_?)del",
        "NC_000002.12:g.(?_110104900)_(110207160_?)del",
        "NC_000024.10:g.(?_14076802)_(57165209_?)del"
    ]:
        await validator_checks(
            test_tokenizer, test_classifier, q, val_params, validator_instance
        )

    for q in [
        "NC_000023.11:g.(?_156013167)_(155980375_?)del",
        "NC_000024.10:g.(14076805_14076804)_(14076803_14076802)del"

    ]:
        await validator_checks(
            test_tokenizer, test_classifier, q, val_params, validator_instance,
            expect_valid=False
        )


@pytest.mark.asyncio
async def test_protein_insertion(test_tokenizer, test_classifier, val_params):
    """Test that protein insertion validator works correctly"""
    validator_instance = ProteinInsertion
    for q in [
        "NP_005219.2:p.Asp770_Asn771insGlyLeu",
        "BRAF T599_V600insV"
    ]:
        await validator_checks(
            test_tokenizer, test_classifier, q, val_params, validator_instance
        )

    for q in [
        "NP_005219.2:p.Asn770_Gly771insGlyLeu",
        "NP_005219.2:p.Asp770_Gly771insGlyLeu",
        "BRAF E599_V600insV"
    ]:
        await validator_checks(
            test_tokenizer, test_classifier, q, val_params, validator_instance,
            expect_valid=False
        )


@pytest.mark.asyncio
async def test_cdna_insertion(test_tokenizer, test_classifier, val_params):
    """Test that cdna insertion validator works correctly"""
    validator_instance = CdnaInsertion
    for q in [
        "ENST00000000442.11:c.426_500insT",
        "NM_007294.3:c.2902_2903insTC",
        "ENST00000331728.9:c.2049_2050insA",
        "LIMK2 c.2049_2050insA"
    ]:
        await validator_checks(
            test_tokenizer, test_classifier, q, val_params, validator_instance
        )

    for q in [
        "NM_007294.3:c.7224_7225insTC",
        "LIMK2 c.486488_48649545656530insA"
    ]:
        await validator_checks(
            test_tokenizer, test_classifier, q, val_params, validator_instance,
            expect_valid=False
        )


@pytest.mark.asyncio
async def test_genomic_insertion(test_tokenizer, test_classifier, val_params):
    """Test that genomic insertion validator works correctly"""
    validator_instance = GenomicInsertion
    for q in [
        "NC_000022.10:g.30051593_30051594insT",
        "NC_000017.10:g.37880993_37880994insGCTTACGTGATG",
        "ERBB2 g.37880993_37880994insGCTTACGTGATG",
        "chr6-31239170-C-CA"
    ]:
        await validator_checks(
            test_tokenizer, test_classifier, q, val_params, validator_instance
        )

    for q in [
        "NC_000022.10:g.51304566_51304567insT",
        "NC_000022.10:g.51304567_51304568insT"
    ]:
        await validator_checks(
            test_tokenizer, test_classifier, q, val_params, validator_instance,
            expect_valid=False
        )


@pytest.mark.asyncio
async def test_genomic_duplication(test_tokenizer, test_classifier, val_params):
    """Test that genomic duplication validator works correctly"""
    validator_instance = GenomicDuplication
    for q in [
        "NC_000003.12:g.49531262dup",
        "NC_000016.10:g.2087938_2087948dup"
    ]:
        await validator_checks(
            test_tokenizer, test_classifier, q, val_params, validator_instance
        )

    for q in [
        "NC_000003.12:g.495312625165465465465dup",
        "NC_000016.10:g.2087948_2087938dup"
    ]:
        await validator_checks(
            test_tokenizer, test_classifier, q, val_params, validator_instance,
            expect_valid=False
        )


@pytest.mark.asyncio
async def test_genomic_duplication_ambiguous(
    test_tokenizer, test_classifier, val_params
):
    """Test that genomic duplication ambiguous validator works correctly"""
    validator_instance = GenomicDuplicationAmbiguous
    for q in [
        "NC_000020.11:g.(?_30417576)_(31394018_?)dup",
        "NC_000023.11:g.(?_154021812)_154092209dup",
        "NC_000023.11:g.154021812_(154092209_?)dup",
        "NC_000023.11:g.(31060227_31100351)_(33274278_33417151)dup",
        "NC_000023.10:g.(31078344_31118468)_(33292395_33435268)dup"
    ]:
        await validator_checks(
            test_tokenizer, test_classifier, q, val_params, validator_instance
        )

    await validator_checks(
        test_tokenizer, test_classifier, "NC_000023.11:g.(?_154092209)_154021812dup",
        val_params, validator_instance, expect_valid=False
    )


@pytest.mark.asyncio
async def test_amplification(test_tokenizer, test_classifier, val_params):
    """Test that amplification validator works correctly"""
    validator_instance = Amplification
    await validator_checks(
        test_tokenizer, test_classifier, "BRAF Amplification", val_params,
        validator_instance
    )
