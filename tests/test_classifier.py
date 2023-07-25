"""Module for testing classifiers"""
from variation.schemas.classification_response_schema import (
    AmplificationClassification, ProteinSubstitutionClassification,
    ProteinStopGainClassification, ProteinReferenceAgreeClassification,
    CdnaSubstitutionClassification, GenomicSubstitutionClassification,
    CdnaReferenceAgreeClassification, GenomicReferenceAgreeClassification,
    ProteinDelInsClassification, CdnaDelInsClassification, GenomicDelInsClassification,
    ProteinDeletionClassification, CdnaDeletionClassification,
    GenomicDeletionClassification, GenomicDeletionAmbiguousClassification,
    ProteinInsertionClassification, CdnaInsertionClassification,
    GenomicInsertionClassification, GenomicDuplicationClassification,
    GenomicDuplicationAmbiguousClassification
)


def _get_classification(test_tokenizer, test_classifier, q):
    """Get classification for a query"""
    assert isinstance(q, str)
    tokens = test_tokenizer.perform(q, [])
    return test_classifier.perform(tokens)


def assert_match(
    test_tokenizer, test_classifier, q, expected_classification_instance
):
    """Check that a single classification match is working correctly"""
    classification = _get_classification(test_tokenizer, test_classifier, q)
    assert isinstance(
        classification,
        expected_classification_instance
    ), f"{q} is instance {type(classification)}"


def assert_no_match(test_tokenizer, test_classifier, q):
    """Check that no classification match is working correctly"""
    classification = _get_classification(test_tokenizer, test_classifier, q)
    assert classification is None, q


def test_amplification(test_tokenizer, test_classifier):
    """Test that amplification classifier works"""
    for q in ["BRAF amplification", "braf AMPLIFICATION"]:
        assert_match(test_tokenizer, test_classifier, q, AmplificationClassification)

    for q in ["gene amplification"]:
        assert_no_match(test_tokenizer, test_classifier, q)


def test_protein_substitution(test_tokenizer, test_classifier):
    """Test that protein substitution classifier works"""
    for q in [
        "BRAF V600E",
        "braf V600E",
        "NRAS G13V",
        "NP_004324.2:p.Val600Glu",
        "NP_065681.1:p.Met918Thr"
    ]:
        assert_match(
            test_tokenizer, test_classifier, q, ProteinSubstitutionClassification
        )

    for q in [
        "braf v600e",
        "BRAFV600E",
        "v600z",
        "V600E BRAF",
        "BRAF V600E foo",
        "BRAF",
        "V600E",
        "(V600E)",
        "NP_065681.1:c.Met918Thr",
    ]:
        assert_no_match(test_tokenizer, test_classifier, q)


def test_cdna_substitution(test_tokenizer, test_classifier):
    """Test that cdna substitution classifier works"""
    for q in [
        "NM_000551.3:c.292T>C",
        "BRAF V600E c.23T>A"
    ]:
        assert_match(test_tokenizer, test_classifier, q, CdnaSubstitutionClassification)

    for q in [
        "V170 (c.509F>A)",
        "RX_:g.292TC",
        "V170D (c.509T>A)",
        "NM_000551.3:c.292TC",
        "foo Y98H (c.292T>C)",
        "LRG_199t1:c.54G>H",
    ]:
        assert_no_match(test_tokenizer, test_classifier, q)


def test_genomic_substitution(test_tokenizer, test_classifier):
    """Test that genomic substitution classifier works"""
    for q in [
        "NC_000017.10:g.292T>C",
        "BRAF V600E g.23T>A",
        "7-292-A-C",
        "chrX-292-A-T",
        "chromosome10-292-G-A"
    ]:
        assert_match(
            test_tokenizer, test_classifier, q, GenomicSubstitutionClassification
        )

    for q in [
        "V170 (g.509F>A)",
        "RX_:c.292TC",
        "V170D (g.509T>A)",
        "NC_000017.10:g.292TC",
        "foo Y98H (g.292T>C)"
    ]:
        assert_no_match(test_tokenizer, test_classifier, q)


def test_protein_stop_gain(test_tokenizer, test_classifier):
    """Test that protein stop gain classifier works"""
    assert_match(
        test_tokenizer, test_classifier, "ENSP00000343204.4:p.Trp690Ter",
        ProteinStopGainClassification
    )

    assert_no_match(test_tokenizer, test_classifier, "ENS00000343204.4:c.Trp690Ter")


def test_protein_reference_agree(test_tokenizer, test_classifier):
    """Test that protein reference agree classifier works"""
    assert_match(
        test_tokenizer, test_classifier, "NP_000213.1:p.Leu862=",
        ProteinReferenceAgreeClassification
    )

    for q in [
        "Leu862==",
        "NP_000213.1:p.Leu862=="
    ]:
        assert_no_match(test_tokenizer, test_classifier, q)


def test_cdna_reference_agree(test_tokenizer, test_classifier):
    """Test that cdna reference agree classifier works"""
    assert_match(
        test_tokenizer, test_classifier, "NM_004006.2:c.123=",
        CdnaReferenceAgreeClassification
    )

    for q in [
        "CODING_DNA_:c.123=",
        "g.123=",
        "foo VHL c.123="
    ]:
        assert_no_match(test_tokenizer, test_classifier, q)


def test_genomic_reference_agree(test_tokenizer, test_classifier):
    """Test that genomic reference agree classifier works"""
    for q in [
        "NC_000017.10:g.123=",
        "chr11-252-t-t"
    ]:
        assert_match(
            test_tokenizer, test_classifier, q, GenomicReferenceAgreeClassification
        )

    for q in [
        "GENOMIC_:g.123=",
        "c.123=",
        "foo VHL g.123="
    ]:
        assert_no_match(test_tokenizer, test_classifier, q)


def test_protein_delins(test_tokenizer, test_classifier):
    """Test that protein delins classifier works"""
    for q in [
        "NP_001333827.1:p.Leu747_Thr751delinsPro",
        "NP_001333827.1:p.Leu747delinsProArg",
        "NP_005219.2:p.Glu746_Thr751delinsValAla",
        "NP_005219.2:p.G776delinsVC"
    ]:
        assert_match(test_tokenizer, test_classifier, q, ProteinDelInsClassification)

    for q in [
        "EGFR E709_G719delins11",
        "N:p.Leu747_Thr751delinsPro"
    ]:
        assert_no_match(test_tokenizer, test_classifier, q)


def test_cdna_delins(test_tokenizer, test_classifier):
    """Test that cdna delins classifier works"""
    for q in [
        "NM_005157.6:c.1423_1424delinsGT",
        "ENST00000277541.6:c.7330delinsACA",
        "NM_000551.3:c.615delinsAA",
        "ENST00000257290.5:c.2524_2525delinsAT"
    ]:
        assert_match(test_tokenizer, test_classifier, q, CdnaDelInsClassification)

    for q in [
        "N_005157.6:g.1423_1424delinsGT",
        "c.1423delinsX",
        "NM_000797.3:c.812_829delins908_925",
        "foo c.131_234delinsA",
        "foo NM_005157.6:c.1423_1424delinsGT",
        "LRG_199t1:c.79_80delinsTT",
        "LRG_199:c.79_80delinsTT",
    ]:
        assert_no_match(test_tokenizer, test_classifier, q)


def test_genomic_delins(test_tokenizer, test_classifier):
    """Test that genomic delins classifier works"""
    for q in [
        "NC_000017.10:g.1423_1424delinsGT",
        "NC_000017.10:g.7330delinsACA",
        "NC_000003.12:g.10149938delinsAA"
    ]:
        assert_match(test_tokenizer, test_classifier, q, GenomicDelInsClassification)

    for q in [
        "N_000017.10:c.1423_1424delinsGT",
        "g.1423delinsX",
        "NC_000017.10:g.812_829delins908_925",
        "foo g.131_234delinsA",
        "foo  NC_000017.10:g.1423_1424delinsGT"

    ]:
        assert_no_match(test_tokenizer, test_classifier, q)


def test_protein_deletion(test_tokenizer, test_classifier):
    """Test that protein deletion classifier works"""
    for q in [
        "NP_004439.2:p.Leu755_Thr759del",
        "NP_000213.1:p.Val560del",
        "NP_000213.1:p.Lys550_Lys558del",
        "KIT D419del",
        "KIT E554_V559del",
        "CTNNB1 Y30_I35del",
        "ENSP00000256474.2:p.Phe76del",
        "EGFR L747_T751delLREAT"
    ]:
        assert_match(test_tokenizer, test_classifier, q, ProteinDeletionClassification)

    for q in [
        "fakegene g.Leu755_Thr759delLeu",
        "GENE c.L755del",
        "NP_004439.2:c.Leu755_Thr759del",
        "LRG_199p1:p.Val7del",
        "LRG_199p1:p.(Val7del)",
        "NP_003997.1:p.(Lys23_Val25del)"
    ]:
        assert_no_match(test_tokenizer, test_classifier, q)


def test_cdna_deletion(test_tokenizer, test_classifier):
    """Test that cdna deletion classifier works"""
    for q in [
        "ENST00000269571.5:c.2263_2277del",
        "NM_004448.3:c.2263_2277delTTGAGGGAAAACACA",
        "NM_000535.6:c.2117delA",
        "ENST00000256474.2:c.163delG",
        "MLH1 c.1852_1854delAAG"
    ]:
        assert_match(test_tokenizer, test_classifier, q, CdnaDeletionClassification)

    for q in [
        "GENE c.1799_1800delTGinsAT",
        "GENE c.2524_2525delinsAT",
        "NM_004333.4:c.1799_1800delTGinsAT"
    ]:
        assert_no_match(test_tokenizer, test_classifier, q)


def test_genomic_deletion(test_tokenizer, test_classifier):
    """Test that genomic deletion classifier works"""
    for q in [
        "NC_000017.10:g.37880219_37880233del",
        "NC_000004.11:g.55593610_55593615delTTGTTG",
        "NC_000003.11:g.10183645del",
        "NC_000003.11:g.10188302delG",
        "Y-1313-ATTGAC-a",
        "chr1-2-ca-C"
    ]:
        assert_match(test_tokenizer, test_classifier, q, GenomicDeletionClassification)

    for q in [
        "GENE g.152419920_152419921delinsAG",
        "GENE g.152419920_152419921delAinsG"
    ]:
        assert_no_match(test_tokenizer, test_classifier, q)


def test_genomic_deletion_ambiguous(test_tokenizer, test_classifier):
    """Test that genomic deletion ambiguous classifier works"""
    for q in [
        "NC_000023.11:g.(?_31120496)_(33339477_?)del",
        "NC_000023.11:g.(?_155980375)_(156013167_?)del",
        "NC_000023.11:g.(31060227_31100351)_(33274278_33417151)del",
        "BRAF g.(31060227_31100351)_(33274278_33417151)del"
    ]:
        assert_match(
            test_tokenizer, test_classifier, q, GenomicDeletionAmbiguousClassification
        )

    for q in [
        "GENE (?_155980375)_(156013167_?)del",
        "accession:g.(?_155980375)_(156013167_?)del",
        "NC_000023.11:g.(?_155980375)_(156013167_?)del foo",
        "GENE (?_31120496)_(33339477_?)del"
    ]:
        assert_no_match(test_tokenizer, test_classifier, q)


def test_protein_insertion(test_tokenizer, test_classifier):
    """Test that protein insertion classifier works"""
    for q in [
        "NP_005219.2:p.Cys770_Gly771insGlyLeu",
        "NP_001333827.1:p.Ala763_Tyr764insPheGlnGluAla",
        "BRAF T599_V600insV",
        "EGFR A763_Y764insFQEA"
    ]:
        assert_match(test_tokenizer, test_classifier, q, ProteinInsertionClassification)

    for q in [
        "GENE p.Lys23insAsp",
        "GENE Lys23insAsp",
        "GENE p.His4_Gln5insAlaG",
        "ACCESSION_23042.2:p.His4_Gln5insAla"
    ]:
        assert_no_match(test_tokenizer, test_classifier, q)


def test_cdna_insertion(test_tokenizer, test_classifier):
    """Test that cdna insertion classifier works"""
    for q in [
        "NM_000551.3:c.230_231insTCT",
        "NM_000551.3:c.358_359insAC"
    ]:
        assert_match(test_tokenizer, test_classifier, q, CdnaInsertionClassification)

    for q in [
        "GENE 358_359insAC",
        "accession:c.358_359insAC",
        "NM_004006.2:c.849_850ins858_895",
        "NM_000551.3:c.358_359insAC foo"
    ]:
        assert_no_match(test_tokenizer, test_classifier, q)


def test_genomic_insertion(test_tokenizer, test_classifier):
    """Test that genomic insertion classifier works"""
    for q in [
        "NC_000023.10:g.32867861_32867862insT",
        "NC_000023.10:g.32862923_32862924insCCT",
        "NC_000009.11:g.5070053_5070054insG",
        "20-14223252-T-TATGCATG",
        "chr17-131543-G-GA"
    ]:
        assert_match(test_tokenizer, test_classifier, q, GenomicInsertionClassification)

    for q in [
        "GENE 32867861_32867862insT",
        "accession:g.32867861_32867862insT",
        "NC_000023.10:g.32867861_32867862insT foo"
    ]:
        assert_no_match(test_tokenizer, test_classifier, q)


def test_genomic_duplication(test_tokenizer, test_classifier):
    """Test that genomic duplication classifier works"""
    for q in [
        "NC_000003.12:g.49531262dup",
        "NC_000016.10:g.2087938_2087948dup",
        "BRAF g.2087938_2087948dup"
    ]:
        assert_match(
            test_tokenizer, test_classifier, q, GenomicDuplicationClassification
        )

    for q in [
        "foo (?_30417576)_(31394018_?)dup",
        "Accession:g.49531262dup"
    ]:
        assert_no_match(test_tokenizer, test_classifier, q)


def test_genomic_duplication_ambiguous(test_tokenizer, test_classifier):
    """Test that genomic duplication ambiguous classifier works"""
    for q in [
        "NC_000020.11:g.(?_30417576)_(31394018_?)dup",
        "NC_000023.11:g.(31060227_31100351)_(33274278_33417151)dup",
        "NC_000023.11:g.(?_154021812)_154092209dup"
    ]:
        assert_match(
            test_tokenizer, test_classifier, q,
            GenomicDuplicationAmbiguousClassification
        )

    for q in [
        "foo (?_30417576)_(31394018_?)dup",
        "Accession:g.49531262dup",
        "NC_000023.11:g.(31060227_33274278)_(31100351_33417151)"
    ]:
        assert_no_match(test_tokenizer, test_classifier, q)
