"""Module containing regex"""
import re
from typing import List, Tuple

from variation.schemas.token_response_schema import TokenType
from variation.schemas.classification_response_schema import ClassificationType

CDNA_GENOMIC_SUBSTITUTION = re.compile(
    r"^(?P<pos>\d+)(?P<ref>[ACTGN])>(?P<alt>[ACTGN])$"
)

CDNA_GENOMIC_REFERENCE_AGREE = re.compile(r"^(?P<pos>\d+)=$")

CNDA_GENOMIC_DELETION = re.compile(
    r"^(?P<pos0>\d+)(_(?P<pos1>\d+))?del(?P<deleted_sequence>[ACTGN]+)?$"
)

GENOMIC_DELETION_AMBIGUOUS_1 = re.compile(
    r"^\((?P<pos0>\?|\d+)_(?P<pos1>\?|\d+)\)_\((?P<pos2>\?|\d+)_(?P<pos3>\?|\d+)\)del$"
)

GENOMIC_DELETION_AMBIGUOUS_2 = re.compile(
    r"^\((?P<pos0>\?|\d+)_(?P<pos1>\?|\d+)\)_(?P<pos3>\d+)del$"
)

GENOMIC_DELETION_AMBIGUOUS_3 = re.compile(
    r"^(?P<pos0>\d+)_\((?P<pos2>\?|\d+)_(?P<pos3>\?|\d+)\)del$"
)

CDNA_GENOMIC_DELINS = re.compile(
    r"^(?P<pos0>\d+)(_(?P<pos1>\d+))?delins(?P<inserted_sequence>[ACTGN]+)$"
)

CDNA_GENOMIC_INSERTION = re.compile(
    r"^(?P<pos0>\d+)_(?P<pos1>\d+)ins(?P<inserted_sequence>[ACTGN]+)$"
)

PROTEIN_SUBSTITUTION = re.compile(
    r"^(?P<ref>[a-zA-z]+)(?P<pos>\d+)(?P<alt>([a-zA-Z]|Ter|\*)+)$"
)

PROTEIN_INSERTION = re.compile(
    r"^(?P<aa0>[a-zA-z]+)(?P<pos0>\d+)_(?P<aa1>[a-zA-z]+)(?P<pos1>\d+)ins(?P<inserted_sequence>[a-zA-z]+)$"  # noqa: E501
)

PROTEIN_DELINS = re.compile(
    r"^(?P<aa0>[a-zA-z]+)(?P<pos0>\d+)(_(?P<aa1>[a-zA-z]+)(?P<pos1>\d+))?delins(?P<inserted_sequence>[a-zA-z]+)$"  # noqa: E501
)

PROTEIN_DELETION = re.compile(
    r"^(?P<aa0>[a-zA-z]+)(?P<pos0>\d+)(_(?P<aa1>[a-zA-z]+)(?P<pos1>\d+))?del(?P<deleted_sequence>[a-zA-z]+)?$"  # noqa: E501
)

PROTEIN_REFERENCE_AGREE = re.compile(
    r"^(?P<ref>[a-zA-z]+)(?P<pos>\d+)=$"
)

GENOMIC_DUPLICATION = re.compile(
    r"^(?P<pos0>\d+)(_(?P<pos1>\d+))?dup$"
)

GENOMIC_DUPLICATION_AMBIGUOUS_1 = re.compile(
    r"^\((?P<pos0>\?|\d+)_(?P<pos1>\?|\d+)\)_\((?P<pos2>\?|\d+)_(?P<pos3>\?|\d+)\)dup$"
)

GENOMIC_DUPLICATION_AMBIGUOUS_2 = re.compile(
    r"^\((?P<pos0>\?|\d+)_(?P<pos1>\?|\d+)\)_(?P<pos3>\d+)dup$"
)

GENOMIC_DUPLICATION_AMBIGUOUS_3 = re.compile(
    r"^(?P<pos0>\d+)_\((?P<pos2>\?|\d+)_(?P<pos3>\?|\d+)\)dup$"
)

PROTEIN_REGEXPRS: List[Tuple[any, TokenType, ClassificationType]] = [
    (
        PROTEIN_SUBSTITUTION,
        TokenType.PROTEIN_SUBSTITUTION,
        ClassificationType.PROTEIN_SUBSTITUTION
    ),
    (
        PROTEIN_REFERENCE_AGREE,
        TokenType.PROTEIN_REFERENCE_AGREE,
        ClassificationType.PROTEIN_REFERENCE_AGREE
    ),
    (
        PROTEIN_DELINS,
        TokenType.PROTEIN_DELINS,
        ClassificationType.PROTEIN_DELINS
    ),
    (
        PROTEIN_DELETION,
        TokenType.PROTEIN_DELETION,
        ClassificationType.PROTEIN_DELETION
    ),
    (
        PROTEIN_INSERTION,
        TokenType.PROTEIN_INSERTION,
        ClassificationType.PROTEIN_INSERTION
    )
]

CDNA_REGEXPRS: List[Tuple[any, TokenType, ClassificationType]] = [
    (
        CDNA_GENOMIC_SUBSTITUTION,
        TokenType.CODING_DNA_SUBSTITUTION,
        ClassificationType.CODING_DNA_SUBSTITUTION
    ),
    (
        CDNA_GENOMIC_DELINS,
        TokenType.CODING_DNA_DELINS,
        ClassificationType.CODING_DNA_DELINS
    ),
    (
        CNDA_GENOMIC_DELETION,
        TokenType.CODING_DNA_DELETION,
        ClassificationType.CODING_DNA_DELETION
    ),
    (
        CDNA_GENOMIC_INSERTION,
        TokenType.CODING_DNA_INSERTION,
        ClassificationType.CODING_DNA_INSERTION
    )
]

GENOMIC_REGEXPRS: List[Tuple[any, TokenType, ClassificationType]] = [
    (
        CDNA_GENOMIC_SUBSTITUTION,
        TokenType.GENOMIC_SUBSTITUTION,
        ClassificationType.GENOMIC_SUBSTITUTION
    ),
    (
        CDNA_GENOMIC_DELINS,
        TokenType.GENOMIC_DELINS,
        ClassificationType.GENOMIC_DELINS
    ),
    (
        CDNA_GENOMIC_INSERTION,
        TokenType.GENOMIC_INSERTION,
        ClassificationType.GENOMIC_INSERTION
    )
]
