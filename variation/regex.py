"""Module containing regex patterns"""
import re
from typing import Any, List, Tuple

from variation.schemas.app_schemas import AmbiguousRegexType
from variation.schemas.classification_response_schema import ClassificationType
from variation.schemas.token_response_schema import TokenType

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

PROTEIN_REFERENCE_AGREE = re.compile(r"^(?P<ref>[a-zA-z]+)(?P<pos>\d+)=$")

GENOMIC_DUPLICATION = re.compile(r"^(?P<pos0>\d+)(_(?P<pos1>\d+))?dup$")

# (#_#)_(#_#) OR (?_#)_(#_?)
GENOMIC_DUPLICATION_AMBIGUOUS_1 = re.compile(
    r"^\((?P<pos0>\?|\d+)_(?P<pos1>\?|\d+)\)_\((?P<pos2>\?|\d+)_(?P<pos3>\?|\d+)\)dup$"
)

# (?_#)_#, (#_?)_# OR (#_#)_#
GENOMIC_DUPLICATION_AMBIGUOUS_2 = re.compile(
    r"^\((?P<pos0>\?|\d+)_(?P<pos1>\?|\d+)\)_(?P<pos2>\d+)dup$"
)

# #_(#_?) OR #_(#_#)
GENOMIC_DUPLICATION_AMBIGUOUS_3 = re.compile(
    r"^(?P<pos0>\d+)_\((?P<pos2>\?|\d+)_(?P<pos3>\?|\d+)\)dup$"
)

# (#_#)_(#_#) OR (?_#)_(#_?)
GENOMIC_DELETION_AMBIGUOUS_1 = re.compile(
    r"^\((?P<pos0>\?|\d+)_(?P<pos1>\?|\d+)\)_\((?P<pos2>\?|\d+)_(?P<pos3>\?|\d+)\)del$"
)

# (?_#)_#, (#_?)_# OR (#_#)_#
GENOMIC_DELETION_AMBIGUOUS_2 = re.compile(
    r"^\((?P<pos0>\?|\d+)_(?P<pos1>\?|\d+)\)_(?P<pos2>\d+)del$"
)

# #_(#_?) OR #_(#_#)
GENOMIC_DELETION_AMBIGUOUS_3 = re.compile(
    r"^(?P<pos0>\d+)_\((?P<pos2>\?|\d+)_(?P<pos3>\?|\d+)\)del$"
)

# _REGEXPRS are used to help group the regex pattern and associated token type and
# classification type

# Note: Order matters for regexprs
PROTEIN_REGEXPRS: List[Tuple[Any, TokenType, ClassificationType]] = [
    (PROTEIN_DELINS, TokenType.PROTEIN_DELINS, ClassificationType.PROTEIN_DELINS),
    (PROTEIN_DELETION, TokenType.PROTEIN_DELETION, ClassificationType.PROTEIN_DELETION),
    (
        PROTEIN_SUBSTITUTION,
        TokenType.PROTEIN_SUBSTITUTION,
        ClassificationType.PROTEIN_SUBSTITUTION,
    ),
    (
        PROTEIN_REFERENCE_AGREE,
        TokenType.PROTEIN_REFERENCE_AGREE,
        ClassificationType.PROTEIN_REFERENCE_AGREE,
    ),
    (
        PROTEIN_INSERTION,
        TokenType.PROTEIN_INSERTION,
        ClassificationType.PROTEIN_INSERTION,
    ),
]

# Note: Order matters for regexprs
CDNA_REGEXPRS: List[Tuple[Any, TokenType, ClassificationType]] = [
    (CDNA_GENOMIC_DELINS, TokenType.CDNA_DELINS, ClassificationType.CDNA_DELINS),
    (CNDA_GENOMIC_DELETION, TokenType.CDNA_DELETION, ClassificationType.CDNA_DELETION),
    (
        CDNA_GENOMIC_SUBSTITUTION,
        TokenType.CDNA_SUBSTITUTION,
        ClassificationType.CDNA_SUBSTITUTION,
    ),
    (
        CDNA_GENOMIC_REFERENCE_AGREE,
        TokenType.CDNA_REFERENCE_AGREE,
        ClassificationType.CDNA_REFERENCE_AGREE,
    ),
    (
        CDNA_GENOMIC_INSERTION,
        TokenType.CDNA_INSERTION,
        ClassificationType.CDNA_INSERTION,
    ),
]

# Note: Order matters for regexprs
GENOMIC_REGEXPRS: List[Tuple[Any, TokenType, ClassificationType]] = [
    (CDNA_GENOMIC_DELINS, TokenType.GENOMIC_DELINS, ClassificationType.GENOMIC_DELINS),
    (
        CNDA_GENOMIC_DELETION,
        TokenType.GENOMIC_DELETION,
        ClassificationType.GENOMIC_DELETION,
    ),
    (
        CDNA_GENOMIC_SUBSTITUTION,
        TokenType.GENOMIC_SUBSTITUTION,
        ClassificationType.GENOMIC_SUBSTITUTION,
    ),
    (
        CDNA_GENOMIC_REFERENCE_AGREE,
        TokenType.GENOMIC_REFERENCE_AGREE,
        ClassificationType.GENOMIC_REFERENCE_AGREE,
    ),
    (
        CDNA_GENOMIC_INSERTION,
        TokenType.GENOMIC_INSERTION,
        ClassificationType.GENOMIC_INSERTION,
    ),
    (
        GENOMIC_DUPLICATION,
        TokenType.GENOMIC_DUPLICATION,
        ClassificationType.GENOMIC_DUPLICATION,
    ),
]


# Note: Order matters for regexprs
GENOMIC_DUP_AMBIGUOUS_REGEXPRS: List[
    Tuple[Any, TokenType, ClassificationType, AmbiguousRegexType]
] = [
    (
        GENOMIC_DUPLICATION_AMBIGUOUS_1,
        TokenType.GENOMIC_DUPLICATION_AMBIGUOUS,
        ClassificationType.GENOMIC_DUPLICATION_AMBIGUOUS,
        AmbiguousRegexType.REGEX_1,
    ),
    (
        GENOMIC_DUPLICATION_AMBIGUOUS_2,
        TokenType.GENOMIC_DUPLICATION_AMBIGUOUS,
        ClassificationType.GENOMIC_DUPLICATION_AMBIGUOUS,
        AmbiguousRegexType.REGEX_2,
    ),
    (
        GENOMIC_DUPLICATION_AMBIGUOUS_3,
        TokenType.GENOMIC_DUPLICATION_AMBIGUOUS,
        ClassificationType.GENOMIC_DUPLICATION_AMBIGUOUS,
        AmbiguousRegexType.REGEX_3,
    ),
]


# Note: Order matters for regexprs
GENOMIC_DEL_AMBIGUOUS_REGEXPRS: List[
    Tuple[Any, TokenType, ClassificationType, AmbiguousRegexType]
] = [
    (
        GENOMIC_DELETION_AMBIGUOUS_1,
        TokenType.GENOMIC_DELETION_AMBIGUOUS,
        ClassificationType.GENOMIC_DELETION_AMBIGUOUS,
        AmbiguousRegexType.REGEX_1,
    ),
    (
        GENOMIC_DELETION_AMBIGUOUS_2,
        TokenType.GENOMIC_DELETION_AMBIGUOUS,
        ClassificationType.GENOMIC_DELETION_AMBIGUOUS,
        AmbiguousRegexType.REGEX_2,
    ),
    (
        GENOMIC_DELETION_AMBIGUOUS_3,
        TokenType.GENOMIC_DELETION_AMBIGUOUS,
        ClassificationType.GENOMIC_DELETION_AMBIGUOUS,
        AmbiguousRegexType.REGEX_3,
    ),
]
