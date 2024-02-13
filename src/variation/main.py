"""Main application for FastAPI."""
import datetime
import traceback
from enum import Enum
from typing import List, Optional, Union
from urllib.parse import unquote

import pkg_resources
from bioutils.exceptions import BioutilsError
from cool_seq_tool.schemas import Assembly, ResidueMode
from fastapi import FastAPI, Query
from ga4gh.vrs import models
from hgvs.exceptions import HGVSError
from pydantic import ValidationError

from variation import logger
from variation.query import QueryHandler
from variation.schemas import NormalizeService, ServiceMeta, ToVRSService
from variation.schemas.copy_number_schema import (
    AmplificationToCxVarService,
    ParsedToCnVarQuery,
    ParsedToCnVarService,
    ParsedToCxVarQuery,
    ParsedToCxVarService,
)
from variation.schemas.gnomad_vcf_to_protein_schema import GnomadVcfToProteinService
from variation.schemas.hgvs_to_copy_number_schema import (
    HgvsToCopyNumberChangeService,
    HgvsToCopyNumberCountService,
)
from variation.schemas.normalize_response_schema import (
    HGVSDupDelModeOption,
    TranslateIdentifierService,
)
from variation.schemas.service_schema import (
    ToCdnaService,
    ToGenomicService,
)
from variation.schemas.vrs_python_translator_schema import (
    TranslateFromFormat,
    TranslateFromQuery,
    TranslateFromService,
    TranslateToHGVSQuery,
    TranslateToQuery,
    TranslateToService,
    VrsPythonMeta,
)
from variation.version import __version__


class Tag(Enum):
    """Define tag names for endpoints"""

    MAIN = "Main"
    SEQREPO = "SeqRepo"
    TO_PROTEIN_VARIATION = "To Protein Variation"
    VRS_PYTHON = "VRS-Python"
    TO_COPY_NUMBER_VARIATION = "To Copy Number Variation"
    ALIGNMENT_MAPPER = "Alignment Mapper"


query_handler = QueryHandler()


app = FastAPI(
    title="The VICC Variation Normalizer",
    description="Services and guidelines for normalizing variations.",
    version=__version__,
    contact={
        "name": "Alex H. Wagner",
        "email": "Alex.Wagner@nationwidechildrens.org",
        "url": "https://www.nationwidechildrens.org/specialties/institute-for-genomic-medicine/research-labs/wagner-lab",
    },
    license={
        "name": "MIT",
        "url": "https://github.com/cancervariants/variation-normalization/blob/main/LICENSE",
    },
    docs_url="/variation",
    openapi_url="/variation/openapi.json",
    swagger_ui_parameters={"tryItOutEnabled": True},
)

translate_summary = (
    "Translate a HGVS, gnomAD VCF and Free Text descriptions to VRS" " variation(s)."
)
translate_description = (
    "Translate a human readable variation description to "
    "VRS variation(s)."
    " Performs fully-justified allele normalization. "
    " Does not do any liftover operations or make any inferences "
    "about the query."
)
translate_response_description = "A  response to a validly-formed query."
q_description = "HGVS, gnomAD VCF or Free Text description on GRCh37 or GRCh38 assembly"


@app.get(
    "/variation/to_vrs",
    summary=translate_summary,
    response_description=translate_response_description,
    response_model=ToVRSService,
    response_model_exclude_none=True,
    description=translate_description,
    tags=[Tag.MAIN],
)
async def to_vrs(
    q: str = Query(..., description=q_description),
) -> ToVRSService:
    """Translate a HGVS, gnomAD VCF and Free Text descriptions to VRS variation(s).
    Performs fully-justified allele normalization. Does not do any liftover operations
    or make any inferences about the query.

    :param q: HGVS, gnomAD VCF or Free Text description on GRCh37 or GRCh38 assembly
    :return: ToVRSService model for variation
    """
    return await query_handler.to_vrs_handler.to_vrs(unquote(q))


normalize_summary = (
    "Normalizes and translates a HGVS, gnomAD VCF or Free Text description on GRCh37 "
    "or GRCh38 assembly to a single VRS Variation."
)
normalize_response_description = "A response to a validly-formed query."
normalize_description = (
    "Normalizes and translates a human readable variation description to a single VRS "
    "Variation. Performs fully-justified allele normalization. Will liftover to GRCh38 "
    "and aligns to a priority transcript. Will make inferences about the query."
)
q_description = "HGVS, gnomAD VCF or Free Text description on GRCh37 or GRCh38 assembly"
hgvs_dup_del_mode_decsr = (
    "This parameter determines how to interpret HGVS dup/del expressions in VRS."
)


@app.get(
    "/variation/normalize",
    summary=normalize_summary,
    response_description=normalize_response_description,
    response_model=NormalizeService,
    response_model_exclude_none=True,
    description=normalize_description,
    tags=[Tag.MAIN],
)
async def normalize(
    q: str = Query(..., description=q_description),
    hgvs_dup_del_mode: Optional[HGVSDupDelModeOption] = Query(
        HGVSDupDelModeOption.DEFAULT, description=hgvs_dup_del_mode_decsr
    ),
    baseline_copies: Optional[int] = Query(
        None,
        description="Baseline copies for HGVS duplications and deletions represented as Copy Number Count Variation",
    ),
    copy_change: Optional[models.CopyChange] = Query(
        None,
        description="The copy change for HGVS duplications and deletions represented as Copy Number Change Variation.",
    ),
) -> NormalizeService:
    """Normalize and translate a HGVS, gnomAD VCF or Free Text description on GRCh37
    or GRCh38 assembly to a single VRS Variation. Performs fully-justified allele
    normalization. Will liftover to GRCh38 and aligns to a priority transcript. Will
    make inferences about the query.

    :param q: HGVS, gnomAD VCF or Free Text description on GRCh37 or GRCh38 assembly
    :param hgvs_dup_del_mode: This parameter determines how to interpret HGVS dup/del
        expressions in VRS.
    :param baseline_copies: Baseline copies for HGVS duplications and deletions.
        Required when `hgvs_dup_del_mode` is set to `copy_number_count`.
    :param copy_change: The copy change for HGVS duplications and deletions represented
        as Copy Number Change Variation. If not set, will use default `copy_change` for
        query.
    :return: NormalizeService for variation
    """
    return await query_handler.normalize_handler.normalize(
        unquote(q),
        hgvs_dup_del_mode=hgvs_dup_del_mode,
        baseline_copies=baseline_copies,
        copy_change=copy_change,
    )


@app.get(
    "/variation/translate_identifier",
    summary="Given an identifier, use SeqRepo to return a list of aliases.",
    response_description="A response to a validly-formed query.",
    response_model=TranslateIdentifierService,
    response_model_exclude_none=True,
    description="Return list of aliases for an identifier",
    tags=[Tag.SEQREPO],
)
def translate_identifier(
    identifier: str = Query(..., description="The identifier to find aliases for"),
    target_namespaces: Optional[str] = Query(
        None, description="The namespaces of the aliases, separated by commas"
    ),
) -> TranslateIdentifierService:
    """Return data containing identifier aliases.

    :param str identifier: The identifier to find aliases for
    :param Optional[str] target_namespaces: The namespaces of the aliases,
        separated by commas
    :return: TranslateIdentifierService data
    """
    aliases = []
    warnings = []
    identifier = identifier.strip()
    try:
        aliases = query_handler.seqrepo_access.sr.translate_identifier(
            identifier, target_namespaces=target_namespaces
        )
    except KeyError:
        warnings = [f"Identifier, {identifier}, does not exist in SeqRepo"]
    except Exception as e:
        warnings = [f"SeqRepo could not translate identifier, {identifier}:" f" {e}"]

    return TranslateIdentifierService(
        identifier_query=identifier,
        warnings=warnings,
        aliases=aliases,
        service_meta_=ServiceMeta(
            version=__version__,
            response_datetime=datetime.datetime.now(tz=datetime.timezone.utc),
        ),
    )


from_fmt_descr = (
    "Format of input variation to translate. Must be one of `beacon`, "
    "`gnomad`, `hgvs`, or `spdi`"
)


@app.get(
    "/variation/translate_from",
    summary="Given variation as beacon, gnomad, hgvs or spdi representation, "
    "return VRS Allele object using vrs-python's translator class",
    response_description="A response to a validly-formed query.",
    response_model_exclude_none=True,
    description="Return VRS Allele object",
    response_model=TranslateFromService,
    tags=[Tag.VRS_PYTHON],
)
def vrs_python_translate_from(
    variation: str = Query(
        ...,
        description="Variation to translate to VRS object."
        " Must be represented as either beacon, "
        "gnomad, hgvs, or spdi.",
    ),
    fmt: Optional[TranslateFromFormat] = Query(None, description=from_fmt_descr),
) -> TranslateFromService:
    """Given variation query, return VRS Allele object using vrs-python"s translator
        class

    :param str variation: Variation to translate to VRS object. Must be represented
        as either beacon, gnomad, hgvs, or spdi
    :param Optional[TranslateFromFormat] fmt: Format of variation. If not supplied,
        vrs-python will infer its format.
    :return: TranslateFromService containing VRS Allele object
    """
    variation_query = unquote(variation.strip())
    warnings = []
    vrs_variation = None
    try:
        resp = query_handler.vrs_python_tlr.translate_from(variation_query, fmt)
    except (
        KeyError,
        ValueError,
        ValidationError,
    ) as e:
        warnings.append(f"vrs-python translator raised {type(e).__name__}: {e}")
    except HGVSError as e:
        warnings.append(f"hgvs raised {type(e).__name__}: {e}")
    except BioutilsError as e:
        warnings.append(f"bioutils raised {type(e).__name__}: {e}")
    else:
        vrs_variation = resp.model_dump(exclude_none=True)

    return TranslateFromService(
        query=TranslateFromQuery(variation=variation_query, fmt=fmt),
        warnings=warnings,
        variation=vrs_variation,
        service_meta_=ServiceMeta(
            version=__version__,
            response_datetime=datetime.datetime.now(tz=datetime.timezone.utc),
        ),
        vrs_python_meta_=VrsPythonMeta(
            version=pkg_resources.get_distribution("ga4gh.vrs").version
        ),
    )


g_to_p_summary = (
    "Given GRCh38 gnomAD VCF, return VRS Variation object on MANE protein coordinate."
)
g_to_p_response_description = "A response to a validly-formed query."
g_to_p_description = (
    "Return VRS Variation object on protein coordinate for variation provided."
)
q_description = (
    "GRCh38 gnomAD VCF (chr-pos-ref-alt) to normalize to MANE protein variation."
)


@app.get(
    "/variation/gnomad_vcf_to_protein",
    summary=g_to_p_summary,
    response_description=g_to_p_response_description,
    response_model_exclude_none=True,
    description=g_to_p_description,
    response_model=GnomadVcfToProteinService,
    tags=[Tag.TO_PROTEIN_VARIATION],
)
async def gnomad_vcf_to_protein(
    q: str = Query(..., description=q_description),
) -> GnomadVcfToProteinService:
    """Return VRS representation for variation on protein coordinate.

    :param q: gnomad VCF to normalize to protein variation.
    :return: GnomadVcfToProteinService for variation
    """
    q = unquote(q.strip())
    return await query_handler.gnomad_vcf_to_protein_handler.gnomad_vcf_to_protein(q)


hgvs_dup_del_mode_decsr = (
    "This parameter determines how to interpret HGVS dup/del expressions in VRS."
)


def _get_allele(
    request_body: Union[TranslateToQuery, TranslateToHGVSQuery], warnings: List
) -> Optional[models.Allele]:
    """Return VRS allele object from request body. `warnings` will get updated if
    exceptions are raised

    :param Union[TranslateToQuery, TranslateToHGVSQuery] request_body: Request body
        containing `variation`
    :param List warnings: List of warnings
    :return: VRS Allele object if valid
    """
    allele = None
    try:
        allele = models.Allele(**request_body["variation"])
    except ValidationError as e:
        warnings.append(f"`allele` is not a valid VRS Allele: {e}")
    return allele


@app.post(
    "/variation/translate_to",
    summary="Given VRS Allele object as a dict, return variation expressed as "
    "queried format using vrs-python's translator class",
    response_description="A response to a validly-formed query.",
    description="Return variation in queried format representation. "
    "Request body must contain `variation` and `fmt`. `variation` is"
    " a VRS Allele object represented as a dict. `fmt` must be either"
    " `spdi` or `hgvs`",
    response_model=TranslateToService,
    response_model_exclude_none=True,
    tags=[Tag.VRS_PYTHON],
)
async def vrs_python_translate_to(request_body: TranslateToQuery) -> TranslateToService:
    """Given VRS Allele object as a dict, return variation expressed as queried
    format using vrs-python's translator class

    :param TranslateToQuery request_body: Request body. `variation` is a VRS Allele
        object represented as a dict. `fmt` must be either `spdi` or `hgvs`
    :return: TranslateToService containing variation represented as fmt representation
        if valid VRS Allele, and warnings if found
    """
    query = request_body
    request_body = request_body.model_dump(by_alias=True)
    warnings = []

    allele = _get_allele(request_body, warnings)

    variations = []
    if allele:
        try:
            variations = query_handler.vrs_python_tlr.translate_to(
                allele, request_body["fmt"]
            )
        except ValueError as e:
            warnings.append(f"vrs-python translator raised {type(e).__name__}: {e}")

    return TranslateToService(
        query=query,
        warnings=warnings,
        variations=variations,
        service_meta_=ServiceMeta(
            version=__version__,
            response_datetime=datetime.datetime.now(tz=datetime.timezone.utc),
        ),
        vrs_python_meta_=VrsPythonMeta(
            version=pkg_resources.get_distribution("ga4gh.vrs").version
        ),
    )


to_hgvs_descr = (
    "Return variation as HGVS expressions. Request body must"
    " contain `variation`, a VRS Allele object represented as a dict. "
    "Can include optional parameter `namespace`. If `namespace` is not"
    " None, returns HGVS strings for the specified namespace. If "
    "`namespace` is None, returns HGVS strings for all alias translations."
)


@app.post(
    "/variation/vrs_allele_to_hgvs",
    summary="Given VRS Allele object as a dict, return HGVS expression(s)",
    response_description="A response to a validly-formed query.",
    description=to_hgvs_descr,
    response_model=TranslateToService,
    response_model_exclude_none=True,
    tags=[Tag.VRS_PYTHON],
)
async def vrs_python_to_hgvs(request_body: TranslateToHGVSQuery) -> TranslateToService:
    """Given VRS Allele object as a dict, return variation expressed as HGVS
        expression(s)

    :param TranslateToHGVSQuery request_body: Request body. `variation` is a VRS Allele
        object represented as a dict. Can provide optional parameter `namespace`.
        If `namespace` is not None, returns HGVS strings for the specified namespace.
        If `namespace` is None, returns HGVS strings for all alias translations.
    :return: TranslateToService containing variation represented as HGVS representation
        if valid VRS Allele, and warnings if found
    """
    query = request_body
    request_body = request_body.model_dump(by_alias=True)
    warnings = []

    allele = _get_allele(request_body, warnings)

    variations = []
    if allele:
        try:
            variations = query_handler.vrs_python_tlr._to_hgvs(
                allele, namespace=request_body.get("namespace") or "refseq"
            )
        except ValueError as e:
            warnings.append(f"vrs-python translator raised {type(e).__name__}: {e}")

    return TranslateToService(
        query=query,
        warnings=warnings,
        variations=variations,
        service_meta_=ServiceMeta(
            version=__version__,
            response_datetime=datetime.datetime.now(tz=datetime.timezone.utc),
        ),
        vrs_python_meta_=VrsPythonMeta(
            version=pkg_resources.get_distribution("ga4gh.vrs").version
        ),
    )


@app.get(
    "/variation/hgvs_to_copy_number_count",
    summary="Given HGVS expression, return VRS Copy Number Count Variation",
    response_description="A response to a validly-formed query.",
    description="Return VRS Copy Number Count Variation",
    response_model=HgvsToCopyNumberCountService,
    response_model_exclude_none=True,
    tags=[Tag.TO_COPY_NUMBER_VARIATION],
)
async def hgvs_to_copy_number_count(
    hgvs_expr: str = Query(..., description="Variation query"),
    baseline_copies: Optional[int] = Query(
        ..., description="Baseline copies for duplication"
    ),
    do_liftover: bool = Query(
        False, description="Whether or not to liftover " "to GRCh38 assembly."
    ),
) -> HgvsToCopyNumberCountService:
    """Given hgvs expression, return copy number count variation

    :param hgvs_expr: HGVS expression
    :param baseline_copies: Baseline copies number
    :param do_liftover: Whether or not to liftover to GRCh38 assembly
    :return: HgvsToCopyNumberCountService
    """
    return await query_handler.to_copy_number_handler.hgvs_to_copy_number_count(
        unquote(hgvs_expr.strip()),
        baseline_copies,
        do_liftover,
    )


@app.get(
    "/variation/hgvs_to_copy_number_change",
    summary="Given HGVS expression, return VRS Copy Number Change Variation",
    response_description="A response to a validly-formed query.",
    description="Return VRS Copy Number Change Variation",
    response_model=HgvsToCopyNumberChangeService,
    response_model_exclude_none=True,
    tags=[Tag.TO_COPY_NUMBER_VARIATION],
)
async def hgvs_to_copy_number_change(
    hgvs_expr: str = Query(..., description="Variation query"),
    copy_change: models.CopyChange = Query(..., description="The copy change"),
    do_liftover: bool = Query(
        False, description="Whether or not to liftover " "to GRCh38 assembly."
    ),
) -> HgvsToCopyNumberChangeService:
    """Given hgvs expression, return copy number change variation

    :param hgvs_expr: HGVS expression
    :param copy_change: copy change
    :param do_liftover: Whether or not to liftover to GRCh38 assembly
    :return: HgvsToCopyNumberChangeService
    """
    return await query_handler.to_copy_number_handler.hgvs_to_copy_number_change(
        unquote(hgvs_expr.strip()),
        copy_change,
        do_liftover,
    )


@app.post(
    "/variation/parsed_to_cn_var",
    summary="Given parsed genomic components, return VRS Copy Number Count "
    "Variation",
    response_description="A response to a validly-formed query.",
    description="Return VRS Copy Number Count Variation",
    response_model=ParsedToCnVarService,
    response_model_exclude_none=True,
    tags=[Tag.TO_COPY_NUMBER_VARIATION],
)
def parsed_to_cn_var(request_body: ParsedToCnVarQuery) -> ParsedToCnVarService:
    """Given parsed genomic components, return Copy Number Count Variation.

    :param request_body: Request body
    :return: ParsedToCnVarService containing Copy Number Count variation and list of
        warnings
    """
    try:
        resp = query_handler.to_copy_number_handler.parsed_to_copy_number(request_body)
    except Exception:
        traceback_resp = traceback.format_exc().splitlines()
        logger.exception(traceback_resp)

        return ParsedToCnVarService(
            copy_number_count=None,
            warnings=["Unhandled exception. See logs for more details."],
            service_meta_=ServiceMeta(
                version=__version__,
                response_datetime=datetime.datetime.now(tz=datetime.timezone.utc),
            ),
        )
    else:
        return resp


@app.post(
    "/variation/parsed_to_cx_var",
    summary="Given parsed genomic components, return VRS Copy Number Change "
    "Variation",
    response_description="A response to a validly-formed query.",
    description="Return VRS Copy Number Change Variation",
    response_model=ParsedToCxVarService,
    response_model_exclude_none=True,
    tags=[Tag.TO_COPY_NUMBER_VARIATION],
)
def parsed_to_cx_var(request_body: ParsedToCxVarQuery) -> ParsedToCxVarService:
    """Given parsed genomic components, return Copy Number Change Variation

    :param request_body: Request body
    :return: ParsedToCxVarService containing Copy Number Change variation and list of
        warnings
    """
    try:
        resp = query_handler.to_copy_number_handler.parsed_to_copy_number(request_body)
    except Exception:
        traceback_resp = traceback.format_exc().splitlines()
        logger.exception(traceback_resp)

        return ParsedToCxVarService(
            copy_number_count=None,
            warnings=["Unhandled exception. See logs for more details."],
            service_meta_=ServiceMeta(
                version=__version__,
                response_datetime=datetime.datetime.now(tz=datetime.timezone.utc),
            ),
        )
    else:
        return resp


amplification_to_cx_var_descr = (
    "Translate amplification to VRS Copy Number Change "
    "Variation. If `sequence`, `start`, and `end` are "
    "all provided, will return a SequenceLocation with "
    "those properties. Else, gene-normalizer will be "
    "used to retrieve the SequenceLocation."
)


@app.get(
    "/variation/amplification_to_cx_var",
    summary="Given amplification query, return VRS Copy Number Change Variation",
    response_description="A response to a validly-formed query.",
    description=amplification_to_cx_var_descr,
    response_model=AmplificationToCxVarService,
    response_model_exclude_none=True,
    tags=[Tag.TO_COPY_NUMBER_VARIATION],
)
def amplification_to_cx_var(
    gene: str = Query(..., description="Gene query"),
    sequence_id: Optional[str] = Query(None, description="Sequence identifier"),
    start: Optional[int] = Query(
        None, description="Start position as residue coordinate"
    ),
    end: Optional[int] = Query(None, description="End position as residue coordinate"),
) -> AmplificationToCxVarService:
    """Given amplification query, return Copy Number Change Variation
    Parameter priority:
        1. sequence, start, end (must provide ALL)
        2. use the gene-normalizer to get the SequenceLocation

    :param gene: Gene query
    :param sequence_id: Sequence ID for the location. If set, must also provide `start`
        and `end`
    :param start: Start position as residue coordinate for the sequence location.
        If set, must also provide `sequence` and `end`
    :param end: End position as residue coordinate for the sequence location. If set,
        must also provide `sequence` and `start`
    :return: AmplificationToCxVarService containing Copy Number Change and
        list of warnings
    """
    return query_handler.to_copy_number_handler.amplification_to_cx_var(
        gene=gene,
        sequence_id=sequence_id,
        start=start,
        end=end,
    )


@app.get(
    "/variation/alignment_mapper/p_to_c",
    summary="Translate protein representation to cDNA representation",
    response_description="A response to a validly-formed query.",
    description="Given protein accession and positions, return associated cDNA "
    "accession and positions to codon(s)",
    response_model=ToCdnaService,
    response_model_exclude_none=True,
    tags=[Tag.ALIGNMENT_MAPPER],
)
async def p_to_c(
    p_ac: str = Query(..., description="Protein RefSeq accession"),
    p_start_pos: int = Query(..., description="Protein start position"),
    p_end_pos: int = Query(..., description="Protein end position"),
    residue_mode: ResidueMode = Query(
        ResidueMode.RESIDUE,
        description="Residue mode for `p_start_pos` and `p_end_pos`",
    ),
) -> ToCdnaService:
    """Translate protein representation to cDNA representation

    :param str p_ac: Protein RefSeq accession
    :param int p_start_pos: Protein start position
    :param int p_end_pos: Protein end position
    :param ResidueMode residue_mode: Residue mode for `p_start_pos` and `p_end_pos`.
    :return: ToCdnaService containing cDNA representation, warnings, and
        service meta
    """
    try:
        c_data, w = await query_handler.alignment_mapper.p_to_c(
            p_ac, p_start_pos, p_end_pos, residue_mode
        )
    except Exception as e:
        logger.error("Unhandled exception: %s", str(e))
        w = "Unhandled exception. See logs for more information."
        c_data = None
    return ToCdnaService(
        c_data=c_data,
        warnings=[w] if w else [],
        service_meta=ServiceMeta(
            version=__version__,
            response_datetime=datetime.datetime.now(tz=datetime.timezone.utc),
        ),
    )


@app.get(
    "/variation/alignment_mapper/c_to_g",
    summary="Translate cDNA representation to genomic representation",
    response_description="A response to a validly-formed query.",
    description="Given cDNA accession and positions for codon(s), return associated genomic"
    " accession and positions for a given target genome assembly",
    response_model=ToGenomicService,
    response_model_exclude_none=True,
    tags=[Tag.ALIGNMENT_MAPPER],
)
async def c_to_g(
    c_ac: str = Query(..., description="cDNA RefSeq accession"),
    c_start_pos: int = Query(..., description="cDNA start position for codon"),
    c_end_pos: int = Query(..., description="cDNA end position for codon"),
    cds_start: Optional[int] = Query(
        None, description="CDS start site. If not provided, this will be computed."
    ),
    residue_mode: ResidueMode = Query(
        ResidueMode.RESIDUE,
        description="Residue mode for `c_start_pos` and `c_end_pos`",
    ),
    target_genome_assembly: Assembly = Query(
        Assembly.GRCH38, description="Genomic assembly to map to"
    ),
) -> ToGenomicService:
    """Translate cDNA representation to genomic representation

    :param str c_ac: cDNA RefSeq accession
    :param int c_start_pos: cDNA start position for codon
    :param int c_end_pos: cDNA end position for codon
    :param Optional[int] cds_start: CDS start site. If not provided, this will be
        computed.
    :param ResidueMode residue_mode: Residue mode for `c_start_pos` and `c_end_pos`.
    :param Assembly target_genome_assembly: Genome assembly to get genomic data for
    :return: ToGenomicService containing genomic representation, warnings, and
        service meta
    """
    try:
        g_data, w = await query_handler.alignment_mapper.c_to_g(
            c_ac,
            c_start_pos,
            c_end_pos,
            cds_start=cds_start,
            residue_mode=residue_mode,
            target_genome_assembly=target_genome_assembly,
        )
    except Exception as e:
        logger.error("Unhandled exception: %s", str(e))
        w = "Unhandled exception. See logs for more information."
        g_data = None
    return ToGenomicService(
        g_data=g_data,
        warnings=[w] if w else [],
        service_meta=ServiceMeta(
            version=__version__,
            response_datetime=datetime.datetime.now(tz=datetime.timezone.utc),
        ),
    )


@app.get(
    "/variation/alignment_mapper/p_to_g",
    summary="Translate protein representation to genomic representation",
    response_description="A response to a validly-formed query.",
    description="Given protein accession and positions, return associated genomic "
    "accession and positions for a given target genome assembly",
    response_model=ToGenomicService,
    response_model_exclude_none=True,
    tags=[Tag.ALIGNMENT_MAPPER],
)
async def p_to_g(
    p_ac: str = Query(..., description="Protein RefSeq accession"),
    p_start_pos: int = Query(..., description="Protein start position"),
    p_end_pos: int = Query(..., description="Protein end position"),
    residue_mode: ResidueMode = Query(
        ResidueMode.RESIDUE,
        description="Residue mode for `p_start_pos` and `p_end_pos`",
    ),
    target_genome_assembly: Assembly = Query(
        Assembly.GRCH38, description="Genomic assembly to map to"
    ),
) -> ToGenomicService:
    """Translate protein representation to genomic representation

    :param str p_ac: Protein RefSeq accession
    :param int p_start_pos: Protein start position
    :param int p_end_pos: Protein end position
    :param ResidueMode residue_mode: Residue mode for `p_start_pos` and `p_end_pos`.
    :param Assembly target_genome_assembly: Genome assembly to get genomic data for
    :return: ToGenomicService containing genomic representation, warnings, and
        service meta
    """
    try:
        g_data, w = await query_handler.alignment_mapper.p_to_g(
            p_ac,
            p_start_pos,
            p_end_pos,
            residue_mode=residue_mode,
            target_genome_assembly=target_genome_assembly,
        )
    except Exception as e:
        logger.error("Unhandled exception: %s", str(e))
        w = "Unhandled exception. See logs for more information."
        g_data = None
    return ToGenomicService(
        g_data=g_data,
        warnings=[w] if w else [],
        service_meta=ServiceMeta(
            version=__version__,
            response_datetime=datetime.datetime.now(tz=datetime.timezone.utc),
        ),
    )