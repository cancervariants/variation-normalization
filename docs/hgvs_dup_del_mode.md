# HGVS Dup Del Mode

This mode helps us interpret deletions and duplications that are represented as HGVS expressions.\
The mode can be set to `default`, `absolute_cnv`, `relative_cnv`, `repeated_seq_expr`, or `literal_seq_expr`.


## Default Characteristics

- if baseline_copies is not set and endpoints are ambiguous:
    - relative_cnv
    - if relative_copy_class not provided:
        - relative_copy_class = `EFO:0030067` (copy number loss) if del, `EFO:0030070` (copy number gain) if dup
- elif baseline_copies is provided:
    - absolute_cnv
    - copies are baseline_copies + 1 for dup, baseline_copies - 1 for del
- elif len del or dup > 100bp: (use outermost coordinates)
    - repeated_seq_expr with a derived_seq_expr subject (Allele)
- else:
    - literal_seq_expr (normalized LiteralSequenceExpression Allele)

# Notes

- Ambiguous ranges are of the form:
    - `(#_#)_(#_#)`
    - `(?_#)_(#_?)`
    - `(?_#)_#`
    - `#_(#_?)`
- We do not normalize any ambiguous ranges
- We do not change the molecular context for ambiguous ranges.
- The `/to_vrs` endpoint uses the default mode for HGVS deletions and duplications.
- The `/normalize` endpoint uses the default mode for HGVS deletions and duplications if a mode is not set.
