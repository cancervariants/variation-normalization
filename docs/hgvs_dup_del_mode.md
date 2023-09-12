# HGVS Dup Del Mode

This mode helps us interpret deletions and duplications that are represented as HGVS expressions.\
The mode can be set to `default`, `copy_number_count`, `copy_number_change`, or `allele`


## Default Characteristics

- if baseline_copies is not set and endpoints are ambiguous:
    - copy_number_change
    - if copy_change not provided:
        - copy_change = `efo:0030067` (loss) if del, `efo:0030070` (gain) if dup
- elif baseline_copies is provided:
    - copy_number_count
    - copies are baseline_copies + 1 for dup, baseline_copies - 1 for del
  else:
    - allele

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
