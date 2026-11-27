# HGVS deletions and duplications

The `hgvs_dup_del_mode` parameter controls how to interpret HGVS deletions and
duplications.

## Default behavior

| Input | VRS representation |
| --- | --- |
| Ambiguous endpoints without `baseline_copies` | Copy Number Change |
| `baseline_copies` provided | Copy Number Count |
| Unambiguous endpoints without `baseline_copies` | Allele |

For Copy Number Change, `copy_change` defaults to loss for a deletion and gain
for a duplication. For Copy Number Count, the count is `baseline_copies - 1`
for a deletion and `baseline_copies + 1` for a duplication.

The `/to_vrs` endpoint uses the default mode.

The `/normalize` endpoint accepts `default`, `copy_number_count`, `copy_number_change`,
and `allele`.

## Ambiguous ranges

In these supported ranges, `#` stands for a position and `?` marks an unknown bound.

- `(#_#)_(#_#)`
- `(?_#)_(#_?)`
- `(?_#)_#`
- `#_(#_?)`

Ambiguous ranges retain their molecular context and are not normalized.
