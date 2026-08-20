"""Features that legitimately do not satisfy the usual arithmetic.

Three unrelated things in real annotation produce a CDS whose length is not a
whole number of codons, or which has no translation to check against:

* **pseudogenes** — 526 of 5,808 loci in *M. aeruginosa*. A damaged copy carries
  no `/translation`.
* **partial features** — RefSeq's `partial=true` with `start_range`/`end_range`,
  and AUGUSTUS predictions that run off the end of a contig. The coding sequence
  is truncated by the assembly, not by an error.
* **`/artificial_location`** — 78 records in *F. graminearum*. NCBI has adjusted
  the coordinates to work around a known assembly or sequencing problem, so the
  span deliberately does not translate cleanly.

Every completeness check that measures a protein against its coding sequence
must skip all three, and it must skip them by calling one predicate. Written as a
single module now rather than three times at milestone 8, because three copies
drift and the third one is always the one nobody updates.
"""

from __future__ import annotations

from dataclasses import dataclass

#: Attribute keys whose presence means "the arithmetic is not expected to hold".
#: Each is the verbatim qualifier or attribute name as it appears in input.
NON_CONFORMING_FLAGS: tuple[str, ...] = (
    "pseudo",  # GenBank /pseudo and /pseudogene; RefSeq GFF3 pseudo=true
    "pseudogene",
    "partial",  # RefSeq partial=true; set by the GenBank tokenizer for <,>
    "artificial_location",  # NCBI has adjusted the coordinates deliberately
    "ribosomal_slippage",  # a base is translated twice, so length/3 is wrong
    "transl_except",  # a codon is translated against the standard table
    "exception",  # /exception="rearrangement required for product"
)


@dataclass(frozen=True)
class Exemption:
    """Why a feature is exempt, so a report can say which reason applied."""

    flag: str
    value: str | None


def exemption_for(attributes: dict[str, list[str]]) -> Exemption | None:
    """Return the reason this feature is exempt from length arithmetic, or None.

    Matching is case-insensitive because the same concept is spelled differently
    across formats, and presence-only because these flags are frequently valueless
    in GenBank (`/pseudo` has no value at all) while GFF3 writes `pseudo=true`.
    """
    lowered = {key.lower(): values for key, values in attributes.items()}
    for flag in NON_CONFORMING_FLAGS:
        if flag in lowered:
            values = lowered[flag]
            value = values[0] if values else None
            # An explicit false is not an exemption. RefSeq never writes it, but
            # a hand-edited file might, and reading `partial=false` as "exempt"
            # would silently disable the check on a feature that should have it.
            if value is not None and value.strip().lower() in ("false", "0", "no"):
                continue
            return Exemption(flag=flag, value=value)
    return None


def is_exempt(attributes: dict[str, list[str]]) -> bool:
    return exemption_for(attributes) is not None
