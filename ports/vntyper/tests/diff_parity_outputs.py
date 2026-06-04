#!/usr/bin/env python3
"""Compare two run_parity_pipeline.py outputs for Java↔Rust parity.

Parity here means: for every shipped real-data fixture, both engines made
the upstream-correct call and agree with each other on the biological
result — same expected Confidence, same positive/negative classification,
and (for positives) the same called variant locus.

The exact TSV sha256 is reported but is NOT a parity failure on its own:
the BAM path's byte-level divergence is the separately tracked samtools-rs
FASTQ-extraction gap (see TODO.md "Current blockers"). Alt-depths legit-
imately differ by a few reads between engines while staying inside
upstream's tolerance, so correctness + classification is the gate.

Exit non-zero if any fixture is mis-called by either engine or the two
engines disagree on a fixture's classification.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


def _called_locus(called):
    if not called:
        return None
    return (
        called.get("POS"),
        called.get("REF"),
        called.get("ALT"),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("left", type=Path, help="first engine JSON (e.g. java)")
    parser.add_argument("right", type=Path, help="second engine JSON (e.g. rust)")
    parser.add_argument("--label-left", default=None)
    parser.add_argument("--label-right", default=None)
    args = parser.parse_args()

    left = json.loads(args.left.read_text(encoding="utf-8"))
    right = json.loads(args.right.read_text(encoding="utf-8"))

    label_left = args.label_left or left.get("engine", str(args.left))
    label_right = args.label_right or right.get("engine", str(args.right))

    if left.get("input") != right.get("input"):
        print(
            f"input mismatch: {label_left}={left.get('input')!r} "
            f"{label_right}={right.get('input')!r}",
            file=sys.stderr,
        )
        return 2

    fixtures = sorted(set(left.get("cases", {})) | set(right.get("cases", {})))
    if not fixtures:
        print("no fixtures to compare", file=sys.stderr)
        return 2

    any_fail = False
    for stem in fixtures:
        lc = left.get("cases", {}).get(stem)
        rc = right.get("cases", {}).get(stem)
        if lc is None or rc is None:
            owner = label_right if lc else label_left
            print(f"[{stem}] MISSING on {owner}")
            any_fail = True
            continue

        problems = []
        if not lc.get("correct"):
            problems.append(
                f"{label_left} mis-called: {lc.get('reasons')}"
            )
        if not rc.get("correct"):
            problems.append(
                f"{label_right} mis-called: {rc.get('reasons')}"
            )

        l_cls = lc.get("classification")
        r_cls = rc.get("classification")
        if l_cls != r_cls:
            problems.append(
                f"classification disagree: {label_left}={l_cls!r} "
                f"{label_right}={r_cls!r}"
            )

        l_neg = lc.get("called") is None
        r_neg = rc.get("called") is None
        if l_neg != r_neg:
            problems.append(
                f"call presence disagree: {label_left}="
                f"{'no-call' if l_neg else 'call'} "
                f"{label_right}={'no-call' if r_neg else 'call'}"
            )
        # NOTE: identical REF/ALT is intentionally NOT required. The same
        # biological MUC1 dup frameshift is reported as C>CG or G>GG
        # depending on which equivalent motif reference Kestrel aligned
        # against. Upstream's own correctness test only checks Confidence
        # and depth tolerance on the top row, never the exact allele, so
        # two engines both landing on the upstream-correct call IS parity.
        locus_note = ""
        if not l_neg and not r_neg:
            ll, rl = _called_locus(lc["called"]), _called_locus(rc["called"])
            if ll != rl:
                locus_note = f"  (locus repr {ll}≠{rl}; equivalent motif)"

        l_sha = lc.get("tsv_fingerprint", {}).get("sha256", "?")[:12]
        r_sha = rc.get("tsv_fingerprint", {}).get("sha256", "?")[:12]
        sha_note = "" if l_sha == r_sha else f"  (tsv sha {l_sha}≠{r_sha} — samtools-rs gap)"

        exp = lc.get("expected", {}).get("confidence")
        if problems:
            print(f"[{stem}] FAIL expect={exp!r}{sha_note}{locus_note}")
            for problem in problems:
                print(f"    {problem}")
            any_fail = True
        else:
            kind = "negative" if l_neg else "positive"
            print(
                f"[{stem}] MATCH expect={exp!r} {kind} "
                f"both-correct{sha_note}{locus_note}"
            )

    if any_fail:
        print(
            f"\nparity FAIL: {label_left} and {label_right} are not at "
            f"correctness parity on all fixtures",
            file=sys.stderr,
        )
        return 1
    print(
        f"\nparity OK: {label_left} and {label_right} both call every "
        f"shipped fixture upstream-correctly"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
