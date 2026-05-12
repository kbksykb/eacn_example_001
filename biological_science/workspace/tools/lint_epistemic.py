#!/usr/bin/env python3
"""
lint_epistemic.py — Pre-PDF-build epistemic lint pass for the REAL/RareShield manuscript.

Checks every *.tex and *.md file under workspace/manuscript/, workspace/handoffs/, and sections/
against a catalog of red-flag phrases, one per trap from the 8-trap checklist
(Supp Note S0, \\label{sec:eightTraps}).

Non-zero hits = CI fail. Output is a list of (file, line, trap, matched phrase, suggested fix).

Each trap has two kinds of patterns:
  - RED_FLAG: phrases that almost always indicate the trap
  - GUARDED_CONTEXT: phrases that neutralize the red flag when co-located (within +/- N lines)

Usage:
    python lint_epistemic.py <path-to-manuscript-root>
    # returns non-zero exit code on any red-flag hit that isn't guarded

Author: Philosophy Agent (agent-mozur8ub) for BioSci CI integration.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path
from dataclasses import dataclass
from typing import Iterable


# Traps per workspace/handoffs/philo_supp_note_S0.tex (Supp Note S0).
# Each rule: (trap_id, human_label, red_flag_regex, guard_regex_or_None, suggested_fix).
RULES: list[tuple[str, str, str, str | None, str]] = [
    # S0.1 Label leakage
    (
        "S0.1",
        "Label leakage",
        r"\b(marker|celltype|cell[- ]type)\s+(gene|label|panel)s?\s+(for|in)\s+(detection|integration|seeding|witnessing)\b",
        r"(out[- ]of[- ]band|verification[- ]only|post[- ]hoc|interpretation[- ]only|not\s+used\s+for\s+detection)",
        "Markers must be out-of-band (verification/interpretation only). Restate: "
        "\"marker panels are used only for post-hoc interpretation of flagged motifs, "
        "not for detection or integration.\"",
    ),
    # S0.2 Smoothness-as-truth
    (
        "S0.2",
        "Smoothness-as-truth",
        r"\b(smooth(ness)?|regulari[sz]ation)\b[^.]{0,80}\b(is|implies|proves|means)\b[^.]{0,80}\b(biolog|correct|preserv)",
        r"(over[- ]smoothing|smoothness\s+is\s+a\s+regulari|smoothness\s+penalty|P_n|penalis(ing|e)\s+over[- ]smoothing)",
        "Smoothness is a regularizer, not ground truth. Rephrase to reference the "
        "over-smoothing penalty ($P_n$) or the DTM-PH / spectral-gap invariants that "
        "detect destroyed small features.",
    ),
    # S0.3 Abundance invariance
    (
        "S0.3",
        "Abundance invariance",
        r"\b(ARI|NMI|ASW|silhouette|scIB(?!-U))\b[^.]{0,60}\b(preserve|retain|capture|detect)\s+(rare|minority)",
        None,
        "Volume-weighted label-based scores do NOT degrade gracefully for rare cells. "
        "Remove the claim or restrict to abundant cell types. Use LRS(w) for rare-cell claims.",
    ),
    # S0.4 Benchmark realism
    (
        "S0.4",
        "Benchmark realism",
        r"\b(pan[- ]cancer|atlas[- ]scale|Kang\s*2024)\b[^.]{0,120}\b(safe|reliable|validated|certif(ied|y)|preserv(ed|ation))\b",
        r"(within\s+the\s+envelope|above[- ]floor|\\Cref\{thm:minimax\}|\\Cref\{thm:S1-minimax\}|operating\s+envelope|\\pi[\\\s]*Delta\^?2|1\.5.{0,10}10\^?\{?-4)",
        "Pan-cancer-scale claims must be qualified by the Theorem 1 minimax envelope "
        "(detectable regime \\pi Delta^2 >= 1.5e-4 at Kang depth). Add the envelope "
        "caveat or restrict the claim to the above-floor stratum.",
    ),
    # S0.5 Annotation completeness
    (
        "S0.5",
        "Annotation completeness",
        r"\b(Kang|HLCA|Tabula\s+Sapiens)\b[^.]{0,80}\b(ground\s+truth|complete(ly)?\s+annotat|comprehensive\s+(annotation|dictionary))",
        r"(do\s+not\s+treat|not\s+as\s+ground\s+truth|label[- ]free|incomplete\s+dictionary|historical\s+.+\s+unknown)",
        "No published atlas is a complete cell-type dictionary (see AS DC / LAMP3+ / "
        "FOLR2+ historical cases). Restate: \"we do not treat the published annotation "
        "as ground truth; REAL detection runs label-free.\"",
    ),
    # S0.6 Post-hoc detection != solution
    (
        "S0.6",
        "Post-hoc detection != solution",
        r"\b(recover(ed|ing|y)|restore[ds]?)\b[^.]{0,80}\b(equal[sl]?|same|equivalent)\b[^.]{0,40}\b(preserv|never\s+erased)",
        r"(post[- ]hoc|CellANOVA|counterfactual\s+is\s+unobserved|not\s+equivalent\s+to\s+a\s+motif\s+never\s+erased)",
        "Post-hoc recovery is not equivalent to prevention. Distinguish REAL (diagnosis) "
        "from RareShield (protection); do not conflate.",
    ),
    # S0.7 Evaluation circularity
    (
        "S0.7",
        "Evaluation circularity",
        r"\b(validat(ed|ion)|benchmark)\b[^.]{0,80}\bon\b[^.]{0,60}\b(integrated|post[- ]integration)\b[^.]{0,80}\b(atlas|annotation)",
        r"(Baron|Travaglini\s*2020|FACS[- ]sort|manually?\s+annotated|pre[- ]integration\s+annotation|independent\s+of\s+batch[- ]integration)",
        "Avoid circularity: at least one dataset in the evaluation panel must be "
        "annotated independently of batch integration (Baron pancreas FACS+manual; "
        "Travaglini 2020 lung).",
    ),
    # S0.8 Failure-to-reject as success
    (
        "S0.8",
        "Failure-to-reject as success",
        r"\b(no|zero|low)\s+(metric\s+drop|score\s+drop|LRS\s*loss|loss)\b[^.]{0,80}\b(safe|proves?|implies|means?|shows?)\b[^.]{0,60}\b(preserv|nothing\s+lost|integrity|discover)",
        r"(within\s+the\s+calibrated\s+sensitivity|above\s+the\s+power\s+floor|\\Cref\{thm:minimax\}|failure[- ]to[- ]reject|not\s+a\s+guarantee|not\s+prove\s+absence)",
        "A PASS (low |W| / low LRS-loss) is evidence of no detectable loss above the "
        "calibrated sensitivity, NOT evidence of preservation. Explicitly cite "
        "\\Cref{thm:minimax} operating envelope.",
    ),
]


@dataclass
class Hit:
    file: Path
    line: int
    trap: str
    label: str
    matched: str
    context: str
    fix: str


def compile_rules() -> list[tuple[str, str, re.Pattern, re.Pattern | None, str]]:
    compiled: list[tuple[str, str, re.Pattern, re.Pattern | None, str]] = []
    for trap_id, label, red_flag, guard, fix in RULES:
        compiled.append(
            (
                trap_id,
                label,
                re.compile(red_flag, re.IGNORECASE),
                re.compile(guard, re.IGNORECASE) if guard else None,
                fix,
            )
        )
    return compiled


def check_file(path: Path, rules, context_lines: int = 4) -> list[Hit]:
    try:
        text = path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return []
    lines = text.splitlines()
    hits: list[Hit] = []
    for i, line in enumerate(lines):
        stripped = line.strip()
        if not stripped or stripped.startswith("%"):
            continue
        for trap_id, label, red, guard, fix in rules:
            m = red.search(line)
            if not m:
                continue
            start = max(0, i - context_lines)
            end = min(len(lines), i + context_lines + 1)
            window = "\n".join(lines[start:end])
            if guard is not None and guard.search(window):
                continue
            hits.append(
                Hit(
                    file=path,
                    line=i + 1,
                    trap=trap_id,
                    label=label,
                    matched=m.group(0).strip(),
                    context=line.strip()[:200],
                    fix=fix,
                )
            )
    return hits


def iter_target_files(root: Path) -> Iterable[Path]:
    extensions = {".tex", ".md"}
    skip_segments = {".git", "node_modules", "__pycache__", "figures"}
    for p in root.rglob("*"):
        if p.is_dir():
            continue
        if p.suffix not in extensions:
            continue
        if any(seg in skip_segments for seg in p.parts):
            continue
        yield p


def render_hits(hits: list[Hit]) -> str:
    if not hits:
        return "lint:epistemic — OK (no ungarded red-flag phrases found)\n"
    out = [
        "lint:epistemic — FAIL",
        f"{len(hits)} red-flag hit(s) without a matching guard. Fix or add a guard "
        "within {+/-4} lines. See workspace/handoffs/philo_supp_note_S0.tex "
        "(sec:eightTraps) for trap definitions.",
        "",
    ]
    for h in hits:
        out.append(f"[{h.trap} {h.label}] {h.file}:{h.line}")
        out.append(f"  matched : {h.matched!r}")
        out.append(f"  context : {h.context}")
        out.append(f"  fix     : {h.fix}")
        out.append("")
    return "\n".join(out)


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print("usage: python lint_epistemic.py <path>", file=sys.stderr)
        return 2
    root = Path(argv[1]).resolve()
    if not root.exists():
        print(f"lint:epistemic — path does not exist: {root}", file=sys.stderr)
        return 2
    rules = compile_rules()
    all_hits: list[Hit] = []
    for f in iter_target_files(root):
        all_hits.extend(check_file(f, rules))
    print(render_hits(all_hits))
    return 0 if not all_hits else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
