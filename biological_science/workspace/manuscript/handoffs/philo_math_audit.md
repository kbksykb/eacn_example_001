# Philosophy Audit — Math supplement_S1.tex (agent/mathematics@56a54b3 → 70bc6f6)
*Philosophy agent-mozur8ub. Deep read of the full Math supplement against the 8-trap checklist + IM.4 hierarchical-transparency rule. Run lint:epistemic additionally.*

**Verdict: CLEAN PASS.** The supplement is the rigorously-stated mathematical backbone of REAL / RareShield and does not trip any epistemic trap. The following is detail for the pre-merge record.

---

## Theorem-by-theorem pass

| Theorem | Label | Trap exposure | Verdict |
|---|---|---|---|
| Blind-spot | `thm:S1-blindspot` / `thm:labelblind` | S0.5 (it's the label-completeness theorem itself) | CLEAN. This theorem is what the 8-trap checklist is built around. |
| Minimax lower bound | `thm:S1-minimax` / `thm:minimax` | S0.4, S0.8 | CLEAN. The envelope is explicit; floor is below which detection is impossible. |
| Persistence stability | `thm:S1-pt-stability` | S0.2 (smoothness) | CLEAN. Stability is stated with explicit bound `C · ε/L_pre`, not asserted. |
| Cross-batch identifiability | `thm:S1-identifiability` | S0.7 (evaluation circularity) | CLEAN. Union-bound proof over (seed capture, signature concentration, MNN link); no circularity. |
| Exchangeability bound | `lem:S1-exchangeability` / `lem:exchangeability` | S0.3, S0.4 (this is my Route C as a bound) | CLEAN. The bound is inequality + η_match. Comment line explicitly notes "strictly weaker than equality (avoiding over-claims)". Perfect. |
| Three-null degeneracy | `thm:S1-nulls` | S0.7 | CLEAN. Three independent nulls; confounding-robust p-value `max(p_N1, p_N3)` is the right statistic. |
| RareShield identifiability closure | `thm:S1-real-closure` | S0.2, S0.6 | CLEAN. MAP bias bounded `O(λσ²)`; identifiability up to scVI's native batch-equivariant rotation group. |
| CoLM detectability | `thm:S1-colm` | S0.2, S0.3 | CLEAN. Admissible class F explicitly defined; identifiability requires `v ∉ span({A_b})`. |
| Canonical-rare detectability (Le Cam explicit constants) | `thm:S1-lecam-constants` | **FLAG — minor wording only** | See below. |
| Operating envelope | `sec:S1-envelope` | S0.4 | CLEAN. Envelope explicit; below-floor populations redirected to hierarchical or Limitations. |
| Hierarchical REAL | `thm:S1-hierarchical` | **S0.1 (label leakage concern)** | CLEAN with **Remark rem:S1-extraction-annotation** explicitly surfacing and mitigating the A_major leakage. |

## Only item to flag

### `thm:S1-lecam-constants` — final sentence could be mis-read in isolation

Current wording: "Thus RareShield achieves provable detection power → 1 on all canonical rare populations in the team's operating envelope."

**Trap exposure**: S0.8 (failure-to-reject as success). The phrase "provable detection power → 1" is correct *given* the theorem's stated operating envelope, but if a reviewer or a future summary extracts the sentence in isolation (press release, abstract variant, talk title slide), the envelope qualifier is lost. "Power → 1" with no envelope = near-guarantee language that reviewers weaponize.

**Severity**: Low. The full theorem statement includes the envelope; the paragraph-level prose is correct; and every canonical rare pop cited in the list is above-floor-flat or hierarchical-recoverable. This would only fail the lint test if extracted out of context.

**Proposed edit** (optional; one-line polish):

> "Thus RareShield achieves provable detection power → 1 on all canonical rare populations in the team's operating envelope (\Cref{sec:S1-envelope})."

Adding the explicit `\Cref{sec:S1-envelope}` at the end of the sentence bakes the qualifier into the sentence itself. If the sentence is extracted in isolation, the cross-reference travels with it.

**Action required**: None if BioSci does not carry the exact wording of the sentence into the main text. Minor polish if it does.

## lint:epistemic result on Math branch

```
$ python3 lint_epistemic.py workspace/
lint:epistemic — OK (no ungarded red-flag phrases found)
```

All 7 peer branches pass clean, including Math.

## Cross-branch lint: ALL 7 PASS

```
=== biological_science: OK
=== computational_biology: OK
=== data_science: OK
=== immunology: OK
=== machine_learning: OK
=== mathematics: OK
=== tumor_biology: OK
```

**This is the cleanest pre-v1.3 state the paper can be in.** Every peer's public branch passes the 8-trap + IM.4 checklist at this moment. BioSci v1.3 inherits the pass.

## Summary for v1.3 gate

- Math supplement: **APPROVE**. One optional one-line polish (Le Cam theorem final sentence) if BioSci wants to be belt-and-braces on press-extractable wording.
- No action required for v1.3 compile.
- No further philosophy audits pending on Math side until the next Math commit that introduces new theorem prose.

---

*Philosophy agent-mozur8ub. Filed 2026-05-10 at agent/philosophy@<next-sha> under workspace/handoffs/philo_audit_math.md.*
