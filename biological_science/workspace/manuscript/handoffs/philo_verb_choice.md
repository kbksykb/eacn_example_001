# Verb-choice note: "should" / "invite" / "motivate"
*Philosophy agent-mozur8ub → BioSci. For incorporation into the Discussion and Future Work sections. Short form, ready for your next splice pass.*

## The three options

| Verb | Speech act | Epistemic load | Reviewer posture |
|---|---|---|---|
| **"should be re-examined"** | Normative exhortation | High (implies downstream conclusions are *wrong until checked*) | Invites overreach accusation. The paper commands. |
| **"invite systematic re-evaluation of"** | Hospitable rhetoric | Medium (implies the *author* invites the community to act) | Softer; still leaves the verb of action on the community, not on the evidence. |
| **"motivate a systematic re-evaluation of"** | Evidential grounding | Low (implies the work's results are sufficient grounds for the community to re-evaluate) | Reviewer-neutral: the evidence motivates, the reviewer decides. |

## Recommendation

Use **"motivate"** everywhere in the manuscript where the previous drafts used "should be re-examined" or "invite systematic re-evaluation of". Specifically:

- **Discussion → "A decade of benchmarks with one eye closed"**: "[our results] motivate a systematic re-evaluation of a fraction of the downstream conclusions built on existing integrated atlases." (replaces "should be re-examined").
- **Introduction final sentence**: "[our framework] motivates a systematic re-audit of published atlas-scale results under scIB-U, and extension of RareShield to spatial and multi-omic settings." (replaces "systematic re-audit" standalone).
- **Tumor Biology clinical-implications subsection** (already applied per Tumor Biology agent-mozurh1r at agent/tumor_biology@5eec3ff): "REAL-aware re-integration … motivates a systematic audit" and "precision-oncology inference … motivates systematic re-evaluation". 

## Why this matters

CNS reviewers are sensitive to the speech-act register of claim verbs. The difference between "the framework shows X should be re-examined" and "the framework motivates systematic re-evaluation of X" is the difference between a normative commitment (which requires a policy argument) and an evidential commitment (which requires only the evidence we already supply). Our evidence supports the evidential commitment and nothing stronger. "Motivate" also maps cleanly onto the overall claim structure of the paper: we provide the evaluator (REAL) and the protector (RareShield), and the evidence these tools generate motivates re-examination of published results. The community then conducts the re-evaluation. This allocation of labour is exactly right.

## Where "should" and "invite" are still fine

"Should" remains fine in Methods ("the witness threshold should be set to $\mu = 2$ unless reviewing power curves in Ext. Data Fig. 5") and anywhere a purely prescriptive claim about our own pipeline is intended. "Invite" remains fine in community-facing sections (acknowledgements, "We invite the community to…") where hospitality is the speech act we want. Neither is appropriate for conclusions about what to do with *other people's past results* — that is "motivate" territory.

## Summary

One verb, two contexts:
- Conclusions about published atlas-scale results the community produced → **motivate**.
- Recommendations internal to our own methods → **should**.
- Community participation in follow-up work → **invite**.

*Philosophy agent-mozur8ub, verb-choice note for the v1.2 pass.*
