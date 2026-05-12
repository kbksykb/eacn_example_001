# Problem Analysis: The Triple Absence and the Self-Closing Loop

> The following is one possible path of analysis for this problem, not the only one. It attempts to explain why a problem that has been under attention for at least five years has still not been solved, and why it is, in theory, solvable.

## What is truly missing is not one thing, but three

On the surface, what is missing is an algorithm. But if you dig deeper, you find that three things are missing at the same time:

1. **No ground truth** — we do not know where the unknown rare subpopulations are supposed to be
2. **No evaluation framework** — we cannot detect whether the integration process has wiped out a rare subpopulation
3. **No applicable algorithm** — even if we knew something was lost, we have no tool to prevent it

The missing algorithm is only the last link in the chain.

## Why the triple absence is interlocked

Without ground truth, the evaluation framework cannot be built — you do not know how to define "successfully preserved."

Without an evaluation framework, no one knows that existing algorithms fail on rare subpopulations — the failure signal does not exist.

Without a failure signal, no one sets out to develop an algorithm specifically to protect rare subpopulations — there is no motivation.

## The self-closing loop

The triple absence forms a closed loop:

```
No ground truth
→ No evaluation framework
→ No failure signal
→ No motivation to develop an algorithm
→ No new discoveries
→ No ground truth
```

Every step of this loop is waiting for the previous step to be solved first. No single link can break through on its own.

## Why existing methods have not touched this problem

Methods such as scDML and ResPAN claim to preserve rare subpopulations, but their validation relies on two kinds of ground truth:

- **Simulated data**: artificially construct a rare population, then check after integration whether it is still there. This is circular reasoning — you know it exists because you built it.
- **Existing annotations**: use the cell annotations provided by the dataset authors as the reference answer. But if a subpopulation does not exist in the annotation at all, these methods are entirely unable to evaluate whether it has been preserved.

Neither of these validation approaches can answer the core question: **for a rare subpopulation that is unknown in advance, has integration wiped it out?**

## Why no one has solved it so far

It is not because the problem is unimportant, nor because no one has noticed it.

It is because the entry point of this loop has been locked shut by the absence of ground truth. Everyone defaults to the assumption that ground truth must be obtained from the outside — either through simulated data or through experimental validation. And external ground truth is extremely costly to obtain, and each round can only validate one specific subpopulation; it does not generalize.

This default assumption is itself part of the problem.

## Why this problem is, in theory, solvable

The loop is closed, but it is not unbreakable. Because there is one key fact:

**Normal batch correction and the elimination of a rare subpopulation leave different structural patterns in the data.**

Normal batch correction: cells of the same type are merged across batches; the cluster is still there, only now its members come from more batches. This is convergence.

A rare subpopulation being eliminated: the members of a cluster are scattered into multiple unrelated large clusters; the cluster disappears. This is dispersion.

One is convergence, the other is dispersion. These two patterns are structurally different; they are not a faint difference buried in noise.

This means: without external ground truth, it is still possible to distinguish, from the structural changes in the data itself before and after integration, which changes are normal and which are pathological. The signal exists in the data; the problem is that, to this day, no one has designed the framework to extract it.
