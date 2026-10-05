# Reflective projective billiards: manuscripts and verification code

This repository collects three standalone manuscripts by Jorge Lucas González, together with the computational checkers associated with their arguments.

## Manuscripts

| Manuscript | Status | Files |
|---|---|---|
| [Reflective projective billiards in the plane: odd periods, period five, and metric rigidity](03-projective-billiards-planar/README.md) | Preprint; [arXiv:2609.28817](https://arxiv.org/abs/2609.28817) | [PDF](03-projective-billiards-planar/main.pdf) · [LaTeX](03-projective-billiards-planar/main.tex) |
| [Reflective projective billiards in higher dimensions: monodromy, threshold phenomena, and explicit polytopes](04-projective-billiards-higher-dimensional/README.md) | Manuscript; not peer-reviewed | [PDF](04-projective-billiards-higher-dimensional/main.pdf) · [LaTeX](04-projective-billiards-higher-dimensional/main.tex) |
| [The period spectrum of reflective projective billiards: joins, suspensions, and metric obstructions](05-projective-billiards-period-spectrum/README.md) | Manuscript; not peer-reviewed | [PDF](05-projective-billiards-period-spectrum/main.pdf) · [LaTeX](05-projective-billiards-period-spectrum/main.tex) |

The planar directory contains an updated, self-contained manuscript snapshot compiled on 3 October 2026. It may differ from the version in the arXiv record; use arXiv for that record's version history.

## Reproducibility

The `verificadores/` folders contain scripts for specific algebraic and combinatorial checks discussed in the manuscripts. The scripts supplement the written proofs; their scope is described in each paper's README and source comments.

Install the external Python dependencies with:

```bash
python -m pip install -r requirements.txt
```

The other imports use Python's standard library. No journal publication or peer-review status is implied for the two manuscripts marked as such above.
