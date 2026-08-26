# Working in this directory

This is the LaTeX source for the ModelSEED 2026 update manuscript, targeting
**Nucleic Acids Research, Database Issue**.

**Before changing anything about formatting, length, figures, references, or the
abstract, read [`NAR_REQUIREMENTS.md`](NAR_REQUIREMENTS.md).** It carries NAR's
and OUP's actual rules with citations and verbatim quotes — page budget, figure
formats and dpi, the mandatory graphical abstract, reference style, and the
submission checklist. It is the authority; do not restate requirements from
memory, and do not add one to that file without a source.

Structure and build instructions are in [`README.md`](README.md). In short:
`main.tex` holds all formatting and `\input`s the prose; every section lives in
`sections/`. Do not put prose in `main.tex`.
