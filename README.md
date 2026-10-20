# Fréchet mean set inference

Reproducibility code, fixed numerical references and figure/table generators
for confidence sets for Fréchet means.

| Output | Data in `reference/` | Generated file |
|---|---|---|
| Coverage comparisons | `coverage.csv` | `figure-coverage.tex` |
| Finite-ranking set cardinalities | `ranking.csv` | `table-ranking.tex` |
| Circle mean angular lengths | `circle.csv` | `table-circle.tex` |

Keep the complete repository layout and run commands from its root.
The numerical implementations are in `euclid/`, `ranking/`, `circle/` and
`common/`; the fixed environment is in `environment/`. The CSVs and plotting
templates in `reference/` are sufficient to generate the bundled figure and
tables. Their mapping and hashes are recorded in
[`reference/manifest.json`](reference/manifest.json).

## Generate the figure and tables

Use Python 3.9 or later:

```sh
python3 verify_release.py
python3 export_results.py export --output build/results
```

The three files above and `plot-style.tex` are written to `build/results/`.
For LaTeX rendering, load `amssymb` and `plot-style.tex` in the preamble, then
load the generated figure and table fragments in the document body:

```latex
% Preamble
\usepackage{amssymb}
\input{build/results/plot-style.tex}
% Document body
\input{build/results/figure-coverage.tex}
\input{build/results/table-ranking.tex}
\input{build/results/table-circle.tex}
```

Each fragment includes its figure or table environment. Circle angles are in
degrees and include all connected components for Global, Local and EL.
Existing outputs are reused only if identical; otherwise choose a new
directory under `build/`.

## Recompute the numerical results

Use the recorded Linux/WSL environment in
[`environment/runtime.dcf`](environment/runtime.dcf). `setup.R install` restores
the packages fixed in [`environment/packages.dcf`](environment/packages.dcf).
Full recomputation produces large certificate archives; allow ample disk space.

```sh
Rscript --vanilla setup.R install
Rscript --vanilla run.R run runs/reproduction 1000 32
python3 verify_release.py --recomputed runs/reproduction
```

This is one run: 1,000 replications per condition and at most 32 workers.
It generates observations and 1,999 multiplier draws per dataset with
`settings.dcf` unchanged, then writes the figure, tables and style to
`runs/reproduction/results/`. Keep the automatically generated sibling directory
`runs/reproduction-ranking-candidates/` with the run for reference verification.

Repeat the same R command to resume an interrupted run. Changed execution source,
settings or worker count require a new output directory. Unresolved results at
the fixed budget exit nonzero as `CERTIFICATION_INCOMPLETE`; valid bounds are
retained.

## Verification

`verify_release.py` checks the bundled files and generated outputs.
Use `--recomputed RUN` to compare a completed run with the fixed reference.

The [numerical contract](OPTIMIZATION_AND_CERTIFICATION.md) specifies certification,
rounding and statistical limits. Original observations, certificate archives
and private execution records are not included. No new full numerical run or
unconditional cross-platform bit equality is claimed. Retain the release
commit or lock hash independently. License: [MIT](LICENSE).
