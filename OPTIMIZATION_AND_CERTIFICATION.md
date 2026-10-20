# Numerical contract and coverage decisions

This contract covers G/L and comparator methods for Euclidean data, the three finite-ranking metrics, the continuous ranking graph, and the Circle. The reference-output mapping is fixed in `reference/manifest.json`; execution instructions are in `README.md`.

## Fixed inputs

The numerical problem conditions on the observations, observed tuning, and a finite set of multiplier draws. The RNG uses L'Ecuyer-CMRG with Inversion and master seed 202610074. Streams for all conditions are reserved in advance, so pilot runs and parallel scheduling do not change samples. Ranking and graph cases share observations and the same 24-by-B Gaussian array.

Observed tuning uses the minimum empirical risk:

$$
\varepsilon_n=F_n^\star\sqrt{\log\log n/n},\qquad \tau_n=n^{1/5}.
$$

Population means are used only to evaluate coverage and decide when that evaluation can stop. Population answers are not supplied to minimization, feasible candidates, interval upper bounds, or observed tuning.

## Inversion and zero variance

$$
F_n(q)=n^{-1}\sum_i d^2(Y_i,q),\qquad
\Delta_n(q,r)=F_n(q)-F_n(r),
$$

$$
\widehat\sigma_n^2(q,r)=n^{-1}\sum_i
\{d^2(Y_i,q)-d^2(Y_i,r)-\Delta_n(q,r)\}^2.
$$

The variance denominator is n. Over comparisons with positive variance, the inverted statistic is

$$
T_n(q)=\max\{0,\sup_r\sqrt n\Delta_n(q,r)/\widehat\sigma_n(q,r)\}.
$$

A zero-variance comparison with a positive risk difference excludes the candidate at every finite cutoff. Inclusion allows equality; exclusion requires a strict inequality. For a nonnegative cutoff c, each comparison is checked through

$$
\Delta_n(q,r)\le0\quad\text{or}\quad
c^2\widehat\sigma_n^2(q,r)-n\Delta_n(q,r)^2\ge0.
$$

Continuous domains require a cover of all comparisons. A finite grid alone does not certify inclusion. Small positive variances are not set to zero by an arbitrary threshold, and no ridge is added.

## Calibration domains and bounds

G uses all ordered comparisons with positive variance. L additionally requires

$$
F_n(q)-F_n^\star\le\varepsilon_n,\qquad
\sqrt n\Delta_n(q,r)\ge-\tau_n\widehat\sigma_n(q,r).
$$

For each draw, feasible candidates satisfying the original constraints provide a lower bound, and a cover of the full domain provides an upper bound:

$$
\underline S_{n,b}^{\#}\le S_{n,b}^{\#}\le\overline S_{n,b}^{\#}
\quad\Longrightarrow\quad
\underline c_n:=\underline S_{n,(1900)}^{\#}
\le c_n\le\overline S_{n,(1900)}^{\#}=:\overline c_n.
$$

B=1999, and calibration uses the 1900th order statistic without interpolation. Set monotonicity gives

$$
\mathcal C_n(\underline c_n)\subseteq\mathcal C_n(c_n)
\subseteq\mathcal C_n(\overline c_n).
$$

The decision is IN if the lower-cutoff set contains every population mean, OUT if at least one population mean lies outside the upper-cutoff set, and UNRESOLVED otherwise. Coverage concerns simultaneous inclusion of all population means; no single mean is selected. A certified IN/OUT decision agrees with the ideal fixed-B procedure even when the cutoff enclosure remains wide. Attaining cutoff width .01 is not required to stop a resolved coverage evaluation. Numerical stability or termination of a local optimizer alone does not prove a global error bound.

## Geometry-specific treatment

| Geometry | Calibration and membership |
|---|---|
| Euclidean | Projection norm onto the centered sample's column space; verified Gram calculations and an exact rational fallback |
| Finite ranking | Exhaustive enumeration of 552 nontrivial ordered comparisons among 24 candidates; checks on integer ranges and ambiguous selector conditions |
| Continuous graph | All comparisons across 36 edges; piecewise quadratic risks, cone upper bounds, and interval subdivision; evaluation of comparator edges for the full population mean set |
| Circle | Charts determined by sample antipodes, prefix moments, verified profiling, and all-chart comparisons |

Circle population means and mathematical pi are enclosed by intervals. A single floating-point representative does not replace a membership certificate for a population mean that is not representable in binary64. Upper bounds are retained near uncertain roots and zero-variance configurations. Exact integer calculations are distinguished from calculations requiring outward rounding.

For Euclidean delta=0, the sample is represented in planar N2(0,I2) coordinates. All three original latent Gaussian columns are generated before selecting two columns, preserving the RNG sequence. This avoids introducing a spurious residual direction through a three-dimensional floating-point projection. The generator for delta>0 is unchanged. Candidates outside the plane are excluded by the support condition.

With covariance denominator n, the delta=0 Gaussian model satisfies

$$
n\bar Y_n^{\mathrm T}\widehat\Sigma_n^+\bar Y_n
\sim\frac{2n}{n-2}F_{2,n-2}.
$$

The chi-square Wald method is an asymptotic comparator; the Gaussian F method is a finite-sample reference. Circle EL retains all components of ordinary score inversion and has its own calibration, separate from G/L. Circle KDE Wald is a floating-point comparator, distinct from G/L numerical certification.

## Reporting and reproducibility

All R replications are retained. If I replications have certified inclusion and U remain unresolved, the numerical coverage range is [I/R,(I+U)/R], with width U/R. The Monte Carlo interval uses the Clopper--Pearson lower endpoint for I successes and upper endpoint for I+U successes. These are individual 95% intervals, not simultaneous 95% intervals over all conditions. Numerical uncertainty and Monte Carlo sampling error are reported separately.

The decimal .01 audit converts saved binary64 endpoints to exact rational numbers and compares their difference with 1/100. Cutoff precision, numerical resolution of coverage, and nominal statistical coverage of 95% are separate properties. The fixed-B rank argument 1900/2000 requires joint convergence of the evaluation statistic and draws to independent limit variables with a common continuous distribution. It does not imply fixed-B cutoff consistency or finite-sample 95% coverage for every model.

The execution environment, inputs, source, budgets, and output hashes are recorded. A validation PASS does not replace formal verification of the complete upper-bound algorithm or a statistical theorem. Runtime and byte-identical PDFs are not substitutes for scientific reproducibility.

Compiler-setting overrides and custom Makevars are rejected. The permitted compiler configuration and build-rule hashes are included in the frozen runtime and native-cache identity.

## Finite-ranking counts

Calibration checks 24-by-23=552 ordered comparison pairs. The inverted confidence set is a subset of 24 rankings. These two candidate counts are reported separately. G/L invert the same full-comparison statistic; only the Local calibration domain is smaller. Consequently,

$$
c_L\le c_G,\qquad \mathcal C_L\subseteq\mathcal C_G.
$$

For certified candidate indicator bounds g_q^-,g_q^+,l_q^-,l_q^+,

$$
\sum_{q=1}^{24}g_q^-\le\#\mathcal C_G\le\sum_{q=1}^{24}g_q^+,
\qquad
\sum_{q=1}^{24}l_q^-\le\#\mathcal C_L\le\sum_{q=1}^{24}l_q^+.
$$

For each candidate, the implementation checks possible states (G,L) in {(0,0),(1,0),(1,1)} against these bounds. It also encloses counts for the intersection, G only, L only, and the symmetric difference. The L-only count is zero. When every candidate's state is determined, the counts are exact integers. Unresolved count bounds are not claimed to be sharp after accounting for the common cutoffs.

A `*_possible` report means the available bounds do not exclude that outcome; it does not prove realizability by one common pair of cutoffs for all candidates. A `*_certified` report proves set equality or strict shrinkage. `ranking-counts` is separate postprocessing of saved certificates; it changes neither calibration nor the original coverage results.

## Deterministic certification of total Circle angle

`circle/length.hpp` shares full-circle partitions, fixed data, and convexity certificates across G/L/EL. `circle/length.R` refines Global/Local cutoffs using the same 1,999 draws. Global uses an independent ROOT cone envelope over the full domain; Local uses subdivision directed at the relevant order statistic. A maximum over feasible points alone is not an upper bound. Inner cells at the lower cutoff and outer cells at the upper cutoff enclose the ideal-cutoff set. Remaining cutoff-band uncertainty triggers further calibration.

The full circumference is represented by M=2^24 equal-angle cells. If I cells are certified inside and O outside, the reported total-angle interval is given by the exact rational bounds

$$
360I/M\ \le\ \operatorname{angle}(\mathcal C)\ \le\ 360(M-O)/M.
$$

All disconnected components contribute to the total. Unresolved cells and cells enclosing antipodes remain in the upper bound. Cutoff and membership-statistic monotonicity enforce L as a subset of G; contradictory certificates are execution errors.

The direct stopping condition is an angle width no greater than the binary64 setting 0.57 degrees for each sample and method. This is separate from the decimal 1/100 cutoff audit. Means and paired differences G-L, G-EL, and L-EL are aggregated using Fraction arithmetic, with displayed decimal endpoints rounded outward. Public tables and summaries use degrees only. Internal hexadecimal radian bounds are retained for independent checks using true pi.

The generated angle table displays midpoints of certified mean-angle intervals, rounded to one decimal place with ties to even. `export_results.py` computes midpoints and rounding using exact rationals and checks all 36 displayed values: each lies less than 0.3 degrees from both interval endpoints. The certified intervals themselves are retained.

This is a numerical enclosure guarantee conditional on fixed observations, tuning, and draws. It is distinct from population coverage, Monte Carlo error in expected angles, and the accuracy of individual boundary endpoints. If a budget limit leaves a wide interval, its valid bounds and counts are saved.
