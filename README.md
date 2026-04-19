# multicomp_test

MATLAB permutation and false-discovery-rate tools for multiple-comparison testing on 1-D and 2-D data — computing significance while properly controlling family-wise error or false-discovery rate across many bins.

Part of the Prerau Lab MATLAB codebase. Can also be used standalone.

## Two families of tests

| Family | What it controls | When to use |
|---|---|---|
| **Permutation (max-stat)** | Family-wise error rate (FWER) | You want strong control: "no false positive anywhere in the dataset" |
| **FDR (Benjamini-Hochberg)** | False-discovery rate | You can tolerate some false positives in exchange for higher power |

## Functions

| Function | Purpose |
|---|---|
| `permtest` | Permutation test on 1-D grouped data; returns significant bins via max t-stat threshold |
| `permtest2` | Thin 2-D wrapper around `permtest` — reshapes to 2-D, calls the 1-D core, reshapes the result |
| `gpermtest` | "Global" acceptance-bounds permutation test for 1-D data with a custom statistic function |
| `gpermtest2` | 2-D wrapper around `gpermtest` |
| `FDR_1D` | Benjamini-Hochberg FDR control on 1-D t-test results with optional `prctile` / `std` thresholds |
| `FDR_2D` | 2-D wrapper around `FDR_1D` |

Each function runs its own demo with no inputs — e.g., `permtest()` — which is the fastest way to see the output format.

## Quick examples

### 1-D permutation test (max-stat FWER control)

```matlab
% Two groups, 100 bins, 20 trials each
group1 = randn(100, 20);
group2 = randn(100, 20) + [zeros(40,20); ones(20,20)*0.7; zeros(40,20)];

[sigbins, tstat, thresh, perm_tmax] = permtest(group1, group2, 0.05, 10000, true);
% sigbins    : 1x100 logical, true where a significant group difference survives
%              the max-stat permutation threshold
% tstat      : 1x100 observed t-statistic
% thresh     : scalar, the (1-alpha)-percentile of permuted max |t|
% perm_tmax  : 1x10000 null distribution of max |t|
```

### 2-D permutation test

```matlab
% Two groups, 50x50 images, 20 / 33 trials
g1 = randn(50, 50, 20);
g2 = randn(50, 50, 33) + randn(50, 50);

[sigbins_2D, tstat_2D] = permtest2(g1, g2);   % plots a 2-D heat-map with significant regions contoured
```

### FDR control on 1-D data

```matlab
% Returns p-values, thresholded logical mask, and the BH-cutoff
[sigbins, pvals, cutoff] = FDR_1D(group1, group2, 'FDR', 0.05, 'method', 't', 'paired', false);
```

## Entry-point signatures

All functions accept `varargin` and use `inputParser`. Common name-value pairs:

| Name | Default | Description |
|---|---|---|
| `alpha_level` (perm tests) | `0.05` | Significance threshold α |
| `iterations` (perm tests) | `10000` | Number of permutation shuffles |
| `FDR` (FDR tests) | `0.05` | Target false-discovery rate |
| `method` (FDR tests) | `'t'` | `'t'` for t-test, `'sign'` for sign test |
| `paired` | `false` | Paired vs independent two-sample |
| `nonparam` (FDR tests) | `false` | Use nonparametric rank statistic |
| `ploton` | `true` | Plot results |
| `statfcn` (gpermtest) | `@(x) mean(x, 2, 'omitnan')` | Test statistic function |

See `help <function>` at the MATLAB prompt for the full docstring of each.

## Method notes

### Max-stat permutation (`permtest`, `permtest2`)

Computes the observed t-statistic per bin. Under the null, shuffles group labels `iterations` times; for each shuffle, records `max(|t|)` across all bins. The (1-α)-percentile of that null distribution is the threshold. A bin is significant if its observed |t| exceeds the threshold.

This is strong FWER control — the probability of any false positive across the entire bin set is ≤ α.

### Global acceptance bounds (`gpermtest`)

More flexible than `permtest`: you pass a custom statistic function (default: `mean`). Under the null, computes the statistic at each bin for each permutation; the empirical `[α/2, 1-α/2]` percentiles at every bin form the "global acceptance bounds." Observed values outside those bounds are significant.

Use this when you want bounds on, say, the mean difference or a phase statistic, rather than a t-distribution-based test.

### FDR (`FDR_1D`, `FDR_2D`)

Runs a t-test (or sign test) at each bin, collects p-values, applies Benjamini-Hochberg to get the step-up cutoff. Less conservative than max-stat permutation: you lose strong FWER control but get substantially more power for dense-signal datasets.

## Install

```matlab
addpath('/path/to/multicomp_test');
```

## Dependencies

- **MATLAB R2020a or later**
- **Statistics and Machine Learning Toolbox** for `ttest` / `ttest2` inside `FDR_1D`. Permutation tests use `ttest2` too but can be swapped for any pairwise statistic.
- **Parallel Computing Toolbox** (optional) — the permutation tests use `parfor`; without the toolbox they just run serially.

## Design principles

- **R2020a-compatible.** No `arguments` blocks, no name=value call syntax.
- **2-D is a thin wrapper.** `permtest2` / `gpermtest2` / `FDR_2D` each reshape to 2-D, delegate to the 1-D core, reshape the result. Small, auditable wrappers.
- **Plotting is optional.** All functions accept `'ploton', false` for programmatic use.

## License

BSD 3-Clause. See [`LICENSE`](LICENSE).
