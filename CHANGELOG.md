# Changelog

Version numbers have the form MAJOR.MINOR.PATCH. MAJOR changes when existing calls may stop working.

## Not released yet

Requires the current development versions of matlab-utils (`multiple_testing_fdr`, new `verbose` argument of `convex_sampling_hit_and_run`) and of the Metabolic Network Toolbox (`stm_options.mu_sampling_method`). See "Erratum" in `README.md` for changes that affect results.

- Thermodynamic forces are sampled uniformly from the feasible polytope (`stm_options.mu_sampling_method = 'hit-and-run'`); the previous selection of polytope vertices is available as `'extreme_points'`
- `multiple_testing_fdr` moved to matlab-utils and corrected; new function `fdr_two_sided`
- Cooperativities `h` are used when computing turnover rates; fixed sign error with `mu_eqconstraint`; no clipping of small mass-action ratios; fluxes given by `v_min` above zero are respected
- A seed sets the random number generator only once, so different model variables are no longer built from identical random numbers; `interaction_statistics` no longer draws identical samples
- Inactive reactions: `no_equilibrium` is applied correctly, and elasticities of reactions at equilibrium are computed from the one-way rates
- Multi-reactions: fixed crash with second-order elasticities and overwriting of other elasticities
- Numerical rank in the large-network branch of `stm_sample_elasticities` with a relative tolerance
- `spectral_density_to_smooth_variance`: variances were too large by a factor of pi
- `solve_quadratic2`: corrected derivation (factor 1/2) and the state of the random number generator is restored
- Fixed `stm_simulate_perturbation_sample` (means over the runs), `stm_discover_perturbation_sample` (default priors), `check_hyperantagonistic`, `pvalue_from_countnumber` (nan values), `stm_compare_ensembles` (scalar outputs), `segre_epistasis` (lethal mutations), `energetics_control_and_noise` (given concentrations), `matrix_circle_scale`, likelihood functions (kinetic law of the model)
- `synergy_count_cycles`: help text states which cycles are counted
- Unstable states: `stm_sample_model` and `stm_sample_elasticities` return `control.stable = 0` and nan control coefficients; new function `stm_sample_stable_model` draws models until one is stable (with an error after `stm_options.n_max_unstable` unsuccessful attempts), and is used for ensembles, perturbation samples, interaction statistics and MCMC priors; `stm_sample_ensemble` reports the number of discarded models
- `stm_sample_steady_state`: forces above `dmu_limit` that follow from given data give a warning instead of an error
- `influence_anova` and `interaction_anova` warn if there are too few permutations for the false discovery rate (`fdr_check_permutations`)
- New option `stm_options.flag_control`: skip control analysis
- New tests in `tests/test_regressions.m`

## 1.0.0 (2026-10-31)

First numbered version.

- HTML documentation removed from the repository (it can be regenerated with m2html)
- Version function `stm_version`; `stm_dependencies` checks the required versions
- Fixed demo scripts and `stm_default_options`, tests and README
- `stm_default_options` removed: the identical function in the Metabolic Network Toolbox is used
- `stm_sample_model` and `stm_sample_elasticities` check the dependencies automatically, once per Matlab session
- All functions and folders starting with `es_` (and the folder `es-utils`) renamed to `stm_` (`stm-utils`); the structs `es_options` and `es_constraints` are now `stm_options` and `stm_constraints`

The state before the cleanup of October 2026 is available under the git tag `v0.9-legacy`.
