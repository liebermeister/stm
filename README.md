Structural thermokinetic modelling for Matlab
=============================================

[Structural thermokinetic modelling](https://www.metabolic-economics.de/stm/) (STM) is a computational method to determine consistent states of kinetic metabolic models. It follows the paradigm of Structural Kinetic Modelling, but accounts for thermodynamic laws to yield a consistent description of reversible biochemical reaction kinetics.

This repository contains Matlab functions for Structural thermokinetic modelling that extend the [Metabolic Network Toolbox](https://github.com/liebermeister/metabolic-network-toolbox).

## Installation

Download or clone the repository and add the `stm` folder, with its subfolders, to your Matlab path:

```matlab
addpath(genpath('/path/to/stm/stm'));
savepath;   % optional: keep the folders on the path in future sessions
```

Then install the required packages (see below). To check which of them are installed, run `stm_dependencies`.

For getting started, run one of the demo scripts in the folder `stm/demo`:

- `demo_stm_glycolysis`: interactive demo (a model of glycolysis; press a key to continue after each step)
- `demo_stm_noor_2016_modified`: model of E. coli central metabolism
- `demo_stm_workflow`: the complete workflow, with possible variants (commented out by `if 0`); a good starting point for applying STM to your own models

and see the instructions at http://www.metabolic-economics.de/stm/workflow-matlab.html

## Requirements

- Matlab (tested with R2025b)
- [Matlab utility functions](https://github.com/liebermeister/matlab-utils) (required)
- [Metabolic Network Toolbox](https://github.com/liebermeister/metabolic-network-toolbox) (required)
- Optimization Toolbox (for linear programs and `fsolve`). Optionally, Gurobi can be used
  for linear programs instead of `linprog` (see `help optim_solver` in matlab-utils).

Some functions require further packages:

| Package | Needed for |
|---------|------------|
| [SBtab for Matlab](https://github.com/liebermeister/sbtab-matlab) | reading models and data from SBtab files (demos, E. coli model) |
| [libSBML](https://sbml.org/software/libsbml/) (Matlab interface) | SBML import and export (`demo_stm_glycolysis`, `stm_sbml_to_model`) |
| [Tensor Toolbox](https://www.tensortoolbox.org) | second-order elasticities and response coefficients (`stm_options.flag_second_order = 1`, the default) |
| [efmtool](https://csb.ethz.ch/tools/software/efmtool.html) | removing thermodynamically infeasible cycles from flux distributions |
| [Parameter balancing](https://github.com/liebermeister/parameter-balancing-matlab) | `demo_stm_workflow` and the scripts for building the E. coli model |
| Statistics and Machine Learning Toolbox | synergy graphics and clustering, beta-distributed saturation values (`stm_options.set_alpha_nonuniform = 1`) |

## Contents

- **`stm`**: the main functions
  - `stm_default_options`: default settings (structs `stm_options` and `stm_constraints`; type `help stm_default_options`)
  - `stm_sample_model`: sample a model (steady state and elasticities) for a given network
  - `stm_sample_steady_state`, `stm_sample_elasticities`: the two steps of `stm_sample_model`
  - `stm_sample_ensemble`, `stm_compare_ensembles`: create model ensembles and compare them statistically
- **`stm/stm`**: more specific functions: graphics, prediction of perturbation effects,
  enzyme synergies, metabolic fluctuations, Bayesian (MCMC) sampling
- **`stm/demo`**: demo scripts
- **`resources/models/e_coli_ccm`**: scripts for building the E. coli model from the STM article
- **`resources`**: models and data (see the README files there)
- **`tests`**: tests

Each folder in `stm` contains a `Contents.m` file with an overview (type, e.g., `help stm`).

## E. coli model

The E. coli central metabolism model, the running example in the article, can be found in the folder `resources/models-article/Escherichia_coli_ccm/`. Details about the model and model construction are given there.

## Tests

From the repository folder, run

```matlab
results = runtests('tests')
```

## Erratum

Errors found in a code review in October 2026 have been corrected in this repository, and no erratum was sent to the journal. The corrections can change results obtained with earlier versions, including results in published work. Earlier versions are available as git tags.

- **Thermodynamic forces were not sampled uniformly.** With the sampling method 'v and mu', reaction Gibbs free energies were picked from a few vertices of the feasible polytope. They are now sampled uniformly (`stm_options.mu_sampling_method = 'hit-and-run'`); the old procedure is available as `'extreme_points'`. The default limit `dmu_limit` was lowered from 100 to 50 kJ/mol. Uniform sampling gives typical forces of about half this limit. To use another limit, set `stm_constraints.dmu_min` and `dmu_max` as well, because they are initialised from the default. Ensemble statistics change.
- **Smoothed variances were too large by a factor of pi** (`spectral_density_to_smooth_variance`). This scales the variances of metabolite fluctuations; their relative values do not change.
- **The false discovery rate was about twice the nominal value** (`multiple_testing_fdr`), and the two one-sided tests in `influence_anova` and `interaction_anova` doubled the number of tests. Lists of significant interactions may have been too long.
- **Cooperativities were ignored in turnover rates.** Models with `h ~= 1` (used automatically for large stoichiometric coefficients) were thermodynamically inconsistent when evaluated by `modular_velocities`; the rate laws in `result.kinetics` did not vanish at equilibrium.
- **Fixed random seeds correlated model variables.** With a numerical `stm_options.seed`, variables sampled in different functions were built from identical random numbers. `interaction_statistics` used a fixed seed for every model, so the samples of an ensemble were identical.
- **Inactive reactions** (v = 0): with `no_equilibrium`, the enzyme was not switched off correctly, and elasticities of reactions at equilibrium were wrong.
- **Control coefficients of unstable states** were previously returned with values around 10^10. `control_coefficients` now gives an error or, in the STM functions, nan values with `result.control.stable = 0`. Ensembles contain stable models only.

## Versions

To see the installed version, type `stm_version`. Older versions can be downloaded from the [releases page](https://github.com/liebermeister/stm/releases) of this repository, and changes between versions are listed in `CHANGELOG.md`.

## License
This package is released under the [GNU General Public License](LICENSE), version 2 or later.

## Contact
Please contact [Wolfram Liebermeister](mailto:wolfram.liebermeister@gmail.com) with any questions or comments.

## References
Liebermeister W. (2013), *Elasticity sampling links thermodynamics to metabolic control*,
[arXiv:1309.0267](https://arxiv.org/abs/1309.0267)

Liebermeister W. (2022), *Structural thermokinetic modelling*,
[Metabolites 12(5), 434](https://www.mdpi.com/2218-1989/12/5/434)
