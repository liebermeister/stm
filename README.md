Structural thermokinetic modelling for Matlab
=============================================

[Structural thermokinetic modelling](https://www.metabolic-economics.de/stm/) (STM) is a computational method to determine consistent states of kinetic metabolic models. It follows the paradigm of Structural Kinetic Modelling, but accounts for thermodynamic laws to yield a consistent description of reversible biochemical reaction kinetics.

This repository contains Matlab functions for Structural thermokinetic modelling that extend the [Metabolic Network Toolbox](https://github.com/liebermeister/metabolic-network-toolbox).

## Installation

Download or clone the repository and add the `matlab` folder, with its subfolders, to your Matlab path:

```matlab
addpath(genpath('/path/to/stm/matlab'));
savepath;   % optional: keep the folders on the path in future sessions
```

Then install the required packages (see below). To check which of them are installed, run `es_dependencies`.

For getting started, run one of the demo scripts in the folder `matlab/demo`:

- `demo_stm_glycolysis`: interactive demo (a model of glycolysis; press a key to continue after each step)
- `demo_stm_noor_2016_modified`: model of E. coli central metabolism
- `demo_stm_workflow`: the complete workflow, with possible variants (commented out by `if 0`); a good starting point for applying STM to your own models

and see the instructions at http://www.metabolic-economics.de/stm/workflow-matlab.html

## Requirements

- Matlab (tested with R2025b)
- [Matlab utility functions](https://github.com/liebermeister/matlab-utils) (required)
- [Metabolic Network Toolbox](https://github.com/liebermeister/metabolic-network-toolbox) (required)
- Optimization Toolbox (for linear programs and `fsolve`). If IBM CPLEX is on the Matlab path, it is used
  for linear programs instead of `linprog` (see `help lp_solve` in matlab-utils).

Some functions require further packages:

| Package | Needed for |
|---------|------------|
| [SBtab for Matlab](https://github.com/liebermeister/sbtab-matlab) | reading models and data from SBtab files (demos, E. coli model) |
| [libSBML](https://sbml.org/software/libsbml/) (Matlab interface) | SBML import and export (`demo_stm_glycolysis`, `es_sbml_to_model`) |
| [Tensor Toolbox](https://www.tensortoolbox.org) | second-order elasticities and response coefficients (`es_options.flag_second_order = 1`) |
| [efmtool](https://csb.ethz.ch/tools/software/efmtool.html) | removing thermodynamically infeasible cycles from flux distributions |
| [Parameter balancing](https://github.com/liebermeister/parameter-balancing-matlab) | `demo_stm_workflow` and the scripts for building the E. coli model |
| Statistics and Machine Learning Toolbox | synergy graphics and clustering, beta-distributed saturation values (`es_options.set_alpha_nonuniform = 1`) |

## Contents

- **`matlab`**: the main functions
  - `es_default_options`: default settings (structs `es_options` and `es_constraints`; type `help es_default_options`)
  - `es_sample_model`: sample a model (steady state and elasticities) for a given network
  - `es_sample_steady_state`, `es_sample_elasticities`: the two steps of `es_sample_model`
  - `es_sample_ensemble`, `es_compare_ensembles`: create model ensembles and compare them statistically
- **`matlab/stm`**: more specific functions: graphics, prediction of perturbation effects,
  enzyme synergies, metabolic fluctuations, Bayesian (MCMC) sampling
- **`matlab/demo`**: demo scripts
- **`matlab/models/e_coli_ccm`**: scripts for building the E. coli model from the STM article
- **`resource-data`**: models and data (see the README files there)
- **`tests`**: tests
- **`doc`**: HTML documentation of all functions

Each folder in `matlab` contains a `Contents.m` file with an overview (type, e.g., `help stm`).

## E. coli model

The E. coli central metabolism model, the running example in the article, can be found in the folder `resource-data/models-article/Escherichia_coli_ccm/`. Details about the model and model construction are given there.

## Tests

From the repository folder, run

```matlab
results = runtests('tests')
```

## License
This package is released under the [GNU General Public License](LICENSE), version 2 or later.

## Contact
Please contact [Wolfram Liebermeister](mailto:wolfram.liebermeister@gmail.com) with any questions or comments.

## References
Liebermeister W. (2013), *Elasticity sampling links thermodynamics to metabolic control*,
[arXiv:1309.0267](https://arxiv.org/abs/1309.0267)

Liebermeister W. (2022), *Structural thermokinetic modelling*,
[Metabolites 12(5), 434](https://www.mdpi.com/2218-1989/12/5/434)
