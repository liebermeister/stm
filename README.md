Structural thermokinetic modelling for Matlab
=============================================

[Structural thermokinetic modelling](https://www.metabolic-economics.de/stm/) is a computational method to determine consistent states of kinetic metabolic models. It follows the paradigm of Structural Kinetic Modelling, but accounts for thermodynamic laws to yield a consistent description of reversible biochemical reaction kinetics.

This repository contains Matlab functions for Structural thermokinetic modelling that extend the [Metabolic Network Toolbox](https://github.com/liebermeister/metabolic-network-toolbox).

For demo scripts, see the subdirectory ```matlab/demo```. To see the commands used, please have a look at the scripts.

Relevant models and data can be found in the folder ```resources-data```.

## Dependencies

The functions in this toolbox depend on the

  * Matlab utility functions (https://github.com/liebermeister/matlab-utils)

  * Metabolic Network Toolbox (https://github.com/liebermeister/metabolic-network-toolbox).

Some functions depend on the following Matlab toolboxes:

  * SBML Toolbox               (http://sbml.org/Software/SBMLToolbox)

  * SBtab functions            (https://github.com/liebermeister/sbtab-matlab)

  * Tensor Toolbox             (http://www.sandia.gov/~tgkolda/TensorToolbox/index-2.5.html)

  * efmtool                    (http://www.csb.ethz.ch/tools/efmtool)

Please make sure that these matlab packages are installed in your system and that all these directories and subdirectories are included in your matlab path. If CPLEX for matlab is installed in your system, this will speed up some calculations.

For getting started, please see the instructions at

http://www.metabolic-economics.de/stm/workflow-matlab.html

## E. coli model

The E. coli central metabolism model, which serves as a running example in the article, canbe found in the folder ``resource-data/models-article/Escherichia_coli_ccm/``.

## License
This package is released under the [GNU General Public License](LICENSE).

## Contact
Please contact [Wolfram Liebermeister](mailto:wolfram.liebermeister@gmail.com) with any questions or comments.

## References
Liebermeister W. (2013), *Elasticity sampling links thermodynamics to metabolic control*
[arXiv:1309.0267](https://arxiv.org/abs/1309.0267)
