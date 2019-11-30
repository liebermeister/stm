Elasticity Sampling functions for Matlab
========================================

[Elasticity sampling](https://www.metabolic-economics.de/elasticity-sampling/index.html) is a computational method to determine consistent states of kinetic metabolic models. It follows the paradigm of Structural Kinetic Modelling, but accounts for thermodynamic laws to yield a consistent description of reversible biochemical reaction kinetics.

This repository contains Matlab functions for Elasticity Sampling that extend the [Metabolic Network Toolbox](https://github.com/liebermeister/metabolic-network-toolbox).

## Dependencies: 

The functions for elasticity sampling depend on the Metabolic Network Toolbox (https://github.com/liebermeister/metabolic-network-toolbox). Some functions depend on the following Matlab toolboxes:

  o SBML Toolbox               (http://sbml.org/Software/SBMLToolbox)

  o SBtab functions            (https://github.com/liebermeister/sbtab-matlab)

  o Tensor Toolbox             (http://www.sandia.gov/~tgkolda/TensorToolbox/index-2.5.html)

  o efmtool                    (http://www.csb.ethz.ch/tools/efmtool)

Please make sure that these matlab packages are installed in your system and that all these directories and subdirectories are included in your matlab path.

If CPLEX for matlab is installed in your system, this will speed up some calculations,

For getting started, please see the instructions at

http://www.metabolic-economics.de/elasticity-sampling/workflow-matlab.html

## License
This package is released under the [GNU General Public License](LICENSE).

## Contact
Please contact [Wolfram Liebermeister](wolfram.liebermeister@gmail.com) with any questions or comments.

## References
Elasticity sampling links thermodynamics to metabolic control
Liebermeister W. (2013), [arXiv:1309.0267](https://arxiv.org/abs/1309.0267)