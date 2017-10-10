% ELASTICITY_SAMPLING - Matlab toolbox for thermodynamically feasible elasticity sampling
%
% Main functions
%   es_default_options     - Create option structures 'es_options' and 'es_constraints' with default parameters
%   es_sample_model        - Sample all model parameters for a network; 
%                            wrapper around 'es_sample_steady_state' and 'es_sample_elasticities'
%   es_sample_multiple     - Create model ensemble and collect values for an output function
%   es_compare_ensembles   - Compare the output function for two model ensembles
%
% Functions for individal phases of model sampling (called by 'es_sample_model')
%   es_sample_steady_state - Sample feasible stationary states for a network
%   es_sample_elasticities - Given a stationary state, sample elasticities, run MCA, and reconstruct model
%
% Options for all functions are stored in matlab structs 'es_options' and 'es_constraints'.
%
% Demos
%   See directory 'demo'
%
% MATLAB toolboxes required
%   Metabolic Network Toolbox              (https://github.com/wolframliebermeister/mnt)
%   SBMLtoolbox    - SBML import / export  (http://sbml.org/Software/SBMLToolbox)
%   SBtab toolbox  - SBtab format          (https://github.com/wolframliebermeister/sbtab-matlab)
%   Tensor toolbox -                       (http://www.sandia.gov/~tgkolda/TensorToolbox/index-2.5.html)
%   efmtool        - Elementary flux modes (http://www.csb.ethz.ch/tools/efmtool)
%
% Copyright (C) 2011-2013
% Wolfram Liebermeister  <wolfram.liebermeister@gmail.com>