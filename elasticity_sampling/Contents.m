% ELASTICITY_SAMPLING - Matlab toolbox for thermodynamically feasible elasticity sampling
%
% Options for all functions are stored in matlab structs 'es_options' and 'es_constraints'.
%
% Main functions
%   es_default_options         - Set default values in 'es_options' and 'es_constraints'
%   es_sample_model            - Sample all model parameters for a network
%   es_sample_multiple         - Create model ensemble and collect values for an output function
%   es_compare_ensembles - Compare the output function for two model ensembles
%
% Functions for phases of model sampling (called by 'es_sample_model')
%   es_sample_steady_state - Sample feasible stationary state for a network
%   es_sample_elasticities     - Sample elasticities, run MCA, and reconstruct model
%
% Demos
%   See directory 'demo'
%
% Other MATLAB Toolboxes required:
%   Metabolic Network Toolbox 
%   Tensor Toolbox, Sandia National Labs
%     see http://www.sandia.gov/~tgkolda/TensorToolbox/index-2.5.html
%
% Copyright (C) 2011-2013
% Wolfram Liebermeister  <wolfram.liebermeister@gmail.com>