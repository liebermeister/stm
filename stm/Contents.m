% STRUCTURAL THERMOKINETIC MODELLING - Matlab toolbox for thermodynamically feasible elasticity sampling
%
% Demos
%   For example applications of STM, see the scripts in the 'demo' folder
%
% Main functions
%   stm_default_options     - Create structures 'stm_options' and 'stm_constraints' with default parameters
%   stm_sample_steady_state - Generate or sample a reference state for a network
%   stm_sample_elasticities - Sample elasticities and control properties for a network (with given reference state)
%   stm_sample_model        - Sample all model parameters (wrapper for 'stm_sample_steady_state' and 'stm_sample_elasticities')
%   stm_sample_ensemble     - Create model ensemble and collect values for an output variable
%   stm_compare_ensembles   - Compare the output variables for two model ensembles
%
%   Options for all functions are stored in matlab structs 'stm_options' and 'stm_constraints'.
%
% Other STM functions
%   More specific functions, including functions for graphical output, can be found in the subfolder 'stm'
%
% Models
%   Matlab code for generating the E. coli ccm model (from the STM article) can be found in the repository folder 'resources/models-article/Escherichia_coli_ccm/scripts'
%
% Required packages (see README.md for details)
%   Matlab utility functions  (https://github.com/liebermeister/matlab-utils)
%   Metabolic Network Toolbox (https://github.com/liebermeister/metabolic-network-toolbox)
%
% To check which packages are installed, run 'stm_dependencies'
%
% Copyright (C) 2011-2026
% Wolfram Liebermeister  <wolfram.liebermeister@gmail.com>
