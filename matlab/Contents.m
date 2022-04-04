% STRUCTURAL THERMODYNAMIC MODELLING - Matlab toolbox for thermodynamically feasible elasticity sampling
%
% Demos
%   For example applications of STM, see the scripts in the 'demo' folder
%
% Main functions
%   es_default_options     - Create structures 'es_options' and 'es_constraints' with default parameters
%   es_sample_steady_state - Generate or sample reference state(s) for a network 
%   es_sample_elasticities - Sample elasticities and control properties for a network (with given reference state)
%   es_sample_model        - Sample all model parameters (wrapper for 'es_sample_steady_state' and 'es_sample_elasticities')
%   es_sample_ensemble     - Create model ensemble and collect values for an output variable
%   es_compare_ensembles   - Compare the output variables for two model ensembles
%
% Extra convenience function 
%   sample_elasticities    -  Wrapper function around 'es_sample_elasticities' with different function arguments
%
%   Options for all functions are stored in matlab structs 'es_options' and 'es_constraints'.
%
% Other STM functions
%   More specific functions, including functions for graphical output, can be found in the subfolder 'stm'    
%   
% Models 
%   Matlab code for generating the E. coli ccm model (from the STM article) can be found in the folder '/models/e_coli_ccm'
%   
% MATLAB toolboxes required
%   Metabolic Network Toolbox              (https://github.com/wolframliebermeister/mnt)
%   SBMLtoolbox    - SBML import / export  (http://sbml.org/Software/SBMLToolbox)
%   SBtab toolbox  - SBtab format          (https://github.com/wolframliebermeister/sbtab-matlab)
%   Tensor toolbox -                       (http://www.sandia.gov/~tgkolda/TensorToolbox/index-2.5.html)
%   efmtool        - Elementary flux modes (http://www.csb.ethz.ch/tools/efmtool)
%
% Copyright (C) 2011-2022
% Wolfram Liebermeister  <wolfram.liebermeister@gmail.com>