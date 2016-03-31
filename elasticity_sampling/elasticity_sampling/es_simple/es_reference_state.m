function [result, fluxes, es_options, es_constraints] = es_reference_state(network, es_options, es_constraints)

% [result, fluxes, es_options, es_constraints] = es_reference_state(network, es_options, es_constraints)
%
% Elasticity sampling
%
% Construct a reference state and compute the elasticities, 
% assuming half-saturation (i.e., saturation values alpha = 0.5)
%
% This function is a wrapper around 'es_sample_model.m'

eval(default('epsilon_stationary','10^-5','method_flux_sampling','''sample_and_discard'''));

es_options_default = struct('ind_ignore',[]);
es_options         = join_struct(es_options_default, es_options);


% ------------------------------------------------------------
% initialise variables

N       = network.N;
W       = network.regulation_matrix;
[nm,nr] = size(N);
ind_ext = find(network.external);


% --------------------------------------------------------------
% set standard es_options and sample model

es_options.verbose             = 1;
es_options.set_alpha_to_half   = 1;
es_options.seed                = 0;
es_options.epsilon_stationary  = epsilon_stationary;

result = es_sample_model(N, W, ind_ext, es_constraints, es_options);

% correct = es_check_flux_constraints(result.v,N,ind_ext,es_constraints,1,epsilon_stationary,es_options.cycles)

% --------------------------------------------------------------

fluxes.v_mean   = es_constraints.v_mean;
fluxes.sample_v = result.v;
fluxes.best_v   = result.v;
