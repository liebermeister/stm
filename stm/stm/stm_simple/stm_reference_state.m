function [result, fluxes, stm_options, stm_constraints] = stm_reference_state(network, stm_options, stm_constraints)

%STM_REFERENCE_STATE - Elasticity sampling
%
% [result, fluxes, stm_options, stm_constraints] = stm_reference_state(network, stm_options, stm_constraints)
%
% Construct a reference state and compute the elasticities,
% assuming half-saturation (i.e., saturation values alpha = 0.5)
%
% This function is a wrapper around 'stm_sample_model.m'

eval(default('epsilon_stationary','10^-5','method_flux_sampling','''sample_and_discard'''));

es_options_default = struct('ind_ignore',[]);
stm_options         = join_struct(es_options_default, stm_options);


% ------------------------------------------------------------
% initialise variables

N       = network.N;
W       = network.regulation_matrix;
[nm,nr] = size(N);
ind_ext = find(network.external);


% --------------------------------------------------------------
% set standard stm_options and sample model

stm_options.verbose             = 1;
stm_options.set_alpha_to_half   = 1;
stm_options.seed                = 0;
stm_options.epsilon_stationary  = epsilon_stationary;

result = stm_sample_model(N, W, ind_ext, stm_constraints, stm_options);

% correct = mnt_check_flux_constraints(result.v,N,ind_ext,stm_constraints,1,epsilon_stationary,stm_options.cycles)

% --------------------------------------------------------------

fluxes.v_mean   = stm_constraints.v_mean;
fluxes.sample_v = result.v;
fluxes.best_v   = result.v;
