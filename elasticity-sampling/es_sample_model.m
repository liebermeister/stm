function [result, es_options, es_constraints] = es_sample_model(N, W, ind_ext, es_constraints, es_options)

% ES_SAMPLE_MODEL - Sample all model parameters for a network
%
% result = es_sample_model(N, W, ind_ext, es_constraints, es_options)
% 
% This script reconstructs one model instance, including the sampled elasticities.
% For samplig multiple model instances, see 'es_sample_multiple'.
%
%
% Inputs (with nm: # metabolites; nr: # reactions)
%   N         - Stoichiometric matrix (nm x nr)
%   W         - Allosteric regulation matrix (nr x nm)
%   ind_ext   - indices of external metabolites
%
%   For the inputs 'es_constraints' and 'es_options', see 'es_default_options'
%
%
% Outputs
%   All output data are stored in a data structure 'result' (see 'es_sample_elasticities')
%   es_options and es_constraints are returned with the (possible) changes made.
%
%
% Sampling consists of two stages:
%   1. Steady state          (calling 'es_sample_steady_state')
%   2. Elasticities and MCA  (calling 'es_sample_elasticities')


% ----------------------------------------------------------------
% Initialise

es_options.flag_test = 0;

if ~isnan(es_options.seed), randn('state',es_options.seed); end


% ----------------------------------------------------------------
% Check reaction sum formulae for large stoichoimetric coefficients
% they are compensated by setting smaller cooperativities
% such that m+ and m- reach maximally the threshold value
% (in es_options.limit_cooperativity)
% Currently, the settings for reaction affinities are NOT changed

M_max = es_options.limit_cooperativity;

if  max(max(abs(N))) > M_max,
  display(' * Large stoichiometric coefficients encountered. Compensating them by small');
  display('   reaction cooperativities h; predefined reaction affinities and equilibrium'); 
  display('   constants are adjusted');
  ind = find(max(abs(N)) > M_max);
  es_options.h(ind) = M_max./max(abs(N(:,ind)));
  %es_constraints.dmu_fix(ind) = es_constraints.dmu_fix(ind) ./ es_options.h(ind);
  %es_constraints.Keq_fix(ind) = es_constraints.Keq_fix(ind) .^ es_options.h(ind);
end


% ----------------------------------------------------------------
% Steady-state phase: sample feasible steady quantities

[c0, c, u, J, Keq, mu, K, Kma, zeta, A] = es_sample_steady_state(N, W, ind_ext, es_constraints, es_options);

% ----------------------------------------------------------------
% Kinetics phase: sample saturation values compute all the rest

result = es_sample_elasticities(N, W, ind_ext, es_constraints, es_options, c0, c, u, J, Keq, mu, K, Kma, zeta, A);
