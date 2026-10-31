function [result, stm_options, stm_constraints] = stm_sample_model(N, W, ind_ext, stm_constraints, stm_options)

% STM_SAMPLE_MODEL - Sample all model parameters for a network
%
% result = stm_sample_model(N, W, ind_ext, stm_constraints, stm_options)
%
% This script reconstructs one model instance, including the sampled elasticities.
% For samplig multiple model instances, see 'stm_sample_ensemble'.
%
%
% Inputs (with nm: # metabolites; nr: # reactions)
%   N         - Stoichiometric matrix (nm x nr)
%   W         - Allosteric regulation matrix (nr x nm)
%   ind_ext   - indices of external metabolites
%
%   For the inputs 'stm_constraints' and 'stm_options', see 'stm_default_options'
%
%
% Outputs
%   All output data are stored in a data structure 'result' (see 'stm_sample_elasticities')
%   stm_options and stm_constraints are returned with the (possible) changes made.
%
%
% Sampling consists of two stages:
%   1. Steady state          (calling 'stm_sample_steady_state')
%   2. Elasticities and MCA  (calling 'stm_sample_elasticities')


% ----------------------------------------------------------------
% Initialise

if exist('dependencies_check_once','file')
  dependencies_check_once('stm_dependencies', true);
end

stm_options.flag_test = 0;

if ~isnan(stm_options.seed), rng(stm_options.seed); end

% the following functions must not seed the generator again: otherwise, different
% model variables would be built from identical random numbers (and thus be correlated)
stm_options_no_seed = stm_options;
stm_options_no_seed.seed = nan;


% ----------------------------------------------------------------
% Check reaction sum formulae for large stoichoimetric coefficients
% they are compensated by setting smaller cooperativities
% such that m+ and m- reach maximally the threshold value
% (in stm_options.limit_cooperativity)
% Currently, the settings for reaction affinities are NOT changed

M_max = stm_options.limit_cooperativity;

if  max(max(abs(N))) > M_max
  display(' * Large stoichiometric coefficients encountered. Compensating them by small');
  display('   reaction cooperativities h; predefined reaction affinities and equilibrium');
  display('   constants are not changed');
  ind = find(max(abs(N)) > M_max);
  stm_options.h(ind) = M_max./max(abs(N(:,ind)));
  stm_options_no_seed.h = stm_options.h;
  %stm_constraints.dmu_fix(ind) = stm_constraints.dmu_fix(ind) ./ stm_options.h(ind);
  %stm_constraints.Keq_fix(ind) = stm_constraints.Keq_fix(ind) .^ stm_options.h(ind);
end


% ----------------------------------------------------------------
% Steady-state phase: sample feasible steady quantities

[c0, c, u, J, Keq, mu, K, Kma, zeta, A] = stm_sample_steady_state(N, W, ind_ext, stm_constraints, stm_options_no_seed);


% ----------------------------------------------------------------
% Kinetics phase: sample saturation values compute all the rest

result = stm_sample_elasticities(N, W, ind_ext, stm_constraints, stm_options_no_seed, c0, c, u, J, Keq, mu, K, Kma, zeta, A);
