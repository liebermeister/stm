function [result, stm_options, stm_constraints] = stm_sample_model_read_flux_task(N,W,ind_ext,stm_constraints,stm_options,filename)

%STM_SAMPLE_MODEL_READ_FLUX_TASK - Sample model with thermodynamic calculations
%
% result = stm_sample_model_read_flux_task(N, W, ind_ext, stm_constraints, stm_options, filename)
%
% Same as stm_sample_model, but the steady state (concentrations, fluxes, Gibbs free energies)
% is read from a .mat file (containing a struct 'model'), computed by external software
% from a task file written by stm_sample_model_write_flux_task.
% The field 'model.predicted_rxn_gibbs' is assumed to contain reaction Gibbs free energies
% (Delta G = -A, negative in flux direction)


% ----------------------------------------------------------------
% Initialise

stm_options.flag_test = 0;

if ~isnan(stm_options.seed), rng(stm_options.seed); end

if  max(max(abs(N))) > stm_options.limit_cooperativity
  display('Large stoichiometric coefficients encountered. Compensating them by small reaction cooperativities');
  ind_compensate = find(max(abs(N)) > stm_options.limit_cooperativity);
  stm_options.h(ind_compensate) = stm_options.limit_cooperativity./max(abs(N(:,ind_compensate)));
end

% ----------------------------------------------------------------
% Thermodynamics phase: read feasible steady quantities (Naama's solution)

load(filename); % contains struct 'model'

[nm,nr] = size(N);
K    = sparse(null(full(N(setdiff(1:nm,ind_ext),:)),'r'));
u    = exp(stm_constraints.log_u_mean + stm_constraints.log_u_std .* randn(nr,1));
c    = model.predicted_conc;
A    = -model.predicted_rxn_gibbs;
mu0  = model.predicted_met_gibbs;
mu   = mu0 + RT * log(c);
Kma  = exp(model.S' * log(c));
c0   = exp(-mu0/RT);  % likely to be numerically intractable;
J    = model.predicted_flux;
Keq  = exp(A/RT) .* Kma;
zeta = exp(stm_options.h.*A/RT);

% ----------------------------------------------------------------
% Kinetics phase: sample alpha, values and compute beta and gamma values

result = stm_sample_elasticities(N, W, ind_ext, stm_constraints, stm_options, c0, c, u, J, Keq, mu, K, Kma, zeta, A);
