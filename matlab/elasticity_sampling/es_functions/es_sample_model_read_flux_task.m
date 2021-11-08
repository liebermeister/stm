function [result, es_options, es_constraints] = es_sample_model_read_flux_task(N,W,ind_ext,es_constraints,es_options,filename)

%ES_SAMPLE_MODEL_READ_FLUX_TASK - Sample model with thermodynamic calculations
% 
% result = es_sample_model_read_flux_task(N, W, ind_ext, es_constraints, es_options, filename)
% 
% Same as es_sample_model, but using thermodynamics calculations


% ----------------------------------------------------------------
% Initialise

es_options.flag_test = 0;

if ~isnan(es_options.seed), randn('state',es_options.seed); end

if  max(max(abs(N))) > es_options.limit_cooperativity,
  display('Large stoichiometric coefficients encountered. Compensating them by small reaction cooperativities');
  ind_compensate = find(max(abs(N)) > es_options.limit_cooperativity);
  es_options.h(ind_compensate) = 1./max(abs(N(:,ind_compensate))); 
end

% ----------------------------------------------------------------
% Thermodynamics phase: read feasible steady quantities (Naama's solution)

load(filename);

K    = sparse(null(network.N(find(network.external==0),:),'r'));
u    = exp(es_constraints.log_u_mean + es_constraints.log_u_std .* randn(nr,1));
c    = model.predicted_conc;
A    = model.predicted_rxn_gibbs;
mu0  = model.predicted_met_gibbs;
mu   = mu0 + RT * log(c);
Kma  = exp(model.S' * log(c));
c0   = exp(-mu0/RT);  % likely to be numerically intractable; 
J    = model.predicted_flux;
Keq  = exp(A/RT) .* Kma;
zeta = exp(-es_options.h.*A/RT);

% ----------------------------------------------------------------
% Kinetics phase: sample alpha, values and compute beta and gamma values

result = es_sample_elasticities(N, W, ind_ext, es_constraints, es_options, c0, c, u, J, Keq, mu, K, Kma, zeta, A);
