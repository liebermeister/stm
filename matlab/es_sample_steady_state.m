function [c0, c, u, v, Keq, mu, K, Kma, zeta, A] = es_sample_steady_state_state(N, W, ind_ext, es_constraints, es_options)

% ES_SAMPLE_STEADY_STATE - Sample feasible stationary state for a network
%
% [c0, c, u, v, Keq, mu, K, Kma, zeta, A] = es_sample_steady_state_state(N, W, ind_ext, es_constraints, es_options)
%
% Sample thermodynamically consistent steady states, including concentrations,
%   fluxes, and thermodynamic quantities. Enzyme levels are freely sampled.
%
% Inputs (with nm: # metabolites; nr: # reactions)
%   N                            - Stoichiometric matrix (nm x nr)
%   W                            - Allosteric regulation matrix (nr x nm)
%   ind_ext                      - indices of external metabolites
%   es_constraints, es_options:  - explained below (also see 'es_default_options')
%
% Outputs (with nm: # metabolites; nr: # reactions)
%   c0        - Vector of chemical equilibrium concentrations (optional)
%   c         - Vector of concentrations
%   u         - Vector of enzyme levels
%   J         - Flux vector
%   Keq       - Vector of equilibrium constants
%   mu        - Vector of chemical potentials
%   K         - Kernel matrix 
%   Kma       - Vector of mass-action ratios
%   zeta      - zeta vector computed from reaction affinities
%   A         - Vector of reaction affinities (thermodynamic driving forces)
%
%
% Sampling methods used:
%
%   Fluxes: 
%     If a flux vector is predefined in es_constraints.v_fix (no missing values), 
%     this flux vector is directly used. Otherwise, a method to determine the
%     fluxes must be given in 'es_options.sampling_method':
%     
%     'accept_flux':         use values given in  es_constraints.v_mean
%     'sample_and_discard':  sample flux values from the distribution defined by
%                            es_constraints.v_mean, es_constraints.v_std, and es_constraints.v_sign
%                            and omit unfeasible cycles
%     'convex_optimisation': sample fluxes based on convex FBA from the distribution defined by 
%                            es_constraints.v_mean, es_constraints.v_std, and es_constraints.dmu_fix
%                            and omit unfeasible cycles
%     'v and mu',            sample flux values using sample_feasible_v.m
%     'v from data':         sample flux values from the distribution defined by
%                            es_constraints.v_mean, es_constraints.v_std, and es_constraints.v_sign
%   
%   Concentrations: 
%     Concentrations given in log_c_fix will be directly used. The remaining 
%     concentrations are randomly sampled based on log_c_mean and log_c_std.
%   
%   Chemical potentials
%     The method chosen to determine concentrations and flux depends on the 
%     information provided in the "[..]_fix" fields of es_constraints:
%     mu_fix, dmu_fix, Keq_fix, mu0_fix.
%
%
%  Fields in 'es_constraints' and 'es_options' relevant to this function:
%    es_options.seed               - Random seed used
%    es_options.sampling_methods   - Alternative sampling procedures
%                                    'v from data'
%                                    'v and mu'
%                                    'c0 and c'
%                                    'sample_and_discard'
%                                    'convex_optimisation'  
%    es_constraints.v_fix          - Predefined fluxes (vector, NaNs are ignored)
%    es_constraints.v_min          - Lower flux bounds (vector, NaNs are ignored)
%    es_constraints.v_max          - Upper flux bounds (vector, NaNs are ignored)
%    es_constraints.v_sign         - Flux signs (vector, NaNs are ignored)
%    es_constraints.log_u_mean     - Mean enzyme level (for sampling)
%    es_constraints.log_u_std      - Std dev for enzyme levels (for sampling)
%    es_constraints.log_c_mean     - Mean metabolite concentration (for sampling)
%    es_constraints.log_c_std      - Std dev for metabolite concentrations (for sampling)
%  
%   Additional fields required, depending on 'es_options.sampling_method':
%     'v from data'          es_constraints.v_mean
%                            es_constraints.v_std
%   
%     'v and mu'             es_constraints.mu_fix
%   
%     'c0 and c'             es_constraints.log_c0
%                            es_constraints.log_c0_std
%   
%     'sample_and_discard'   es_options.cycle_correction
%
%     'convex_optimisation'  es_constraints.dmu_fix
%                            es_options.cycle_correction
%  
%   For some methods, additionally required:
%      es_options.ind_ignore  reactions to be ignored in thermodynamic loops 
%                             (only needed for flux correction by loop substraction)
%                             (only needed for flux correction by convex FBA)
%

% ---------------------------------------------------------------------------------
% initialise

if ~isnan(es_options.seed), 
  randn('state',es_options.seed); 
  rand('state',es_options.seed); 
end

[nm,nr] = size(N);

all_v_given   = sum(isnan(es_constraints.v_fix))     == 0;
all_c_given   = sum(isnan(es_constraints.log_c_fix)) == 0;
all_mu0_given = sum(isnan(es_constraints.mu0_fix))   == 0;

if all_mu0_given,
  es_constraints.Keq_fix = exp(-1/RT * N' * es_constraints.mu0_fix);
end

all_keq_given = sum(isnan(es_constraints.Keq_fix)) == 0;

if all_c_given * all_mu0_given,
  es_constraints.mu_fix = es_constraints.mu0_fix + RT * es_constraints.log_c_fix;
end

all_mu_given  = sum(isnan(es_constraints.mu_fix)) == 0;

if all_mu_given,
  es_constraints.dmu_fix = N' * es_constraints.mu_fix;
end

if all_keq_given * all_c_given,
  es_constraints.dmu_fix = -RT * [es_constraints.Keq_fix - N' * es_constraints.log_c_fix];
end

all_dmu_given = sum(isnan(es_constraints.dmu_fix)) == 0;

% check 

if ~sum(isnan([es_constraints.v_fix; es_constraints.Keq_fix; es_constraints.log_c_fix])),
  my_v = es_constraints.v_fix;
  my_theta = log(es_constraints.Keq_fix) - N'*es_constraints.log_c_fix;
  if sum([sign(my_v.*my_theta)==-1].*[my_v~=0]), 
    sign([my_v,my_theta])
    error('Predefined metabolic state is thermodynamically infeasible')
  end
end


% ---------------------------------------------------------------------------------
% determine fluxes v

if all_v_given,
  
  display(' Using given flux distribution');
  v = es_constraints.v_fix; 

else

  v = [];
  es_constraints.v_min(es_constraints.v_sign ==  1) = 0;
  es_constraints.v_max(es_constraints.v_sign == -1) = 0;

  switch es_options.sampling_method,
    
    case 'accept_flux',
      v = es_constraints.v_mean;      
    
    case 'sample_and_discard',
      es_options.sampling_cycle_correction = 1;    
      [sample_v, best_v] = sample_fluxes_given_data(N, ind_ext, es_constraints.v_mean, es_constraints.v_std, 1, 1, 1, es_options.ind_ignore, es_constraints.v_sign);
      v = best_v.v;
    
    case 'convex_optimisation',
      es_options.sampling_cycle_correction = 1;
      [sample_v, best_v, res] = sample_fluxes_by_convex_fba(N,ind_ext,es_constraints.v_mean,es_constraints.v_std,struct,-es_constraints.dmu_fix);
      v                   = best_v.v;
      es_constraints.dmu_fix = -res.A;    % 
    
    case 'v and mu',
      v = sample_feasible_v(N,ind_ext,es_constraints,es_options);  
    
    case 'v from data',
      sample_v = sample_fluxes_given_data(N, ind_ext, es_constraints.v_mean, es_constraints.v_std, 1, 1, 1, es_options.ind_ignore, es_constraints.v_sign);
      v        = sample_v.v{1};
      
    otherwise,
    error(sprintf('unknown method %s',  es_options.sampling_method));
  end  

  if length(v),
    if es_options.sampling_cycle_correction,
      %% correct        = es_check_flux_es_constraints(es_constraints.v_mean,N,ind_ext,es_constraints,1,epsilon_stationary,cycles);
      %% correct mean fluxes by Nathan's method (requires variable 'ind_ignore')
      [v_feas, C] = eba_make_feasible(v, N, 'loose', nan, es_options.ind_ignore, 'efmtool');
      v = v_feas;  % if not given -> fluxes are sampled again
    end
  end

end

% ---------------------------------------------------------------------------------
% given the fluxes, determine c, mu, Kma, Keq

%% sample c freely, determine Keq, A, zeta, depending on given information

c0       = [];
mu       = [];
K        = [];

if all_c_given,
  c        = exp(es_constraints.log_c_fix);
  display(' Using given metabolite profile');
else
  c        = exp(es_constraints.log_c_mean + es_constraints.log_c_std .* randn(nm,1));
  ind_c    = find(isfinite(es_constraints.log_c_fix));
  c(ind_c) = exp(es_constraints.log_c_fix(ind_c));
end

% mass action ratios
Kma = exp(N' * log(c));

% safety fix
Kma(Kma<10^-10) = 10^-10;

if all_keq_given,

  display(' Using given equilibrium constants');
  Keq  = es_constraints.Keq_fix;
  A    = RT * log(Keq ./ Kma);
  zeta = exp(es_options.h .* A/RT);

elseif all_mu_given,

  display(' Using given chemical potentials');
  mu   = es_constraints.mu_fix;
  A    = - N' * mu;
  zeta = exp(es_options.h .* A/RT);
  Keq  = exp(A/RT) .* Kma;

elseif all_dmu_given,

  display(' Using given chemical potential differences');
  A    = - es_constraints.dmu_fix;
  zeta = exp(es_options.h .* A/RT);
  Keq  = exp(A/RT) .* Kma;
  
  %% safety fix
  %% Keq(Keq>10^50) = 10^50;

  %% JUST FOR SAFETY ..REMOVE THIS LATER:
  if ~isfield(es_constraints,'mu_eqconstraint'),
    es_constraints.mu_eqconstraint = [];
  end
  
  if isempty(es_constraints.mu_eqconstraint),
    mu = -pinv(full(N'))*A;
  else,
    mu = pinv([-full(N'); es_constraints.mu_eqconstraint.matrix'])*[-A; es_constraints.mu_eqconstraint.vector];
  end
  
  display(' * Computing mu from A by pseudoinverse (in es_sample_steady_state)');

else,

  %% sample mu
  if es_options.verbose, display(' Computing extreme feasible mu vectors ..'); end
  [mu_list, feasible] = sample_feasible_mu(N, ind_ext, v, es_constraints, es_options);
  if es_options.verbose, display(' .. done'); end
  if ~feasible, error('Unfeasible flux distribution'); end
  display('  Choosing chemical potentials mu randomly from extreme points of the allowed region');
  mu     = mu_list(:,ceil(rand*size(mu_list,2)));
  log_c  = log(c);
  log_c0 = log_c - mu/RT;
  Keq    = exp(N'*log_c0);
  Kma    = exp(N' * log(c));
  zeta   = [Keq./Kma].^es_options.h;
  A      = RT * log(Keq./Kma);
end

if isempty(mu),
  mu   = pinv(full(N')) * es_constraints.dmu_fix;
  c0   = exp(mu/RT) ./ c;
end

if find([v~=0] .* [abs(A) < es_constraints.dmu_limit_min]),
  ind = find([v~=0] .* [abs(A) < es_constraints.dmu_limit_min]);
  error(sprintf('Overly small reaction affinity %f encountered', A(ind) .* sign(v(ind)))); 
end

if find([v~=0].*[abs(A) > [1+10^-5]*es_constraints.dmu_limit]), 
  error(sprintf('Overly large reaction affinity %f encountered', max(abs(A)))); 
end


% ----------------------------------------------------------------------
% output quantities

u = exp(es_constraints.log_u_mean + es_constraints.log_u_std .* randn(nr,1));


% ----------------------------------------------------------------------
% safety checks

if sum([v .* A]<0), 
  [v, A, sign(v)~=sign(A)]
  error('Signs of reaction rates and affinities disagree!!!'); 
end 

if sum([log(Keq./Kma) .* A] <0), error('Problem with equilibrium constants'); end 

if sum(Keq==0),     error('Vanishing equilibrium constant'); end
if sum(isinf(Keq)), error('Infinite equilibrium constant'); end
if sum(isnan(Keq)), error('Unknown equilibrium constant'); end
