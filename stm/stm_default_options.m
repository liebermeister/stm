function [stm_options, stm_constraints] = stm_default_options(N)

% STM_DEFAULT_OPTIONS - Set default values in 'stm_options' and 'stm_constraints'
%
% [stm_options, stm_constraints] = stm_default_options(N)
%
% The default values are read from two SBtab Config files. To change the defaults, edit these files.
%  - General flux and thermodynamic constraints and options (fields marked with * below) are taken from
%    'flux_constraints_default_options' (Metabolic Network Toolbox; Config file 'flux-constraints_Config.tsv')
%  - The STM-specific options are read from 'structural-thermokinetic-modelling_Config.tsv'
%    (in the same folder as this function; TableID 'ConfigureSTM'; option IDs 'stm_options.<name>'
%    and 'stm_constraints.<name>')
% For vectors and matrices, the Config files contain a single value that is used for all elements;
% the default reaction Gibbs free energy bounds dmu_min and dmu_max (empty in the Config file)
% are given by -dmu_limit and dmu_limit.
%
% Initialise structures 'stm_options' and 'stm_constraints' used to define tasks in elasticity sampling
%
% Function arguments:
% N can be the (complete) stoichiometric matrix or a 'network' data structure from the Metabolic Network Toolbox)
%
% Fields marked with * : general flux and thermodynamic settings (see 'flux_constraints_default_options')
%
% Function outputs:
%
% Structure 'stm_constraints': directives for defining the constraints in model sampling
%
% Field name                        Size    Meaning
% ----------                        ----    -------
%*stm_constraints.v_fix              nr x 1  fluxes predefined values (overrides vmin and vmax)
%*stm_constraints.v_min              nr x 1  fluxes lower bounds (overrides -vmax)
%*stm_constraints.v_max              nr x 1  fluxes upper bounds
%*stm_constraints.v_sign             nr x 1  fluxes signs (overrides vmin and vmax)
%*stm_constraints.v_mean             nr x 1  fluxes data (can also be a matrix)
%*stm_constraints.v_std              nr x 1  fluxes data error bars (can also be a matrix)
%*stm_constraints.Keq_fix            nr x 1  equilibrium constants predefined values
%*stm_constraints.ext_signs          nm x 1  signs of production rates of external metabolites (nan for internal metabolites)
%*stm_constraints.log_c_fix          nm x 1  ln metabolite levels predefined values
%*stm_constraints.log_c_min          nm x 1  ln metabolite levels lower bounds
%*stm_constraints.log_c_max          nm x 1  ln metabolite levels upper bounds
%*stm_constraints.log_c_mean         nm x 1  ln metabolite levels data
% stm_constraints.log_c0_mean        nm x 1  ln equilibrium metabolite levels data
% stm_constraints.log_c0_std         nm x 1  ln equilibrium metabolite levels data error bars
%*stm_constraints.mu0_fix            nm x 1  std chemical potentials predefined values
%*stm_constraints.mu0_min            nm x 1  std chemical potentials lower bounds
%*stm_constraints.mu0_max            nm x 1  std chemical potentials upper bounds
%*stm_constraints.mu_fix             nm x 1  chemical potentials predefined values
%*stm_constraints.mu_min             nm x 1  chemical potentials lower bounds
%*stm_constraints.mu_max             nm x 1  chemical potentials upper bounds
%*stm_constraints.dmu_fix            nr x 1  Reaction Gibbs free energies predefined values
%*stm_constraints.dmu_min            nr x 1  Reaction Gibbs free energies lower bounds
%*stm_constraints.dmu_max            nr x 1  Reaction Gibbs free energies upper bounds
%*stm_constraints.dmu_limit           1 x 1  Reaction Gibbs free energies upper limit for absolute values
%                                            (when raising it, raise dmu_min and dmu_max as well)
%*stm_constraints.dmu_limit_min       1 x 1  Reaction Gibbs free energies lower limit for absolute values
% stm_constraints.log_u_mean         nr x 1  ln enzyme levels data
% stm_constraints.log_u_std          nr x 1  ln enzyme levels data error bars
% stm_constraints.beta_M_fix         nr x nm saturation values beta_M predefined values
% stm_constraints.beta_A_fix         nr x nm saturation values beta_A predefined values
% stm_constraints.beta_I_fix         nr x nm saturation values beta_I predefined values
% stm_constraints.alpha_min          1 x 1  minimum value for uniform sampling of alpha values
% stm_constraints.alpha_max          1 x 1  maximum value for uniform sampling of alpha values
% stm_constraints.alpha_A_mean       nr x nm mean values for alpha_A    (only for beta distribution)
% stm_constraints.alpha_A_std        nr x nm std dev values for alpha_A (only for beta distribution)
% stm_constraints.alpha_I_mean       nr x nm mean values for alpha_I    (only for beta distribution)
% stm_constraints.alpha_I_std        nr x nm std dev values for alpha_I (only for beta distribution)
% stm_constraints.alpha_M_mean       nr x nm mean values for alpha_M    (only for beta distribution)
% stm_constraints.alpha_M_std        nr x nm std dev values for alpha_M (only for beta distribution)
%*stm_constraints.rho                 1 x 1  upper bound on v+/v or v-/v
%*stm_constraints.dmu_eqconstraint  struct  Fields .matrix, .vector define equality constraints
%*stm_constraints.mu_eqconstraint   struct  Fields .matrix, .vector define equality constraints
%
% with nr = #reactions; nm = #metabolites; nx = #external metabolites
% For defining only some entries in a vector, replace the others by 'nan'
%
%
% 2. Structure 'stm_options': general directives for the algorithm
%
% Name                        Type           Meaning
% ----                        ----           -------
%*stm_options.seed                            Random seed
% stm_options.sampling_method                 Strategy for steady state sampling
%                                              'v from data'
%                                              'v and mu'
%                                              'sample_and_discard'
%                                              'accept_flux'
% stm_options.mu_sampling_method   string     How to choose reaction Gibbs free energies (with 'v and mu'):
%                                              'hit-and-run': uniform sampling from the feasible polytope
%                                              'extreme_points': vertices of the polytope along the coordinate axes
% stm_options.flux_units                      'dimensionless', 'mM/s', 'mol/s'
% stm_options.n_samples            integer    Number of samples in multiple sampling (model ensemble)
% stm_options.n_max_unstable       integer    Maximal number of unsuccessful attempts to sample a stable model (then an error is thrown)
% stm_options.flag_control         boolean    Compute control and response coefficients? (0: faster; result.control.stable is empty)
%*stm_options.n_dmu_samples        integer    For each of these: number of reaction affinity samples
% stm_options.kinetic_law          string     Type of modular rate law (see modular rate laws, e.g., 'cs', 'ms')
% stm_options.h                    nr x 1     Predefined reaction cooperativities
% stm_options.set_alpha_to_half    boolean    Set all saturation values to 0.5 instead of random sampling?
% stm_options.set_alpha_nonuniform boolean    Draw alpha values from beta distribution?
%                                              Use distribution parameters from the sparse matrices
%                                              stm_constraints.alpha_A_mean etc  (only nonzero entries are used):
% stm_options.KV_prior_mean        float      kV value to be substituted in inactive reactions
% stm_options.limit_cooperativity  float      threshold value for reaction cooperativity
% stm_options.no_equilibrium       boolean    If rates vanish, assume that zero enzyme levels are the reason
%*stm_options.ind_ignore           vector     Indices of reactions to be ignored when computing thermodynamic
%                                              loops (only in for flux correction by loop substraction)
% stm_options.zc                   vector     Metabolite derivative of objective function
% stm_options.zv                   vector     Flux derivative of objective function
% stm_options.limit_elasticities_in_multireactions (boolean) normalise scaled elasticities in multireactions (eg biomass reaction) to an sum of absolute values = 1
%*stm_options.verbose              (boolean)  Verbose output?
% stm_options.graphics_flag        (boolean)  Show graphics?
% stm_options.print_graphics       (boolean)  Save graphics to file?
% stm_options.flag_test            (boolean)  Run tests?
% stm_options.sampling_cycle_correction (boolean) Correct sampled fluxes for thermodynamic cycles
%*stm_options.epsilon_stationary   float      Tolerance for stationarity
%*stm_options.cycles                          Thermodynamic cycles (NaN: computed when needed)
% stm_options.flag_second_order    boolean    Compute second-order elasticities and response coefficients?

% TO BE IMPLEMENTED (matrices should rather be in 'constraints')
% stm_options.sample_K     Sample KM, KA, KI values instead if saturation values
%                         Use distribution parameters from the sparse matrices:
%                         (only nonzero entries are used):
% stm_constraints.KA_mean  mean values for alpha_A
% stm_constraints.KA_std   std dev values for alpha_A
% stm_constraints.KI_mean  mean values for alpha_I
% stm_constraints.KI_std   std dev values for alpha_I
% stm_constraints.KM_mean  mean values for alpha_M
% stm_constraints.KM_std   std dev values for alpha_M

if isstruct(N), N = N.N; end

[nm,nr] = size(N);

% general flux and thermodynamic constraints and options (Metabolic Network Toolbox)

[stm_options, stm_constraints] = flux_constraints_default_options(N);

% STM-specific options

config_file = fullfile(fileparts(mfilename('fullpath')), 'structural-thermokinetic-modelling_Config.tsv');
config      = load_config_table(config_file);
stm_options     = join_struct(stm_options,     config.stm_options);
stm_constraints = join_struct(stm_constraints, config.stm_constraints);

% expand the values from the Config file to vectors and matrices

for f = {'log_u_mean','log_u_std'}
  stm_constraints.(f{1}) = stm_constraints.(f{1}) * ones(nr,1);
end

for f = {'log_c0_mean','log_c0_std'}
  stm_constraints.(f{1}) = stm_constraints.(f{1}) * ones(nm,1);
end

for f = {'beta_M_fix','beta_A_fix','beta_I_fix'}
  stm_constraints.(f{1}) = stm_constraints.(f{1}) * ones(nr,nm);
end

stm_options.h = stm_options.h * ones(nr,1); % reaction cooperativity
