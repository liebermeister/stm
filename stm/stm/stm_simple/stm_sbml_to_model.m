function [network, result, options, constraints] = stm_sbml_to_model(sbml_filename, external_metabolites, cofactors, v_sign)

%STM_SBML_TO_MODEL Wrapper function around stm_sample_model
%
% [network, result, options, constraints] = stm_sbml_to_model(sbml_filename, external_metabolites, cofactors, v_sign)
%
% Wrapper function around stm_sample_model
%
% Translates sbml file into matlab model structure with modular rate law kinetics and
% samples elasticities with some default choice of the sampling constraints and parameters
%
% Function arguments:
%  sbml_filename         sbml file
%  external_metabolites  will be set external (in addition to "dead ends" detected automatically)
%  cofactors             will also be set external
%  v_sign                vector of predefined signs for flux distribution

eval(default('external_metabolites','{}','cofactors','{}','v_sign','[]'));

network   = network_sbml_import(sbml_filename);

% remove metabolites that do not take part in any reaction
% (the kinetics from the SBML file is removed; it is replaced by the sampled kinetics below)
if isfield(network,'kinetics'), network = rmfield(network,'kinetics'); end
network   = network_choose(network, [], 1:length(network.actions), 0);

deadends      = network.metabolites(network_find_ext_metabolites(network));
met_external  = unique([column(cofactors); column(deadends); column(external_metabolites)]);
network       = network_set_external(network,1,met_external);
network       = netgraph_make_graph(network);


% --------------------------------------------
% set parameters for elasticity sampling

[nm,nr] = size(network.N);

c = ones(nm,1);  % metabolite concentrations
u = ones(nr,1);  % enzyme concentrations

[options, constraints] = stm_default_options(network);

options.kinetic_law       = 'cs';
options.sampling_method   = 'v and mu';
options.set_alpha_to_half = 1;
options.verbose           = 1;
options.seed              = 0;

N       = network.N;
W       = zeros(size(N'));
ext_ind = find(network.external);

if ~isempty(v_sign)
  constraints.v_fix(v_sign==0)  = 0;
  constraints.v_sign            = v_sign;
  constraints.v_sign(v_sign==0) = nan;
end

constraints.v_min             = -1*ones(nr,1);
constraints.v_max             =  1*ones(nr,1);
constraints.log_c_mean        = log(c);
constraints.log_c_std         = zeros(nm,1);
constraints.log_u_mean        = log(u);
constraints.log_u_std         = zeros(nr,1);


% --------------------------------------------
% translate network into kinetic model

result           = stm_sample_model(N,W,ext_ind,constraints,options);
network.kinetics = result.kinetics;
