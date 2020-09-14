% Elasticity sampling: analysis of hepatonet ATP regeneration network
%
% See Article Figure 3 


% -----------------------------------------------------
% Load model (network structure)
% -----------------------------------------------------

model_dir  = [es_BASEDIR '/../resources/models-article/Hepatonet_CCM_Fluxes'];
model_name = 'Hepatonet_Aerobic_ATP_Regeneration'; 

cd(model_dir); load(model_name);


% -----------------------------------------------------
% Determine thermodynamic state variables by parameter balancing
% -----------------------------------------------------

% load thermodynamic data
% (since human data are not available, E coli data are used as proxies)

state_data_dir = [es_BASEDIR '/../resources/data/kinetic-data/'];

data_files     = {[state_data_dir 'E_coli_glycolysis_GFE.tsv'], ... 
		  [state_data_dir 'GFE.tsv'], ...
                  [state_data_dir 'E_coli_glucose_concentration.tsv'], ...
                  [state_data_dir 'ATPase_affinity.tsv']}';

% predefine some of the concentrations

ind_water     = label_names({'H2O [c]','H2O [m]'},  network.metabolites);
ind_phosphate = label_names({'Pi [c]', 'Pi [m]'},   network.metabolites); 
ind_ATP       = label_names({'ATP [c]', 'ATP [m]'}, network.metabolites);

ind_water     = ind_water(find(ind_water));
ind_phosphate = ind_phosphate(find(ind_phosphate));
ind_ATP       = ind_ATP(find(ind_ATP));

c_fix                = nan * ones(size(network.metabolites));
c_fix(ind_phosphate) = 40 * ones(size(ind_phosphate)); 
c_fix(ind_ATP)       = 10 * ones(size(ind_ATP)); 

pb_options = struct('A_min', 1, 'A_max', 10, 'conc_min', 0.0001, 'conc_max', 100, 'sigma_mu0',3, 'data_refer_to_molar',0, 'ind_water',ind_water,'c_fix',c_fix);

pb_options.virtual_reactions = {};

% add ATP hydrolysis as a virtual reaction (to control its driving force)
% actually balance in mitochondria would also be relevant ...

if find(strcmp('ATP [c]',network.metabolites)),
  pb_options.virtual_reactions(1).metabolites     = {'ATP [c]', 'H2O [c]', 'ADP [c]', 'Pi [c]'};
  pb_options.virtual_reactions(1).stoichiometries = [-1,-1,1,1];
end

if find(strcmp('ATP [m]',network.metabolites)),
  pb_options.virtual_reactions(2).metabolites     = {'ATP [m]', 'H2O [m]', 'ADP [m]', 'Pi [m]'};
  pb_options.virtual_reactions(2).stoichiometries = [-1,-1,1,1];
end

pb_options.virtual_A_lower = 5;
pb_options.virtual_A_upper = 60;

[c, mu0, Keq, A, kinetic_data] = parameter_balancing_thermodynamic(network, v, data_files, pb_options);


% --------------------------------------------
% Run elasticity sampling
% --------------------------------------------

[es_options, es_constraints] = es_default_options(network);

es_constraints.v_fix         = v;
es_constraints.log_c_fix     = log(c); 
es_constraints.mu_fix        = mu0 + RT * log(c); 
es_constraints.dmu_fix       = network.N' * es_constraints.mu_fix;

es_options.ind_ignore        = [];
es_options.flag_second_order = 0;
es_options.sampling_method   = 'v from data';
es_options.kinetic_law       = 'cs'
es_options.flag_second_order = 1;
es_options.print_graphics    = 0;

[result, fluxes, es_options, es_constraints] = es_reference_state(network, es_options, es_constraints);


% ---------------------------------------
% Graphics
% ---------------------------------------

% For clarity, omit NAD compounds from graphics 
network.graphics_par.omitmetabolites       = {'NAD+ [c]', 'NAD+ [m]', 'NADH [c]', 'NADH [m]'};

% Target reaction in graphics: ATP production
% 'ADP(c) + ATP(m) --> ADP(m) + ATP(c)' (from 'hepatonet1_reference_states') 

target_reaction    = 'r2418'; 

es_reference_state_graphics(network, es_options, result, target_reaction);
