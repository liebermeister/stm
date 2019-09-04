% elasticity_sampling_workflow - Demo script
%
% Demo script for thermodynamically feasible elasticity sampling 
% 
% In this script, we load a metabolic network model, run elasticity sampling 
% with two types of rate laws, and compare the results of the two model ensembles 


%-------------------------------------------------------------------------------
% STEP 0: Preparations
%-------------------------------------------------------------------------------

clear

% Check whether all necessary m-files have been installed:

es_dependencies;

% Resource directory of the elasticity sampling package, in case it s needed

resource_dir = [es_BASEDIR '/../resources/'];


%-------------------------------------------------------------------------------
% STEP 1: Define the model's network structure
%-------------------------------------------------------------------------------

model_dir = [resource_dir '/models/demo_e_coli_glycolysis/'];

%-------------------------------------------------------------------------------
% Define data structure 'network' for network structure

% Load network model from an SBtab file

if 1,
  network = sbtab_to_network([model_dir 'demo_e_coli_glycolysis_model.tsv'],struct('one_sbtab_file', 1));
end

% OR: Load network structure from SBML file

if 0,
  %% NOTE that KEGG annotations are not read from SBML file
  %% therefore, further below, none of the imported
  %% data will be mapped to the model!
  network = network_sbml_import([model_dir 'demo_e_coli_glycolysis_model.xml']);
end

% OR: Load a previously constructed network structure from a .mat file

if 0,
  M = load([model_dir 'demo_e_coli_glycolysis_model.mat'],'network'); network=M.network;
end

% OR: Directly construct a network structure within matlab by using the function 'network_construct'
% (see "help network_construct")

%-------------------------------------------------------------------------------
% Read network graphics information from a file (SBtab format)
% Network graphics information can be edited using functions from the mnt toolbox

network = netgraph_read_positions(network,[model_dir 'demo_e_coli_glycolysis_Positions.tsv']);
 
%-------------------------------------------------------------------------------
% Some of the following commands requires the stoichiometric matrix N, the allosteric 
% regulation matrix W, and the indices of external metabolites as direct inputs

% We extract these variables from the network structure.
 
N       = network.N;                  % stoichiometric matrix
W       = network.regulation_matrix;  % allosteric regulation matrix 
ind_ext = find(network.external);     % indices of external metabolites

%-------------------------------------------------------------------------------
% Set a target reaction for the following analyses

target_reaction    = 'ATPase'; 

%-------------------------------------------------------------------------------
% STEP 2: Define a flux distribution
%-------------------------------------------------------------------------------

% ------------------------------------------------------------------------------
% Compute a flux distribution (which need not be stationary, but thermodynamically feasible

% Use some (approximately known) fluxes, as well as standard deviations and possibly signs, and determine
% stationary fluxes resembling them; In this example, we only give a flux data value for the very first reaction
% Instead of using 'guess_flux_std', standard deviations from flux measurements can be used

if 1,
  v_mean = [1; nan;  nan;  nan;  nan;  nan;  nan;  nan;  nan;  nan;  nan;  nan;  nan];
  v_sign = [1; nan;  nan;  nan;  nan;  nan;  nan;  nan;  nan;  nan;  nan;  nan;  nan];
  v_std  = guess_flux_std(v_mean,0.1,0.1); 
  v      = project_fluxes(network.N, find(network.external), v_mean, v_std, v_sign);
end

% OR: Use Flux balance analysis (FBA), possibly followed by FBA with minimal sum of fluxes

if 0,
  fba_constraints    = fba_default_options(network);
  fba_constraints.zv = zeros(13,1);
  fba_constraints.zv(label_names(target_reaction,network.actions)) = 1;
  
  v = fba(network,fba_constraints);  
  %% DOES NOT WORK HERE !!!v     = pmf(network,fba_constraints,1,v_fba);
end

% OR: Compute a stationary flux distribution from the right kernel matrix
%   This is a pragmatic, and usually not very realistic choice
%   Compute a right-kernel matrix K and, for this example, choose its
%   first column as our flux mode v (the first flux is supposed to be positive)

if 0,
  K = full(network_analyse(network));
  v = K(:,1) * sign(K(1));              
end


% ------------------------------------------------------------------------------
% If necessary, modify the flux distribution such that it becomes thermodynamically feasible ("loopless")

v = eba_make_feasible(v,network.N);
 
%-------------------------------------------------------------------------------
% STEP 3: Define thermodynamically feasible concentrations, equilibrium constants,
%         and chemical potentials
%-------------------------------------------------------------------------------

% For later use, we load some model-related data and map them onto the model

resource_dir = [es_BASEDIR '../resources/data/'];

data_files = {[resource_dir 'data-thermodynamic/E_coli_glycolysis_GFE.tsv'], ...
              [resource_dir 'data-thermodynamic/GFE.tsv'], ...
              [resource_dir 'data-thermodynamic/ATPase_affinity.tsv'], ...
              [resource_dir 'data-metabolite/E_coli_glucose_concentration.tsv']}';

data_quantities = {'standard chemical potential','equilibrium constant', 'Michaelis constant', 'concentration','reaction affinity'}';
model_data      = kinetic_data_load(data_quantities, [], network, data_files, 0, 1, 1, 0);
% kinetic_data_print(model_data,network,1);

% There are various possibilities to do this: here, we see an estimation based on existing
% experimental data (concentrations, Gibbs free energies of formation, and reaction affinities)
% The estimation is performed by parmater balancing

if 1,
  pb_options = struct('A_min', 3, 'A_max', 20, 'conc_min', 0.00001, 'conc_max', 100, 'sigma_mu0',3, 'data_refer_to_molar',0);
  [c, mu0, Keq, A, kinetic_data] = parameter_balancing_thermodynamic(network, v, model_data, pb_options);
end

% An alternative matlab function for the same purpose is called thermo_pb.m
% It does not rely on the data structure, but allows for providing all function arguments directly

if 0,  
  my_c               = model_data.c.median;
  my_c(isnan(my_c))  = 1; % Guess of typical concentration 1 mM of non-measured metabolites
  mydG0              = -log(model_data.Keq.median);
  mopt.c_min_default = 0.0000001; % mM
  mopt.c_max_default = 200;       % mM
  thermo_pb_options  = mdf_default_options(network,mopt);
  thermo_pb_options.c_min(thermo_pb_options.c_min<0.0001) = 0.0001;
  thermo_pb_options.c_max(thermo_pb_options.c_max>10)     = 10;
  thermo_pb_options.delta_G0                              = mydG0;
  thermo_pb_options.delta_G0_std                          = ones(size(mydG0));
  %% thresholding is necessary for quadratic problem in thermo_pb
  thermo_pb_options.delta_G0_std(~isfinite(thermo_pb_options.delta_G0_std)) = 5;
  thermo_pb_options.delta_G0_std(thermo_pb_options.delta_G0_std<1)  = 1;
  thermo_pb_options.delta_G0_std(thermo_pb_options.delta_G0_std>10) = 10;
  thermo_pb_options.log_c                = log(my_c);
  thermo_pb_options.log_c_std            = log(1.1)*ones(size(network.metabolites));
  thermo_pb_options.dG_threshold         = 0.5*RT;  

  [c, dG0, A] = thermo_pb(network.N, v, thermo_pb_options);
end

% Another possibility is to use the Max-Min driving force method, a heuristics 
% that tries to avoid small driving forces in any reaction

if 0,
  mopt.c_min_default = 0.0000001; % mM
  mopt.c_max_default = 200;       % mM
  mdf_options = mdf_default_options(network,mopt);

 [c, A] = mdf(network.N, v, [], mdf_options);
end

%-------------------------------------------------------------------------------
% STEP 4: Set options for elasticity sampling, and run it
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% User options for elasticity sampling are stored in two matlab structs,
% 'es_options' and 'es_constraints'.  
% 'es_options' contains general settings for the calculation workflow, 
% 'es_constraints' contains numerical details about the model.
% To learn more about these structs, try 'help es_default_options'

% Default settings can be obtained as follows:
 
[es_options, es_constraints] = es_default_options(N);
 
%-------------------------------------------------------------------------------
% Here we assume that the reference state is already defined by the previously 
% compute variables v, c, and dmu, and we insert them into the 'es_options' struct.
% The functions for elasticity sampling also provide options for computing v, c, and dmu
% directly, based on other types of information. 
% For further information, please type "help es_sample_steady_state"

 %-------------------------------------------------------------------------------
% We insert the flux mode into 'es_constraints' and run elasticity sampling 
% All output is stored in the data structure 'result'
 
es_options.sampling_method   = 'accept_flux';
es_constraints.v_fix         = v;
es_constraints.log_c_fix     = log(c); 
es_constraints.dmu_fix       = -A; 
es_constraints.dmu_limit_min = pb_options.A_min;

% If desired: use known (or assumed) KM values to obtain saturation values
if 1,
  estimated_KM_values = full(abs(network.N')) * diag(c);
  estimated_KM_values(find(isfinite(model_data.KM.median))) = model_data.KM.median(find(isfinite(model_data.KM.median)));
  [alpha_M,beta_M] = k_to_alpha(estimated_KM_values,c);
  if 0,
    es_constraints.beta_M_fix = beta_M;
  else
    %% Alternatively, if desired: sample beta-distributed saturation values
    es_options.set_alpha_nonuniform = 1;
    es_constraints.alpha_M_mean = beta_M;
    es_constraints.alpha_M_std  = 0.1 * beta_M; 
  end
end

% Second order effects:
%  - either all second order coefficients are computed (which can take long for large models:

if 0, 
  es_options.flag_second_order = 1;
end

%  - or second order coefficients are computed only for one output function (a linear function 
%    of concentrations and fluxes); In this example, we define only the linear coefficients 
%    for fluxes, namely, a coefficient of 1 for our target reaction (previously defined), and 0 otherwise

if 1, 
  es_options.zv                = zeros(size(network.actions)); % define target function for second order coefficients
  es_options.zv(label_names(target_reaction, network.actions))  = 1; 
end

% ----------------------------------------------------------
% Now, run elasticity sampling

result = es_sample_model(N, W, ind_ext, es_constraints, es_options);

% Alternatively, the function can also be directly applied to the fields of 'network' 
% result = es_sample_model(network.N, network.regulation_matrix, find(network.external), es_constraints, es_options);


%-------------------------------------------------------------------------------
% STEP 5: Display results of elasticity smapling
%-------------------------------------------------------------------------------

% Show graphics without saving .eps files
if 1, 
  es_options.print_graphics = 0;
  psfile_dir = [];
end

% Show graphics with saving .eps files
if 0,
  es_options.print_graphics = 0;
  psfile_dir = MY_PSFILE_DIRECTORY; % Please insert your directory path
end

% The mnt toolbox allows for defining a second network structure 'network_CoHid', 
% which contains the graphics information for network graphics with some metabolites
% (e.g. cofactors) and reactions hidden. In this example, we simply use 'network' again

network_CoHid = network; 
es_filenames  = struct; 

% Show graphics for reference state

es_reference_state_graphics(network, es_options, result, target_reaction, es_filenames, network_CoHid, psfile_dir);

% Show graphics for interaction effects (second order control analysis)
% Note that this is not possible after es_sample_model has been run with the option es_options.flag_second_order = 0; 

interaction_analysis_graphics(network, network_CoHid, result, es_filenames, target_reaction, psfile_dir);


%-------------------------------------------------------------------------------
% STEP 6: Save results of elasticity smapling and the reconstructed kinetic model
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% The matlab struct 'result' contains, among other types of information,  a field 
% 'result.kinetics'. containing the kinetic rate laws obtained from elasticity sampling
% By inserting this field into 'network', the network becomes a kinetic model.
 
network.kinetics = result.kinetics;
 
% The model can now be simulated with functions from the Metabolic Network Toolbox
 
% --------------------------------------------------------------------
% Exporting the model to files

% Export the model (including rate laws) to SBML
%   Note that SBMLmodel yields a matlab data structure; 
%   If a filename (4th argument)is provided, an SBML file is written to disk

if 0,
  SBMLmodel = network_sbml_export(network,1, 'MY MODEL NAME', MY_SBML_EXPORT_FILE_NAME);
end

% OR: Export the model (including rate laws) to SBtab

if 0,
  network_to_sbtab(network,struct('filename', MY_SBTAB_EXPORT_FILE_NAME));
end

% OR: Export the model (including rate laws) to a .mat file

if 0,
  save(MY_MAT_EXPORT_FILE_NAME,'network','result','fluxes','es_options','es_constraints', 'target_reaction');
end

%-------------------------------------------------------------------------------
% STEP 7: Create model ensembles and test their results 
%         for statistically significant differences 
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Now we try multiple elasticity sampling runs. We first modify the settings
 
es_options.n_samples         = 10;
es_options.set_alpha_to_half = 0;
 
% We will sample the interactive effects of enzyme pairs on a specific output function.
% As an output function, we choose the sum of squares of all flux control coefficients.
 
output_function  = inline('sum(result.control.CJ .^2)','result','other');
 
% Now we will create two ensembles, one with cs kinetics and one with ms kinetics
 
%-------------------------------------------------------------------------------
% First ensemble: common saturable modular rate law ('cs')
 
es_options.kinetic_law = 'cs';
 
output1 = es_sample_multiple(N, W, ind_ext, es_constraints, es_options, output_function);
 
% This was the first ensemble
 
%-------------------------------------------------------------------------------
% Second ensemble:  simultaneous binding modular rate law ('ms')
 
es_options.kinetic_law = 'ms';
 
output2 = es_sample_multiple(N, W, ind_ext, es_constraints, es_options, output_function);

% This was the second ensemble
 
%-------------------------------------------------------------------------------
% Now we try to determine significant differences between both ensembles
% We choose a false discovery rate of 5 percent.
 
false_discovery_rate = 0.05;
 
% The result of the comparison are stored in a data structure result
 
result = es_compare_ensembles(output1, output2, false_discovery_rate);
