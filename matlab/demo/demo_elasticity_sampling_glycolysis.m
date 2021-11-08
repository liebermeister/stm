% demo_elasticity_sampling_glycolysis - Demo script for thermodynamically feasible elasticity sampling 

clear; clearvars -global

es_dependencies;

model_dir = [es_BASEDIR '/../resources/models/demo_e_coli_glycolysis'];

cd(model_dir);

echo on;
clc
%---------------------------------------------------------
% Elasticity sampling demo
%
% In this script, we load a metabolic network model, 
% run elasticity sampling with two types of rate laws, 
% and compare the results of the two model ensembles 
%---------------------------------------------------------
 
% Press key to continue

pause
clc
%-------------------------------------------------------------------------------
% Load network model from SBML file

network = network_sbml_import('demo_e_coli_glycolysis_model.xml')

% Possible alternatives: 
% 1. Load network model from SBtab file
%      nn = sbtab_load_network_model_and_data('demo_e_coli_glycolysis_model.tsv'); 
%      network = nn.network;
% 2. Define network structure within matlab by using the function 'network_construct'
 
% Press key to continue
 
pause
clc
%-------------------------------------------------------------------------------
% Read graphics information for the network model from file
%
network = netgraph_read_positions(network,'demo_e_coli_glycolysis_Positions.tsv');
 
% Press key to continue
 
pause
clc
%-------------------------------------------------------------------------------
% The matlab command for elasticity sampling requires the stoichiometric matrix N 
% and the allosteric regulation matrix W.

% We extract these matrices from the network structure.
 
N       = network.N;                  % stoichiometric matrix
 
W       = network.regulation_matrix;  % allosteric regulation matrix
 
ind_ext = find(network.external);     % indices of external metabolites
 
% Then, we compute a right-kernel matrix K and, for this example, choose its
% first column as our flux mode v (the first flux is supposed to be positive)
 
K = full(network_analyse(network));   % stationary fluxes
 
v = K(:,1) * sign(K(1));              
 
% Press key to continue
pause
clc
%-------------------------------------------------------------------------------
% User options for elasticity sampling are stored in two matlab structs,
% 'es_options' and 'es_constraints'. 
 
% 'es_options' contains general settings for the calculation workflow, 
% 'es_constraints' contains numerical details about the model.
 
% Default settings can be obtained as follows:
 
[es_options, es_constraints] = es_default_options(N);
 
% To learn more about these structs, try 'help es_default_options'
 
% Now press key to continue
pause
clc
%-------------------------------------------------------------------------------
% After inserting the flux mode into 'es_constraints', we run elasticity sampling 
% All output will be stored in a data structure 'result'
 
es_constraints.v_fix = v;
 
result = es_sample_model(N, W, ind_ext, es_constraints, es_options);
 
% Press key to continue
pause
clc
%-------------------------------------------------------------------------------
% The matlab struct 'result' contains, among other types of information,  a field 
% 'kinetics'. containing the kinetic rate laws obtained from elasticity sampling
 
result.kinetics
 
% With these rate laws, the network becomes a kinetic model.
 
network.kinetics = result.kinetics;
 
% It can be simulated with functions from the Metabolic Network Toolbox
 
% Press any key to continue
pause
clc
% --------------------------------------------------------------------
% We can also export the model, including the rate laws, to SBML format

SBMLmodel = network_sbml_export(network,1);
 
% Press key to continue
pause
clc
%-------------------------------------------------------------------------------
% Now we try multiple elasticity sampling runs. We first modify the settings
 
es_options.n_samples         = 10;
 
es_options.set_alpha_to_half = 0;
 
% We will sample the interactive effects of enzyme pairs on a specific output function.
% As an output function, we choose the sum of squares of all flux control coefficients.
 
output_function  = inline('sum(result.control.CJ .^2)','result','other');
 
% Now we will create two ensembles, one with cs kinetics and one with ms kinetics
 
% Press key to continue
pause
clc
%-------------------------------------------------------------------------------
% First ensemble: common saturable modular rate law ('cs')
 
es_options.kinetic_law = 'cs';
 
output1 = es_sample_ensemble(N, W, ind_ext, es_constraints, es_options, output_function);
 
% This was the first ensemble
 
% Press key to continue
pause
clc
%-------------------------------------------------------------------------------
% Second ensemble:  simultaneous binding modular rate law ('ms')
 
es_options.kinetic_law = 'ms';
 
output2 = es_sample_ensemble(N, W, ind_ext, es_constraints, es_options, output_function);

% This was the second ensemble
 
% Press key to continue
pause
clc
%-------------------------------------------------------------------------------
% Now we try to determine significant differences between both ensembles
% We choose a false discovery rate of 5 percent.
 
false_discovery_rate = 0.05;
 
% The result of the comparison are stored in a data structure res
 
res = es_compare_ensembles(output1, output2, false_discovery_rate);
 
% All variables are still in your workspace
 
% Enjoy working with elasticity sampling!
 
% Press key to finish
pause
return