function res = interaction_statistics(es_filenames, target_reaction, condition1, condition2, input_file1, input_file2, n_mc, n_per, fdr, fluxes_fixed)

%INTERACTION_STATISTICS - Run elasticity sampling repeatedly and compute statistics over the influence values
%
% res = interaction_statistics(es_filenames, target_reaction, condition1, condition2, input_file1, input_file2, n_mc, n_per, fdr, fluxes_fixed)
%
% Run elasticity sampling repeatedly and compute statistics 
% over the resulting influence values (1st order control coefficients) 
% and interaction values (i.e. 2nd order control coefficients) 
% 
% Output:
% res.p_influence       p-values for influences
% res.p_interaction     p-values for interactions
% res.influence_sig     
% res.influence_a_sig   
% res.influence_b_sig   
% res.interaction_sig   
% res.interaction_a_sig 
% res.interaction_b_sig


% ---------------------------------------------
% result files contain data structures 'es_options', 'es_constraints', 'result'
% predefined flux distribution is contained in 'es_constraints.vfix'

load(es_filenames.network_file);

N          = network.N;
W          = network.regulation_matrix;
ind_ext    = find(network.external);
[nm,nr]    = size(N);
ind_target = find(strcmp(target_reaction, network.actions));


% ---------------------------------------------

cd(es_filenames.es_dir); load(input_file1);
es_options.set_alpha_to_half  = 0;
es_options.seed               = 0;
es_options.epsilon_stationary = 0.001;

if fluxes_fixed,
  es_options.method_flux_sampling = 'accept_flux';
  es_constraints.v_fix = result.v;
end

% ---------------------------------------------

for it = 1:n_mc,
  fprintf('MC sample %d/%d:',it,n_mc);
  result         = es_sample_model(N,W,ind_ext,es_constraints,es_options);
  R_target_u_sc  = result.control.RJu_sc(ind_target,:);
  R_target_uu_sc = squeeze(result.control.RJuu_sc(ind_target,:,:));
  influence_list(:,it,1)     = R_target_u_sc;
  interaction_list(:,:,it,1) = R_target_uu_sc - diag(diag(R_target_uu_sc));
end

cd(es_filenames.es_dir); 
load(input_file2); % es_options, result; from script 'ycm_elasticity_sampling'
es_options.set_alpha_to_half  = 0;
es_options.seed               = 0;
es_options.epsilon_stationary = 0.001;

if fluxes_fixed,
  es_options.method_flux_sampling = 'accept_flux';
  es_constraints.v_fix = result.v;
end

for it = 1:n_mc,
  fprintf('MC sample %d/%d:',it,n_mc);
  result         = es_sample_model(N,W,ind_ext,es_constraints,es_options);
  R_target_u_sc  = result.control.RJu_sc(ind_target,:);
  R_target_uu_sc = squeeze(result.control.RJuu_sc(ind_target,:,:));
  influence_list(:,it,2)     = R_target_u_sc;
  interaction_list(:,:,it,2) = R_target_uu_sc - diag(diag(R_target_uu_sc));
end


% --------------------------------------------------------------
% Detect significant values and differences in the interaction matrices 
% by a permutation test

ind_intervention = 1:nr;

p_influence   = influence_anova(    influence_list, n_per, fdr, ind_intervention);
p_interaction = interaction_anova(interaction_list, n_per, fdr, ind_intervention);


% --------------------------------------------------------------
% use significance for thresholding

influence_sig = nan*ones(nr,1);
influence_sig(ind_intervention) = abs(p_influence.mean_total_significant).*p_influence.mean_total;
interaction_sig = nan*ones(nr);
dummi = abs(p_interaction.mean_total_significant).*p_interaction.mean_total; dummi(isnan(dummi)) = 0; dummi = dummi + dummi';
interaction_sig(ind_intervention,ind_intervention) = dummi;

influence_a_sig = nan*ones(nr,1);
influence_a_sig(ind_intervention) =abs(p_influence.mean_a_significant).*p_influence.mean_a;
interaction_a_sig = nan*ones(nr);
dummi = abs(p_interaction.mean_a_significant).*p_interaction.mean_a; dummi(isnan(dummi)) = 0; dummi = dummi + dummi';
interaction_a_sig(ind_intervention,ind_intervention) = dummi;

influence_b_sig = nan*ones(nr,1);
influence_b_sig(ind_intervention) =abs(p_influence.mean_b_significant).*p_influence.mean_b;
interaction_b_sig = nan*ones(nr);
dummi = abs(p_interaction.mean_b_significant).*p_interaction.mean_b; dummi(isnan(dummi)) = 0; dummi = dummi + dummi';
interaction_b_sig(ind_intervention,ind_intervention) = dummi;

res.p_influence       = p_influence;
res.p_interaction     = p_interaction;
res.influence_sig     = influence_sig;
res.influence_a_sig   = influence_a_sig;
res.influence_b_sig   = influence_b_sig;
res.interaction_sig   = interaction_sig;
res.interaction_a_sig = interaction_a_sig;
res.interaction_b_sig = interaction_b_sig;

res.ind_intervention = ind_intervention;
