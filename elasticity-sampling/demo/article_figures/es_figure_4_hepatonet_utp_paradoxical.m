% Hepatonet elasticity sampling // analysis for UTP rephosphorylation
% Analysis of effective self-inhibition
% Example anaerobic UTP rephosphorylation with saturation values alpha = 0.5 in reference state

% ----------------------------------------------------------------------------------
% Load model 
% ----------------------------------------------------------------------------------

model_dir  = [es_BASEDIR '/../resources/models-article/Hepatonet_CCM_Fluxes'];
model_name = 'Hepatonet_Anaerobic_UTP_Regeneration'; 
ref_state  = 'Hepatonet_Anaerobic_UTP_Regeneration_reference_state'; 

cd(model_dir); load(model_name); load(ref_state);

network.kinetics = reference_state.kinetics;

fd_name      = flux_distribution_name;
N            = network.N;
W            = network.regulation_matrix;
ind_ext      = find(network.external);
target       = 'UTP [c]'; % TARGET VARIABLE: Production rate of (external) UTP
it_target    = label_names({target},network.metabolites);
it_v_target  = find(N(it_target,:));


% ----------------------------------------------------------------------------------
% Enzyme increase in the target reaction can be counterproductive!
% Test model with alpha = 0.5 
% ----------------------------------------------------------------------------------

nn          = network; 
nn.kinetics = reference_state.kinetics;

nn.kinetics.u(1) = 1;
[c_ss,v_ss] = network_steady_state(nn,nn.kinetics.c);
[t1, C1,dum,dum2,V1] = network_integrate(nn, c_ss, 10);
nn.kinetics.u(1) = 1.5;
[t2, C2,dum,dum2,V2] = network_integrate(nn, c_ss, 50);

figure(100); clf; subplot(2,1,1); set(gca,'fontsize',18);
plot([t1 - t1(end); t2],[V1(1,:), V2(1,:)],'Linewidth',2); hold on; axis tight
plot([0 0],[0,1.2*V1(1,1)],'k'); hold off; 
xlabel('Time (a.u.)'); ylabel('Rate (a.u.)'); 
title('Temporal response to downregulation');

uu_list = 0.5:0.05:2; clear vv_list
for itt = 1:length(uu_list),
  nn.kinetics.u(it_v_target) = uu_list(itt) .* network.kinetics.u(it_v_target);  
  [c,v] = network_steady_state(nn,c_ss);
  vv_list(itt) = nn.N(it_target,:) * v;
end

figure(101); clf; subplot(2,1,1); set(gca,'fontsize',18);
plot([1 1],[min(vv_list) max(vv_list)],'k-'); hold on 
plot(uu_list,vv_list,'Linewidth',2); hold on; plot(1,nn.N(it_target,:)*v_ss,'b.','Markersize',20); hold off; axis tight;  
xlabel('Relative enzyme level'); ylabel('Stationary rate'); 
title('Stationary response to activity change');


% ----------------------------------------------------------------------------------
% Again, with sampled alpha values
% ----------------------------------------------------------------------------------

zz           = network.N(it_target,:);
readout      = inline('zz * reference_state.control.RJu_un','reference_state','zz');

my_options                   = es_options;
my_options.n_samples         = 1000;
my_options.set_alpha_to_half = 0;
my_options.flag_second_order = 0;

output = es_sample_multiple(N, W, ind_ext, es_constraints, my_options, readout,[],[],zz);

R_prod_joint_upregulation = output(it_v_target,:)' * ones(length(it_v_target),1);

n_neg = sum(R_prod_joint_upregulation(:,1)<0);
n_pos = sum(R_prod_joint_upregulation(:,1)>0);

n_models    = my_options.n_samples
n_found     = n_neg;
p_neg_naive = n_neg/n_models;
p_neg_mean  = [n_neg+1] / [n_neg+n_pos+2];
p_neg_sigma = sqrt([ [n_neg+1]*[n_models - n_neg + 1]] / [ [n_models+2]^2*[n_models+3]]);

figure(102); clf;
subplot(2,1,1); set(gca,'fontsize',18);
plot([0 0],[0 1],'k-'); hold on; 
plot(sort(R_prod_joint_upregulation),[1:1000]/1000,'Linewidth',2);
xlabel(sprintf( 'Scaled control coefficient.  Prob(C<0)= %1.3f',p_neg_mean));
ylabel('Cumulative probability'); 
title('Role of elasticities');
axis([-2 1 0 1]);


% ----------------------------------------------------------------------------------
% Varying thermodynamic forces
% ----------------------------------------------------------------------------------

dmu             = es_constraints.dmu_fix;
dmu_factor_list = 10.^[-1.9:0.1:1.2];

RJu_A_ref = N(it_target,:) * reference_state.control.RJu_un(:,it_v_target) * ones(length(it_v_target),1);

clear RJu_A_list
for it = 1:length(dmu_factor_list),
  my_constraints = es_constraints;
  my_constraints.dmu_limit = inf;
  my_constraints.dmu_limit_min = 0;
  my_constraints.mu_fix = nan *   my_constraints.mu_fix;
  my_constraints.dmu_fix = dmu * dmu_factor_list(it);
  my_options = es_options;
  my_options.flag_second_order = 0;
  my_reference_state      = es_sample_model(N,W,ind_ext,my_constraints,my_options);
  RJu_A_list(it) = [N(it_target,:) * my_reference_state.control.RJu_un(:,it_v_target) * ones(length(it_v_target),1)]';
end

figure(103); clf;
subplot(2,1,1); set(gca,'fontsize',18); 
plot([min(dmu_factor_list), max(dmu_factor_list)],[0,0],'k-'); hold on; 
plot(dmu_factor_list, RJu_A_list,'Linewidth',2); axis tight;
plot(1,RJu_A_ref,'b.','Markersize',20); hold off
xlabel('Scaling of reaction affinities'); ylabel('Scaled control coefficient'); set(gca, 'XScale','log');
title(sprintf('Role of reversibility'));

figure(104); clf;
netgraph_concentrations(network_CoHid,[],[],1,struct('arrowstyle','fluxes','arrowvalues',v,'arrowsize',0.03,'omitmetabolites',{'H+(PG) [c]'}))
