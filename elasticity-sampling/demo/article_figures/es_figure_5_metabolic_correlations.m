% Prediction of static metabolite fluctuations for different types of 
% variations (variable 'scenario_name') in a single reference state 
%
% See Article Figure 6

% -----------------------------------------------------------------------------
% load variables network, network_CoHid, v, mu, reference_state
% ------------------------------------------------------------

model_dir  = [es_BASEDIR '/../resources/models-article/Hepatonet_CCM_Fluxes'];
model_name = 'Hepatonet_Aerobic_ATP_Regeneration_model'; 
ref_state  = 'Hepatonet_Aerobic_ATP_Regeneration_reference_state'; 

cd(model_dir); load(model_name); load(ref_state);

% For clarity, omit NAD compounds from graphics 
network_CoHid.graphics_par.omitmetabolites       = {'NAD+ [c]', 'NAD+ [m]', 'NADH [c]', 'NADH [m]'};

% ------------------------------------------------------------
% Choose model and scenario for variations
% ------------------------------------------------------------

model_name    = 'hepatonet1_atp_glucose_o2';
scenario_name = 'only_enzymes';

% Other possible choices: 
%  scenario_name = 'all_external';
%  scenario_name = 'one_external'; my_ext_metab = 'Glucose [s]';
%  (Variable 'my_ext_metab' must be set for scenario 'one_external')

condition_name = 'cs';

% ------------------------------------------------------------
% load model and precalculated reference state

[nm,nr] = size(network.N);
n_ext   = sum(network.external);

%% check stable steady state
if ~reference_state.control.stable, error('Unstable steady state'); end

c       = reference_state.c;
RSu_sc  = reference_state.control.RSu_sc;
RSuu_sc = reference_state.control.RSuu_sc;
RSx_sc  = reference_state.control.RSs_sc;
RSxx_sc = reference_state.control.RSss_sc;


% ------------------------------------------------------------
% initialise covariances

switch scenario_name
  
  case 'only_enzymes',
    Cov_log_u     = log(2)^2 * eye(nr); % enzyme levels
    Cov_log_x     = zeros(n_ext);     % external concentrations

  case 'all_external';
    Cov_log_u     = zeros(nr);    % enzyme levels
    Cov_log_x     = log(2)^2 * eye(n_ext); % external concentrations

  case 'one_external';
    ext_metab     = network.metabolites(find(network.external));
    ind_ext       = label_names({my_ext_metab}, ext_metab);
    Cov_log_u     = zeros(nr);    % enzyme levels
    Cov_log_x     = zeros(n_ext); % external concentrations
    Cov_log_x(ind_ext,ind_ext)= log(2)^2;

end


% ------------------------------------------------------------
% compute covariance matrices

[log_c_mean_1, log_c_std_1, log_c_cov_1,  log_c_mean_2, log_c_std_2, log_c_cov_2] = network_analytical_distribution(c, RSu_sc, RSuu_sc, Cov_log_u, RSx_sc, RSxx_sc, Cov_log_x);

% the variance of the external metabolites should be shown!


% ------------------------------------------------------------
% graphics

% show concentration covariances as matrix

figure(1); 
gp = struct('squaresize',0.03,'arrowsize',0.02,'actprintnames',0,'flag_edges',1,'FontSize',6,'linewidth',1);
netgraph_concentrations(network_CoHid,network.external,[],1,gp);

figure(2); clf
im(log_c_cov_1,[]); axis square; colormap(rb_colors); colorbar

% show concentration covariances on network (curved arows)
% omit values lower than threshold (5 percent of maximum value)

figure(3); clf
log_c_cov_1_thr = log_c_cov_1 / max(abs(log_c_cov_1(:)));
log_c_cov_1_thr(abs(log_c_cov_1_thr)<0.05) = 0;
netgraph_metabolite_interactions(network_CoHid, log_c_std_1, log_c_cov_1_thr, rb_colors, gp);

% same, for 2nd-order approximation

figure(4); clf;
log_c_cov_2_thr = log_c_cov_2 / max(abs(log_c_cov_2(:)));
log_c_cov_2_thr(abs(log_c_cov_2_thr)<0.05) = 0;
netgraph_metabolite_interactions(network_CoHid, log_c_std_2, log_c_cov_2_thr, rb_colors, gp);

% More positive than negative covariances?

figure(5); clf
log_c_cov_1_stats = triu(log_c_cov_1,1);
values = sort(log_c_cov_1_stats(log_c_cov_1_stats~=0));
plot([0,0],[1,length(values)],'k-'); hold on; 
plot(values,1:length(values)); hold off; axis tight
xlabel('Covariances for ln c, 1st order approximation')
ylabel('Cumulative distribution (for metabolite pairs)'); 

% --------------------------------------------------

figure(1); title('External metabolites');
figure(2); title('Log concentration covariance matrix');
figure(3); title(sprintf('Standard deviations and covariances of log. concentrations\n (1st order approximation)'))
figure(4); title(sprintf('Standard deviations and covariances of log. concentrations\n (2nd order approximation)'))
figure(5); title(sprintf('Log concentration covariances - cumulative distribution'))
