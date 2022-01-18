function r = es_sample_model_write_flux_task(network,es_constraints,es_options,filename)

% ES_SAMPLE_MODEL_WRITE_FLUX_TASK
%
% result = es_sample_model_write_flux_task(N, W, ind_ext, es_constraints, es_options)
%
% analogous to es_sample_model; write task file for naama to compute feasible flux distribution

if 0,
  addpath ~/projekte/stm/matlab/models/ 
  model_name = 'eccm'; condition_name = 'ishii_growth_07';  
  filenames       = stm_model_filenames(model_name,condition_name); 
  cd([filenames.es_dir]); load(filenames.model_data_file);   
  filename = '~/projekte/stm/matlab/models_for_naama/eccm_for_naama';
  r = es_sample_model_write_flux_task(network,es_constraints,es_options);
end

eval(default('filename','[]'));
nx = sum(network.external);
nn = network_add_transport(network, network.metabolites(find(network.external)),1);  
[nm,nr] = size(nn.N);
external_prod_sign = sign(network.N(find(network.external),:)*es_constraints.v_mean);

r.rxns                             = nn.actions;
r.mets                             = nn.metabolites;
r.S                                = sparse(nn.N);
r.lb                               = nan * ones(nr,1); 
r.ub                               = nan * ones(nr,1);
r.lb(1:end-nx)                     = es_constraints.v_min;
r.ub(1:end-nx)                     = es_constraints.v_max;
r.lb(end-nx+find(external_prod_sign ==-1)) = -inf;
r.ub(end-nx+find(external_prod_sign ==-1)) = inf;
r.lb(end-nx+find(external_prod_sign ~=-1)) = -inf;
r.ub(end-nx+find(external_prod_sign ~=-1)) = 0;
r.lb(es_constraints.v_sign>0)         = 0;
r.ub(es_constraints.v_sign<0)         = 0;
r.lb(isfinite(es_constraints.v_fix))  = es_constraints.v_fix(isfinite(es_constraints.v_fix)); 
r.ub(isfinite(es_constraints.v_fix))  = es_constraints.v_fix(isfinite(es_constraints.v_fix)); 
r.c_lb                             = exp(es_constraints.log_c_min);
r.c_ub                             = exp(es_constraints.log_c_max);
%r.gfe_formation                    = es_constraints.mu0_fix;
r.gfe_formation_lb                 = es_constraints.mu0_min;
r.gfe_formation_ub                 = es_constraints.mu0_max;
r.gfe_reaction_lb                  = es_constraints.dmu_min;
r.gfe_reaction_ub                  = es_constraints.dmu_max;
r.flux_data_mean                   = nan * ones(nr,1); 
r.flux_data_std                    = nan * ones(nr,1); 
r.flux_data_mean(1:end-nx)         = es_constraints.v_mean;
r.flux_data_std(1:end-nx)          = es_constraints.v_std;

r.gibbs_threshold                  = 0.1; % kJ/mol, necessary for driving a reaction

if length(filename),
  save(filename, 'r');
end
