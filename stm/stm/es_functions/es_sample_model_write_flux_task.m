function r = es_sample_model_write_flux_task(network,es_constraints,es_options,filename)

% ES_SAMPLE_MODEL_WRITE_FLUX_TASK - Write task file for computing a feasible flux distribution with external software
%
% r = es_sample_model_write_flux_task(network, es_constraints, es_options, filename)
%
% Analogous to es_sample_model; the struct r (with the flux and thermodynamic constraints)
% is saved to the .mat file 'filename' (optional). The solution can be read with
% es_sample_model_read_flux_task

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

if length(filename)
  save(filename, 'r');
end
