function r = stm_sample_model_write_flux_task(network,stm_constraints,stm_options,filename)

% STM_SAMPLE_MODEL_WRITE_FLUX_TASK - Write task file for computing a feasible flux distribution with external software
%
% r = stm_sample_model_write_flux_task(network, stm_constraints, stm_options, filename)
%
% Analogous to stm_sample_model; the struct r (with the flux and thermodynamic constraints)
% is saved to the .mat file 'filename' (optional). The solution can be read with
% stm_sample_model_read_flux_task

eval(default('filename','[]'));
nx = sum(network.external);
nn = network_add_transport(network, network.metabolites(find(network.external)),1);
[nm,nr] = size(nn.N);
external_prod_sign = sign(network.N(find(network.external),:)*stm_constraints.v_mean);

r.rxns                             = nn.actions;
r.mets                             = nn.metabolites;
r.S                                = sparse(nn.N);
r.lb                               = nan * ones(nr,1);
r.ub                               = nan * ones(nr,1);
r.lb(1:end-nx)                     = stm_constraints.v_min;
r.ub(1:end-nx)                     = stm_constraints.v_max;
r.lb(end-nx+find(external_prod_sign ==-1)) = -inf;
r.ub(end-nx+find(external_prod_sign ==-1)) = inf;
r.lb(end-nx+find(external_prod_sign ~=-1)) = -inf;
r.ub(end-nx+find(external_prod_sign ~=-1)) = 0;
r.lb(stm_constraints.v_sign>0)         = 0;
r.ub(stm_constraints.v_sign<0)         = 0;
r.lb(isfinite(stm_constraints.v_fix))  = stm_constraints.v_fix(isfinite(stm_constraints.v_fix));
r.ub(isfinite(stm_constraints.v_fix))  = stm_constraints.v_fix(isfinite(stm_constraints.v_fix));
r.c_lb                             = exp(stm_constraints.log_c_min);
r.c_ub                             = exp(stm_constraints.log_c_max);
%r.gfe_formation                    = stm_constraints.mu0_fix;
r.gfe_formation_lb                 = stm_constraints.mu0_min;
r.gfe_formation_ub                 = stm_constraints.mu0_max;
r.gfe_reaction_lb                  = stm_constraints.dmu_min;
r.gfe_reaction_ub                  = stm_constraints.dmu_max;
r.flux_data_mean                   = nan * ones(nr,1);
r.flux_data_std                    = nan * ones(nr,1);
r.flux_data_mean(1:end-nx)         = stm_constraints.v_mean;
r.flux_data_std(1:end-nx)          = stm_constraints.v_std;

r.gibbs_threshold                  = 0.1; % kJ/mol, necessary for driving a reaction

if length(filename)
  save(filename, 'r');
end
