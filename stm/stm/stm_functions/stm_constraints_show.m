function stm_constraints_show(stm_constraints,network,pars)

% STM_CONSTRAINTS_SHOW  Display options (from struct stm_constraints)
%
% stm_constraints_show(stm_constraints,network,pars)

eval(default('pars','struct'));
[nm,nr] = size(network.N);

stm_constraints.v_sign(find(isfinite(stm_constraints.v_fix))) = ...
    sign(stm_constraints.v_fix(find(isfinite(stm_constraints.v_fix))));

pars.arrowstyle  = 'fluxes';
pars.arrowvalues = stm_constraints.v_sign;
pars.arrowvalues(isnan(pars.arrowvalues)) = 0;

pars.actvalues   = nan*stm_constraints.v_fix;
pars.actvalues(find(isfinite(stm_constraints.v_fix))) = 1;
pars.actvalues(find(stm_constraints.v_fix==0))        = -1;

if isfield(stm_constraints,'ext_signs')
  pars.metvalues = stm_constraints.ext_signs;
end

pars.colorbar = 0;
clf; netgraph_concentrations(network,[],[],1,pars);
