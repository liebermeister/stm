function es_constraints_show(es_constraints,network,pars)

% es_constraints_show(es_constraints,network,pars)

eval(default('pars','struct'));
[nm,nr] = size(network.N);

es_constraints.v_sign(find(isfinite(es_constraints.v_fix))) = ...
    sign(es_constraints.v_fix(find(isfinite(es_constraints.v_fix))));

pars.arrowstyle  = 'fluxes';
pars.arrowvalues = es_constraints.v_sign;
pars.arrowvalues(isnan(pars.arrowvalues)) = 0;

pars.actvalues   = nan*es_constraints.v_fix;
pars.actvalues(find(isfinite(es_constraints.v_fix))) = 1;
pars.actvalues(find(es_constraints.v_fix==0))        = -1;

pars.metvalues = es_constraints.ext_sign;
pars.colorbar=0;
clf; netgraph_concentrations(network,[],[],1,pars); 
